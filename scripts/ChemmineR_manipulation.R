#Use ChemmineR package to develop atom pair sets for two sets of food compounds and 
#develop a set of all food compounds with Tanimoto similarity score >0.85

library(ChemmineR, parallel)

#Read in data
combo_sdfset <- read.SDFset('~/DirectedStudy/Objective3/SDF_retrieval/combo_SDF.txt')
list_A <- read.table("~/DirectedStudy/Objective3/setA_info.txt", sep="\t", header=T, comment.char="")
list_B <- read.table("~/DirectedStudy/Objective3/setB_info.txt", sep='\t', header=T, quote="")
colnames(list_B) <- c("B_name", "dbB_ID", "B_structure")
list_B$index <- row.names(list_B) #Create index for list B

#Remove invalid SDFs from food compounds list
valid <- validSDF(combo_sdfset)
combo_sdfset <- combo_sdfset[valid]

#Create AP set for combo set
combo_apset <- sdf2ap(combo_sdfset)

#Function to take in a cutoff value, write a file containing similar cmpds, and return the resulting datasets 
similarity_search <- function (T_cutoff) {
  similarity <- data.frame()
  for (i in 1:nrow(list_A)) { #Loop through list A cmpds looking for list B cmpds above cutoff
    search <- cmp.search(combo_apset, combo_apset[i], cutoff=T_cutoff, type=3, quiet=T)
    search$Origin_cmpd <- list_A$A_name[i]
    search$dbA_ID <- list_A$dbA_ID[i]
    search$Gene_target <- list_A$Gene[i]
    similarity <- rbind(similarity, search) #Append "origin compound" info to each result
  }
  
  similarity <- subset(similarity, similarity$index > nrow(list_A), select=-cid) #Remove "self-hits" (a.k.a. CMP1-CMP272)
  similarity$index <- similarity$index - nrow(list_A) #Reset indices so they correspond to the list_A data frame
  colnames(similarity)[2] <- "Tanimoto_score"
  
  results <- merge(similarity, list_B, by="index") #Merge similarity and list_B tables to associate results with compound names
  results <- results[order(results$Origin_cmpd),c(1,3:5,2,6:8)] #Rearrangement
  results <- results[results$Origin_cmpd!=toupper(results$B_name),] #Remove more "self-hits"
  
  write.table(results[,-1], paste0('~/DirectedStudy/Objective3/ChemmineR/Results_', Tan_cutoff, '.txt'), col.names=T, row.names=F, sep="\t", quote=F)
  
  return(results)
}

#Create a list of Tanimoto cutoffs and send them through the similarity function
cutoff_list <- list(0.65, 0.75, 0.85, 0.95)
results_tables <- mclapply(cutoff_list, similarity_search) #Parallelizes lapply
names(results_tables) <- cutoff_list

#Diagnostics
results_65 <- results_tables[[1]]; results_65$Foodcmpd_name <- as.character(results_65$B_name); table(results65$B_name) #Freq. table of food cmpds
num_failedSDFs <- length(valid) - length(combo_sdfset)
failedAPs <- which(sapply(as(combo_apset, "list"), length)==1)