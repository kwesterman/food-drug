#Use ChemmineR package to develop atom pair sets for a combined set of drug 
#compounds and FooDB food compounds and develop a set of all food compounds
#with Tanimoto similarity score >0.85 with respect to one of the drugs

#Read in combination SDF set for 272 drugs followed by all FooDB compounds
combo_sdfset <- read.SDFset('~/DirectedStudy/Objective3/SDF_retrieval/combo_SDF.txt')
#Read in list of drugs, Chembl IDs, SMILES structures, and gene targets
drug_info <- read.table("~/DirectedStudy/Objective3/chembl_20_mysql/Drug_info.txt", sep="\t", header=T, comment.char="")
#Read in list of food compound info and change formatting
food_info <- read.table("~/DirectedStudy/Objective3/foodb_mysql/foodcmpds_info.txt", sep='\t', header=T, quote="")
colnames(food_info) <- c("Foodcmpd_name", "FooDB_ID", "Cmpd_smiles")
food_info$index <- row.names(food_info)

#Remove invalid SDFs from food compounds list
valid <- validSDF(combo_sdfset)
combo_sdfset <- combo_sdfset[valid]

#Create AP set for combo set
combo_apset <- sdf2ap(combo_sdfset)
failedAPs <- which(sapply(as(combo_apset, "list"), length)==1)

#Perform similarity search and append "drug of origin" info to each result
similarity <- data.frame()
for (i in 1:nrow(drug_info)) {
  search <- cmp.search(combo_apset, combo_apset[i], cutoff=0.85, type=3, quiet=T)
  search$Origin_drug <- drug_info$Drug_name[i]
  search$CHEMBL_ID <- drug_info$Chembl_ID[i]
  search$Gene_target <- drug_info$Gene[i]
  similarity <- rbind(similarity, search)
}

#Trim results to remove any results that are drugs (a.k.a. CMP1-CMP272), and reset indices so they correspond to the food_info data frame
similarity <- subset(similarity, similarity$index > nrow(drug_info), select=-cid)
similarity$index <- similarity$index - nrow(drug_info)
colnames(similarity)[2] <- "Tanimoto_score"

#Merge similarity results with food info to match results with compound names
results <- merge(similarity, food_info, by="index")
results <- results[order(results$Origin_drug),c(1,3:5,2,6:8)]

#Write final table
write.table(results[,-1], "~/DirectedStudy/Objective3/ChemmineR/Results.txt", col.names=T, row.names=F, sep="\t", quote=F)

#Diagnostics
self_return <- results[as.character(results$Origin_drug)==toupper(as.character(results$Foodcmpd_name)),]