#Set up annotated PPI network

library(plyr)
library(igraph)
library(org.Hs.eg.db)
library(GOstats)
library(biomaRt)

ppi_sources <- list(Bioplex='mmc3.txt', NCBI_data='human_ppiset.txt')
choice <- 1 #1 for Bioplex (~20,000, all TAP-MS), 2 for GeneRIF interactions set (~200,000, many sources)

setwd(paste0('~/Project/PPI/', names(ppi_sources)[choice]))

ppi <- read.delim(ppi_sources[[choice]]) #List of human PPIs
colnames(ppi)[1:2] <- c('geneA_id', 'geneB_id')
TG_genes <- read.delim('../../TG_genes.txt') #Larry's list of 416 genes related to TG, with annotations

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org") #Set up biomart object
TG <- getBM(attributes=c('hgnc_symbol','entrezgene'), filters='hgnc_symbol', values=TG_genes$Gene, mart=ensembl) #Get EntrezIDs
TG_notinBM <- setdiff(TG_genes$Gene, TG$hgnc_symbol) #Gene symbols not found in biomart


##### GRAPH CREATION #####

TG_links <- ppi[ppi$geneA_id %in% TG$entrezgene | ppi$geneB_id %in% TG$entrezgene,] #Filter PPI list for TG genes
TG_withPartners <- unique(append(TG_links$geneA_id, TG_links$geneB_id)) #Vector of TG genes and partners
TG_notInDB <- setdiff(TG$entrezgene, TG_withPartners) #Genes not found in dataset (EntrezIDs)
TG_notInDB <- TG[TG$entrezgene %in% TG_notInDB, 'hgnc_symbol'] #Convert to symbols

TG_network <- ppi[ppi$geneA_id %in% TG_withPartners | ppi$geneB_id %in% TG_withPartners,] #Get tertiary nodes

edgelist <- as.matrix(TG_network[,c('geneA_id','geneB_id')]) #Edgelist must be a matrix
class(edgelist) <- 'character' #Graph construction is wonky with integers as nodes

TG_graph <- graph_from_edgelist(edgelist, directed=F) #Create undirected graph from edge list
TG_graph <- igraph::simplify(TG_graph) #Removes multiple edges ('igraph::' clarifies which 'simplify' method to use)


##### CLUSTERING #####

communities = list()
communities$fast_greedy <- cluster_fast_greedy(TG_graph)
communities$walktrap <- cluster_walktrap(TG_graph)
communities$leading_eigenvector <- cluster_leading_eigen(TG_graph)

com <- communities[['walktrap']] #Choice of clustering algorithm


##### PLOT #####

#Series of layouts
E(TG_graph)$weight <- ifelse(crossing(com, TG_graph), 0.01, 10) ###BETTER WAY TO WEIGHT EDGES HIGHER WITHIN CLUSTERS?

coords_circle <- layout_(TG_graph, in_circle())
coords_kk <- layout_(TG_graph, with_kk(kkconst=2000)) ##kk may be inappropriate for this size graph
coords_fr <- layout_(TG_graph, with_fr(dim=3))#(area = 10*vcount(TG_graph)^2))

#Set attributes
V(TG_graph)$color <- membership(com)
V(TG_graph)$size <- map(degree(TG_graph), c(1,4)) ##Use another measure for size reasons? Somehow include a prior towards a medium size?
V(TG_graph)$label.cex <- ifelse(V(TG_graph)$size > 30, 1, 0.01)
V(TG_graph)$label <- NA

plot(TG_graph, layout = coords_fr)
#legend("bottomright", legend=c("TG","non-TG"), col=c("blue", "red"), lty=1, bty='n')


##### ANALYSIS #####

#Basic stats
print(paste0('Number of genes: ', length(V(TG_graph)))) #Number of genes (nodes)
print(paste0('Number of interactions: ', length(E(TG_graph)))) #Number of interactions (edges)
print(paste0('Number of clusters: ', length(unique(membership(com))))) #Number of clusters

#Centrality measures for each node
degrees <- as.data.frame(igraph::degree(TG_graph))
closenesses <- as.data.frame(closeness(TG_graph))
betweennesses <- as.data.frame(betweenness(TG_graph))

cluster_mem <- data.frame(gene = names(membership(com)), cluster = as.vector(membership(com)))
#cluster_mem <- cluster_mem[order(cluster_mem$cluster),c('cluster','gene')]

#Cluster sizes
cluster_sizes <- table(cluster_mem$cluster) #Table of cluster sizes
cluster_sizes <- sort(cluster_sizes, decreasing=T) #Sort table
plot(cluster_sizes) #Look at distribution of cluster sizes

#Create list of clusters
clusters <- lapply(names(cluster_sizes), function(x) cluster_mem[cluster_mem$cluster == x,]) #Better way to accomplish this?
threshold <- 10 #Threshold for acceptable cluster size for further analysis
clusters <- clusters[sapply(clusters, nrow) >= 10] #Filter for clusters with size above threshold


##### GO ENRICHMENT ANALYSIS #####

entrez_object <- org.Hs.egGO #Mapping between EntrezID and GO annotation
mapped_genes <- mappedkeys(entrez_object) #Get all genes mapped to GO IDs (to become the "universe")

GO_enrichment <- function(clust) {
  params <- new('GOHyperGParams', geneIds=as.character(clust$gene), universeGeneIds=mapped_genes, ontology='BP',
                pvalueCutoff=0.001, conditional=F, testDirection='over', annotation="org.Hs.eg.db")
  return(hyperGTest(params))
}

GOs <- lapply(clusters, GO_enrichment)


#### FURTHER ANALYSIS -- GETTING TOP (DRUGGABLE) GENE FROM EACH CLUSTER ####
druggable <- read.delim("../druggable_genes_1679_ProteinAtlas_20151113.txt")
druggable <- getBM(attributes='entrezgene', filters='hgnc_symbol', values=druggable$Gene, mart=ensembl)
druggable <- druggable$entrezgene

#df_mem <- data.frame(gene=names(membership(com)), cluster_no=as.vector(membership(com)))

best_druggable <- function(clust, deg) {
  deg <- deg[deg$gene %in% clust & deg$gene %in% druggable & deg$gene %in% TG$entrezgene,] #Filter to only druggable TG proteins (& deg$gene %in% TG$entrezgene)
  return (ifelse(nrow(deg) != 0, deg$gene[1], 'no_druggable'))
}

degrees$gene <- rownames(degrees) #For ease of use in best_druggable function
degrees <- degrees[order(degrees[1], decreasing=T),] #For ease of use in best_druggable function
genes_of_interest <- sapply(clusters, function(x) best_druggable(x$gene, degrees)) #Retrieve "best" gene for each cluster
genes_of_interest <- data.frame(cluster_id=seq(1,length(genes_of_interest)), best_gene=genes_of_interest)


##### WRITE DATA OF INTEREST #####

write.graph(TG_graph, "TG_graph.gml", format="gml")
write.table(TG_notInDB, "notInBioplex.txt", sep='\t', row.names=F, quote=F) ###Get symbols here, why are there NA??
write.table(degrees[degrees[1] > 50,], "highDegree_nodes.txt", sep='\t', row.names=F, quote=F)
write.table(genes_of_interest, 'genes_of_interest.txt', row.names=F, quote=F)

# Statistics:
# Number of genes: 4826
# Number of edges: 10806
# Number of clusters (using walktrap algorithm): 26
# Cluster size: [1,670], median = 8