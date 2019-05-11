#Cluster enrichment analysis with hypergeometric method
library("org.Hs.eg.db")
library("clusterProfiler")
set.seed(20)

data1<-read.table(file.choose(),header = TRUE, sep = "\t")
data1$Gene<-rownames(grph_data)

#convert IDs
Clustlist<- bitr(data1$Gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop=TRUE)
Clustlist = Clustlist$ENTREZID


geneList = data1$Sample # any numeric values, this is necessary for one of the tests

names(geneList) = as.character(Clustlist)
geneList = sort(geneList, decreasing = TRUE)

#GO enrichment
#gsecc <- gseGO(geneList=geneList, ont="BP", OrgDb=org.Hs.eg.db, verbose=F, pvalueCutoff = 1)
#gsecc <- gseKEGG(geneList=geneList, organism = "hsa", verbose=F, pvalueCutoff = 1)
gsecc <- enrichKEGG(gene=Clustlist, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)
gsecc1 <- enrichGO(gene=Clustlist, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5)

head(summary(gsecc))
write.table(summary(gsecc), file = "enrichKEGG.txt")
write.table(summary(gsecc1), file = "enrichGOBP.txt")