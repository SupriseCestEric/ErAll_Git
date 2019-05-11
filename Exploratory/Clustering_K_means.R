#Script for estimating optimal number of clusters for K-means analysis
#Plotting of within-group sum of squares, gap statistic and shilouette. 
#Based on blog: 2-bit bio


library("limma")
library("affy")
library("Biobase")
library("clusterProfiler")
library("cluster")
library("ComplexHeatmap")
library("reshape")
library("ggplot2")

#read in data
df1<-read.table(file.choose(), header = TRUE, sep = ";", row.names = 1)
df1<- as.matrix(df1)

#Differential expression analysis to filter out non-variable genes. 

###### FOR MICROARRAY######
#cut data by median expression of gene of interest
Triage<- df1[grep("Gene_of_interest", rownames(df1)), ]
median(Triage)
Triagedf <- as.data.frame(Triage)
res<- apply(Triagedf,2,function(x) ifelse(x>median(Triage),1,0))

###assigned the new dichotomized variable for UGT high or low to the phenodata
Pheno<- as.matrix(res)
APheno<- new("AnnotatedDataFrame", data = as.data.frame(Pheno))
## Build new Expressionset for analysis
StatsSet<- ExpressionSet(assayData = df1, phenoData = APheno)
Design = model.matrix(~StatsSet$Triage)

## fit the model
fit <- lmFit(StatsSet, Design)
## Improve Variance estimation with eBayes
eBfit <- eBayes(fit)
results <- topTable(eBfit, n=2500, adjust = "BH")
df1<-df1[rownames(df1) %in% rownames(results),]

###### FOR RNASEQ #########

library(DESeq2)
library(coseq)
Triage<- df1[grep("UGT2B17_n3", rownames(df1)), ]
median(Triage)
Triagedf <- as.data.frame(Triage)
res<- apply(Triagedf,2,function(x) ifelse(x>median(Triage),1,0))
res<-as.data.frame(res)
res$Triage<-as.factor(res$Triage)
df1<-round(df1,0)

colData <- res

dds <- DESeqDataSetFromMatrix(countData = df1,
                              colData = colData,
                              design = ~ Triage)

dds <- DESeq(dds, test="LRT", reduced = ~1)
dds_res <- results(dds,alpha=0.05)
summary(dds_res)

#view pval disrt
hist(dds_res$pvalue)

#subset data for clustering analysis, choose only padj < 0.05 and mean count > 50
alpha=0.5
df1_counts <- df1[which(dds_res$padj < alpha & 
                              rowMeans(counts(dds, norm=TRUE)>50)),]


###############################################################

#cutoff for clustering was 2500 genes to alleviate computations. 
df1<-df1[rownames(df1) %in% rownames(results[1:2500,]),]


#mean-center data
df1<-df1-rowMeans(df1)

#within group sum of squares
wss <- (nrow(df1)-1)*sum(apply(df1,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(df1,
                                     centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")



#gap statistic
gap <- clusGap(df1, kmeans, 20, B = 100, verbose = interactive())
plot(gap, main = "Gap statistic")
abline(v=which.max(gap$Tab[,3]), lty = 2)


#silhouette
sil <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(df1, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(df1))
  sil[i] <- mean(ss[, 3])
}

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)


#Perform K-means clustering. 

set.seed(20)
kClust <- kmeans(df1, centers=4, nstart = 1000, iter.max = 20) # replace k as needed
kClusters <- kClust$cluster



findclust<- kClusters[grep("UGT2B17", names(kClusters))]
findclust<- kClusters[kClusters == 2]
# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid,df1, kClusters)

#A posteriori verification of orthogonality of cluster centroids
cor(kClustcentroids)




#Subset the cores molten dataframe so we can plot the core
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')

#plot
p1 <- ggplot(Kmolten, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Sample") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Sample",color = "Cluster")
p1


ht = Heatmap(df1, km = 3, cluster_columns = TRUE) # replace km with actual clust number
rcl.list <- row_order(ht)





