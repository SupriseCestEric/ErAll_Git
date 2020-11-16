

#Dependencies
library("FactoMineR")
library("edgeR")
library("ggplot2")
library("Biobase")
library("limma")
library("gplots")
library("dplyr")
library("dendextend")
library("plotly")
library(dplyr)
library(viper)
library(dorothea)
library(gprofiler2)
library(tidyr)
library("genefilter")


set.seed(234)
#Read in data
count_tab <- read.table("counts_allsamples.txt", header = T, sep = "\t",stringsAsFactors = F)

#format ENSG IDS
df <- count_tab %>% separate(Geneid, c("Gene","version"),"[.]")
df <- df[!(df$Gene %in% df$Gene[duplicated(df$Gene)]),]
rownames(df) <- df$Gene

#This is for ICGC data
# count_tab <- read.table("icgc_wide.txt", header = T, sep = "\t",stringsAsFactors = F)
# #format ENSG IDS
# df <- count_tab %>% separate(gene_id, c("Gene","version"),"[.]")
# df <- df[!(df$Gene %in% df$Gene[duplicated(df$Gene)]),]
# rownames(df) <- df$Gene
# 
# #get only counts
# raw_counts <- as.matrix(df[,3:296])

#get only counts
raw_counts <- as.matrix(df[,8:107]) # select only numeric rows

#low-count filter
f_counts <- raw_counts[rowMeans(raw_counts)>10,] # low-expression filter

#get Gene of interest
UGT2B17 <- f_counts[grep("ENSG00000197888", rownames(f_counts)),] # get vector of gene expression for UGT2B17

#convert IDs to HGNC
gc1 <- gconvert(query = rownames(f_counts), organism = "hsapiens",
         target="HGNC", mthreshold = Inf, filter_na = TRUE) # use g:profiler API to convert IDs

#format rownames to HGNC
length(gc1$target[duplicated(gc1$target)]) #see if there are any duplicates
gc1 <- gc1[!(gc1$target %in% duplicated(gc1$target)),]
gc1 <- gc1[!(gc1$input %in% duplicated(gc1$input)),]

f_counts <- f_counts[rownames(f_counts) %in% gc1$input,]
f_counts <- f_counts[order(rownames(f_counts)),]
gc1 <- gc1[gc1$input %in% rownames(f_counts),]
gc1 <- gc1[order(gc1$input),] # check for other dupplicated IDs using duplicated(), some may still remain... if so then remove
findT<- duplicated(gc1$input)
which(findT %in% "TRUE") # find index of duplicated IDs
gc1 <- gc1[-c(18047,21736),] # remove them

rownames(f_counts) <- gc1$target #change rownmames to HGNC


# Get high confidence regulons from Dorothea DB 
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

#calculate TF activities from expression matrix
tf_activities <- run_viper(f_counts, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))

write.table(tf_activities, file = "LL100_dorothea_TF_activity.txt", sep = "\t", quote = F)

subj <- tf_activities[rownames(tf_activities)=="IRF1",]
cor.test(UGT2B17, subj ,method = "spearman")

#visualize PCA
svd1<-svd(t(tf_activities))
plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",xlab="1st PC", main="clustering of LL100 cell lines")


#Plot TF activities in UGT2B17-high cells for TFs of interest
goi <- c("STAT3","BCL6","IRF1","PAX5", "NFKB1","NFKB2","RELA","RELB") #some TFs of interest
#Plot TF activities in UGT2B17-high cells for significant results
res <- rowTtest(tf_activities[,names(UGT2B17[UGT2B17>55])], tf_activities[,names(UGT2B17[UGT2B17<55])]) #do a row T-Test based on the median-split UGT2B17 expression vector, or other relevant threshold
FDR <- p.adjust(res$p.value, method = "BH")
names(FDR) <- rownames(res$p.value) # fix P-value to FDR


#Melt data for ggplot
mtf <- melt(tf_activities)
names(mtf)<- c("TF","CL","Activity_score")
m_long <- ifelse(mtf$CL %in% names(UGT2B17[UGT2B17>55]),"UGT2B17-High","UGT2B17-Low") # create dummy variable for High / Low expression
mtf$UGT2B17 <- m_long
mtf$UGT2B17 <- as.factor(mtf$UGT2B17)
mtf$CL <- as.factor(mtf$CL)
mtf$TF <- as.factor(mtf$TF)

#subset the data for each of the desired plots
subset <- mtf[mtf$TF %in% goi,]
subset1 <- mtf[mtf$TF %in% names(FDR[FDR<0.05]),]

# plot graph
ggplot(subset,aes(x=TF,y=Activity_score, fill = UGT2B17)) + geom_boxplot(position="dodge") + xlab("Transcription Factor") + ylab("TF activity score") + theme_bw()
ggplot(subset1,aes(x=TF,y=Activity_score, fill = UGT2B17)) + geom_boxplot(position="dodge") + xlab("Transcription Factor") + ylab("TF activity score") + theme_bw()

#for Heatmaps (This also shows how to access TF target info)
IRF1 <- regulons$target[regulons$tf == "IRF1"]
STAT3 <- regulons$target[regulons$tf == "STAT3"]
NFKB1 <- regulons$target[regulons$tf == "NFKB1"]
ATF4 <- regulons$target[regulons$tf == "ATF4"]

#Normalize expression data to make it comparable between-sample
dge <- DGEList(f_counts)
dge <- calcNormFactors(dge, method = "TMM")
pseudo_TMM <- log2(scale(dge$counts,center=FALSE,scale=dge$samples$norm.factors)+1)
pseudo_counts_cen <-pseudo_TMM - rowMeans(pseudo_TMM)



#plot heatmaps
HM1<-pseudo_counts_cen[rownames(pseudo_counts_cen) %in% ATF4,] # modify for whatever TF
HM2 <- cbind(rowMeans(HM1[,colnames(HM1) %in% names(UGT2B17[UGT2B17>55])]), rowMeans(HM1[,colnames(HM1) %in% names(UGT2B17[UGT2B17<55])]))
colnames(HM2) <- c("UGT2B17-high", "UGT2B17-low")
Heatmap(HM2, heatmap_legend_param = list(title = "Centered log2 Expression"))


#adjust signage on target genes for each TF based on TF activities for that gene (-1/1)
signs<- regulons[regulons$tf == "ATF4",]
signs<- signs[order(signs$target),]
signs<- signs[signs$target %in% rownames(HM1),]
HM1 <- HM1[order(rownames(HM1)),]
graph_pattern<-HM1*signs$mor

### Some additionnal exploratory plots
#Formatting to create a plot showing:
# mean +- SD of target gene expression per cell line and activity scores. Then overlay target gene expression
HM1<-pseudo_counts_cen[rownames(pseudo_counts_cen) %in% ATF4,] #ATF4 can be swapped out for any TF, or could be looped over
signs<- regulons[regulons$tf == "ATF4",] #swap STF here as well, and at line 138
signs<- signs[order(signs$target),]
signs<- signs[signs$target %in% rownames(HM1),]

HM1 <- HM1[order(rownames(HM1)),]
graph_pattern<-HM1*signs$mor
gr <- rbind(colMeans(graph_pattern), apply(graph_pattern, 2, sd), pseudo_counts_cen[rownames(pseudo_counts_cen) == "UGT2B17",], tf_activities[rownames(tf_activities) == "ATF4",])
rownames(gr) <- c("Means", "SD", "UGT2B17","TFa")

# take only cells expressing UGT2B17, makes the plot smaller
gr <- gr[,gr[3,]>-2.37] #-2.37 here is the minima, so min(gr[3,]) could be used instead of -2.37. These are essentially samples with 0 counts
grt<- t(gr)
grt <- as.data.frame(grt)
grt$Cells_index <- seq(1:37)

#more graphs
# ON TOP: plot Mean expression +- SD of target genes with blue line indicating UGT2B17 expression 
# ON BOTTOM: plot TF activity of IRF1 calculated with Viper with blue line for UGT2B17 expression 
target_distro<- ggplot(grt, aes(x=Cells_index, y=Means))+geom_ribbon(aes(ymin=Means-SD, ymax=Means+SD), alpha = 0.1)+geom_line(aes(Cells_index, UGT2B17), color = "blue")+geom_line()+theme_bw()+xlab("")+ylab("Centered log-expression") + ylim(-12,12)
activity_distro<- ggplot(grt, aes(x=Cells_index, y=TFa))+geom_line(aes(Cells_index, UGT2B17), color = "blue")+geom_line()+ylim(-12,12)+theme_bw() + ylab("ATF4 Activity") + xlab("Sample index")
ggarrange(target_distro, activity_distro, nrow=2)

# Co-expression analysis of mean-centered data with UGt2B17by simple gene-wise K--means clustering
# DGE analysis to reduce the number of genes.....
# Just to lighten computations. 
factor1 <- factor(ifelse(UGT2B17>55,1,0))
factor1 <- factor(factor1, levels = c(0,1))
design_matrix <- model.matrix(~ factor1)

dge <- estimateGLMCommonDisp(dge, design_matrix)
dge <- estimateGLMTrendedDisp(dge, design_matrix)
dge <- estimateGLMTagwiseDisp(dge, design_matrix)

par(mfrow=c(1,1))
plotBCV(dge, main = paste0("BCV plot"))

fit<-glmFit(dge,design_matrix)
res <- glmLRT(fit, coef = 2)


pvals<- data.frame(raw.pvalue = res$table$PValue,
                   FDR=p.adjust(res$table$PValue, "BH"))
hist(pvals$raw.pvalue,150)

res<-topTags(res,nrow(dge$counts))
alpha=0.1

#Quick check to see how many genes are flagged as DE and if there are any imbalances. 
DEGs<-res$table[res$table$FDR<=alpha,]
table(sign(DEGs$logFC))

df1<- pseudo_counts_cen[rownames(pseudo_counts_cen) %in% rownames(DEGs),]

#Kmeans
#Because the number of clusters is unknown, we investigate a number varying from 2 to 50.
K_choice <- 2:50
km <- vector("list", length(K_choice))
names(km) <- paste("K=", K_choice, sep="")
for(K in K_choice) {
  cat("K =", K, "...\n")
  km[[paste("K=",K,sep="")]] <- kmeans(df1, centers=K,
                                       nstart=5)
}

#To choose the best number of clusters, we examine the within-sum of squares, Alternatively, BIC or other criterion could be better
withinss <- unlist(lapply(km, function(x) x$tot.withinss))
betweenss <- unlist(lapply(km, function(x) x$betweenss))
totss <- unlist(lapply(km, function(x) x$totss))

par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(K_choice, withinss, type="b", xlab="Number of Clusters",
     ylab="Within-group SS")


#Looking at plots, 50 clusters chosen, no clear minimum
km_final <- km[["K=50"]] # look for elbow, not easy with this set.....
km_labels <- km_final$cluster
km_centers <- km_final$centers
km_size <- km_final$size

#Visualize
par(mfrow=c(1,1))
matplot(t(km_centers), type="l", col="black", lty=1, lwd=2,
        ylab="Kmeans cluster centers", xlab="sample index")

#find 2B17 cluster
head(km_final$cluster[names(km_final$cluster) %in% "UGT2B17"])
UGT <- names(km_final$cluster[km_final$cluster == 49])

intersect(STAT3, UGT) 
intersect(UGT, IRF1) 
intersect(NFKB1, UGT) 

# I usually also print out the list of co-expressed genes for Chantal, and run some motif / pathway enrichment on the results





# For survival analyses

cd <- read.table(file.choose(), header = T, sep = "\t")
cd <- cd[cd$icgc_donor_id %in% colnames(tf_activities),]
cd <- cd[order(cd$icgc_donor_id),]
cd$IRF <- tf_activities[rownames(tf_activities)=="IRF1",]
cd$STAT3 <- tf_activities[rownames(tf_activities)=="STAT3",]
cd$NFKB1 <- tf_activities[rownames(tf_activities)=="NFKB1",]


library('survival')
tmp_df$donor_vital_status[tmp_df$donor_vital_status == "alive"] <- 0
tmp_df$donor_vital_status[tmp_df$donor_vital_status == "deceased"] <- 1
tmp_df$donor_vital_status <- as.numeric(tmp_df$donor_vital_status)
tmp_df$donor_survival_time <- tmp_df$donor_survival_time/365


model1<-survfit(formula=Surv(tmp_df$donor_survival_time,tmp_df$donor_vital_status)~tmp_df$STAT3, data=tmp_df)#replace STAT3 by whichever TF
survdiff(formula=Surv(tmp_df$donor_survival_time,tmp_df$donor_vital_status)~tmp_df$STAT3, data=tmp_df)
summary(model1)
#survival univariate
a<-coxph(formula=Surv(tmp_df$donor_survival_time,tmp_df$donor_vital_status)~tmp_df$STAT3, data=tmp_df)
summary(a)