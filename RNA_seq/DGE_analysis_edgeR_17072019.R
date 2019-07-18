#!/usr/bin/env/Rscript

"

Title:  Rmarkdown script for easier visualization of differential gene expression.
Author: Eric Allain



"
#Setting up your data
#This analysis will require two files, a raw count matrix in tab-seperated format, and a design matrix, also in .tsv format. 
#The count matrix should be organized with sample names as column headers and feature names as row labels. 
#The design matrix should be inverted, with sample labels as row names and phenotype data (Treatment, Sex, Age, Cell line, Tissue, etc.) as column headers.
#The first column in the design matrix should always be the sample groups which aggregate all technical replicates. Other factors will follow. 

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

#Set your working directory ()

setwd("~/Path/to/your/wd")

raw_counts <- read.table("../data/counts.txt", header = TRUE, row.names = 1, stringsAsFactors = F, sep = "\t")
raw_counts <- as.matrix(raw_counts)

#loose filtering for removal of low-abundance transcripts
raw_counts_filt <- raw_counts[rowSums(raw_counts) > 0, ]

#Load up phenoData in tab-delimited format, samples on rows, features on columns. 
Pheno1<- read.table("~/Path/to/PhenoData.tsv", header=TRUE, row.names = 1, sep = "\t")

#store list of all phenotype data names for exploratory analysis
covariables<-colnames(Pheno1)

Pheno2<- new("AnnotatedDataFrame", data = as.data.frame(Pheno1))

#Build new Expressionset for analysis
ESet<- ExpressionSet(assayData = raw_counts_filt, phenoData = Pheno2)


#create DGE object
dge <- DGEList(raw_counts_wn)


#Set some graphical parameters
##add the number of colours which corresponds to the number of experimental conditions. The end argument (each = 3) corresponds to the number of replicates.
##This could be set to colour graphs based on any variable. 

colors=rep(c("dodgerblue","firebrick1", "MediumVioletRed","SpringGreen","orange","lightskyblue","mediumpurple1","darkolivegreen3"),each=3)


#raw library sizes
barplot(dge$samples$lib.size,col=colors,main="Raw library size")

#density plots of data distributions

#transform to log-scale for easier visualization
pseudo_counts <- log2(dge$counts+1)

#plot the density distribution for each sample triplicate (by=3), fetch graph titles from the main factor column in the phenotype data file (Pheno1)
#mfrow may need to be changes depending on number of samples. Default can accomodate 9 samples in triplicate (28 columns).
#this allows us to quickly check for any possible outliers or large changes in data distributions which may bias results.

par(mfrow=c(3,3))
for (i in seq(1,ncol(pseudo_counts),by=3))
{
  plot(density(pseudo_counts[,i]),col=colors[i],main=Pheno1$group[i])
  lines(density(pseudo_counts[,i+1]),lwd=2,col=colors[i+1])
  lines(density(pseudo_counts[,i+2]),lwd=2,col=colors[i+2])
}

#We can also plot all experimental conditions on one single plot. 
i=1
par (mfrow=c(1,1))
plot(density(pseudo_counts[,i]),main="Count distribution - all samples",col=colors[i])
for (i in 2:24)
  lines(density(pseudo_counts[,i]),lwd=2,col=colors[i])

#Eucledian distance-based methods such as singular vector decomposition (SVD) can then be used to view reduced-dimention data.
#We can then colour-code datapoints to show trends
#Here we can see the variance explained by each principal component. 

#first, center the data
pseudo_counts_cen <-pseudo_counts - rowMeans(pseudo_counts)

#Now lets see how much variance is explained by each principal component after dimentionnality reduction.
svd1<-svd(pseudo_counts_cen)
svd_percent<-(svd1$d^2/sum(svd1$d^2))
plot(svd_percent,main="% variance explained",col=2)

#Then we can colour datapoints by factor levels to observe trends.

par(mfrow=c(3,3))
for(item in covariables){
  plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",xlab="1st PC",col=as.factor(Pheno1[[item]]), main=item)
}


#between sample normalization. In this case the Trimmed mean of M-values (TMM) method is used, which is a robust normalization method for RNAseq data.
dge <- calcNormFactors(dge, method = "TMM")
pseudo_TMM <- log2(scale(dge$counts,center=FALSE,scale=dge$samples$norm.factors)+1)

i=1
par (mfrow=c(1,1))
plot(density(pseudo_TMM[,i]),main="Count distribution - all samples",col=colors[i])
for (i in 2:24)
  lines(density(pseudo_TMM[,i]),lwd=2,col=colors[i])

#Now we can generate our previous plots with normalized data to check that nothing has changed.

pseudo_counts_cen <-pseudo_TMM - rowMeans(pseudo_TMM)
svd1<-svd(pseudo_counts_cen)
svd_percent<-(svd1$d^2/sum(svd1$d^2))
plot(svd_percent,main="% variance explained",col=2)
par(mfrow=c(3,3))
for(item in covariables){
  plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",xlab="1st PC",col=as.factor(Pheno1[[item]]), main=item)
}

#we can even view this interactively by altering the threshold.

#We can even view this interactively wih shiny.

pca<-as.data.frame(svd1$v)

butt<-list()
for (item in covariables) {
  but1<-list(method = "restyle", args = list("marker.color", list(Pheno1[[item]]), "text", list(factor(Pheno1[[item]]))),label = item)
  butt<-c(butt,but1)
}

p <- plot_ly(pca, 
             x = ~V1, 
             y = ~V2) %>%
  add_markers(marker = list(line = list(color = "black", width = 1))) %>%
  
  layout(
    title = "Principal component analysis",
    updatemenus = list(
      list(
        y = 0.75,
        buttons = butt
      )
    )
  )
p



#Differential gene expression analysis and statistics
#Setting the references is done to constrain our statistical model. 
#We are loking to estimate the parameters which best fit the model, through the maimum likelihood function. The maximum of this function is the best-fit parameter.
#Multiple maxima for the best-fit function can be found if we don't fix a reference condition, so applying this constraint allows for a reference condition and replicability.
#Untreated control is the reference value for gene expression in most models.

factor1 <- relevel(as.factor(design$factor1), ref = "WT")
factor2 <- relevel(as.factor(design$factor2), ref = "control")

#Make sure to go through each factor in this way to be sure the control condition is appropriately set.

#build the linear model. This is where careful thought should be applied to determine the variables to include. 
#models can be compared, but avoid overfitting. Contrary to popular belief, subjectivity is key to building the simplest model which best explains the data.

design_matrix <- model.matrix(~ design$replicat + factor1 + factor2 +
                                factor1:factor2)


#Estimate dispertions, a necessary step to assess gene-level scatter in our data.
dge <- estimateGLMCommonDisp(dge, design_matrix)
dge <- estimateGLMTrendedDisp(dge, design_matrix)
dge <- estimateGLMTagwiseDisp(dge, design_matrix)

#plotting a BCV graph shows us if dispertion is adequately controlled with low-abundance genes. 
#If dispertion is too high, consider using a more stringent cutoff for low-abundance genes. 

par(mfrow=c(1,1))
plotBCV(dge, main = paste0("BCV plot"))




#completely arbitraty filtering, anything is good as long as we are increasing the stringency of our threshold. 
ridx <- rowSums(cpm(dge) > 1) >= nrow(design)/2 
dge.f <- dge[ridx,]
cat("Number of genes after filtering:", nrow(dge.f$counts),"\n")

#estimate dispertion again with new filter
dge.f$samples$lib.size<-colSums(dge.f$counts)
dge.f <- calcNormFactors(dge.f, method = "TMM")
dge.f <- estimateGLMCommonDisp(dge.f, design_matrix)
dge.f <- estimateGLMTrendedDisp(dge.f, design_matrix)
dge.f<- estimateGLMTagwiseDisp(dge.f, design_matrix)

plotBCV(dge.f, main = paste0("BCV plot"))

#By looking at the design matrix colnames(design_matrix), 
#we deduce that we want to test the fourth coefficient.
design_matrix

#fit the model
fit.f<-glmFit(dge.f,design_matrix)

#1st coef is always intercept, set coefficients to whichever contrast is needed. 
res <- glmLRT(fit.f, coef = 4)

#Adjust for multiple comparisons
pvals<- data.frame(raw.pvalue = res$table$PValue,
                   FDR=p.adjust(res$table$PValue, "BH"))
hist(pvals$raw.pvalue,150)

#if p-val graph is normal for all contrasts, good, if not, it means that your N is perhaps too weak to include interaction terms in the model
#Prudence in interpretation can be warranted. 

res<-topTags(res,nrow(dge.f$counts))
alpha=0.05

#Quick check to see how many genes are flagges as DE and if there are any imbalances. 
DEGs<-res$table[res$table$FDR<=alpha,]
table(sign(DEGs$logFC))

#volcano
res$table$threshold<-ifelse(res$table$FDR < alpha, "FDR < 0.05", "N.S.")
res$table$threshold<-as.factor(res$table$threshold)
volcano_dat<-res$table
volcano_dat$name<-rownames(res$table)
g = ggplot(data=volcano_dat, aes(x=logFC, y=-log10(FDR), colour=threshold))+geom_point(alpha=0.4,size=1.75)+xlab("Log2 Fold-Change") + ylab("-Log10 (adjusted P value)")+theme_bw()+theme(axis.text=element_text(size=20),axis.title=element_text(size=20))+theme(legend.position="none")
ggplotly(g+geom_text(aes(label=name)))

p <- plot_ly(volcano_dat, x = ~logFC, y = ~-log10(FDR), type = "scatter", color = volcano_dat$threshold, text = volcano_dat$name)

sessionInfo()
