
###needs some form of link to local lib... unable to load packages. 

#dependencies

library("ggplot2")
library("Biobase")
library("limma")
library("gplots")
library("dplyr")
library("dendextend")

makeHierarch<-function(file1, prefix){
  #this makes a hierarchecal clustering dendogram of samples
  name1<-paste(prefix,"clust",Sys.Date(),".pdf")
  pdf(name1, onefile = TRUE)
  dist1=dist(t(as.matrix(file1)))
  hm<-heatmap(as.matrix(dist1),main="Heatmap of eucledian distance between samples")
  hc<-hclust(dist(t(as.matrix(file1))))
  hc1<-as.dendrogram(hc)
  print(hm)
  plot(hc,hang=-1)
  dev.off()
}

makeDistrPlot<-function(file1, prefix){
  #this makes a histogram showing the data distributions
  name1<-paste(prefix,"hist",Sys.Date(),".pdf")
  pdf(name1, onefile = TRUE)
  bPlot<-boxplot(as.matrix(file1),col=2,main="Raw Data")
  lbPlot<-boxplot(as.matrix(log2(file1)+1),col=2,main="Log2-transformed data")
  print(bPlot)
  print(lbPlot)
  dev.off()
}

makeSVDplot<-function(file1,prefix,pheno){
  #this makes an singular vector decomposition of data and plots
  #the first 2 principal components in 3D with colour-coded covariates
  name1<-paste(prefix,"SVD",Sys.Date(),".pdf")
  pdf(name1, onefile = TRUE)
  file1_center = file1 - rowMeans(file1)
  svd1<-svd(file1_center)
  svd_percent<-(svd1$d^2/sum(svd1$d^2))
  plot(svd_percent,main="% variance explained",col=2)
  for(item in colnames(pheno)){
    plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC",xlab="1st PC",col=as.factor(pheno[[item]]), main=item)
  }
  dev.off()
}


main<-function(){
  #
  #Usage: Rscript ExpressionSetExplore.R /path/to/wdir/ assayData.tsv phenoData.tsv graphfile_prefix_string 
  #
  
  
  
  args<-commandArgs(trailingOnly = TRUE)
  #Load up a tab-delim file with standard eSet format - features on rows and samples on columns
  folder=args[1]
  df=args[2]
  pd=args[3]
  pref=args[4]
  setwd(folder)
  Assay1<-read.table(df, header=TRUE, stringsAsFactors = F)
  filt_Assay1<-filter(Assay1,rowMeans(Assay1)>1)
  
  #Load up phenoData in tab-delim format, samples on rows, features on columns. 
  Pheno1<- read.table(pd, header=TRUE)
  
  #store list of all phenotype data names for exploratory analysis
  covariables<-colnames(Pheno1)
  
  Pheno2<- new("AnnotatedDataFrame", data = as.data.frame(Pheno1))
  ## Build new Expressionset for analysis
  StatsSet<- ExpressionSet(assayData = as.matrix(Assay1), phenoData = Pheno2)
  
  makeHierarch(filt_Assay1,pref)
  makeDistrPlot(filt_Assay1,pref)
  makeSVDplot(filt_Assay1,pref,Pheno1)
  
}

main()













