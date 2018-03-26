#!/usr/bin/env Rscript
#PUT THIS SCRIPT IN THE SAME DIR AS UGT ANNOTATION GFF FILE 
#INVOKE AT THE COMMAND LINE WITH EACH UGT TO BE PLOTTED AS ONE ARGUMENT
#GRAPH WILL DISPLAY FROM TOP TO BOTTOM< AS ENTERED AT COMMAND LINE
#WILL SAVE THE PNG UNDER THE NAME OF FIRST ARG, 
args = commandArgs(trailingOnly = TRUE)
# cat("Please enter your path/to/R/libraries (this may be retrieved with .libPaths in an R session):   ")
# input=readLines("stdin", n=1)
# print(input)
# .libPaths(input)
if(!require(Biobase)){source("https://bioconductor.org/biocLite.R")
  biocLite("Biobase")
  }
if(!require(GenomeGraphs)){source("https://bioconductor.org/biocLite.R")
  biocLite("GenomeGraphs")
  }
if(!require(biomaRt)){source("https://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
  }
if(!require(grid)){install.packages("grid")}

library("Biobase")
library("GenomeGraphs")
library("biomaRt")
library("grid")

data1<- read.table("all_UGT_v3.gff", fill=TRUE)
data1<- na.omit(data1)
colnames(data1)<-c("chrom","tool","type","start","end","dot", "strand","frame","info", "gene", "name", "split", "transcript")
#add gene at top? maybe with biomart annotation? The ideogram needs to come after loop to check chromosome...then append to moelList


SpliceGraphs<-function(){
  
  #SpliceGraphs is intended to plot variants from a custom .gff or .gtf annotation
  #swap out the filename above to your own custom annotation in the same directory and you can plot exon-intron plots with ideogram and legend below
  #
  #The generated filename will be the string of the 1st argument
  #arguments should correspond to transcript IDs in the .gff file 
  
  
  modelList = list()
  #loop over any rows which have the term to search from args and append section start-stop coords to the empty list
  for (item in args) local({
    ar=new("GeneModel", exonStart = data1[grep(toString(item), data1$transcript),4], exonEnd = data1[grep(toString(item), data1$transcript),5])
    modelList <<- c(modelList, ar)
  })
  print(modelList)
  #check if same chromosome and if so make ideogram, else print error
  for (item in args){chroms<<-data1[grep(toString(item),data1$transcript),1]}
  print(chroms)
  if (length(droplevels(chroms,chroms[1])) != length(chroms)) {
    stop("You must input at least one UGT isoform, with all queries on the same chromosome")
  } else
    if (length(chroms)==length(droplevels(chroms,chroms[1]))){
      ideo<<-new("Ideogram", chromosome = strsplit(substr(gsub("\\s+", "", data1[,1]), 4, 4), '')[[1]])
    }

  # this gets the range for the x axis, but right now only grabs the min and max from 1st arg
  mincoords<-min(modelList[[1]]@exonStart)
  maxcoords<-max(modelList[[1]]@exonEnd)
  #for (item in args){mincoords<<-as.vector(rbind(data1[grep(toString(item),data1$transcript),4],mincoords))}
  #for (item in args){maxcoords<<-as.vector(rbind(data1[grep(toString(item),data1$transcript),4],maxcoords))}
  print(mincoords)
  print(maxcoords)
  
  #cat1<- as.vector(rbind(mincoords, maxcoords))
  chrom.axis <- makeGenomeAxis(add53 = TRUE, add35=TRUE)

  modelList = c(ideo, modelList, chrom.axis)
  #geneModel=new("GeneModel", exonStart = data[,4], exonEnd = data[,5])
  #Maybe add a title at the top? perhaps the command args should be gene, transcript1, transcript2...transcript n.
  #define exon box width? may look better
  #Need to grab chromosome from one of the input transcripts and check to see if same locus, if not return error. 
  png(filename=toString(args[1]),width=1000,height=800)
  gdPlot(modelList, minBase = mincoords-5000, maxBase=maxcoords+5000)
  dev.off()
}
#chromosome = as.numeric(strsplit(substr(gsub("\\s+", "", data1[,1]), 4, 4), '')[[1]])


SpliceGraphs()