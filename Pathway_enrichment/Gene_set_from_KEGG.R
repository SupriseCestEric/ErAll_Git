
#getting gene sets from kegg pathway
library("gage")
library(rtracklayer)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
org <- org.Hs.eg.db

kg<-kegg.gsets(species = "hsa", id.type = "kegg", check.new=FALSE)

AMPK<-kg$kg.sets$`hsa04152 AMPK signaling pathway`
STER<-kg$kg.sets$`hsa00100 Steroid biosynthesis`
XENO<-kg$kg.sets$`hsa00980 Metabolism of xenobiotics by cytochrome P450`
CCAR<-kg$kg.sets$`hsa05204 Chemical carcinogenesis`


AMPK1 <- mapIds(org, keys = AMPK, keytype = "ENTREZID", column = "SYMBOL")
STER1 <- mapIds(org, keys = STER, keytype = "ENTREZID", column = "SYMBOL")
XENO1 <- mapIds(org, keys = XENO, keytype = "ENTREZID", column = "SYMBOL")
CCAR1 <- mapIds(org, keys = CCAR, keytype = "ENTREZID", column = "SYMBOL")
