## Big Thanks to Charles Joly-Beauparlant

library(rtracklayer)
library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
org <- org.Hs.eg.db

# Import gene_set
####
gene_set <- gr_universe[sample(1:length(gr_universe), 150)] # Ã€ changer
gene_set_1 <- mapIds(org, keys = gene_set, keytype = "SYMBOL", column = "ENSEMBL")
gene_set_2<-gr_universe[gr_universe$ensembl %in% gene_set_1,]
####

# Importer les peaks
####
filenames <- list.files() %>%
  set_names(c("sample1", "sample2", "sample3")) # Set sample names for each file
####

#if narrowpeaks
extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

#if bedgraph
extraCols1<- c(scores = "numeric")


import_peak <- function(x) {
    if (str_detect(x, "\\.bedgraph$")) {
        res <- import(x, format = "BED", extraCols = extraCols1) %>%
            keepStandardChromosomes(pruning.mode = "coarse")
    } else if (str_detect(x, "\\.narrowPeak$")) {
        res <- import(x, format = "BED", extraCols = extraCols) %>%
            keepStandardChromosomes(pruning.mode = "coarse")
    } else {
        stop("Invalid format")
    }
    seqlevelsStyle(res) <- "UCSC"
    annotatePeak(res, TxDb = txdb, annoDb = "org.Hs.eg.db") %>%
        as.GRanges
}

peaks <- map(filenames, import_peak)

#grab peaks within a 10kb window centered on TSS
peaks_tss <- map(peaks, ~ .x[abs(.x$distanceToTSS) <= 10000])

# Grab all known genes from txdb
gr_universe <- genes(txdb)
gr_universe$ensembl <- mapIds(org, keys = gr_universe$gene_id, keytype = "ENTREZID", column = "ENSEMBL")

# Analysis by bootstrap
bootstrap_analysis <- function(x, nsample = 10000) { # x: narrowpeak to analyze
    size <- length(gene_set_2)
    # 0. Calculate gene_set_value, the observed enrichment
    gene_set_value <- length(findOverlaps(gene_set_2, x)) / size

    # 1. Precalculer dans l'univers, les gÃ¨nes qui se trouvent dans les narrowPeak
    precalculated <- gr_universe$ensembl %in% x$ENSEMBL

    # 2. Vecteur d'index
    index_vector <- sample(seq_along(precalculated), size*nsample, replace = TRUE)

    # 3. Boolean matrix with TRUE / FALSE. True if peaks overlap a gene in the dummy set
    m_boolean <- matrix(precalculated[index_vector], nrow = size)

    # 4. Calculate enrichment for each column in matrix
    bootstrap_result <- colSums(m_boolean)/size

    # 5. See where the observed value falls within the re-samlpled distribution
    #enlever cette section pour générer les graphiques
    sum(bootstrap_result < gene_set_value) / length(bootstrap_result)
}
set.seed(1234)

#set iterations
res <- map(peaks_tss, bootstrap_analysis, nsample = 1000000)



#Plor outcome
library("ggplot2")
g<-qplot(res[[3]], geom = "histogram")
g+theme_bw()+xlab("Fraction of genes targeted by NFkB in sample") + ylab("Count")
