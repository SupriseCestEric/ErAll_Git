#dependencies
library(limma)
library(affy)

MyData <- ReadAffy()
eset <- rma(MyData)
expr_data<-exprs(eset)

## read phenotype data
pheno_data<- read.table(file.choose(), header = TRUE, sep= "\t")
pheno_data<- as.matrix(pheno_data)

## Build new Expressionset for analysis
APheno<- new("AnnotatedDataFrame", data = as.data.frame(pheno_data))
StatsSet<- ExpressionSet(assayData = expr_data, phenoData = APheno)
Design = model.matrix(~StatsSet$X) #X = primary variable


## fit the model
fit <- lmFit(StatsSet, Design)
## Improve Variance estimation with eBayes
eBfit <- eBayes(fit)
toptable(eBfit)


library("hgu133a.db", lib.loc="~/R/win-library/3.2")
res <- topTable(eBfit, n=nrow(eBfit), adjust = "BH")
Probes <- as.character(rownames(res))
NewNames<- select(hgu133a.db, Probes, c("SYMBOL", "GENENAME"))
NewNames<- NewNames[!duplicated(NewNames$PROBEID),]
res <- cbind(res,NewNames)
res_agg <-aggregate(VF[,1:(ncol(res)-1)],list(VF[,-1]),FUN=mean) # aggregated results for gene-level expression
write.table(res_agg,file="file_output_name.txt")