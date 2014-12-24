library(affy)
library(arrayQualityMetrics)
library(sva)

setwd("~/Desktop/science/Lab Research/Placenta/Affy_Control")

# Version 1. Normalize all and do combat
pd = read.AnnotatedDataFrame(filename="pdata_affy_control.txt")
affyData = ReadAffy(phenoData=pd, sampleNames=pd$SampleAccessionNumber)
#arrayQualityMetrics(expressionset = affyData, outdir = "QC_Report_Trimester_log", force = TRUE, do.logtransform = TRUE, intgroup = "Trimester")

eset = rma(affyData)
#arrayQualityMetrics(expressionset = eset, outdir = "QC_Report_RMA_Trimester", force = TRUE, do.logtransform = FALSE, intgroup = "Trimester")

# Remove batch-effect caused by different DataSets
batch = pd$DataSetAccesionNumber
mod = model.matrix(~as.factor(Trimester), data=pData(pd))
combat_edata = ComBat(dat=exprs(eset), batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
exprs(eset) <- combat_edata
#arrayQualityMetrics(expressionset = eset, outdir = "QC_Report_RMA_ComBat_Trimester", force = TRUE, do.logtransform = FALSE, intgroup = "Trimester")

# Version 2. Normalize within each DataSet, do combat

pd = read.AnnotatedDataFrame(filename="pdata_affy_control.txt")

pd12767 <- pd[pd$DataSetAccesionNumber=="GSE12767"]
affyData12767 = ReadAffy(filenames=row.names(pd12767), phenoData=pd12767, sampleNames=pd12767$SampleAccessionNumber)
eset12767 = rma(affyData12767)

pd9984 <- pd[pd$DataSetAccesionNumber=="GSE9984"]
affyData9984 = ReadAffy(filenames=row.names(pd9984), phenoData=pd9984, sampleNames=pd9984$SampleAccessionNumber)
eset9984 = rma(affyData9984)

pd18809 <- pd[pd$DataSetAccesionNumber=="GSE18809"]
affyData18809 = ReadAffy(filenames=row.names(pd18809), phenoData=pd18809, sampleNames=pd18809$SampleAccessionNumber)
eset18809 = rma(affyData18809)

allGenes <- combine(eset12767, eset9984, eset18809)

batch = pd$DataSetAccesionNumber
mod = model.matrix(~as.factor(Trimester), data=pData(pd))
combat_edata = ComBat(dat=exprs(allGenes), batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
exprs(allGenes) <- combat_edata

#arrayQualityMetrics(expressionset = allGenes, outdir = "QC_Report_RMA_ComBat_Trimester_V2", force = TRUE, do.logtransform = FALSE, intgroup = "Trimester")

# DiffExp genes

library(Biobase)
library(limma)

# First and Third
ind13 <- colnames(exprs(allGenes)) %in% pd[pd$Trimester!="Second"]$SampleAccessionNumber
exprs13 <- exprs(allGenes)[,ind13]
types = factor(pd[ind13]$Trimester)
design = model.matrix(~ 0+types)
colnames(design) = levels(types)
contrast.matrix = makeContrasts(First-Third, levels=design)
fit = lmFit(exprs13,design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
topTable(fit2, number=100, sort.by="p")

# Plot p-value distrib
tT <- topTable(fit2, adjust.method="fdr", number=nrow(fit2))
x <- -log10(tT$adj.P.Val)
plot(x, type="l")
sigline <- c(.05, .01, .005, .001,.0005, .0001, .00005)
sigline <- -log10(sigline)
sigcolors <- c("red", "blue", "green", "yellow","pink","purple", "black")
sapply(1:length(sigline), function(x){abline(h=sigline[x], col=sigcolors[x])})

# Annotation

library(hgu133plus2.db)
probes = fit2$genes$ID
fit2$genes$Symbol = unlist(mget(probes,hgu133plus2SYMBOL,ifnotfound=NA))
fit2$genes$EntrezID = unlist(mget(probes,hgu133plus2ENTREZID,ifnotfound=NA))
fit2$genes$Description = unlist(mget(probes,hgu133plus2GENENAME,ifnotfound=NA))

fullList = topTable(fit2, number=100, sort.by="logFC", resort.by="p")
write.table(fullList, file="top100_13trimester_results.txt", row.names=FALSE, sep="\t", quote=FALSE)

# First and Second
ind12 <- colnames(exprs(allGenes)) %in% pd[pd$Trimester!="Third"]$SampleAccessionNumber
exprs12 <- exprs(allGenes)[,ind12]
types = factor(pd[ind12]$Trimester)
design = model.matrix(~ 0+types)
colnames(design) = levels(types)
contrast.matrix = makeContrasts(First-Second, levels=design)
fit = lmFit(exprs12,design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
topTable(fit2, number=100, sort.by="p")

# Plot p-value distrib
tT <- topTable(fit2, adjust.method="fdr", number=nrow(fit2))
x <- -log10(tT$adj.P.Val)
plot(x, type="l")
sigline <- c(.05, .01, .005, .001,.0005, .0001, .00005)
sigline <- -log10(sigline)
sigcolors <- c("red", "blue", "green", "yellow","pink","purple", "black")
sapply(1:length(sigline), function(x){abline(h=sigline[x], col=sigcolors[x])})

# Annotation

library(hgu133plus2.db)
probes = fit2$genes$ID
fit2$genes$Symbol = unlist(mget(probes,hgu133plus2SYMBOL,ifnotfound=NA))
fit2$genes$EntrezID = unlist(mget(probes,hgu133plus2ENTREZID,ifnotfound=NA))
fit2$genes$Description = unlist(mget(probes,hgu133plus2GENENAME,ifnotfound=NA))

fullList = topTable(fit2, number=100, sort.by="logFC", resort.by="p")
write.table(fullList, file="top100_12trimester_results.txt", row.names=FALSE, sep="\t", quote=FALSE)

