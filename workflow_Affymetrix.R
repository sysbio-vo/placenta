############ AFFYMETRIX DATA ANALYSIS
### 
##################################################################################################################################################
########################################## AFFYMETRIX PLATFORM READING FROM CEL ##################################################################
##################################################################################################################################################
### load ALL packages 

#library(hgu133plus2.db)
#library(annotate)
library(affy)
#library(pvclust)
library(vsn)
#library(gplots)
library(Biobase)
#library(GEOquery)
library(limma)
library(arrayQualityMetrics)
#library(hash)
# Download raw data from GEO and unpack

Affymetrix_annotations_GPL570<-c("GSE18809","GSE12767","GSE9984")

### This for-loop executes script to automaticly read (and preprocessing) all GPL570 data 

for(affy_a in Affymetrix_annotations_GPL570){
  setwd(paste("/home/vlad/gene_regulatory network project/data_norma/",affy_a, "_RAW/",sep="",collapse = ""))#9984_RAW/")
  # Create phenodata in each directory
  pd = read.AnnotatedDataFrame(filename="pheno_norma.txt")
  pData(pd)
  pdd<-pData(pd)
  affy_Data = ReadAffy(filenames = pdd$FileName)
  affy_Data
  arrayQualityMetrics(expressionset = affy_Data, outdir = "/home/vlad/gene_regulatory network project/QC_Report_with_SampleType", force = TRUE, do.logtransform = TRUE)
  #Normalization with RMA
  eset = rma(affy_Data)
  #QualityControl after normalization
  arrayQualityMetrics(expressionset = eset, outdir = "/home/vlad/gene_regulatory network project/QC_Report_with_SampleType_after_RMA", force = TRUE, do.logtransform = F)
  #eset = normalize.AffyBatch.quantiles(eset2, type=c("together"))
  eset
  #Remove bath samples if needed
  sampleNames(eset)<-unlist(lapply(X = sampleNames(eset),function(x) paste(x,affy_a,sep="",collapse="")))
  if(affy_a=="GSE24129"){
    all_data1<-eset}
  if(affy_a=="GSE53291"){
    all_data2<-eset}
  if(affy_a=="GSE18809"){
    all_data3<-eset}
  if(affy_a=="GSE12767"){
    all_data4<-eset}
  if(affy_a=="GSE9984"){
    all_data5<-eset}
  ### For each platform:
  #all_data1 - GSE24129 GPL6244
  #all_data2 - GSE53291 GPL18060
  #all_data3 - GSE18809 GPL570
  #all_data4 - GSE12767 GPL570
  #all_data5 - GSE9984  GPL570
}

# combine all GPL570 data
allGPL570<-as.data.frame(cbind(exprs(all_data3),exprs(all_data4),exprs(all_data5)[,c(1,3:8)]),row.names = rownames(all_data3))
#eset = normalize.AffyBatch.quantiles(eset2, type=c("together"))
eset0<-new('ExpressionSet',exprs=as.matrix(allGPL570))#,phenoData=pheno)
arrayQualityMetrics(expressionset = eset0, outdir = "/home/vlad/gene_regulatory network project/QC_Report_with_SampleType_after_RMA", force = TRUE, do.logtransform = F)

######################## CROSS-PLATFORM NORMALIZATION 

# Create phenodata for all 26(5) samples
pd = read.AnnotatedDataFrame(filename="/home/vlad/gene_regulatory network project/data_norma/pdata_GPL570norma")
pData(pd)
pheno<-pData(pd)
batch = as.factor(pheno$Batch)
mod = model.matrix(~1, data=pheno)
edata<-exprs(new('ExpressionSet',exprs=as.matrix(allGPL570),phenoData=pd))
combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
#qValuesComBat = p.adjust(pValuesComBat,method="BH")
combat_edata<-new('ExpressionSet',exprs=as.matrix(combat_edata))
arrayQualityMetrics(expressionset = combat_edata, outdir = "/home/vlad/gene_regulatory network project/QC_Report_with_SampleType_after_RMA", force = TRUE, do.logtransform = F)


######################## DIFFERENTIAL EXPRESSION 
### First, create design matrix

pd = read.AnnotatedDataFrame(filename="/home/vlad/gene_regulatory network project/data_norma/pdata_GPL570norma")
pData(pd)
Group<-factor(pData(pd)[,1], levels = unique(pData(pd)[,1]))
design <- model.matrix(~0 + Group)
#design<-model.matrix(~Group)
colnames(design) <- c("c3rd_trim","c1st_trim","c2nd_trim")
contrast.matrix <- makeContrasts( c1_3_trim = c1st_trim - c3rd_trim, c1_2_trim = c1st_trim - c2nd_trim, c3_2_trim <- c3rd_trim - c2nd_trim, levels=design)
fit <-lmFit(combat_edata,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2<-eBayes(fit2)
topTable(fit2, number=10, coef=4,adjust="fdr")
selected<-fit2$p.value[, 2] <0.05
esetSel<-combat_edata[selected, ]
esults <- decideTests(fit2)
vennDiagram(esults)
o <- order(fit2$F.p.value)
fit2$coefficients[o[1:30], ]
heatmap(exprs(esetSel))
