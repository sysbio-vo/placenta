library(affy)
library(arrayQualityMetrics)
library(sva)
library(Biobase)
library(limma)
# Load for Illumina analysis
library(lumi)
# Load hgu133plus2 custom cdf
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
library(hgu133plus2.db)
# Load hugene10st custom cdf
library(hugene10sthsentrezgprobe)
library(hugene10sthsentrezgcdf)
library(hugene10sthsentrezg.db)
# Load nugohs1a520180hs custom cdf
library(nugohs1a520180hsentrezgprobe)
library(nugohs1a520180hsentrezgcdf) 
library(nugohs1a520180hsentrezg.db) 

setwd("/home/vlad/gene_regulatory network project/NORMA/")


############## Normalization/BRAINARRAY 

#arrayQualityMetrics(expressionset = eset, outdir = "QC_Report_RMA_ComBat_Trimester", force = TRUE, do.logtransform = FALSE, intgroup = "Trimester")
# Version 2. Normalize within each DataSet, do combat

allpData<-list.files()[which(grepl("^pdata_",list.files())==TRUE)]
for(pdata_filename in allpData){
pd = read.AnnotatedDataFrame(filename=pdata_filename)

pdata<-pData(pd)
if("GSE12767" %in% pd$DataSetAccesionNumber){ 
    pd12767 <- pd[pd$DataSetAccesionNumber=="GSE12767"]
    affyData12767 = ReadAffy(filenames=row.names(pd12767), phenoData=pd12767, sampleNames=pd12767$SampleAccessionNumber)
    affyData12767@cdfName<-"hgu133plus2hsentrezgcdf"
    eset12767 = rma(affyData12767)}
  
if("GSE9984" %in% pd$DataSetAccesionNumber){ 
    pd9984 <- pd[pd$DataSetAccesionNumber=="GSE9984"]
    affyData9984 = ReadAffy(filenames=row.names(pd9984), phenoData=pd9984, sampleNames=pd9984$SampleAccessionNumber)
    affyData9984@cdfName<-"hgu133plus2hsentrezgcdf"
    eset9984 = rma(affyData9984)}
  
if("GSE18809" %in% pd$DataSetAccesionNumber){
    pd18809 <- pd[pd$DataSetAccesionNumber=="GSE18809"]
    affyData18809 = ReadAffy(filenames=row.names(pd18809), phenoData=pd18809, sampleNames=pd18809$SampleAccessionNumber)
    affyData18809@cdfName<-"hgu133plus2hsentrezgcdf"
    eset18809 = rma(affyData18809)}
  
if("GSE53291" %in% pd$DataSetAccesionNumber){ 
    pdGSE53291 <- pd[pd$DataSetAccesionNumber=="GSE53291"]
    affyDataGSE53291 = ReadAffy(filenames=row.names(pdGSE53291), phenoData=pdGSE53291, sampleNames=pdGSE53291$SampleAccessionNumber)
    affyDataGSE53291@cdfName<-"nugohs1a520180hsentrezgcdf"
    esetGSE53291 = rma(affyDataGSE53291)}
  
  
if("GSE24129" %in% pd$DataSetAccesionNumber){ 
    pdGSE24129 <- pd[pd$DataSetAccesionNumber=="GSE24129"]
    affyDataGSE24129 = ReadAffy(filenames=row.names(pdGSE24129), phenoData=pdGSE24129, sampleNames=pdGSE24129$SampleAccessionNumber)
    affyDataGSE24129@cdfName<-"hugene10sthsentrezgcdf"
    esetGSE24129 = rma(affyDataGSE24129)}
  
################## ADD ILLUMINA DATA
  
# For GSE30186 ILLUMINA_HumanHT-12 V4.0
annotation<-read.table("/home/vlad/gene_regulatory network project/BLAST/83.3minOverallIdenILLUMINA_HumanHT-12 V4.0 expression beadchip_IDs.txt",sep=" ")
# Select one probe-one EntrezID:
EntrezID<-unique(na.exclude(annotation$V2))
probeIDs<-vector()
for(ID in EntrezID){
  subprobe<-annotation[which(annotation$V2==ID),]
  subprobe<-subprobe[which(subprobe$V4==max(subprobe$V4)),] #Filter with Specifity score
  subprobe<-subprobe[which(subprobe$V3==max(subprobe$V3)),] #Filter with Coverage score
  if(length(as.character(subprobe$V1))>1){
    num<-sample(1:length(subprobe$V1),size = 1)
    probeID<-as.character(subprobe$V1[num])}
    if(length(subprobe$V1)==1){
      probeID<-as.character(subprobe$V1)}
      probeIDs<-append(probeIDs,probeID)}
      to_sub<-7:12
      file_name<-"/home/vlad/gene_regulatory network project/data_norma/GSE30186_non_normalized.txt"
      annotation0<-read.table(file_name,sep="\t",header = T)
      file_name<-"/home/vlad/gene_regulatory network project/data_norma/GSE30186_non_normalized_filtered.txt"
      write.table(annotation0[which(annotation0$ID_REF %in% probeIDs),],file_name,sep="\t",quote = F,col.names=T,row.names=F)
      x.lumi <- lumiR.batch(file_name,transform = "log")#, sampleInfoFile='sampleInfo.txt')
      x.lumi<-x.lumi[,to_sub]
      # Quantile normalization
      lumi.N<-lumiN(x.lumi,method = "quantile")
      lumi.N.Q_GSE30186<-data.frame(exprs(lumi.N))
      colnames(lumi.N.Q_GSE30186)<-c("GSE30186_1","GSE30186_2","GSE30186_3","GSE30186_4","GSE30186_5","GSE30186_6")
      rownames(lumi.N.Q_GSE30186)<-EntrezID[which(probeIDs %in% rownames(lumi.N.Q_GSE30186))]
    
if("GSE25906" %in% pd$DataSetAccesionNumber){ 
    # For GSE25906 ILLUMINA_Human-6 V2.0
    file_name<-"/home/vlad/gene_regulatory network project/data_norma/GSE25906_non-normalized_filtered.txt"
    to_sub<-1:37
    annotation<-read.table("/home/vlad/gene_regulatory network project/BLAST/83.3minOverallIdenILLUMINA_HumanHT-6 V2.0 expression beadchip_IDs.txt",sep=" ")
    # Select one probe-one EntrezID:
    EntrezID<-unique(na.exclude(annotation$V2))
    probeIDs<-vector()
    for(ID in EntrezID){
      subprobe<-annotation[which(annotation$V2==ID),]
      subprobe<-subprobe[which(subprobe$V4==max(subprobe$V4)),] #Filter with Specifity score
      subprobe<-subprobe[which(subprobe$V3==max(subprobe$V3)),] #Filter with Coverage score
      if(length(as.character(subprobe$V1))>1){
        num<-sample(1:length(subprobe$V1),size = 1)
        probeID<-as.character(subprobe$V1[num])}
      if(length(subprobe$V1)==1){
        probeID<-as.character(subprobe$V1)}
      probeIDs<-append(probeIDs,probeID)}
    file_name<-"/home/vlad/gene_regulatory network project/data_norma/GSE25906_non-normalized.txt"
    annotation0<-read.table(file_name,sep="\t",header = T)
    file_name<-"/home/vlad/gene_regulatory network project/data_norma/GSE25906_non-normalized_filtered.txt"
    write.table(annotation0[which(annotation0$ID_REF %in% probeIDs),],file_name,sep="\t",quote = F,col.names=T,row.names=F)
    
    x.lumi <- lumiR.batch(file_name,transform = "log")#, sampleInfoFile='sampleInfo.txt')
    x.lumi<-x.lumi[,to_sub]
    # Quantile normalization
    lumi.N<-lumiN(x.lumi,method = "quantile")
    lumi.N.Q_GSE25906<-data.frame(exprs(lumi.N))
    colnames(lumi.N.Q_GSE25906)<-paste("GSE25906_",to_sub,sep="")
    rownames(lumi.N.Q_GSE25906)<-EntrezID[which(probeIDs %in% rownames(lumi.N.Q_GSE25906))]}

################################
hgu133plus2_probeIDs<-unique(hgu133plus2hsentrezgprobe$Probe.Set.Name)
hgu133plus2_geneEntrezIDs<-as.character(unlist(mget(hgu133plus2_probeIDs,hgu133plus2hsentrezgENTREZID,ifnotfound = NA)))
  
hugene10st_probeIDs<-unique(hugene10sthsentrezgprobe$Probe.Set.Name)
hugene10st_geneEntrezIDs<-as.character(unlist(mget(hugene10st_probeIDs,hugene10sthsentrezgENTREZID)))
  
nugohs1a520180_probeIDs<-unique(nugohs1a520180hsentrezgprobe$Probe.Set.Name)
nugohs1a520180_geneEntrezIDs<-as.character(unlist(mget(nugohs1a520180_probeIDs,nugohs1a520180hsentrezgENTREZID,ifnotfound = NA)))

#gene_ENTREZ_all<-intersect(intersect(intersect(intersect(na.exclude(hgu133plus2_geneEntrezIDs),na.exclude(hugene10st_geneEntrezIDs)),na.exclude(nugohs1a520180_geneEntrezIDs)),na.exclude(unique(rownames(lumi.N.Q_GSE25906)))),na.exclude(unique(rownames(lumi.N.Q_GSE30186))))
gene_ENTREZ_all<-intersect(intersect(intersect(na.exclude(hgu133plus2_geneEntrezIDs),na.exclude(hugene10st_geneEntrezIDs)),na.exclude(nugohs1a520180_geneEntrezIDs)),na.exclude(unique(rownames(lumi.N.Q_GSE30186))))

to_sub_hgu133<-which(hgu133plus2_geneEntrezIDs %in% gene_ENTREZ_all)
to_sub_hu1st<-which(hugene10st_geneEntrezIDs %in% gene_ENTREZ_all)
to_sub_nugohs1<-which(nugohs1a520180_geneEntrezIDs %in% gene_ENTREZ_all)
#to_sub_illumina6v2.0<-which(rownames(lumi.N.Q_GSE25906) %in% gene_ENTREZ_all)
  
to_sub_illuminaHT12v2.0<-which(rownames(lumi.N.Q_GSE30186) %in% gene_ENTREZ_all)
if(pdata_filename=="pdata_12430.txt"){
  allGenes <- cbind(exprs(eset12767)[to_sub_hgu133,], exprs(eset9984)[to_sub_hgu133,], exprs(eset18809)[to_sub_hgu133,],exprs(esetGSE24129)[to_sub_hu1st,],exprs(esetGSE53291)[to_sub_nugohs1,],lumi.N.Q_GSE30186[to_sub_illuminaHT12v2.0,])}#,lumi.N.Q_GSE25906[to_sub_illumina6v2.0,])  
if(pdata_filename=="pdata_444.txt"){
  allGenes <- exprs(eset9984)[to_sub_hgu133,]}
if(pdata_filename=="pdata_449.txt"){
  allGenes <- cbind(exprs(eset9984)[to_sub_hgu133,],exprs(eset18809)[to_sub_hgu133,])}
if(pdata_filename=="pdata_1244.txt"){
  allGenes <- cbind(exprs(eset9984)[to_sub_hgu133,],exprs(eset12767)[to_sub_hgu133,])}
if(pdata_filename=="pdata_1249.txt"){
  allGenes <- cbind(exprs(eset9984)[to_sub_hgu133,],exprs(eset12767)[to_sub_hgu133,],exprs(eset18809)[to_sub_hgu133,])}
if(pdata_filename=="pdata_12415.txt"){
  allGenes <- cbind(exprs(eset9984)[to_sub_hgu133,],exprs(eset12767)[to_sub_hgu133,],exprs(eset18809)[to_sub_hgu133,],lumi.N.Q_GSE30186[to_sub_illuminaHT12v2.0,])}
if(pdata_filename=="pdata_12416.txt"){
  allGenes <- cbind(exprs(eset9984)[to_sub_hgu133,],exprs(eset12767)[to_sub_hgu133,],exprs(eset18809)[to_sub_hgu133,],exprs(esetGSE53291)[to_sub_nugohs1,])}
if(pdata_filename=="pdata_12417.txt"){
  allGenes <- cbind(exprs(eset9984)[to_sub_hgu133,],exprs(eset12767)[to_sub_hgu133,],exprs(eset18809)[to_sub_hgu133,],exprs(esetGSE24129)[to_sub_hu1st,])}  
rownames(allGenes)<-gene_ENTREZ_all
pd = read.AnnotatedDataFrame(filename=pdata_filename)
rownames(pd)<-pd$SampleAccessionNumber
allGenes<-allGenes[,rownames(pd)]

##################################################################
# batch
#rownames(pd)<-colnames(combat_edata)

batch = pd$DataSetAccesionNumber
mod = model.matrix(~as.factor(Trimester), data=pData(pd))
combat_edata = ComBat(dat=allGenes, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
  
allGenes1<-ExpressionSet(assayData=as.matrix(combat_edata),phenoData = pd)
#exprs(allGenes) <- combat_edata
arrayQualityMetrics(expressionset = allGenes1, outdir = paste("QC_Report_RMA_ComBat_Trimester_allPlatformsAffyandIllumina_",pdata_filename,collapse="",sep=""), force = TRUE, do.logtransform = FALSE, intgroup = "Trimester")
  
# DiffExp genes Affymetrix
# First and Third

filename_out1<-paste("13trimester_",pdata_filename,sep="",collapse="")
ind13 <- colnames(exprs(allGenes1)) %in% pd[pd$Trimester!="Second"]$SampleAccessionNumber
exprs13 <- exprs(allGenes1)[,ind13]
types = factor(pd[ind13]$Trimester)
design = model.matrix(~ 0+types)
colnames(design) = levels(types)
contrast.matrix = makeContrasts(First-Third, levels=design)
fit = lmFit(exprs13,design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)

probes = rownames(exprs(allGenes1))
fit2$genes$ProbeID<-probes
  
fullList = topTable(fit2, number=dim(fit2)[1], sort.by="logFC", resort.by="p")
fullList<-fullList[which(Mod(fullList$logFC)>=2 & fullList$P.Value<=0.05),]
write.table(fullList, file=filename_out1, row.names=FALSE, sep="\t", quote=FALSE)
# First and Second

filename_out2<-paste("12trimester_",pdata_filename,sep="",collapse="")
ind12 <- colnames(exprs(allGenes1)) %in% pd[pd$Trimester!="Third"]$SampleAccessionNumber
exprs12 <- exprs(allGenes1)[,ind12]
types = factor(pd[ind12]$Trimester)
design = model.matrix(~ 0+types)
colnames(design) = levels(types)
contrast.matrix = makeContrasts(First-Second, levels=design)
fit = lmFit(exprs12,design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
top2<-topTable(fit2, number=100, sort.by="p")

fit2$genes$ProbeID<-probes
fullList = topTable(fit2, number=dim(fit2)[1], sort.by="logFC", resort.by="p")
fullList<-fullList[which(Mod(fullList$logFC)>=2 & fullList$P.Value<=0.05),]

write.table(fullList, file=filename_out2, row.names=FALSE, sep="\t", quote=FALSE)
#write.table(exprs(allGenes1),"expressionset.txt",sep="\t",row.names=T)
# Second and third

filename_out3<-paste("23trimester_",pdata_filename,sep="",collapse="")
ind23 <- colnames(exprs(allGenes1)) %in% pd[pd$Trimester!="First"]$SampleAccessionNumber
exprs23 <- exprs(allGenes1)[,ind23]
types = factor(pd[ind23]$Trimester)
design = model.matrix(~ 0+types)
colnames(design) = levels(types)
contrast.matrix = makeContrasts(Second-Third, levels=design)
fit = lmFit(exprs23,design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
top2<-topTable(fit2, number=100, sort.by="p")
  
probes<-rownames(allGenes)
fit2$genes$ProbeID<-probes
  
fullList = topTable(fit2, number=dim(fit2)[1], sort.by="logFC", resort.by="p")
fullList<-fullList[which(Mod(fullList$logFC)>=2 & fullList$P.Value<=0.05),]

write.table(fullList, file=filename_out3, row.names=FALSE, sep="\t", quote=FALSE)
#write.table(exprs(allGenes1),"expressionset.txt",sep="\t",row.names=T)
}  

############ ANALYSE RESULTS ############

allpData<-c("pdata_444.txt","pdata_449.txt","pdata_1244.txt","pdata_1249.txt","pdata_12415.txt","pdata_12416.txt","pdata_12417.txt","pdata_12430.txt")
par(mfrow=c(3,8))
files_in<-list.files()
files_in<-files_in[which(grepl("trimester_pdata",files_in)==TRUE)]
par(mfrow=c(3,8))
for(p_data in allpData){
  files_in2<-files_in[which(grepl(p_data,files_in)==TRUE)]
  for(file_in in files_in2){
    data<-read.table(file_in,sep="\t",header = T)
    hist(data$adj.P.Val,main=paste(as.character(strsplit(file_in,paste("_",p_data,sep="",collapse="")))," N=",length(unique(data$ProbeID)),sep="",collapse=""),xlab=strsplit(p_data,".txt"))}
}






