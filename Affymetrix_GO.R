library(affy)
library(arrayQualityMetrics)
library(sva)
setwd("/home/vlad/gene_regulatory network project/NORMA/")

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
##################################################################

### select one probe

select_probe<-function(genesdata, probesdata){
  if(dim(genesdata)[1]==dim(probesdata)[1]){
    genesID<-as.character(probesdata$GeneID)
    to_select<-vector()
    for(geneID in na.omit(unique(genesID))){
      sub_scores<-probesdata[which(genesID==geneID),]
      probes<-as.character(sub_scores$ProbeID)
      sub_genesdata<-genesdata[which(rownames(genesdata) %in% probes),]
      Overallscores<-as.numeric(as.character(sub_scores$Overall_score))
      Overallscores[which(is.na(Overallscores)==TRUE)]<-0
      if(all(Overallscores)!=0){
        max_score<-rownames(sub_scores)[which(Overallscores==max(Overallscores))]
        if(length(max_score)>1){
          max_score<-max_score[sample(1:length(max_score),size = 1)]
        }
        to_select<-append(to_select,as.numeric(max_score))
      }
    }}
  else if(dim(genesdata)[1]!=dim(probesdata)[1]){
    cat("incorrect dimensions")
  }
  return(to_select)
}

##############################################################################################################
dir_to_scoresfile<-"/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_Jetset.csv"
library(Biobase)
library(limma)
library(hgu133plus2.db)

probes_scores<-function(dir_to_scoresfile,num){
  if(num==1){
    add<-"1scores"
  }
  if(num==2){
    add<-"2scores"
  }
  if(num==3){
    add<-"jetsetscores"
  }
probes_allGenes<-rownames(allGenes)
scores_file<-read.table(dir_to_scoresfile,sep=",",header=T)
if(num==3){
  colnames(scores_file)[c(1,3,8)]<-c("ProbeID","GeneID","Overall_score")
}
scores_file<-scores_file[!duplicated(scores_file$ProbeID),]
rownames(scores_file)<-1:dim(scores_file)[1]
allGenes1<-allGenes[probes_allGenes %in% scores_file$ProbeID]
if(num==4){
  to_sub<-which(grepl("x_at",rownames(allGenes))==FALSE & grepl("s_at",rownames(allGenes))==FALSE)
}
else if(num!=4){
to_sub<-select_probe(allGenes1,scores_file)
}
allGenes1<-allGenes[to_sub]
return(allGenes1)}

allGenes1<-probes_scores(dir_to_scoresfile,3)

# batch
batch = pd$DataSetAccesionNumber
mod = model.matrix(~as.factor(Trimester), data=pData(pd))
combat_edata = ComBat(dat=exprs(allGenes1), batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
exprs(allGenes1) <- combat_edata
#arrayQualityMetrics(expressionset = allGenes, outdir = paste("QC_Report_RMA_ComBat_Trimester_V2",add,collapse="",sep=""), force = TRUE, do.logtransform = FALSE, intgroup = "Trimester")
# DiffExp genes
# First and Third
ind13 <- colnames(exprs(allGenes1)) %in% pd[pd$Trimester!="Second"]$SampleAccessionNumber
exprs13 <- exprs(allGenes1)[,ind13]
types = factor(pd[ind13]$Trimester)
design = model.matrix(~ 0+types)
colnames(design) = levels(types)
contrast.matrix = makeContrasts(First-Third, levels=design)
fit = lmFit(exprs13,design)
fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)

probes = rownames(allGenes1)
fit2$genes$ProbeID<-probes
fit2$genes$Symbol = unlist(mget(probes,hgu133plus2SYMBOL,ifnotfound=NA))
fit2$genes$EntrezID = unlist(mget(probes,hgu133plus2ENTREZID,ifnotfound=NA))
fit2$genes$Description = unlist(mget(probes,hgu133plus2GENENAME,ifnotfound=NA))
fullList = topTable(fit2, number=dim(fit2)[1], sort.by="logFC", resort.by="p")
write.table(fullList, file="all_13trimester_results_minus_x_s_at", row.names=FALSE, sep="\t", quote=FALSE)
# First and Second
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
fit2$genes$Symbol = unlist(mget(probes,hgu133plus2SYMBOL,ifnotfound=NA))
fit2$genes$EntrezID = unlist(mget(probes,hgu133plus2ENTREZID,ifnotfound=NA))
fit2$genes$Description = unlist(mget(probes,hgu133plus2GENENAME,ifnotfound=NA))
fullList = topTable(fit2, number=dim(fit2)[1], sort.by="logFC", resort.by="p")
write.table(fullList, file="all_12trimester_results_minus_x_s_at.txt", row.names=FALSE, sep="\t", quote=FALSE)
#write.table(exprs(allGenes1),"expressionset.txt",sep="\t",row.names=T)
# Second and third

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

probes<-rownames(allGenes1)
fit2$genes$ProbeID<-probes
fit2$genes$Symbol = unlist(mget(probes,hgu133plus2SYMBOL,ifnotfound=NA))
fit2$genes$EntrezID = unlist(mget(probes,hgu133plus2ENTREZID,ifnotfound=NA))
fit2$genes$Description = unlist(mget(probes,hgu133plus2GENENAME,ifnotfound=NA))
fullList = topTable(fit2, number=dim(fit2)[1], sort.by="logFC", resort.by="p")
write.table(fullList, file="all_23trimester_results_minus_x_s_at.txt", row.names=FALSE, sep="\t", quote=FALSE)
#write.table(exprs(allGenes1),"expressionset.txt",sep="\t",row.names=T)

#############################################################################################################################################
###########################################################GO################################################################################
#######Load packages

library(biomaRt)
library(topGO)
#library(ALL)
library(XML)
library(Rgraphviz)

############################################ ********************************************* #####################################
# analyze GOterms enrichment

#function needed to create topGOdata object
geneSelect<-function(data){
  return(data < 0.01)
}

#Some add function
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# make a function for simplicity

GO_analysis<-function(file_name){
#load differential expression data
fit_file<-read.table(file_name,sep="\t",header=T)
#remove not annotated probes
fit_file<-fit_file[complete.cases(fit_file$ID.Symbol, fit_file$adj.P.Val),]
#do a list of genes
geneList<-as.double(as.numeric(as.character(fit_file$adj.P.Val)))
names(geneList)<-fit_file$ID.ProbeID #annotation is hgu13a3plus2.db
sampleGOdata <- new("topGOdata",description = "Simple session", ontology = "BP",allGenes = geneList, geneSel = geneSelect,nodeSize = 10,annot = annFUN.db, affyLib = "hgu133plus2.db")
#GOdata<-new("topGOdata",ontology = "BP", allGenes = e = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")  
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 50)
View(allRes)
#Results 
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol)
#outliers of classic/elim method
sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),elim = pValue.elim[sel.go],classic = pValue.classic[sel.go])
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

#Save data
write.table(allRes,paste(file_name,"_GO.txt",sep="",collapse=""),sep="\t",row.names=F)
pdf(paste(file_name,"_GO.pdf",sep="",collapse = ""))

dev.off()
}

#Run pipeline for each trimester

GO_analysis("all_13trimester_result.txt")
GO_analysis("all_13trimester_results_1st.txt")
GO_analysis("all_13trimester_results_2nd.txt")
GO_analysis("all_13trimester_results_jetset.txt")
GO_analysis("all_13trimester_results_minus_x_s_at")

GO_analysis("all_12trimester_results.txt")
GO_analysis("all_12trimester_results_1st.txt")
GO_analysis("all_12trimester_results_2nd.txt")
GO_analysis("all_12trimester_results_jetset.txt")
GO_analysis("all_12trimester_results_minus_x_s_at.txt")

GO_analysis("all_23trimester_results.txt")
GO_analysis("all_23trimester_results_1st.txt")
GO_analysis("all_23trimester_results_2nd.txt")
GO_analysis("all_23trimester_results_minus_x_s_at.txt")

####################PLOT DISTRIBUTION OF COEFFITIENTS AMONG SELECTED GENES 
#
all_files<-list.files()
expression_files<-all_files[which(grepl("expressionset",all_files)==TRUE)]

plotcoefval<-function(num){
  if(num==1){
    scores_file<-read.table("/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_1.txt",sep=",",header=T)
    allfiles<-all_files[which(grepl("all_[0-9][0-9]trimester_results_1st.txt$",all_files)==TRUE)]
    e<-"1st"}
  if(num==2){
    scores_file<-read.table("/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_2.txt",sep=",",header=T)
    allfiles<-all_files[which(grepl("all_[0-9][0-9]trimester_results_2nd.txt$",all_files)==TRUE)]
    e<-"2nd"}
  if(num==3){
    scores_file<-read.table("/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_Jetset.csv",sep=",",header=T)
    allfiles<-all_files[which(grepl("all_[0-9][0-9]trimester_results_jetset.txt$",all_files)==TRUE)]
    colnames(scores_file)[c(1,3,5,6,8)]<-c("ProbeID","GeneID","Spec_score","Cov_score","Overall_score")
    e<-"jetset"}
  pdf(paste("scoresexpr_",e,"_.pdf",sep="",collapse=""))
  for(filename in allfiles){
    if(grepl("12",filename)==TRUE){
      trimester<-"1st_2nd"}
    if(grepl("13",filename)==TRUE){
      trimester<-"1st_3rd"}
    if(grepl("23",filename)==TRUE){
      trimester<-"2nd_3rd"}
      file_expr<-read.table(filename,sep="\t",header=T)
      allprobeIDs<-as.character(file_expr$ID.ProbeID)
      topprobeIDs<-as.character(file_expr[which(file_expr$adj.P.Val<0.01),1])
      otherprobeIDs<-as.character(file_expr[which((file_expr$ID.ProbeID %in% topprobeIDs)==FALSE),1])
      par(mfrow=c(3,3))
      #Specifity Score
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% allprobeIDs,"Spec_score"])),xlab="Spec_score",main=paste(trimester,"all"))
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% topprobeIDs,"Spec_score"])),xlab="Spec_score",ylab=NULL,main=paste(trimester,"top"))
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% otherprobeIDs,"Spec_score"])),xlab="Spec_score",ylab=NULL,main=paste(trimester,"other"))
      #Coverage score
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% allprobeIDs,"Cov_score"])),xlab="Cov_score",main=paste(trimester,"all"))
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% topprobeIDs,"Cov_score"])),xlab="Cov_score",ylab=NULL,main=paste(trimester,"top"))
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% otherprobeIDs,"Cov_score"])),xlab="Cov_score",ylab=NULL,main=paste(trimester,"other"))
      #Overall Score
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% allprobeIDs,"Overall_score"])),xlab="Overall_score",main=paste(trimester,"all"))
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% topprobeIDs,"Overall_score"])),xlab="Overall_score",ylab=NULL,main=paste(trimester,"top"))
      hist(as.numeric(as.character(scores_file[scores_file$ProbeID %in% otherprobeIDs,"Overall_score"])),xlab="Overall_score",ylab=NULL,main=paste(trimester,"other")) 
  }
  dev.off()
}

plotcoefval(1)
plotcoefval(2)
plotcoefval(3)


