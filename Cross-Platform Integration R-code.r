#Cross-Platform Integration Data Analysis Illumina (L3) Affymetrix (L2)
#Arran Turnbull - 18/01/2011

#######################################################################################################################################################################
#Pre-processing of Illumina Data
#######################################################################################################################################################################

#Data background corrected in GenomeStudio

#First change directory in R

#Load required programs in R
library(lumi)

#Load data using lumiR function - do not convert IDs
filename="L3_bg_SPPraw.txt"
raw.lumi <- lumiR(filename, convertNuID = FALSE)

#View a box plot of expression for the raw data
boxplot(raw.lumi, main="Amplitude of Expression Plot - L3 - Raw Data")

#Log2 Transformation (Log2)
#For other transformation methods change "log2" for: "vst", "cubicRoot"
lumi.T=lumiT(raw.lumi, method=c("log2"), verbose=TRUE)

#Quantile Normalisation (QN)
#For other normalisation methods change "quantile" for: "rsn", "ssn", "loess", "vsn", "rankinvariant"
lumi.N=lumiN(lumi.T, method=c("quantile"), verbose=TRUE)

#Perform quality control
lumi.N.Q=lumiQ(lumi.N)
summary(lumi.N.Q, 'QC')

#View a box plot of expression for the normalised and filtered data
boxplot(lumi.N.Q, main="Amplitude of Expression Plot - L3 - QN")

#Create an expression matrix of the normalised and filtered data
data.processed <- exprs(lumi.N.Q)

#Write data to directory as a .txt file (make sure that row names and column names are displayed correctly in file)
write.table(data.processed,file="L3_bg_QN.txt", append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

#Create multidimensional scaling (MDS) plot
dataMatrix <- data.processed
d1=dist(t(scale(dataMatrix)))
dN=dimnames(dataMatrix)[[2]]
nS=length(dN)
d1M=as.matrix(d1)
dimnames(d1M)=list(dN, dN)
d10.2.mds=cmdscale(d1M, k=2)
x=d10.2.mds[,1]
y=d10.2.mds[,2]
plot(x,y)
z=cbind(x,y)
write.table(z, "L3_bg_QN_mds.txt", quote=F, sep="\t")

#Load new data with failed samples filtered out
data.QC3=read.delim("L3_bg_QN.txt",header=TRUE, sep="\t", row.names=1)

#Create multidimensional scaling (MDS) plot
dataMatrix <- data.combat
d1=dist(t(scale(dataMatrix)))
dN=dimnames(dataMatrix)[[2]]
nS=length(dN)
d1M=as.matrix(d1)
dimnames(d1M)=list(dN, dN)
d10.2.mds=cmdscale(d1M, k=2)
x=d10.2.mds[,1]
y=d10.2.mds[,2]
plot(x,y)
z=cbind(x,y)
write.table(z, "L3B1Hyb2QC3_bg_QN_mds.txt", quote=F, sep="\t")

#Drawing similarity/dissimilairy heatmap in R using Gplots
#Determine correlation between columns of dataset
data.cor <- cor(data.QC3, method="pearson")
library(gplots)
heatmap.2(data.cor, col=heat.colors(75), dendrogram="both",Rowv=TRUE, Colv=TRUE, symm=TRUE, trace="none")
write.table(data.cor, "L3B1Hyb2QC3_bg_QN_cor.txt", quote=F, sep="\t")

#Data filtered by probe detection p-value

#######################################################################################################################################################################
#Perform reannoation of the Sample_Probe_Profile to Ensembl Genes and combine expression values for multiple probes for a given gene to give a New_Sample_Gene_Profile
#######################################################################################################################################################################

#Load new data with failed samples filtered out
data.fullyprocessed=read.delim("L3_bg_QN_filtered(genesENSM).txt",header=TRUE, sep="\t", row.names=1)

#View a box plot of expression for the normalised, QC and filtered data
boxplot(data.fullyprocessed, main="Amplitude of Expression Plot - L3_bg_QN_filtered(genesENSM)")

#Create multidimensional scaling (MDS) plot
dataMatrix <- data.processed
d1=dist(t(scale(dataMatrix)))
dN=dimnames(dataMatrix)[[2]]
nS=length(dN)
d1M=as.matrix(d1)
dimnames(d1M)=list(dN, dN)
d10.2.mds=cmdscale(d1M, k=2)
x=d10.2.mds[,1]
y=d10.2.mds[,2]
plot(x,y)
z=cbind(x,y)
write.table(z, "L3_bg_QN_filtered(genesENSM)_mds.txt", quote=F, sep="\t")

#######################################################################################################################################################################
#Pre-processing of Affymetrix Data
#######################################################################################################################################################################

#First change directory in R

#Load required programs in R
library(affy)

#Read .CEL files into R and load the required alternative_cdf file
#The alternative_cdf file will automatically download from the Bioconductor Brain Array Repository
#Alternative_cdf version 12 for Ensembl Gene IDs has been used here
data<-ReadAffy()
data@cdfName<-"HGU133A_Hs_ENSG"

#View a box plot of expression for the raw data
boxplot(data, main="Amplitude of Expression Plot - Non-Normalised Raw Data")

#mas5 normlisation
data.mas5=mas5(data)
write.exprs(data.mas5, file="L2_Affy_mas5.txt")

#Load mas5 Affy data
data.processed=read.delim("L2_Affy_mas5.txt",header=TRUE, sep="\t", row.names=1)

#Log2 transformation
L2_Affy_mas5_Log2<-log(data.processed,base=2)

#View a box plot of expression for the normalised data
boxplot(L2_Affy_mas5_Log2, main="Amplitude of Expression Plot - Normalised and Log2 Data")

#Write mas5 Log2 Affy data to file
write.table(L2_Affy_mas5_Log2, "L2_Affy_mas5_Log2.txt", sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

#mas5 detection calls
calls<-mas5calls(data)
write.exprs(calls, file="L2_Affy_DetCalls.txt")

#Filter based on detection calls
#Present: p < 0.04
#Marginal: 0.04 = p = 0.06
#Absent: p > 0.06

#View a box plot of expression for the normalised and filtered data
data.filtered=read.delim("L2_Affy_mas5_Log2_filtered.txt",header=TRUE, sep="\t", row.names=1)
boxplot(data.filtered, main="Amplitude of Expression Plot - Affymetrix - Filtered, mas5 Normalised")

#######################################################################################################################################################################
#L2 Affymetrix and L3 Illumina Combined Data
#######################################################################################################################################################################

#Load L2 Affy and L3 Illumina combined data
data.combined=read.delim("L2&L3_AffyIllumCombined.txt",header=TRUE, sep="\t", row.names=1)

#View a box plot of expression for the combined dta
boxplot(data.combined, main="Amplitude of Expression Plot - Combined Data")

#Create multidimensional scaling (MDS) plot
dataMatrix <- data.combined
d1=dist(t(scale(dataMatrix)))
dN=dimnames(dataMatrix)[[2]]
nS=length(dN)
d1M=as.matrix(d1)
dimnames(d1M)=list(dN, dN)
d10.2.mds=cmdscale(d1M, k=2)
x=d10.2.mds[,1]
y=d10.2.mds[,2]
plot(x,y)
z=cbind(x,y)
write.table(z, "L2&L3_AffyIllumCombined_mds.txt", quote=F, sep="\t")

################################################################
#XPN correction 
################################################################

#Data correction by mean centering (Cluster 3)

################################################################
#XPN correction 
################################################################

#Data correction by EB ComBat (ArrayMiner.net)

#Load L2 Affy and L3 Illumina combined and XPN corrected data
data.xpn=read.delim("L2&L3_AffyIllumCombined_XPN.txt",header=TRUE, sep="\t", row.names=1)

#View a box plot of expression for the combined dta
boxplot(data.xpn, main="Amplitude of Expression Plot - Combined and XPN Data")

#Create multidimensional scaling (MDS) plot
dataMatrix <- data.xpn
d1=dist(t(scale(dataMatrix)))
dN=dimnames(dataMatrix)[[2]]
nS=length(dN)
d1M=as.matrix(d1)
dimnames(d1M)=list(dN, dN)
d10.2.mds=cmdscale(d1M, k=2)
x=d10.2.mds[,1]
y=d10.2.mds[,2]
plo

#Drawing similarity/dissimilairy heatmap in R using Gplots
#Determine correlation between columns of dataset
data.cor <- cor(data.xpn, method="pearson")
library(gplots)
heatmap.2(data.cor, col=heat.colors(75), dendrogram="both",Rowv=TRUE, Colv=TRUE, symm=TRUE, trace="none")
write.table(data.cor, "L2&L3_AffyIllumCombined_XPN_cor.txt", quote=F, sep="\t")

################################################################
#Optional EB ComBat correction 
################################################################

#Data correction by EB ComBat (ArrayMiner.net) 

#Load L2 Affy and Illumina combined and XPN corrected data
data.eb=read.delim("L2&L3_AffyIllumCombined_ComBat.txt",header=TRUE, sep="\t", row.names=1)

#View a box plot of expression for the combined dta
boxplot(data.eb, main="Amplitude of Expression Plot - Combined and ComBat Data")

#Create multidimensional scaling (MDS) plot
dataMatrix <- data.eb
d1=dist(t(scale(dataMatrix)))
dN=dimnames(dataMatrix)[[2]]
nS=length(dN)
d1M=as.matrix(d1)
dimnames(d1M)=list(dN, dN)
d10.2.mds=cmdscale(d1M, k=2)
x=d10.2.mds[,1]
y=d10.2.mds[,2]
plot(x,y)
z=cbind(x,y)
write.table(z, "L2-(1&3-QSR)-rmaENSM-filtered_MC_mds.txt", quote=F, sep="\t")

#Drawing similarity/dissimilairy heatmap in R using Gplots
#Determine correlation between columns of dataset
data.cor <- cor(data.eb, method="pearson")
library(gplots)
heatmap.2(data.cor, col=heat.colors(75), dendrogram="both",Rowv=TRUE, Colv=TRUE, symm=TRUE, trace="none")

################################################################
#Optional DWD-standard correction 
################################################################

#Data correction by DWD-standard 

#Load L2 Affy and Illumina combined and XPN corrected data
data.dwd=read.delim("L2&L3_AffyIllumCombined_DWDstandard.txt",header=TRUE, sep="\t", row.names=1)

#View a box plot of expression for the combined dta
boxplot(data.dwd, main="Amplitude of Expression Plot - Combined and DWD-Standard Data")

#Create multidimensional scaling (MDS) plot
dataMatrix <- data.dwd
d1=dist(t(scale(dataMatrix)))
dN=dimnames(dataMatrix)[[2]]
nS=length(dN)
d1M=as.matrix(d1)
dimnames(d1M)=list(dN, dN)
d10.2.mds=cmdscale(d1M, k=2)
x=d10.2.mds[,1]
y=d10.2.mds[,2]
plot(x,y)
z=cbind(x,y)
write.table(z, "L2&L3_AffyIllumCombined_DWD-Standard_mds.txt", quote=F, sep="\t")

#Drawing similarity/dissimilairy heatmap in R using Gplots
#Determine correlation between columns of dataset
data.cor <- cor(data.dwd, method="pearson")
library(gplots)
heatmap.2(data.cor, col=heat.colors(75), dendrogram="both",Rowv=TRUE, Colv=TRUE, symm=TRUE, trace="none")
write.table(data.cor, "L2&L3_AffyIllumCombined_DWD-Standard_cor.txt", quote=F, sep="\t")

################################################################
#Correlation Analysis
################################################################

#Setup Heatmap and Colours

library(gplots)

breaks <- seq(from = 0.85, to = 1, length = 20)
colors = "cm.colors"                                

#Load Data	

data.combat=read.table("L2-IllumAffySubsetCombined(probesENSM)AffyCor-filtered-QN-ComBat(No40)_coradapted.txt")
data.dwd=read.table("L2-IllumAffySubsetCombined(probesENSM)AffyCor-filtered-QN-DWD(No40)_coradapted.txt")
data.mc=read.table("L2-IllumAffySubsetCombined(probesENSM)AffyCor-filtered-QN-MC(No40)_coradapted.txt")
data.xpn=read.table("L2-IllumAffySubsetCombined(probesENSM)AffyCor-filtered-QN-XPN(No40)_coradapted.txt")
  
data.mc.cor<-as.matrix(data.mc)
data.combat.cor<-as.matrix(data.combat)
data.dwd.cor<-as.matrix(data.dwd)
data.xpn.cor<-as.matrix(data.xpn)
  
#Correlation Analysis and Heatmap
  
heatmap.2(data.mc.cor, dendrogram="none", col=colors, Rowv=F, Colv=F, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1, breaks=breaks) 
heatmap.2(data.combat.cor, dendrogram="none", col=colors, Rowv=F, Colv=F, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1, breaks=breaks)
heatmap.2(data.dwd.cor, dendrogram="none", col=colors, Rowv=F, Colv=F, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1, breaks=breaks)
heatmap.2(data.xpn.cor, dendrogram="none", col=colors, Rowv=F, Colv=F, scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=1, cexCol=1, breaks=breaks)

