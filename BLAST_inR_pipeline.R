### RUN BLAST AND SCORE AFFYMETRIX PROBES
### 1. Download annotation data file "HG-U133_Plus2.full.flat" from http://affymetrix2.bioinf.fbb.msu.ru/files.html
###  
file_ann<-read.delim("/home/vlad/gene_regulatory network project/BLAST/HG-U133_Plus2.full.flat",sep="\t",header=T)
sequences<-as.character(file_ann$affy_probe_seq)
### BLAST runs very long with large input files. To speed it up, we devide input sequence data into 60 
### files with ~ 10000 sequences in one
vec<-vector() 
for(i in sequences){
  vec<-append(vec,strsplit(i,""))
}
max<-length(vec)
for(i in 0:60){
  from<-i*10000+1
  to<-(i+1)*10000
  if(i==60){
    to<-max
  }
  write.fasta(sequences = vec[from:to],names = as.character(file_ann$affy_probe_id)[from:to],paste("/home/vlad/gene_regulatory network project/BLAST/probe_seq.fasta_",from,"_",to,sep="",collapse=""))
}

### Install BLAST. Short manual: http://www.ncbi.nlm.nih.gov/books/NBK1763/
### Download all RefSeq data from NCBI server
### ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/
### data is splited into all these files in directory
### Download rna-only
### Create local database for BLAST (see manual above)
### then run script_BLAST.sh

########################################### BLAST PIPELINE #####################################################
############# First script to run script_BLAST.sh
### LIST_FILES="$( ls -l input_files/probe_seq.fasta_ | awk -F '/' '{print $NF}' )"
### for file_name in $LIST_FILES;do
### bash script_BLAST.sh $file_name
### done 
###
###
###
###  Second script: run BLAST
###  run script_BLAST.sh
###  cd ~/gene_regulatory\ network\ project/BLAST/
###  blastall -i ./input_files/$1 -d human.rna.fna -p blastn -o out1.txt
###  Very useful. convert output data to csv
###  sudo python /root/ngs-scripts/blast/blast-to-csv.py out1.txt > ./output_files/out1_$1.csv
###  Next, we work with these tables only
####################################### END OF BLAST PIPELINE ################################################# 

### Two ways, first we can use biomart for transcript/gene ID translation. But note, its server is almost
### always unavailable. That is why, second (boring) way: save transcript IDs to different files. Then load
### it manually online and get a result
### First way:

library(biomaRt)
listMarts()

### run a proxy
#options(RCurlOptions = list(proxy="uscache.kcc.com:80",proxyuserpwd="------:-------"))
#ensembl=useMart("ensembl")
setwd("/home/vlad/gene_regulatory network project/BLAST/output_files/")
files_input<-list.files("/home/vlad/gene_regulatory network project/BLAST/output_files/")
 
# file_input<-read.table(files_input[1],sep=",",header=F)
# function, need it further


### Second way: 

to_write<-function(DD,dir){
  write.table(DD,paste(dir,file_name,"_",deparse(substitute(DD)) ,"_",".txt",sep="",collapse=""),sep = "\n",quote = F,row.names = F,col.names = F)
}

for(file_name in files_input[length(files_input)]){
  print(file_name)
  # load file
  file_input<-read.table(file_name,sep=",",header=F)
  # separate transcript IDs
  transcript_ids<-unlist(strsplit(paste(as.character(file_input$V2),sep="",collapse=""),"[|]"))
  # modify a little
  transcript_ids_dot_minus<-unlist(strsplit(transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)],"[.]"))
  transcript_ids_dot_minus<-transcript_ids_dot_minus[1:length(transcript_ids_dot_minus)/2!=1:length(transcript_ids_dot_minus)%/%2]
  transcript_ids<-transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)]
  #add it to data
  file_input$transcript_ids<-transcript_ids
  file_input$transcript_ids_min<-transcript_ids_dot_minus
  #separate 6 classes of IDs
  NM<-transcript_ids_dot_minus[which(grepl("NM",transcript_ids_dot_minus)==TRUE)]
  XM<-transcript_ids_dot_minus[which(grepl("XM",transcript_ids_dot_minus)==TRUE)]
  NR<-transcript_ids_dot_minus[which(grepl("NR",transcript_ids_dot_minus)==TRUE)]
  XR<-transcript_ids_dot_minus[which(grepl("XR",transcript_ids_dot_minus)==TRUE)]
  NP<-transcript_ids_dot_minus[which(grepl("NP",transcript_ids_dot_minus)==TRUE)]
  XP<-transcript_ids_dot_minus[which(grepl("XP",transcript_ids_dot_minus)==TRUE)]
  #save this data
  to_write(NM,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(XM,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(NR,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(XR,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(NP,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(XP,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  #for(i in 1:dim(file_input)[1]){
  # file_input[i,"transcript_id"]<-unlist(strsplit(as.character(file_input[i,2]),"[|]"))[4]
  #}
}
# save IDs to one file
for(DD in c("NM","XM","NR","XR","NP","XP")){
  system(paste("cat ","/home/vlad/gene_regulatory\\ network\\ project/BLAST/ID_files/*",DD,"* > ","/home/vlad/gene_regulatory\\ network\\ project/BLAST/ID_files/",DD,"_data.txt",sep = "",collapse = ""))
}

# files are too large for biomart http://central.biomart.org/converter/#!/ID_converter/gene_ensembl_config_2 , so 
# we need to write only unique values

setwd("/home/vlad/gene_regulatory network project/BLAST/ID_files//")
for(DD in c("NM","XM","NR","XR")){
  file<-read.table(paste(DD,"_data.txt",sep="",collapse=""),sep = "\t",header = F)
  file_ID_lists<-unique(file$V1)
  write.table(file_ID_lists,paste(DD,"_data.txt",sep="",collapse=""),sep = "\n",quote = F,row.names = F,col.names = F)
}

# Then convert ID to Entrez with biomart http://central.biomart.org/converter/#!/ID_converter/
# As the outcome we have 4 files with converted names. merge them in one file for simplicity.
# Can do it by "cat .. " in shell
### OUTPUT FILE here is Transcript_ID_to_Entrez.txt

####                                              END OF DATA PREPARATION PART                       ###

################################################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###############################
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@######################@@@@@@@@@@@@@@@################@@@@@@@@@@@@@@@@@@@@@@
  
# Load package

library(seqinr)

setwd("/home/vlad/gene_regulatory network project/BLAST/output_files/")
# load oligo sequences 
file_ann<-read.delim("/home/vlad/gene_regulatory network project/BLAST/HG-U133_Plus2.full.flat",sep="\t",header=T)
# load annotation
# for Affy
annotation<-read.delim("/home/vlad/gene_regulatory network project/annotation_data/GPL570-13270.csv",sep="\t",header=T)
annotation<-annotation[,c("ID","ENTREZ_GENE_ID")]
annotation$ID<-as.character(annotation$ID)
annotation$ENTREZ_GENE_ID<-as.character(annotation$ENTREZ_GENE_ID)
# PREPARE ANNOTATION (possible, but not needed)
# remove non-specific IDs
# annotation<-annotation[which(grepl("_a_at",annotation$ID)==TRUE | grepl("[0-9]_at",annotation$ID)==TRUE),]
# remove empty values
annotation<-annotation[complete.cases(annotation),]
# ANNOTATION IS READY

Tr_ID_to_Entrez<-read.table("/home/vlad/gene_regulatory network project/BLAST/ID_files/Transcript_ID_to_Entrez.txt",sep="\t",header=T) 
# convert to unique values
e<-vector()
for(id in as.character(unique(Tr_ID_to_Entrez$RefSeq.mRNA..e.g..NM_001195597.))){
  sub<-paste("_",Tr_ID_to_Entrez[which(Tr_ID_to_Entrez$RefSeq.mRNA..e.g..NM_001195597.==id),2],"_",sep="",collapse="")
  e<-append(e,sub)
}
Tr_ID_to_Entrez1<-as.data.frame(cbind(as.character(unique(Tr_ID_to_Entrez[,1])),e))
colnames(Tr_ID_to_Entrez1)<-c("Tr_ID","Entrez_ID")
### ***************** ###
      # MAIN PART #

out_file<-"/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_new_cov_score_from_sub.txt"
for(file_name in files_input[1:length(files_input)]){
  print(file_name)
  # load file
  file_input<-read.table(file_name,sep=",",header=F)
  # all probenames present 
  all_probenames0<-as.character(unique(file_input$V1))
  all_probenames<-vector()
  for(i in all_probenames0){
    i0<-paste(unlist(strsplit(i,"at"))[1] ,"at",sep="",collapse = "")
    all_probenames<-append(all_probenames,i0)
  }
  ###### data modification part
  # filter data
  file_input<-file_input[which(file_input$V3>32),]
  # separate transcript IDs from other trash
  transcript_ids<-unlist(strsplit(paste(as.character(file_input$V2),sep="",collapse=""),"[|]"))
  ### for Affy platform only:
  # modify a little
  transcript_ids_dot_minus<-unlist(strsplit(transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)],"[.]"))
  transcript_ids_dot_minus<-transcript_ids_dot_minus[1:length(transcript_ids_dot_minus)/2!=1:length(transcript_ids_dot_minus)%/%2]
  transcript_ids<-transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)]
  #add it to data
  file_input$transcript_ids<-transcript_ids
  file_input$transcript_ids_min<-transcript_ids_dot_minus
  # add Entrez ID data
  for(i in 1:dim(Tr_ID_to_Entrez1)[1]){
    tr_id<-as.character(Tr_ID_to_Entrez1[i,1])
    entrez_id<-as.character(Tr_ID_to_Entrez1[i,2])
    to_sub<-which(file_input$transcript_ids_min==tr_id)
    file_input[to_sub,"gene_ID"]<-rep(entrez_id,length(to_sub))
  }
  probeset_names<-vector()
  for(i in 1:dim(file_input)[1]){
    id1<-as.character(file_input[i,1])
    id2<-paste(unlist(strsplit(id1,"at"))[1],"at",sep="",collapse = "")
    probeset_names<-append(probeset_names,id2)
  }
  file_input$probeset_names<-probeset_names
  ##### end of data modification part
  # select probeset IDs
  # are present
  file_input<-file_input[complete.cases(file_input),]
  list_of_probes<-as.character(file_input$probeset_names)
  # all
  probeset_IDs<-annotation[which(annotation$ID %in% list_of_probes),1]
  ### 
  for(probeset_ID in probeset_IDs){ #for each probeset, not a probe
    sub<-file_input[which(file_input$probeset_names==probeset_ID),]
    sub_probes<-as.character(unique(sub$V1))
    num_of_sub_probes<-length(sub_probes)
    num_all<-length(which(all_probenames==probeset_ID))
    if(num_of_sub_probes/num_all>=0.5){ ### Probeset should be represented with at least 50% of mid-best
      probe_specific<-vector()          ### alignments
      genes_specific<-vector()
      for(probe in sub_probes){
        sub1<-sub[which(sub$V1==probe),]
        genes<-unique(sub1$gene_ID)
        if(length(genes)==1){
          probe_specific<-append(probe,probe_specific) ### All specifically detected genes
          genes_specific<-append(genes_specific,genes)
      }
    }
    targeted_genes<-summary(as.factor(genes_specific)) 
    max<-targeted_genes[which.max(targeted_genes)] ### Select a targeted gene with the highest frequency  
    if(length(targeted_genes[which(as.numeric(targeted_genes)==as.numeric(max))])==1){ ### if can be detected
      targeted_gene<-unlist(strsplit(names(max),"_"))[2]
      sub_tar_gene<-sub[which(sub$gene_ID==names(max)),]
      sub_tar_gene<-sub_tar_gene[which(as.character(sub_tar_gene$V1) %in% probe_specific),]
      specifity_score<-0
      pr_IDs<-unique(as.character(sub_tar_gene[,1]))
      for(pr_ID in pr_IDs){
        specifity_score<-specifity_score+max(sub_tar_gene[which(sub_tar_gene[,1]==pr_ID),3])/52
      }
      #tar_transcripts<-as.character(unique(sub_tar_gene$transcript_ids_min)) from previous
      tar_transcripts<-as.character(unique(sub_tar_gene$transcript_ids_min)) #new one
      #tr_freq<-vector() from previous
      all_transcripts<-as.character(unique(Tr_ID_to_Entrez1[which(grepl(targeted_gene,Tr_ID_to_Entrez1$Entrez_ID)==TRUE),1]))
      tar_transcripts<-tar_transcripts[which(tar_transcripts %in% all_transcripts==TRUE)]
    #for(abla in tar_transcripts){
    #  e<-length(unique(sub[which(sub$transcript_ids_min==abla),1]))/num_all
    #  tr_freq<-append(tr_freq,e)
    #}
   #tar_transcripts<-tar_transcripts[which(tr_freq>0.5)]
   specifity_score<-specifity_score/num_all
   #coverage_score<-sum(tr_freq)/length(all_transcripts)
   coverage_score<-length(tar_transcripts)/length(all_transcripts) #new one
   if(length(coverage_score)==0){
     coverage_score<-0
   }
   if(length(all_transcripts)==0){
     coverage_score<-"cannot define"
   }
   #print(specifity_score)
  # save results
  line_to_save<-c(probeset_ID,targeted_gene,coverage_score,specifity_score,coverage_score*specifity_score)
  if(length(line_to_save)==5){
    write(paste(line_to_save,sep="\t",collapse = " "),out_file,sep = "\n",append = T)
  }
  }}}
}


######################### To check with ready annotation data

library(hgu133plus2.db)
library(annotate)

Affy_data<-read.table("/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_new_cov_score_from_sub.txt",sep=" ",header=F)
RefTranscript_ID<-Affy_data$V1
Entrez_probe_ID<-select(hgu133plus2.db, as.character(RefTranscript_ID), "ENTREZID", "PROBEID")
RefTranscript_ID<-as.character(RefTranscript_ID)
e<-vector()
for(probe_ID in RefTranscript_ID){
  a<-as.character(Affy_data[which(RefTranscript_ID==probe_ID),2])
  b<-as.character(Entrez_probe_ID[which(Entrez_probe_ID$PROBEID==probe_ID),2])
  if(is.na(a)==FALSE){
    if(is.na(b)==FALSE){
  if(a==b){
    e<-append(e,"TRUE")
  }
}}
}
}



################################################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###############################
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@######################@@@@@@@@@@@@@@@################@@@@@@@@@@@@@@@@@@@@@@
### ANALYZE DATA ###

Affy_data<-read.table("/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_new_cov_score_from_sub.txt",sep=" ",header=F)
#Illumina_data<-read.table("/home/vlad/gene_regulatory network project/BLAST/ILLUMINA_HumanHT-12 V4.0 expression beadchip_IDs.txt",sep=" ",header=F)
x_at<-Affy_data[which(grepl("x_at",Affy_data$V1)==TRUE),]
a_at<-Affy_data[which(grepl("a_at",Affy_data$V1)==TRUE),]
at<-Affy_data[which(grepl("[0-9]_at",Affy_data$V1)==TRUE),]
s_at<-Affy_data[which(grepl("s_at",Affy_data$V1)==TRUE),]

# Plot cumulative distrubutions
plot_cum_dist<-function(e_at,num,name1,add){
  name="score"
  if(num==3){
    coloor=colors()[436]
  }
  if(num==4){
    coloor=colors()[369]
  }
  if(num==5){
    coloor=colors()[371]
  }
  print(summary(fivenum(e_at[,num])))
  sorted<-sort(as.numeric(as.character(e_at[,num])))
  n<-length(sorted)
  if(add==TRUE){
    lines(sorted, (1:n)/n, type = 's', ylim = c(0, 1),xlab="Sample Quantiles",ylab=name,main=name1,col = coloor,lwd = 3)
  }
  if(add==FALSE){
    plot(sorted, (1:n)/n, type = 's', ylim = c(0, 1),xlab="Sample Quantiles",ylab=name,main=name1,col=coloor,lwd = 3)
  }
  legend(0.65,0.4,c("cov","spec","cov*spec"),col=c(colors()[436],colors()[369],colors()[371]),lty=c(1,1,1),lwd=c(2.5,2.5))
}


plot_cum_dist(x_at,3,"x_at",FALSE)
plot_cum_dist(x_at,4,"x_at",TRUE)
plot_cum_dist(x_at,5,"x_at",TRUE)

plot_cum_dist(a_at,3,"a_at",FALSE)
plot_cum_dist(a_at,4,"a_at",TRUE)
plot_cum_dist(a_at,5,"a_at",TRUE)

plot_cum_dist(at,3,"at",FALSE)
plot_cum_dist(at,4,"at",TRUE)
plot_cum_dist(at,5,"at",TRUE)

plot_cum_dist(s_at,3,"s_at",FALSE)
plot_cum_dist(s_at,4,"s_at",TRUE)
plot_cum_dist(s_at,5,"s_at",TRUE)

plot_cum_dist(Affy_data,3,"all",FALSE)
plot_cum_dist(Affy_data,4,"all",TRUE)
plot_cum_dist(Affy_data,5,"all",TRUE)





#####################                                                               ##################### 
##################                                  KEEP CALM                          ##################
##############                  $$$$$$$$$$$$$      THAT IS ALL    $$$$$$$$$$$$             ##############
##################                                                                     ##################
#####################                                                               #####################










### ALMOST. NEXT COMES SCRIPT FOR ILLUMINA. BUT THATS ANOTHER STORY








### ******************************************************************************************* ###
############################################## ILLUMINA ###########################################

file_ann<-read.delim("/home/vlad/gene_regulatory network project/annotation_data/GPL10558-11219.csv",sep="\t",header=T)

### DATA 
library(seqinr)

sequences<-as.character(file_ann$SEQUENCE)
vec<-vector()
for(i in sequences){
  vec<-append(vec,strsplit(i,""))
}
max<-length(vec)
for(i in 0:4){
  from<-i*10000+1
  to<-(i+1)*10000
  if(i==4){
    to<-max
  }
  write.fasta(sequences = vec[from:to],names = as.character(file_ann$ID)[from:to],paste("/home/vlad/gene_regulatory network project/BLAST/ILLUMINA_probe_seq.fasta_",from,"_",to,sep="",collapse=""))
}

### then run script_BLAST.sh

library(biomaRt)
listMarts()

setwd("/home/vlad/gene_regulatory network project/BLAST/output_files/")
files_input<-list.files("/home/vlad/gene_regulatory network project/BLAST/output_files/")
files_input<-files_input[which(grepl("ILLUMINA",files_input)==TRUE)]

# file_input<-read.table(files_input[1],sep=",",header=F)
# function, need it further
to_write<-function(DD,dir){
  write.table(DD,paste(dir,file_name,"_",deparse(substitute(DD)) ,"_",".txt",sep="",collapse=""),sep = "\n",quote = F,row.names = F,col.names = F)
}
for(file_name in files_input){
  print(file_name)
  # load file
  file_input<-read.table(file_name,sep=",",header=F)
  # separate transcript IDs
  transcript_ids<-unlist(strsplit(paste(as.character(file_input$V2),sep="",collapse=""),"[|]"))
  # modify a little
  transcript_ids_dot_minus<-unlist(strsplit(transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)],"[.]"))
  transcript_ids_dot_minus<-transcript_ids_dot_minus[1:length(transcript_ids_dot_minus)/2!=1:length(transcript_ids_dot_minus)%/%2]
  transcript_ids<-transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)]
  #add it to data
  file_input$transcript_ids<-transcript_ids
  file_input$transcript_ids_min<-transcript_ids_dot_minus
  #separate 6 classes of IDs
  NM<-transcript_ids_dot_minus[which(grepl("NM",transcript_ids_dot_minus)==TRUE)]
  XM<-transcript_ids_dot_minus[which(grepl("XM",transcript_ids_dot_minus)==TRUE)]
  NR<-transcript_ids_dot_minus[which(grepl("NR",transcript_ids_dot_minus)==TRUE)]
  XR<-transcript_ids_dot_minus[which(grepl("XR",transcript_ids_dot_minus)==TRUE)]
  NP<-transcript_ids_dot_minus[which(grepl("NP",transcript_ids_dot_minus)==TRUE)]
  XP<-transcript_ids_dot_minus[which(grepl("XP",transcript_ids_dot_minus)==TRUE)]
  #save this data
  to_write(NM,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(XM,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(NR,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(XR,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(NP,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  to_write(XP,"/home/vlad/gene_regulatory network project/BLAST/ID_files/")
  #for(i in 1:dim(file_input)[1]){
  # file_input[i,"transcript_id"]<-unlist(strsplit(as.character(file_input[i,2]),"[|]"))[4]
  #}
}
# save IDs to one file
for(DD in c("NM","XM","NR","XR","NP","XP")){
  system(paste("cat ","/home/vlad/gene_regulatory\\ network\\ project/BLAST/ID_files/out1_ILLUMINA*",DD,"* > ","/home/vlad/gene_regulatory\\ network\\ project/BLAST/ID_files/",DD,"_data.txt",sep = "",collapse = ""))
}

# files are too large for biomart http://central.biomart.org/converter/#!/ID_converter/gene_ensembl_config_2 , so 
# we need to write only unique values

setwd("/home/vlad/gene_regulatory network project/BLAST/ID_files//")
for(DD in c("NM","XM","NR","XR")){
  file<-read.table(paste(DD,"_data.txt",sep="",collapse=""),sep = "\t",header = F)
  file_ID_lists<-unique(file$V1)
  write.table(file_ID_lists,paste(DD,"_data.txt",sep="",collapse=""),sep = "\n",quote = F,row.names = F,col.names = F)
}

# Then convert ID to Entrez with biomart http://central.biomart.org/converter/#!/ID_converter/
# As the outcome we have 4 files with converted names. merge them in one file for simplicity.
# Can do it by "cat .. " in shell

####                                              END OF DATA PREPARATION PART                       ###
# Load 

#package
library(seqinr)
setwd("/home/vlad/gene_regulatory network project/BLAST/output_files/")
# load oligo sequences 
file_ann<-read.delim("/home/vlad/gene_regulatory network project/BLAST/GPL10558-11219.txt",sep="\t",header=T)
# load annotation
# for Affy
annotation<-read.delim("/home/vlad/gene_regulatory network project/annotation_data/GPL10558-11219.csv",sep="\t",header=T)
annotation<-annotation[,c("ID","Entrez_Gene_ID")]
annotation$ID<-as.character(annotation$ID)
#annotation$ENTREZ_GENE_ID<-as.character(annotation$Entrez_Gene_ID)
# PREPARE ANNOTATION
# remove non-specific IDs
#annotation<-annotation[which(grepl("_a_at",annotation$ID)==TRUE | grepl("[0-9]_at",annotation$ID)==TRUE),]
# remove empty values
annotation<-annotation[complete.cases(annotation),]
# ANNOTATION IS READY

# 
Tr_ID_to_Entrez<-read.table("/home/vlad/gene_regulatory network project/BLAST/ID_files/Transcript_ID_to_Entrez_ILLUMINA.txt",sep="\t",header=T) 
# convert to unique values
e<-vector()
for(id in as.character(unique(Tr_ID_to_Entrez$RefSeq.mRNA..e.g..NM_001195597.))){
  sub<-paste("_",Tr_ID_to_Entrez[which(Tr_ID_to_Entrez$RefSeq.mRNA..e.g..NM_001195597.==id),2],"_",sep="",collapse="")
  e<-append(e,sub)
}
Tr_ID_to_Entrez1<-as.data.frame(cbind(as.character(unique(Tr_ID_to_Entrez$RefSeq.mRNA..e.g..NM_001195597.)),e))
colnames(Tr_ID_to_Entrez1)<-c("Tr_ID","Entrez_ID")
### ***************** ###
# MAIN PART #

out_file<-"/home/vlad/gene_regulatory network project/BLAST/ILLUMINA_HumanHT-12 V4.0 expression beadchip_IDs.txt"
for(file_name in files_input[2:5]){
  print(file_name)
  # load file
  file_input<-read.table(file_name,sep=",",header=F)
  # all probenames present 
  all_probenames0<-as.character(unique(file_input$V1))
  all_probenames<-vector()
  # filter data
  file_input<-file_input[which(file_input$V3>64),]
  # separate transcript IDs
  transcript_ids<-unlist(strsplit(paste(as.character(file_input$V2),sep="",collapse=""),"[|]"))
  # modify a little
  transcript_ids_dot_minus<-unlist(strsplit(transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)],"[.]"))
  transcript_ids_dot_minus<-transcript_ids_dot_minus[1:length(transcript_ids_dot_minus)/2!=1:length(transcript_ids_dot_minus)%/%2]
  transcript_ids<-transcript_ids[which(1:length(transcript_ids)/4==1:length(transcript_ids)%/%4)]
  #add it to data
  file_input$transcript_ids<-transcript_ids
  file_input$transcript_ids_min<-transcript_ids_dot_minus
  # add Entrez ID data
  for(i in 1:dim(Tr_ID_to_Entrez1)[1]){
    tr_id<-as.character(Tr_ID_to_Entrez1[i,1])
    entrez_id<-as.character(Tr_ID_to_Entrez1[i,2])
    to_sub<-which(file_input$transcript_ids_min==tr_id)
    file_input[to_sub,"gene_ID"]<-rep(entrez_id,length(to_sub))
  }
  probeset_names<-file_input$V1
  file_input$probeset_names<-probeset_names
  # select probeset IDs
  # are present
  file_input<-file_input[complete.cases(file_input),]
  list_of_probes<-as.character(file_input$probeset_names)
  # all
  probeset_IDs<-annotation[which(annotation$ID %in% list_of_probes),1]
  for(probeset_ID in probeset_IDs){
    sub<-file_input[which(file_input$probeset_names==probeset_ID),]
    genes<-unique(sub$gene_ID)
  if(length(genes)==1){
    targeted_gene_0<-genes
    max<-mean(sub$V3)
    targeted_gene<-unlist(strsplit(targeted_gene_0,"_"))[2]
    sub_tar_gene<-sub[which(sub$gene_ID==targeted_gene_0),]
    tar_transcripts<-as.character(unique(sub_tar_gene$transcript_ids_min))
    tr_freq<-vector()
    all_transcripts<-as.character(unique(Tr_ID_to_Entrez1[which(grepl(targeted_gene_0,Tr_ID_to_Entrez1$Entrez_ID)==TRUE),1]))
    specifity_score<- as.numeric(max)/100
    coverage_score<-length(which(all_transcripts %in% tar_transcripts ==TRUE))/length(all_transcripts)
      if(length(coverage_score)==0){
        coverage_score<-0
      }
      if(length(all_transcripts)==0){          
        coverage_score<-0
        }
        #print(specifity_score)
        # save results
        line_to_save<-c(probeset_ID,targeted_gene,coverage_score,specifity_score,coverage_score*specifity_score)
        if(length(line_to_save)==5){
          write(paste(line_to_save,sep="\t",collapse = " "),out_file,sep = "\n",append = T)
        }
      }}
}









                ### ANALYZE DATA ###

Affy_data<-read.table("/home/vlad/gene_regulatory network project/BLAST/Affymetrix_IDs_new.txt",sep=" ",header=F)
#Illumina_data<-read.table("/home/vlad/gene_regulatory network project/BLAST/ILLUMINA_HumanHT-12 V4.0 expression beadchip_IDs.txt",sep=" ",header=F)

select_mu<-function(Array_data,num){
  gene_names<-as.character(Array_data[,2])
  more_one_genes<-vector()
  diff_vec<-vector()
  w<-vector()
  e<-vector()
  for(gene_name in gene_names){
    sub<-Array_data[which(Array_data[,2]==gene_name),]
    if(dim(sub)[1]==2){
      more_one_genes<-append(more_one_genes,gene_name)
      diff_vec<-append(diff_vec,sub[1,num]-sub[2,num])
      w<-append(w, as.character(sub[1,1]))
      e<-append(e, as.character(sub[2,1]))
    }
  }
  s<-data.frame(cbind(more_one_genes,diff_vec,w,e))
  return(s)
}
Affy_diff<-select_mu(Array_data = Affy_data,3)
#Affy_diff<-Affy_diff[which(Affy_diff$diff_vec!=0),]
hist(as.numeric(as.character(Affy_diff$diff_vec)),breaks = 20,xlab="diff",main = "difference")
x_at<-Affy_data[which(grepl("x_at",Affy_data$V1)==TRUE),]
a_at<-Affy_data[which(grepl("a_at",Affy_data$V1)==TRUE),]
at<-Affy_data[which(grepl("[0-9]_at",Affy_data$V1)==TRUE),]
s_at<-Affy_data[which(grepl("s_at",Affy_data$V1)==TRUE),]


# Plot cumulative distrubutions

plot_cum_dist<-function(e_at,num,name1,add){
  name="score"
  if(num==3){
    coloor=colors()[436]
  }
  if(num==4){
    coloor=colors()[369]
  }
  if(num==5){
    coloor=colors()[371]
  }
  print(summary(fivenum(e_at[,num])))
  sorted<-sort(as.numeric(as.character(e_at[,num])),decreasing = F)
  #sorted<-as.numeric(as.character(e_at[,num]))
  n<-length(sorted)
  if(add==TRUE){
    lines(sorted, (1:n)/n, type = 's', ylim = c(0, 1),xlab="Sample Quantiles",ylab=name,main=name1,col = coloor,lwd = 3)
  }
  if(add==FALSE){
  plot(sorted, (1:n)/n, type = 's', ylim = c(0, 1),xlab="Sample Quantiles",ylab=name,main=name1,col=coloor,lwd = 3)
  }
  legend(0.65,0.4,c("cov","spec","cov*spec"),col=c(colors()[436],colors()[369],colors()[371]),lty=c(1,1,1),lwd=c(2.5,2.5))
  }


plot_cum_dist(x_at,3,"x_at",FALSE)
plot_cum_dist(x_at,4,"x_at",TRUE)
plot_cum_dist(x_at,5,"x_at",TRUE)

plot_cum_dist(a_at,3,"a_at",FALSE)
plot_cum_dist(a_at,4,"a_at",TRUE)
plot_cum_dist(a_at,5,"a_at",TRUE)

plot_cum_dist(at,3,"at",FALSE)
plot_cum_dist(at,4,"at",TRUE)
plot_cum_dist(at,5,"at",TRUE)

plot_cum_dist(s_at,3,"s_at",FALSE)
plot_cum_dist(s_at,4,"s_at",TRUE)
plot_cum_dist(s_at,5,"s_at",TRUE)

plot_cum_dist(Affy_data,3,"all",FALSE)
plot_cum_dist(Affy_data,4,"all",TRUE)
plot_cum_dist(Affy_data,5,"all",TRUE)


plot.new()
plot(ecdf(x_at$V4))


hist(x_at$V3,breaks=10,xlab = "cov score",col="blue",ylim = c(1,30000),add=F,main="Coverage")
#lines(density(x_at$V5,adjust = 0.005,bw = 0.005),xlab = "spec*cov score")
hist(a_at$V3,breaks=10,xlab = "spec*cov score",col="green",fill=NULL,add=T)
hist(at$V3,breaks=10,xlab = "spec*cov score",border="red",fill=NULL,add=T)
hist(s_at$V3,breaks=10,xlab = "spec*cov score",border="black",fill=NULL,add=T)
legend(0.1,10000,c("x_at","a_at","at","s_at"),col=c("blue","green","red","black"), lty = 1, bg = "gray90")




expr_data<-allGPL570[which(rownames(allGPL570) %in% c(as.character(Affy_diff$w),as.character(Affy_diff$e))),]
r_t<-vector()
b<-vector()
for(i in 1:dim(Affy_diff)[1]){
  id1<-as.character(Affy_diff[i,3])
  id2<-as.character(Affy_diff[i,4])
  to_add<-as.numeric(rowMeans(expr_data[which(rownames(expr_data)==id1),1:5]))-as.numeric(rowMeans(expr_data[which(rownames(expr_data)==id2),1:5]))
  if(length(to_add)==1){
    r_t<-append(r_t,to_add)
    b<-append(b,as.numeric(as.character(Affy_diff[i,2])))
}}

a<-(r_t-mean(r_t))/sd(r_t)
b<-(b-mean(b))/sd(b)
plot(a,b)


Illumina_diff<-select_mu(Illumina_data)
hist(as.numeric(as.character(Illumina_diff$diff_vec)),breaks = 20,xlab="diff",main = "difference")















