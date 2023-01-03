# REQUIRED:
# /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1.txt
# /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1_Parameters.txt
# /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1_Intervals.txt
# 
# sed -e 's/,/\./g' /uz/data/Admin/cme_genome_raw/PGD_Exported/PGD043_C1.txt > /uz/data/Admin/cme_genome_raw/PGD043_C1.adj
# /cm/shared/apps/R/3.1.0/bin/R

#--------------------------------------------------------
#-------      Global variables         ------------
rm(list=ls(all=T))
library(limma)
library(signal)
library(plotrix)
library(MASS)

args <- commandArgs(TRUE)
logs <- args[1]
script <- args[2]
analyses <- args[3]
gtype_window <- args[4]
Family <- args[5]
SibPattern <- args[6]

siCHILD_DIR <- paste(script,"/analyses/siChild/",sep="")
FixParam <- read.table(paste0(siCHILD_DIR,"Fixed_Parameters.tsv"),header=T,stringsAsFactors=F,sep="\t")
for(p in FixParam[,1]){ eval(parse(text=paste0(p,"=", FixParam[FixParam[,1]==p,2])))}

parametersFile <- "/storage/projects/ngs/analyses/ResearchProjecten/PGT/TEST/PGD018/PGD018_Parameters.txt"
IntervalFile   <- "/storage/projects/ngs/analyses/ResearchProjecten/PGT/TEST/PGD018/PGD018_Intervals.txt"

SibPattern="PGD018"

outPath <- paste(analyses,Family,"siChild",SibPattern,sep="/")
dir.create(outPath)

if (!grepl('/$',outPath)) {outPath <- paste(outPath,"/",sep="")}

siCHILD_DIR_functions <- paste(siCHILD_DIR,"functions/",sep="")
print("Loading sources...")
#ldh: change paths DONE
load(paste(siCHILD_DIR,"Ideogram_hg19.rda",sep=""))
load(paste(siCHILD_DIR,"REF_24h_QC_illuminaCytoSNP12.rda",sep=""))

siCHILD_functions <- list.files(pattern="[.]R$", path=siCHILD_DIR_functions);
sapply(paste0(siCHILD_DIR_functions,siCHILD_functions), FUN=source)

#Default parameter settings
options(scipen=999)
Func = "mean"
Chroms <- c(1:22,"X")

print("Loading parameters and intervals...")
dataFile       <- "/storage/projects/ngs/analyses/ResearchProjecten/PGT/TEST/PGD018/PGD018.txt"

Params <- read.table(parametersFile,sep="\t",header=T,stringsAsFactors=F)
Parent = Params[Params$Param=="Parent","Value"]
Parent1 = paste(Parent,"_",Family,sep="")

for(i in 1:nrow(Params)) { if(Params[i,"Param"]=="Parent" | Params[i,"Param"]=="Seed"  | Params[i,"Param"]=="GC_File") {eval(parse(text=paste(Params[i,"Param"],"='",Params[i,"Value"],"'",sep="")))} 
			   else {eval(parse(text=paste(Params[i,"Param"],"=",Params[i,"Value"])))}
			 }
if(ExcInt==1){Int <- read.table(IntervalFile,sep="\t",header=T,stringsAsFactors=F)}
Int <- Int[complete.cases(Int),]

print(paste("Reading data file:",dataFile))
data <- read.table(dataFile,sep="\t",header=T,stringsAsFactors=F)
if( grepl("chr",unique(data[,"Chr"])[1]) ) { data$Chr <- sapply(data[,"Chr"],function(w) { unlist(strsplit(w,"chr"))[2]}) }

for(c in c(1,2,grep("GType",colnames(data)))){data[,c]<-as.character(data[,c])}
for(c in c(grep("Log.R.Ratio",colnames(data)))){data[,c]<-as.numeric(as.character(data[,c]))}
for(c in c(grep("B.Allele.Freq",colnames(data)))){data[,c]<-as.numeric(as.character(data[,c]))}

dataraw <- na.omit(data)

dataraw <- dataraw[grep("cnvi",as.character(dataraw$Name),invert=TRUE),]
dataraw <- dataraw[dataraw[,"Position"]!=0,]
dataraw <- dataraw[dataraw[,"Chr"]!="Y",]
dataraw <- dataraw[dataraw[,"Chr"]!="XY",]
dataraw[,"Chr"]<-gsub("chr","",dataraw[,"Chr"])

GC_File <- "/storage/projects/ngs/analyses/ResearchProjecten/PGT/script/siCHILD/SeqStat_CytoSNP12_Window10000.txt"
SibPattern="E*_Bl*"
GC <- read.table(GC_File,header=F,sep="\t",stringsAsFactors=F)
rownames(dataraw) <- as.character(dataraw[,1])
rownames(GC) <- as.character(GC[,4])
rowsTot <- intersect(rownames(dataraw),rownames(GC))
GC <- GC[rowsTot,]
dataraw <- dataraw[rowsTot,]

logRsRaw<- data.frame(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".Log.R.Ratio",colnames(dataraw))],stringsAsFactors=FALSE)
colnames(logRsRaw)[-c(1:3)] <- paste(gsub(".Log.R.Ratio","",colnames(logRsRaw)[-c(1:3)]),"_",Family,sep="")	
write.table(logRsRaw,paste(outPath,Family,"_logRsRaw.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

BAFs <- data.frame(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".B.Allele.Freq",colnames(dataraw))],stringsAsFactors=FALSE)
colnames(BAFs)[-c(1:3)] <- paste(gsub(".B.Allele.Freq","",colnames(BAFs)[-c(1:3)]),"_",Family,sep="")
write.table(BAFs,paste(outPath,Family,"_BAFs.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

Gtypes<- data.frame(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".GType",colnames(dataraw))],stringsAsFactors=FALSE)
colnames(Gtypes) <- colnames(logRsRaw)
write.table(Gtypes,paste(outPath,Family,"_0.75.gtp",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

ChrPos <- Gtypes[,c("Chr","Position")]
QC <- qcgtype(Gtypes,ChrPos,Family,SibPattern,outPath) 
QCbyParents <- qcbyparents(Gtypes,SibPattern)
SnpSpecArts <- artscsnpspar(Gtypes,SibPattern,outPath)
SnpSpecArts[,"Position"] <- as.numeric(SnpSpecArts[,"Position"])

dataPo <- po2(Gtypes,Family,SibPattern,outPath)

ParScore <- patscore2(dataPo,QC,Chroms,Gtypes,Family) 

names(ParScore)<-gsub("E00_Bl000_","",names(ParScore))
save(ParScore,file=paste(outPath,"ParScore_",Family,".rda",sep=""))

window=Win/2
logRs <- meanwindow(logRsRaw,GC,window,Func,Family,ParScore,outPath,SibPattern)
write.table(logRs,paste(outPath,Family,"_logRsAvgWindow.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

gamma=gammaBAF

if(Seed=="Grandparents"){
 Haps <- htypingOpt1(Gtypes,dataPo,ParScore,Parent1,SibPattern,Family,outPath)
 print("-------------------------------------")
 print("Option one haplotyping was applied...")
 print("-------------------------------------")
 PhBAF <- phasebaf2(Gtypes,dataraw,Family,Parent1,gamma,outPath)
}else{
 Haps <- htypingOpt2(Gtypes,Window,Int,dataPo,ParScore,SibPattern,Family,outPath)
 print("-------------------------------------")
 print("Option two haplotyping was applied...")
 print("-------------------------------------")
 PhBAF <- phasebaf4(Gtypes,dataraw,Family,gamma,ParScore,outPath)
}

dataHap <- Haps[["dataHap"]]
dataHapRaw <- Haps[["dataHapRaw"]] 
write.table(dataHap,paste(outPath,Family,"_Itp.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(dataHapRaw,paste(outPath,Family,"_Raw.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

AvgLogRs <- avgint(logRs,Int,Family,outPath)
SegLogRs <- callpcf(logRs,gammaSC,gammaMC,plateau,Family,outPath)

library(MASS)
Intp <- intphappropser3(outPath,Family,ParScore,Haps,SibPattern)

disthap(Family,Int,outPath,SibPattern)
 
save.image(paste(siCHILD_DIR,Family,".rda",sep=""))

outPathChr <- paste(outPath,"ChrSpecPlots/",sep="")
if (!file.exists(outPathChr)){
 dir.create(outPathChr)
}

chrspecbafplot(dataHap,dataHapRaw,dataPo,BAFs,logRs,AvgLogRs,PhBAF,ChrsLength,ideogram,Family,outPathChr,Chroms,Int)
 
outPathGenome <- paste(outPath,"GenomewidePlots/",sep="")
if (!file.exists(outPathGenome)){
 dir.create(outPathGenome)
}

genomebafplot(BAFs,logRs,AvgLogRs ,PhBAF,ChrsLength,ideogram,Family,outPathGenome)

