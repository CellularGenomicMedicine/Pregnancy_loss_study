
rm(list=ls(all=T))
library(limma)
library(signal)
library(plotrix)
library(MASS)
library("data.table")
library("stringr")

args <- commandArgs(TRUE)
Family <- args[1]
SibPattern   <- args[2]
logs         <- args[7]
#Seed       <- args[3]
Parent     <- args[4]
gtype_window <- args[5]
#Family       <- args[6]
#SibPattern   <- args[2]

gtype_window <- 10000
#Family       <- "PL2631"
script       <- "/ifs/data/research/CGM/scripts/Haplarithmisis_RawScripts/"
analyses     <- "/ifs/data/research/CGM/scripts/Haplarithmisis_RawScripts/"
logs         <- "/ifs/data/research/CGM/scripts/Haplarithmisis_RawScripts/logs/"
SibPattern   <- "E*_Bl*"
Seed         <- "Affected"

outPathInit <- "/ifs/data/research/CGM/Projects/Miscarr/1_AllFams/Output"
siCHILD_DIR <- paste(script,"",sep="")
FixParam <- read.table(paste0(script,"Fixed_Parameters.tsv"),header=T,stringsAsFactors=F,sep="\t")
for(p in FixParam[,1]){ eval(parse(text=paste0(p,"=", FixParam[FixParam[,1]==p,2])))}
parametersFile <- paste0("/ifs/data/research/CGM/scripts/Haplarithmisis_RawScripts/Codes/","Generic_Parameters.txt")
IntervalFile <- paste0("/ifs/data/research/CGM/scripts/Haplarithmisis_RawScripts/Codes/","Generic_Intervals.txt")

if (file.exists(paste(outPathInit,"/",Family,sep=""))){
  outPath <- paste(outPathInit,"/" ,Family,"/",sep="")
} else {
  dir.create(file.path(outPathInit, Family))
  outPath <- file.path(outPathInit, "/", Family, "/")}


if (!grepl('/$',outPath)) {outPath <- paste(outPath,"/",sep="")}

siCHILD_DIR_functions <- paste(siCHILD_DIR,"functions/",sep="")
print("Loading sources...")

load(paste(siCHILD_DIR,"Ideogram_hg19.rda",sep=""))
#load(paste(siCHILD_DIR,"REF_24h_QC_illuminaCytoSNP12.rda",sep=""))

siCHILD_functions <- list.files(pattern="[.]R$", path=siCHILD_DIR_functions);
sapply(paste0(siCHILD_DIR_functions,siCHILD_functions), FUN=source)

#Default parameter settings
options(scipen=999)
Func = "mean"
Chroms <- c(1:22,"X")

print("Loading parameters and intervals...")
dataFile <- paste(analyses,Family,"BAF",paste(Family,"txt",sep="."),sep="/")

Params <- read.table(parametersFile,sep="\t",header=T,stringsAsFactors=F)
Parent1 = paste(Parent,"_",Family,sep="")

for(i in 1:nrow(Params)) { if(Params[i,"Param"]=="Parent" | Params[i,"Param"]=="Seed"  | Params[i,"Param"]=="GC_File") {eval(parse(text=paste(Params[i,"Param"],"='",Params[i,"Value"],"'",sep="")))} 
  else {eval(parse(text=paste(Params[i,"Param"],"=",Params[i,"Value"])))}
}
if(ExcInt==1){Int <- read.table(IntervalFile,sep="\t",header=T,stringsAsFactors=F)}
Int <- Int[complete.cases(Int),]

GC_File <- paste0(script,GC_File)

dataPath <- "/ifs/data/research/CGM/Projects/Miscarr/1_AllFams/Families/MedicalAbortions/"
#dataFile <- paste0(dataPath,Family)
dataFile <- paste0(dataPath,Family,".adj")
print(paste("Reading data file:",dataFile))
data <- fread(dataFile,sep="\t",header=T,stringsAsFactors=F)
data <- as.data.frame(data,stringsAsFactors=FALSE)
if(names(data)[1]!="Name") {names(data)[1] <- "Name"}
for(c in c(1,2,grep("GType",colnames(data)))){data[,c]<-as.character(data[,c])}
for(c in c(grep("Log.R.Ratio",colnames(data)))){data[,c]<-as.numeric(as.character(data[,c]))}
for(c in c(grep("B.Allele.Freq",colnames(data)))){data[,c]<-as.numeric(as.character(data[,c]))}

dataraw <- na.omit(data)
dataraw <- dataraw[order(dataraw[,"Chr"],dataraw[,"Position"]),]
for(chr in c(1:22,"X")){datarawchr <- dataraw[dataraw[,"Chr"]==chr,]; datarawchrSort <- datarawchr[order(datarawchr[,"Position"]),];if(chr=="1"){datarawgenome <-datarawchrSort}else{datarawgenome <-rbind(datarawgenome, datarawchrSort)}}
dataraw <- datarawgenome
dataraw <- dataraw[dataraw[,"Position"]!=0,]
if(grepl("chr",dataraw[1,"Chr"])==T) {dataraw[,"Chr"]<-gsub("chr","",dataraw[,"Chr"])}

dataraw <- dataraw[dataraw[,"Chr"]!="Y",]
dataraw <- dataraw[dataraw[,"Chr"]!="XY",]
dataraw <- dataraw[dataraw[,"Chr"]!="MT",]

GC <- fread(GC_File,sep="\t",header=T,colClasses=c("character","integer","integer","character","numeric","numeric","integer","integer","integer","integer","integer","integer","integer"),stringsAsFactors=F)
GC <- as.data.frame(GC,stringsAsFactors=FALSE)
rownames(dataraw) <- as.character(dataraw[,"Name"])
rownames(GC) <- as.character(GC[,4])
for(chr in paste0("chr",c(1:22,"X"))){GCchr <- GC[GC[,1]==chr,]; GCchrSort <- GCchr[order(GCchr[,2]),];if(chr=="chr1"){GCgenome <-GCchrSort}else{GCgenome <-rbind(GCgenome,GCchrSort)}}
GC <- GCgenome
#GC[,1]<-gsub("chr","",GC[,1])
rowsTot <- intersect(rownames(dataraw),rownames(GC))
GC <- GC[rowsTot,]
dataraw <- dataraw[rowsTot,]

logRsRaw<- data.frame(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".Log.R.Ratio",colnames(dataraw))],stringsAsFactors=FALSE)
colnames(logRsRaw)[-c(1:3)] <- paste(gsub(".Log.R.Ratio","",colnames(logRsRaw)[-c(1:3)]),"_",Family,sep="")    
write.table(logRsRaw,paste(outPath,Family,"_logRsRaw.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

BAFs <- data.frame(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".B.Allele.Freq",colnames(dataraw))],stringsAsFactors=FALSE)
colnames(BAFs)[-c(1:3)] <- paste(gsub(".B.Allele.Freq","",colnames(BAFs)[-c(1:3)]),"_",Family,sep="")
write.table(BAFs,paste(outPath,Family,"_BAFs.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

Gtypes <- data.frame(dataraw[,c("Name","Chr","Position")],dataraw[,grep(".GType",colnames(dataraw))],stringsAsFactors=FALSE)
colnames(Gtypes) <- colnames(logRsRaw)
write.table(Gtypes,paste(outPath,Family,"_0.75.gtp",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

ChrPos <- Gtypes[,c("Chr","Position")]
QC <- qcgtype(Gtypes,ChrPos,Family,SibPattern,outPath) 
QCbyParents <- qcbyparents(Gtypes,SibPattern)
SnpSpecArts <- artscsnpspar(Gtypes,SibPattern,outPath)
SnpSpecArts[,"Position"] <- as.numeric(SnpSpecArts[,"Position"])

dataPo <- po2.1(Gtypes,Family,SibPattern,outPath)

ParScore <- patscore2(dataPo,QC,Chroms,Gtypes,Family) 

names(ParScore)<-gsub("E00_Bl000_","",names(ParScore))
save(ParScore,file=paste(outPath,"ParScore_",Family,".rda",sep=""))

window=Win/2
logRs <- meanwindow2(logRsRaw,GC,window,Func,Family,ParScore,outPath,SibPattern)
write.table(logRs,paste(outPath,Family,"_logRsAvgWindow.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

gamma=gammaBAF

if(Seed=="Grandparents"){
  Haps <- htypingOpt1(Gtypes,dataPo,ParScore,Parent1,SibPattern,Family,outPath)
  print("-------------------------------------")
  print("Option one haplotyping was applied...")
  print("-------------------------------------")
  PhBAF <- phasebaf2(Gtypes,dataraw,Family,Parent1,gamma,outPath,SibPattern)
}else{
  Haps <- htypingOpt2(Gtypes,Window,Int,dataPo,ParScore,SibPattern,Family,outPath)
  print("-------------------------------------")
  print("Option two haplotyping was applied...")
  print("-------------------------------------")
  PhBAF <- phasebaf4_mod(Gtypes,dataraw,Family,gamma,ParScore,outPath,SibPattern)
}


write.table(Haps[["dataHap"]],paste(outPath,Family,"_Itp.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
write.table(Haps[["dataHapRaw"]],paste(outPath,Family,"_Raw.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

AvgLogRs <- avgint(logRs,Int,Family,outPath)
SegLogRs <- callpcf(logRs,gammaSC,gammaMC,plateau,Family,outPath)

library(MASS)
Intp <- intphappropser3(outPath,Family,ParScore,Haps,SibPattern)

#disthap(Family,Int,outPath,SibPattern)
#///
dataHap    <- Haps[["dataHap"]]
dataHapRaw <- Haps[["dataHapRaw"]]
#///

outPathGenome <- paste(outPath,"GenomewidePlots/",sep="")
if (!file.exists(outPathGenome)){
  dir.create(outPathGenome)
}
genomebafplot_miscarrFinalRound3(BAFs,logRs,SegLogRs ,PhBAF,ChrsLength,ideogram,Family,outPathGenome)

outPathChr <- paste(outPath,"ChrSpecPlots/",sep="")
if (!file.exists(outPathChr)){
  dir.create(outPathChr)
}
#Chroms <- "5"
chrspecbafplot_miscarrFinal(Haps[["dataHap"]],Haps[["dataHapRaw"]],dataPo,BAFs,logRs,SegLogRs,PhBAF,ChrsLength,ideogram,Family,outPathChr,Chroms,Int)

