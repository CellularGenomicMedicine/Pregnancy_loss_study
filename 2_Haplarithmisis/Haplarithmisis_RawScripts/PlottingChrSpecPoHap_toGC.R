siCHILD_DIR <- "/home/zamani/PGT/script/siCHILD/"

args <- commandArgs(TRUE)
Family <- args[1]
parametersFile <- args[2]
intervalFile <- args[3]
outPath <- args[4]
PGD_EXPORTED_DIR <- paste("/home/zamani/PGT/TEST/",Family,"/",sep="")
if (!grepl('/$',outPath)) {outPath <- paste(outPath,"/",sep="")}

source(paste(siCHILD_DIR,"chrspecbafplot.R",sep=""))
#chrspecbafplot.R
load(paste(siCHILD_DIR,"Ideogram_hg19.rda",sep=""))
Chroms = c(1:22,"X")

Params <- read.table(parametersFile,sep="\t",header=T,stringsAsFactors=F)

for(i in 1:nrow(Params)) { if(Params[i,"Param"]=="Parent" | Params[i,"Param"]=="Seed"  | Params[i,"Param"]=="GC_File") {eval(parse(text=paste(Params[i,"Param"],"='",Params[i,"Value"],"'",sep="")))}
                           else {eval(parse(text=paste(Params[i,"Param"],"=",Params[i,"Value"])))}
                         }
#if(ExcInt==1){Int <- read.table(paste(PGD_EXPORTED_DIR,Family,"_Intervals.txt",sep=""),sep="\t",header=T,stringsAsFactors=F)}
if(ExcInt==1){Int <- read.table(intervalFile,sep="\t",header=T,stringsAsFactors=F)}

dataPath <- outPath
outPath <- paste(outPath,"ChrSpecPlots/",sep="")

if (!file.exists(outPath)){
  dir.create(outPath)
}

P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P1Seg <-  read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
BAFs <- read.table(paste(dataPath,Family,"_BAFs.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRs <- read.table(paste(dataPath,Family,"_logRsAvgWindow.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc300_gammaMc50_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataPo <- read.table(paste(dataPath,Family,".poo",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataHap <- read.table(paste(dataPath,Family,"_Itp2.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)
dataHapRaw <- read.table(paste(dataPath,Family,"_Raw.hap",sep=""),header=T,sep="\t",stringsAsFactors=F)

chrspecbafplot(dataHap,dataHapRaw,dataPo,BAFs,logRs,logRsSeg,P1,P1Seg,P2,P2Seg,M1,M1Seg,M2,M2Seg,ChrsLengths,ideogram,Family,outPath,Chroms,Int,siCHILD_DIR)
