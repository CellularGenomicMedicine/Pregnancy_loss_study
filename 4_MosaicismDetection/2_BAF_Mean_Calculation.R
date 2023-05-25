###########################################################################################################################
# Author: Rick Essers (Adapted from script by Masoud Zamani Esteki)
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To calculate mean BAF values for genomic coordinates of interest. 
#                 These mean BAF values can be used to calculate mosaicism degree (in percentage) of detected aberrations. 
#                 Indicate whether 4 or 6 BAF bands are present in the genomic region of interest. 

# input: Genomic coordinates of interest and number of visible BAF bands are required. 
#        P1, P2, M1, M2, P1Seg, P2Seg, M1Seg, M2Seg, LogRSeg .txt files for a specific pregnancy loss family. 
#        These files are produced by haplarithmisis. 

# output: For both fetal DNA samples (Extraembroynic mesoderm and chorionic villi), one file is produced containing main BAF values
#         for the genomic location of interest. 

###########################################################################################################################

rm(list=ls(all=T))

Family = "PL0001.adj" 

###Insert the genomic coordinates of the aberration. Coordinates are attained by 1_AutomatedSegmentation_PL.R script. 
Int1 <- "Chr7:58814-5001488"

###Insert the number of BAF bands that are visible for the genomic region of interest (valid options are 4 or 6)
NumberofBands <- "4"

dataPath      <- paste0("/Projects/PregnancyLoss/Data/Output/",Family,"/")
outPath       <- paste0("/Projects/PregnancyLoss/Data/Output/",Family,"/")

outPathOutput <- paste(outPath,"/OutputAutSeg/",sep="")
if (!file.exists(outPathOutput)){
  dir.create(outPathOutput)
}

Int <- rbind(c(0,0,0,"Pat"),c(0,0,0,"Mat"))

###Make sure to have the correct gamma value, for the pregnancy loss study gamma value was set to 14.
logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc14_gammaMc14_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

logRs <- read.table(paste(dataPath,Family,"_logRsRaw.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
logRsWind <- read.table(paste(dataPath,Family,"_logRsAvgWindow.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
Gtypes<- read.table(paste(dataPath,Family,"_0.75.gtp",sep=""),header=T,sep="\t",stringsAsFactors=F)
BAFs <- read.table(paste(dataPath,Family,"_BAFs.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P1Seg <- read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
Poo <- read.table(paste(dataPath,Family,".poo",sep=""),header=T,sep="\t",stringsAsFactors=F)

Int <- do.call("rbind",strsplit(Int1,":"))
Chr = gsub("Chr","",Int[1])
Start <- as.numeric(do.call("rbind",strsplit(Int[2],"-"))[1])
Stop <- as.numeric(do.call("rbind",strsplit(Int[2],"-"))[2])


logRsSegInt <- logRsSeg[logRsSeg$Chr==Chr & logRsSeg$Position>=Start & logRsSeg$Position<=Stop,]
logRsInt <- logRs[logRs$Chr==Chr & logRs$Position>=Start & logRs$Position<=Stop,]
logRsWindInt <- logRsWind[logRsWind$Chr==Chr & logRsWind$Position>=Start & logRsWind$Position<=Stop,]
BAFsdInt <- BAFs[BAFs$Chr==Chr & BAFs$Position>=Start & BAFs$Position<=Stop,]
P1Int <- P1[P1$Chr==Chr & P1$Position>=Start & P1$Position<=Stop,]
P2Int <- P2[P2$Chr==Chr & P2$Position>=Start & P2$Position<=Stop,]
M1Int <- M1[M1$Chr==Chr & M1$Position>=Start & M1$Position<=Stop,]
M2Int <- M2[M2$Chr==Chr & M2$Position>=Start & M2$Position<=Stop,]
GtypeInt <- Gtypes[Gtypes$Chr==Chr & Gtypes$Position>=Start & Gtypes$Position<=Stop,]
P1SegInt <- P1Seg[P1Seg$Chr==Chr & P1Seg$Position>=Start & P1Seg$Position<=Stop,]
P2SegInt <- P2Seg[P2Seg$Chr==Chr & P2Seg$Position>=Start & P2Seg$Position<=Stop,]
M1SegInt <- M1Seg[M1Seg$Chr==Chr & M1Seg$Position>=Start & M1Seg$Position<=Stop,]
M2SegInt <- M2Seg[M2Seg$Chr==Chr & M2Seg$Position>=Start & M2Seg$Position<=Stop,]
PooInt <- Poo[Poo$Chr==Chr & Poo$Position>=Start & Poo$Position<=Stop,]


if (NumberofBands == "4"){

#Extract mean values of BAF bands for given coordinates. 
EM_band1    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.90 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<1),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band2    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.50 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.90),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band3    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.50 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.10),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band4    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.10 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0),grep("E*_Bl*",colnames(BAFsdInt))])

CV_band1    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.90 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<1),grep("Affected",colnames(BAFsdInt))])
CV_band2    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.50 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.90),grep("Affected",colnames(BAFsdInt))])
CV_band3    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.50 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.10),grep("Affected",colnames(BAFsdInt))])
CV_band4    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.10 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0),grep("Affected",colnames(BAFsdInt))])

AllBAF <- rbind(EM_band1,EM_band2,EM_band3,EM_band4,"",CV_band1,CV_band2,CV_band3,CV_band4)

Filename1 = paste0("AllBAF_4bands_",Int1)
write.table(AllBAF,paste0(outPathOutput,Filename1,".txt"),sep="\t",col.names = F,row.names = T)

AllBAF
Int


} else if (NumberofBands == "6"){

EM_band1    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.95 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<1),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band2    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.75 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.95),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band3    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.50 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.75),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band4    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.50 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.25),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band5    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.25 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.05),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band6    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.05 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0),grep("E*_Bl*",colnames(BAFsdInt))]) 

CV_band1    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.95 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<1),grep("Affected",colnames(BAFsdInt))])
CV_band2    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.75 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.95),grep("Affected",colnames(BAFsdInt))])
CV_band3    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.50 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.75),grep("Affected",colnames(BAFsdInt))])
CV_band4    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.50 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.25),grep("Affected",colnames(BAFsdInt))])
CV_band5    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.25 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.05),grep("Affected",colnames(BAFsdInt))])
CV_band6    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.05 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0),grep("Affected",colnames(BAFsdInt))])

AllBAF <- rbind(EM_band1,EM_band2,EM_band3,EM_band4,EM_band5,EM_band6,"",CV_band1,CV_band2,CV_band3,CV_band4,CV_band5,CV_band6)

Filename1 = paste0("AllBAF_6bands_",Int1)
write.table(AllBAF,paste0(outPathOutput,Filename1,".txt"),sep="\t",col.names = F,row.names = T)

AllBAF
Int

}

