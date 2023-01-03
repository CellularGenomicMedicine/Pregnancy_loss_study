rm(list=ls(all=T))

Family = "PL2726.adj" 
Indication = "Trisomy_UPD"
Group = "RPL"
Int1 <- "Chr1:58814-249222527"

dataPathInit <- "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/FinalPlots/AllFams/Abnormal/"
outpath <- paste0("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/FinalPlots/AllFams/Abnormal/",Group,"/",Indication,"/",Family,"/OutputAutSeg/")

Int <- rbind(c(0,0,0,"Pat"),c(0,0,0,"Mat"))

dataPath = paste(dataPathInit,"/",Group,"/",Indication,"/",Family,"/",sep="")
dataPath <- gsub("Input","Output",dataPath)

###Check which gamma was used
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

#BAF_MosaicismDetection 
#EM_UpperMedBAF <- median(BAFsdInt[(BAFsdInt[,4]>0.51 & BAFsdInt[,4]<0.95),4])
#EM_LowerMedBAF <- median(BAFsdInt[(BAFsdInt[,4]<0.49 & BAFsdInt[,4]>0.05),4])
EM_band1    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.90 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<1),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band2    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.50 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.90),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band3    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.50 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0.10),grep("E*_Bl*",colnames(BAFsdInt))])
EM_band4    <- mean(BAFsdInt[(BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]<0.10 & BAFsdInt[,grep("E*_Bl*",colnames(BAFsdInt))]>0),grep("E*_Bl*",colnames(BAFsdInt))])

CV_band1    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.90 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<1),grep("Affected",colnames(BAFsdInt))])
CV_band2    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.50 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.90),grep("Affected",colnames(BAFsdInt))])
CV_band3    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.50 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0.10),grep("Affected",colnames(BAFsdInt))])
CV_band4    <- mean(BAFsdInt[(BAFsdInt[,grep("Affected",colnames(BAFsdInt))]<0.10 & BAFsdInt[,grep("Affected",colnames(BAFsdInt))]>0),grep("Affected",colnames(BAFsdInt))])

AllBAF <- rbind(EM_band1,EM_band2,EM_band3,EM_band4,"",CV_band1,CV_band2,CV_band3,CV_band4)
#AllBAF <- rbind(EM_UpperMedBAF,EM_LowerMedBAF,EM_UpperMeanBAF,EM_LowerMeanBAF,CV_UpperMedBAF,CV_LowerMedBAF,CV_UpperMeanBAF,CV_LowerMeanBAF)

Filename1 = paste0("AllBAF_4bands_",Int1)
write.table(AllBAF,paste0(outpath,Filename1,".txt"),sep="\t",col.names = F,row.names = T)

AllBAF
Int


#P1Seg,P2Seg,M1Seg,M2Seg median and average
#EM_medianP1Seg <- median(P1SegInt[,5])
#EM_medianP1Seg
#EM_medianP2Seg <- median(P2SegInt[,5])
#EM_medianP2Seg
#EM_PTotal <- EM_medianP2Seg - EM_medianP1Seg
#EM_PTotal

#EM_medianM1Seg <- median(M1SegInt[,5])
#EM_medianM1Seg
#EM_medianM2Seg <- median(M2SegInt[,5])
#EM_medianM2Seg
#EM_MTotal <- EM_medianM2Seg - EM_medianM1Seg
#EM_MTotal

#EMtotal <- EM_PTotal + EM_MTotal
#EMtotal

#CV_medianP1Seg <- median(P1SegInt[,4])
#CV_medianP2Seg <- median(P2SegInt[,4])
#CV_medianM1Seg <- median(M1SegInt[,4])
#CV_medianM2Seg <- median(M2SegInt[,4])

#ParentalAllSeg <- rbind(EM_medianP1Seg,EM_medianP2Seg,EM_medianM1Seg,EM_medianM2Seg,CV_medianP1Seg,CV_medianP2Seg,CV_medianM1Seg,CV_medianM2Seg)

#Filename1 = paste0("ParentalAllSeg_",Int1,"_",Family)
#write.table(ParentalAllSeg,paste0(outpath,Filename1,".txt"),sep="\t",col.names = F,row.names = T)

#ParentalAllSeg
#Int

#LogR_valueDetection
#MedSeg <-  apply(logRsSegInt[,-c(1:2)],2,median)
#MeanSeg <-  apply(logRsSegInt[,-c(1:2)],2,mean)

#MedRaw <-  apply(logRsInt[,-c(1:2)],2,median)
#MeanRaw <-  apply(logRsInt[,-c(1:2)],2,mean)

#MedWind <-  apply(logRsWindInt[,-c(1:2)],2,median)
#MeanWind <-  apply(logRsWindInt[,-c(1:2)],2,mean)

#AllSeg <- rbind(MedSeg,MeanSeg,MedRaw, MeanRaw, MedWind, MeanWind)

#Filename2 = paste0("AllSeg_",Int1)

#write.table(AllSeg,paste0(outpath,Filename2,".txt"),sep = "\t", col.names = T,row.names = F)

#AllSeg
#Int
