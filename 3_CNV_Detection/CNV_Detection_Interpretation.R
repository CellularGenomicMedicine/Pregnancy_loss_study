###########################################################################################################################
# Author: Rick Essers (Adapted from script by Masoud Zamani Esteki)
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To extract genomic coordinates and aberrant segments based on LogR (copy-number) segments produced by haplarithmisis. 
#                 Segments of logR are extracted by grouping logR segments together that have the same concurrent value. 
#                 The logR segmentation algorithm of haplarithmisis can be adjusted by changes the gamma values for LogR, in the pregnancy loss study gamma of 14 is used. 

# input: P1, P2, M1, M2, P1Seg, P2Seg, M1Seg, M2Seg, LogRSeg .txt files for a specific pregnancy loss family. 
#        These files are produced by haplarithmisis. 

# output: Per DNA sample, one file containing all segments of logR, with Chr, segment length, LogR values, start and stop positions of the segment, 
#         raw parental haplotype values as well as segmented haplotype values, and an interpretation of the raw data. 

###########################################################################################################################

rm(list=ls(all=T))
require(pastecs)

Family      <- "PL0001.adj"

dataPathInit  <- paste0("/Projects/PregnancyLoss/Data/Output","/",Family,"/")
outPath       <- paste0("/Projects/PregnancyLoss/Data/Output/")


autintp <- function(MedLogR,dMat,dPat,Params){
	
	minLogR = -0.3
	maxLogR = 0.15
	
	if((MedLogR < maxLogR & MedLogR > minLogR) & dMat >= Params["Nrms","Par1min"] & dMat <= Params["Nrms","Par1max"] & dPat >= Params["Nrms","Par2min"] & dPat <= Params["Nrms","Par2max"] ){
				autintp = "normal"
	}else if(MedLogR > maxLogR & dMat >= Params["Nrms","Par1min"] & dMat <= Params["Nrms","Par1max"] & dPat >= Params["Nrms","Par2min"] & dPat <= Params["Nrms","Par2max"] ){
				autintp = "TT"
		
	}else if((MedLogR < maxLogR & MedLogR > minLogR) & dMat >= Params["Upds","Par1min"] & dMat <= Params["Upds","Par1max"] & dPat >= Params["Upds","Par2min"] & dPat <= Params["Upds","Par2max"]){
				autintp = "UM"

	}else if((MedLogR < maxLogR & MedLogR > minLogR) & dMat >= Params["Upds","Par2min"] & dMat <= Params["Upds","Par2max"] & dMat >= Params["Upds","Par1min"] & dPat <= Params["Upds","Par1max"]){
				autintp = "UP"
				
	}else if( MedLogR < minLogR & dMat >= Params["Mons","Par2min"] & dMat <= Params["Mons","Par2max"] & dMat >= Params["Mons","Par1min"] & dPat <= Params["Mons","Par1max"]){
				autintp = "MP"
				
	}else if( MedLogR < minLogR & dMat >= Params["Mons","Par1min"] & dMat <= Params["Mons","Par1max"] & dPat >= Params["Mons","Par2min"] & dPat <= Params["Mons","Par2max"]){
				autintp = "MM"
	
	}else if( MedLogR < minLogR & dMat >= Params["Nuls","Par1min"] & dMat <= Params["Nuls","Par1max"] & dPat >= Params["Nuls","Par2min"] & dPat <= Params["Nuls","Par2max"]){
				autintp = "NN"
	
	}else if( MedLogR < minLogR & dMat >= Params["Nuls","Par2min"] & dMat <= Params["Nuls","Par2max"] & dPat >= Params["Nuls","Par1min"] & dPat <= Params["Nuls","Par1max"]){
				autintp = "NN"

	}else if(MedLogR > maxLogR  & dMat >= Params["Tris","Par2min"] & dMat <= Params["Tris","Par2max"] & dMat >= Params["Tris","Par1min"] & dPat <= Params["Tris","Par1max"]){
				autintp = "TP"
				
	}else if( MedLogR > maxLogR & dMat >= Params["Tris","Par1min"] & dMat <= Params["Tris","Par1max"] & dPat >= Params["Tris","Par2min"] & dPat <= Params["Tris","Par2max"]){
				autintp = "TM"
	}else if(MedLogR > maxLogR & dMat >= Params["Upds","Par1min"] & dMat <= Params["Upds","Par1max"] & dPat >= Params["Upds","Par2min"] & dPat <= Params["Upds","Par2max"] ){
				autintp = "TT_m"
				}else{autintp = "ND" }
autintp		
}#end autintp function

source("/Projects/PregnancyLoss/Scripts/Haplarithmisis/CNV_detection/Codes/CNV_Detection_External_Function.R")
Params <- read.table("/Projects/PregnancyLoss/Scripts/Haplarithmisis/CNV_detection/Codes/CNV_Interpretation_Ranges.txt",sep="\t",header=T,stringsAsFactors=F, row.names = 1)

Output_CNV_Det <- paste0(outPath,Family,"/Output_CNV_Det/",sep="")
if (!file.exists(Output_CNV_Det)){
  dir.create(Output_CNV_Det)
}

#Load the P1, P2, M1, M2, P1Seg, P2Seg, M1Seg, M2Seg files for this family created by haplarithmisis	
	dataPath = paste(dataPathInit,Family,"/",sep="")
	P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	P1Seg <-  read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)	
	M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

#Load the LogRseg.txt file for this family created by haplarithmisis	
	logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc14_gammaMc14_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	
	Inds = colnames(logRsSeg)[c(grep("E*_Bl*",colnames(logRsSeg)),grep("Affected*",colnames(logRsSeg)))]

	dMats <- matrix(NA,length(Inds),23)
	colnames(dMats) <- unique(logRsSeg$Chr) 
	dMats[,1] <- Inds 
	dPats <- dMats
	
for(ind in Inds){

	
	for(chr in unique(logRsSeg$Chr)){
		
		P1SegChr <- P1Seg[P1Seg$Chr==chr,c("Name","Chr","Position",ind)]
		P2SegChr <- P2Seg[P2Seg$Chr==chr,c("Name","Chr","Position",ind)]
		M1SegChr <- M1Seg[M1Seg$Chr==chr,c("Name","Chr","Position",ind)]
		M2SegChr <- M2Seg[M2Seg$Chr==chr,c("Name","Chr","Position",ind)]

		RLEseg <- rle(logRsSeg[logRsSeg$Chr==chr,ind])
		RLEseg$cumlengths <- cumsum(RLEseg$lengths)
		RLEseg$dMats <- rep(NA,length(RLEseg$lengths))
		RLEseg$dPats <- rep(NA,length(RLEseg$lengths))
		RLEseg$starts <- rep(NA,length(RLEseg$lengths))
		RLEseg$stops <- rep(NA,length(RLEseg$lengths))
		RLEseg$P1s <- rep(NA,length(RLEseg$lengths))
		RLEseg$P2s <- rep(NA,length(RLEseg$lengths))
		RLEseg$M1s <- rep(NA,length(RLEseg$lengths))
		RLEseg$M2s <- rep(NA,length(RLEseg$lengths))
		RLEseg$Intp <- rep(NA,length(RLEseg$lengths))
				
		for(i in 1:length(RLEseg$values)){
		
		RLEseg$starts[i] <- logRsSeg[RLEseg$cumlengths[i]-(RLEseg$lengths[i]-1),"Position"]
		RLEseg$stops[i] <- logRsSeg[RLEseg$cumlengths[i],"Position"]
		
		P1SegChrC <- P1SegChr[P1SegChr$Position >= logRsSeg[RLEseg$cumlengths[i]-(RLEseg$lengths[i]-1),"Position"]& P1SegChr$Position <= logRsSeg[RLEseg$cumlengths[i],"Position"],]
		
		P2SegChrC <- P2SegChr[P2SegChr$Position >= logRsSeg[RLEseg$cumlengths[i]-(RLEseg$lengths[i]-1),"Position"]& P2SegChr$Position <= logRsSeg[RLEseg$cumlengths[i],"Position"],]

		M1SegChrC <- M1SegChr[M1SegChr$Position >= logRsSeg[RLEseg$cumlengths[i]-(RLEseg$lengths[i]-1),"Position"]& M1SegChr$Position <= logRsSeg[RLEseg$cumlengths[i],"Position"],]
		
		M2SegChrC <- M2SegChr[M2SegChr$Position >= logRsSeg[RLEseg$cumlengths[i]-(RLEseg$lengths[i]-1),"Position"]& M2SegChr$Position <= logRsSeg[RLEseg$cumlengths[i],"Position"],]

		RLEseg$P1s[i] <- nrow(P1SegChrC)
		RLEseg$P2s[i] <- nrow(P2SegChrC)
		RLEseg$M1s[i] <- nrow(M1SegChrC)
		RLEseg$M2s[i] <- nrow(M2SegChrC)
		
		if(nrow(M1SegChrC)!=0 & nrow(P1SegChrC)!=0 & nrow(P2SegChrC)!=0 & nrow(M2SegChrC)!=0){
			
		dPat <- distcompseg(P1SegChrC,P2SegChrC)		
		dMat <- distcompseg(M1SegChrC,M2SegChrC)


	 tryCatch({
		DdMat <- density(na.omit(dMat[,ind]))
		tsMat <- ts(DdMat$y)
		tpMat <- turnpoints(tsMat)
		YPeaksDdMat <- DdMat$y[tpMat$peaks]
		MaxPeakMat <- max(YPeaksDdMat)
		XPeaksDdMat <- DdMat$x[tpMat$peaks]
		XMaxPeakMat <- XPeaksDdMat[YPeaksDdMat==MaxPeakMat]
		RLEseg$dMats[i] <-  XMaxPeakMat
		},error = function(e){print(paste("density couldn't be fitted")); RLEseg$dMats[i]=median(na.omit(dMat[,ind]))})
		
		if(chr != "X"){
	
			tryCatch({
			DdPat <- density(na.omit(dPat[,ind]))
			tsPat <- ts(DdPat$y)
			tpPat <- turnpoints(tsPat)
			YPeaksDdPat <- DdPat$y[tpPat$peaks]
			MaxPeakPat <- max(YPeaksDdPat)
			XPeaksDdPat <- DdPat$x[tpPat$peaks]
			XMaxPeakPat <- XPeaksDdPat[YPeaksDdPat==MaxPeakPat]
			RLEseg$dPats[i]  <-  XMaxPeakPat
		},error = function(e){print(paste("density couldn't be fitted")); RLEseg$dPats[i]=median(na.omit(dPat[,ind]))})

		}else{dPats[dPats[,1]==ind,chr] = 1- as.numeric(dMats[dMats[,1]==ind,chr]) }#end if statement
		
		if(sum(is.na(RLEseg$dMats[i]))==0 & sum(is.na(RLEseg$dPats[i]))==0 ){
		RLEseg$Intp[i] <- autintp(RLEseg$values[i],RLEseg$dMats[i],RLEseg$dPats[i],Params)
		}else{print(paste("Either dPat or dMat for ",chr,":",RLEseg$starts[i],"-",RLEseg$stops[i],"couldn't be determined",sep=""))}
		
		}else{print(paste("no rows for Chr",chr,":",RLEseg$starts[i],"-",RLEseg$stops[i],sep=""))}
		
		}#end i loop
		
		ChrS <- cbind(chr,do.call(cbind,RLEseg))
		if(chr == "1"){ChrSs <- ChrS }else{ChrSs <- rbind(ChrSs,ChrS)}
    
		#Select for aberrant outcome
		Aberrant1 <- ChrSs[which(ChrSs[,"Intp"] == "TT"),]
		Aberrant2 <- ChrSs[which(ChrSs[,"Intp"] == "UM"),]
		Aberrant3 <- ChrSs[which(ChrSs[,"Intp"] == "UP"),]
		Aberrant4 <- ChrSs[which(ChrSs[,"Intp"] == "MP"),]
		Aberrant5 <- ChrSs[which(ChrSs[,"Intp"] == "MM"),]
		Aberrant6 <- ChrSs[which(ChrSs[,"Intp"] == "NN"),]
		Aberrant7 <- ChrSs[which(ChrSs[,"Intp"] == "TP"),]
		Aberrant8 <- ChrSs[which(ChrSs[,"Intp"] == "TM"),]
		Aberrant9 <- ChrSs[which(ChrSs[,"Intp"] == "TT_m"),]
		
		Aberrant <- rbind(Aberrant1,Aberrant2,Aberrant3,Aberrant4,Aberrant5,Aberrant6,Aberrant7,Aberrant8,Aberrant9)
	  Aberrant <- Aberrant[order(Aberrant[,"chr"], Aberrant[,"cumlengths"]),]
		
	}#end chr loop

print(ind)
	write.table(ChrSs,paste(Output_CNV_Det,ind,"_interpretedSegmented.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)

	write.table(Aberrant,paste(Output_CNV_Det,ind,"_Aberrant.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
	
	
}#end ind loop	

