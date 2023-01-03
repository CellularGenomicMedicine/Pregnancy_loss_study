rm(list=ls(all=T))
require(pastecs)

Family = "PL2016.adj" 
Indication = "Polyploid"
Group = "RPL"

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


dataPathInit <- paste0("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/FinalPlots/AllFams/Abnormal/",Group,"/",Indication,"/")
outPath <- paste0("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/FinalPlots/AllFams/Abnormal/",Group,"/",Indication,"/",Family)

outPathOutput <- paste(outPath,"/OutputAutSeg/",sep="")
if (!file.exists(outPathOutput)){
  dir.create(outPathOutput)
}

source("/Users/G10039937/Surfdrive/ClinicalGenetics/MosaicismDetection/distcompmedseg2.R")
Params <- read.table("/Users/G10039937/Surfdrive/ClinicalGenetics/MosaicismDetection/Ranges.txt",sep="\t",header=T,stringsAsFactors=F,row.names = 1)

#
#for(Family in Families){

	dataPath = paste(dataPathInit,Family,"/",sep="")
	P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	P1Seg <-  read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)	
	M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

	logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc14_gammaMc14_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
	
	Inds = colnames(logRsSeg)[c(grep("E*_Bl*",colnames(logRsSeg)),grep("Affected*",colnames(logRsSeg)))]
#	Inds = Inds[-c(grep("ther",Inds))]

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

		MedLogR <- median(logRsSeg[logRsSeg$Chr==chr,ind])
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
			#print(paste("Max peak mat ==>",MaxPeakPat))
			#print(paste("dPat ==>",XMaxPeakPat))

		}else{dPats[dPats[,1]==ind,chr] = 1- as.numeric(dMats[dMats[,1]==ind,chr]) }#end if statement
		
		if(sum(is.na(RLEseg$dMats[i]))==0 & sum(is.na(RLEseg$dPats[i]))==0 ){
		RLEseg$Intp[i] <- autintp(RLEseg$values[i],RLEseg$dMats[i],RLEseg$dPats[i],Params)
		}else{print(paste("Either dPat or dMat for ",chr,":",RLEseg$starts[i],"-",RLEseg$stops[i],"couldn't be determined",sep=""))}
		
		}else{print(paste("no rows for Chr",chr,":",RLEseg$starts[i],"-",RLEseg$stops[i],sep=""))}
		
		}#end i loop
		
		ChrS <- cbind(chr,do.call(cbind,RLEseg))
		if(chr == "1"){ChrSs <- ChrS }else{ChrSs <- rbind(ChrSs,ChrS)}

	}#end chr loop
  
Coord <- paste0("Chr",ChrSs[,"chr"],":",ChrSs[,"starts"],"-",ChrSs[,"stops"])
ChrSs = cbind(Coord,ChrSs)
  
print(ind)
	write.table(ChrSs,paste(outPathOutput,ind,"_interpretedSegmented.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
	
}#end ind loop	
	
print(Family)
#
#} #end Family loop

####

require(pastecs)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/MosaicismDetection/distcompmedseg2.R")
Params <- read.table("/Users/G10039937/Surfdrive/ClinicalGenetics/MosaicismDetection/Ranges.txt",sep="\t",header=T,stringsAsFactors=F,row.names = 1)

###Do all fams 
###for(Family in Families){

dataPath = paste(dataPathInit,Family,"/",sep="")
P1 <- read.table(paste(dataPath,"P1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2 <- read.table(paste(dataPath,"P2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M1 <- read.table(paste(dataPath,"M1_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2 <- read.table(paste(dataPath,"M2_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P1Seg <-  read.table(paste(dataPath,"P1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
P2Seg <- read.table(paste(dataPath,"P2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)	
M1Seg <- read.table(paste(dataPath,"M1Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)
M2Seg <- read.table(paste(dataPath,"M2Seg_",Family,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

logRsSeg <- read.table(paste(dataPath,Family,"_gammaSc14_gammaMc14_logRseg.txt",sep=""),header=T,sep="\t",stringsAsFactors=F)

Inds = colnames(logRsSeg)[c(grep("E*_Bl*",colnames(logRsSeg)),grep("Affected*",colnames(logRsSeg)))]
#	Inds = Inds[-c(grep("ther",Inds))]

dMats <- matrix(NA,length(Inds),23)
colnames(dMats) <- unique(logRsSeg$Chr) 
dMats[,1] <- Inds 
dPats <- dMats

###
#chr     <- "2"
#ind     <- "Affected_PL2074.adj"
#P1_ind  <- "Affected_PL2074.adj"

#ind     <- "E01_Bl01_PL2074.adj"
#P1_ind  <- "E01_Bl01_PL2074.adj"

#chr <- "1"
###

for(ind in Inds){
  
  ###P1
  for(chr in unique(logRsSeg$Chr)){
    
    if(chr != "X"){
      
      P1SegChr <- P1Seg[P1Seg$Chr==chr,c("Name","Chr","Position",ind)]
      
      RLE_P1 <- rle(P1SegChr[P1SegChr$Chr == chr,ind])
      RLE_P1$cumlengths <- cumsum(RLE_P1$lengths)
      RLE_P1$starts <- rep(NA,length(RLE_P1$lengths))
      RLE_P1$stops <- rep(NA,length(RLE_P1$lengths))
      
      for(i in 1:length(RLE_P1$values)){
        
        RLE_P1$starts[i] <- logRsSeg[RLE_P1$cumlengths[i]-(RLE_P1$lengths[i]-1),"Position"]
        RLE_P1$stops[i] <- logRsSeg[RLE_P1$cumlengths[i],"Position"]
        
      }
    }#end i loop
    
    P1_ChrS <- cbind(chr,do.call(cbind,RLE_P1))
    if(chr == "1"){P1_ChrSs <- P1_ChrS }else{P1_ChrSs <- rbind(P1_ChrSs,P1_ChrS)}
    
  }#end chr loop
  
  Coord <- paste0("Chr",P1_ChrSs[,"chr"],":",P1_ChrSs[,"starts"],"-",P1_ChrSs[,"stops"])
  PatMat <- "P1"
  P1_ChrSs = cbind(PatMat,Coord,P1_ChrSs)  
  
  
  ###P2
  for(chr in unique(logRsSeg$Chr)){
    
    if(chr != "X"){
      
      P2SegChr <- P2Seg[P2Seg$Chr==chr,c("Name","Chr","Position",ind)]
      
      RLE_P2 <- rle(P2SegChr[P2SegChr$Chr == chr,ind])
      RLE_P2$cumlengths <- cumsum(RLE_P2$lengths)
      RLE_P2$starts <- rep(NA,length(RLE_P2$lengths))
      RLE_P2$stops <- rep(NA,length(RLE_P2$lengths))
      
      for(i in 1:length(RLE_P2$values)){
        
        RLE_P2$starts[i] <- logRsSeg[RLE_P2$cumlengths[i]-(RLE_P2$lengths[i]-1),"Position"]
        RLE_P2$stops[i] <- logRsSeg[RLE_P2$cumlengths[i],"Position"]
        
      }    
    }#end i loop
    
    P2_ChrS <- cbind(chr,do.call(cbind,RLE_P2))
    if(chr == "1"){P2_ChrSs <- P2_ChrS }else{P2_ChrSs <- rbind(P2_ChrSs,P2_ChrS)}
    
  }#end chr loop
  
  P2_Coord <- paste0("Chr",P2_ChrSs[,"chr"],":",P2_ChrSs[,"starts"],"-",P2_ChrSs[,"stops"])
  P2_PatMat <- "P2"
  P2_ChrSs = cbind(P2_PatMat,P2_Coord,P2_ChrSs)
  
  
  ###M1
  for(chr in unique(logRsSeg$Chr)){
    
    M1SegChr <- M1Seg[M1Seg$Chr==chr,c("Name","Chr","Position",ind)]
    
    RLE_M1 <- rle(M1SegChr[M1SegChr$Chr == chr,ind])
    RLE_M1$cumlengths <- cumsum(RLE_M1$lengths)
    RLE_M1$starts <- rep(NA,length(RLE_M1$lengths))
    RLE_M1$stops <- rep(NA,length(RLE_M1$lengths))
    
    for(i in 1:length(RLE_M1$values)){
      
      RLE_M1$starts[i] <- logRsSeg[RLE_M1$cumlengths[i]-(RLE_M1$lengths[i]-1),"Position"]
      RLE_M1$stops[i] <- logRsSeg[RLE_M1$cumlengths[i],"Position"]
      
    }#end i loop
    
    M1_ChrS <- cbind(chr,do.call(cbind,RLE_M1))
    if(chr == "1"){M1_ChrSs <- M1_ChrS }else{M1_ChrSs <- rbind(M1_ChrSs,M1_ChrS)}
    
  }#end chr loop
  
  M1_Coord <- paste0("Chr",M1_ChrSs[,"chr"],":",M1_ChrSs[,"starts"],"-",M1_ChrSs[,"stops"])
  M1_PatMat <- "M1"
  M1_ChrSs = cbind(M1_PatMat,M1_Coord,M1_ChrSs)  
  
  
  ###M2
  for(chr in unique(logRsSeg$Chr)){
    
    M2SegChr <- M2Seg[M2Seg$Chr==chr,c("Name","Chr","Position",ind)]
    
    RLE_M2 <- rle(M2SegChr[M2SegChr$Chr == chr,ind])
    RLE_M2$cumlengths <- cumsum(RLE_M2$lengths)
    RLE_M2$starts <- rep(NA,length(RLE_M2$lengths))
    RLE_M2$stops <- rep(NA,length(RLE_M2$lengths))
    
    for(i in 1:length(RLE_M2$values)){
      
      RLE_M2$starts[i] <- logRsSeg[RLE_M2$cumlengths[i]-(RLE_M2$lengths[i]-1),"Position"]
      RLE_M2$stops[i] <- logRsSeg[RLE_M2$cumlengths[i],"Position"]
      
    }#end i loop
    
    M2_ChrS <- cbind(chr,do.call(cbind,RLE_M2))
    if(chr == "1"){M2_ChrSs <- M2_ChrS }else{M2_ChrSs <- rbind(M2_ChrSs,M2_ChrS)}
    
  }#end chr loop
  
  M2_Coord <- paste0("Chr",M2_ChrSs[,"chr"],":",M2_ChrSs[,"starts"],"-",M2_ChrSs[,"stops"])
  M2_PatMat <- "M2"
  M2_ChrSs = cbind(M2_PatMat,M2_Coord,M2_ChrSs)  
  
  
  PatMatBAFSeg <- rbind(P1_ChrSs,P2_ChrSs,M1_ChrSs,M2_ChrSs)	
  
  print(ind)
  write.table(PatMatBAFSeg,paste(outPathOutput,ind,"_PatMatBAFSeg.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
  
}#end ind loop	

print(Family)
#
#} #end Family loop



