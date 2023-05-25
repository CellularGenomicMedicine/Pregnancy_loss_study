###########################################################################################################################
# Author: Rick Essers (Adapted from script by Masoud Zamani Esteki)
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To extract genomic coordinates and haplotyping data values of detected aberrations. 
#                 This data is can be for the calculation of mosaicism degree (in percentage) of aberrations. 
#                 The haplotyping data consists of the P1, P2, M1, M2, values that represent the paternal and maternal contribution to the genome. 
#                 and aberration. 

# input: P1, P2, M1, M2, P1Seg, P2Seg, M1Seg, M2Seg, LogRSeg .txt files for a specific pregnancy loss family. 
#        These files are produced by haplarithmisis. 

# output: Per fetal DNA sample (Extraembroynic mesoderm and chorionic villi), a file is produced containing all values of P1, P2, M1, M2, 
#         grouped together by value. In addition, the genomic coordinates, length, cumulative length, start and stop values are provided. 
#         This information can be used to extract the genomic coordinates and/or maternal and paternal haplotyping values that can be used 
#         to calculate the mosaicism degree of detected aberrations. 


###########################################################################################################################

rm(list=ls(all=T))
require(pastecs)

Family = "PL0001.adj" 

dataPathInit  <- paste0("/Projects/PregnancyLoss/Data/Output/",Family,"/")
outPath       <- paste0("/Projects/PregnancyLoss/Data/Output/",Family,"/")

outPathOutput <- paste(outPath,"/OutputAutSeg/",sep="")
if (!file.exists(outPathOutput)){
  dir.create(outPathOutput)
}


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

	dMats <- matrix(NA,length(Inds),23)
	colnames(dMats) <- unique(logRsSeg$Chr) 
	dMats[,1] <- Inds 
	dPats <- dMats
	
	ind <- "E01_Bl01"
	chr <- "1" 
	
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


