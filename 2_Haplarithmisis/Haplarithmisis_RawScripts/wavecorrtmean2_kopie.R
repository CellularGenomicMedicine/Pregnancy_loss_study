#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					   Function: GC correction Median             	   																%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): dataFile outPath
# 
#
#(->): GC-corrected file
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavecorrtmean2 <- function(logRsMed,ParScore,GC2,Func,Family,outPath,SibPattern){

print("Applying wave correction with trimmed mean adjustment ...")


#logRsMed <- logRsMedRaw
#rownames(logRsMed) <- as.character(logRsMed[,"Name"])
#rownames(GC) <- as.character(GC[,"Name"])
rowsTot <- intersect(rownames(logRsMed),rownames(GC2))

data2 <- logRsMed[rowsTot,]


GC2 <-data.frame(GC2[rowsTot,],stringsAsFactors=FALSE)
GC2[,"GC"] <- as.numeric(GC2[,"GC"])

for(ind in colnames(logRsMed)[-c(1:3)]){
	
	logrBe <- data2[,ind]	
	while(sum(is.na(logrBe))>=1){logrBe[which(is.na(logrBe))]<-logrBe[which(is.na(logrBe))-1]}
	FitLogR <- loessFit(logrBe,GC2[,"GC"])
	#logrAfNorm <- ((logrBe - FitLogR$fitted)-mean(na.omit(logrBe - FitLogR$fitted)))/sd(na.omit(logrBe - FitLogR$fitted))
	DipChrs <- na.omit(rownames(ParScore[[ind]])[ParScore[[ind]][,"Par"]==12])
	if(length(DipChrs)>=1){
		
		for(chr in DipChrs){
			DipPo <- which(data2[,"Chr"]==chr)
			if(chr==DipChrs[1]){DipPos<- DipPo}else{DipPos<-c(DipPos,DipPo)}
			}
	logrAfNorm <- (logrBe - FitLogR$fitted)-mean(na.omit((logrBe - FitLogR$fitted)[DipPos]),trim=0.10)
	print(paste("timmed mean value is",mean(na.omit((logrBe - FitLogR$fitted)[DipPos]),trim=0.10)))
        print(paste("timmed mean value after meand and GC-correction is",mean(logrAfNorm[DipPos],trim=0.10)))
	print(paste("mean of GC-corrected logRs is",mean(logrBe - FitLogR$fitted)))

	}else if (sum(grep(SibPattern,ind))==0){
		
	logrAfNorm <- (logrBe - FitLogR$fitted)-mean(na.omit((logrBe - FitLogR$fitted)[-c(which(data2[,2]=="X" | data2[,2]=="Y" | data2[,2]=="XY"))]),trim=0.10)
		print(paste(ind, "is a multi-cell sample"))
	}else{	
		logrAfNorm <- (logrBe - FitLogR$fitted)
		print(paste("No diploid chromosome in",ind))
	}
	
	data2[,ind]<- logrAfNorm
 		
	print(ind)
}#end ind loop

#write.table(data,paste(outPath,Family,"_logRsMed_WaveCorr_Mean.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

print("GC-corrected file was written")

data2

}#end function

