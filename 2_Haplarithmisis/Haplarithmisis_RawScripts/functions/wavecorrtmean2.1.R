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

wavecorrtmean2.1 <- function(logRsMed4,ParScore,GC3,Family,outPath,SibPattern){

print("Applying wave correction with trimmed mean adjustment ...")


for(ind in colnames(logRsMed4)[-c(1:3)]){
	
	logrBe <- logRsMed4[,ind]	
	while(sum(is.na(logrBe))>=1){logrBe[which(is.na(logrBe))]<-logrBe[which(is.na(logrBe))-1]}
	FitLogR <- loessFit(logrBe,GC3[,"GC"])
	#logrAfNorm <- ((logrBe - FitLogR$fitted)-mean(na.omit(logrBe - FitLogR$fitted)))/sd(na.omit(logrBe - FitLogR$fitted))
	DipChrs <- na.omit(rownames(ParScore[[ind]])[ParScore[[ind]][,"Par"]==12])
	if(length(DipChrs)>=1){
		
		for(chr in DipChrs){
			DipPo <- which(logRsMed4[,"Chr"]==chr)
			if(chr==DipChrs[1]){DipPos<- DipPo}else{DipPos<-c(DipPos,DipPo)}
			}
	logrAfNorm <- (logrBe - FitLogR$fitted)-mean(na.omit((logrBe - FitLogR$fitted)[DipPos]),trim=0.10)
	print(paste("timmed mean value is",mean(na.omit((logrBe - FitLogR$fitted)[DipPos]),trim=0.10)))
        print(paste("timmed mean value after meand and GC-correction is",mean(logrAfNorm[DipPos],trim=0.10)))
	print(paste("mean of GC-corrected logRs is",mean(logrBe - FitLogR$fitted)))

	}else if (sum(grep(SibPattern,ind))==0){
		
	logrAfNorm <- (logrBe - FitLogR$fitted)-mean(na.omit((logrBe - FitLogR$fitted)[-c(which(logRsMed4[,2]=="X" | logRsMed4[,2]=="Y" | logRsMed4[,2]=="XY"))]),trim=0.10)
		print(paste(ind, "is a multi-cell sample"))
	}else{	
		logrAfNorm <- (logrBe - FitLogR$fitted)
		print(paste("No diploid chromosome in",ind))
	}
	
	logRsMed4[,ind]<- logrAfNorm
 		
	print(ind)
}#end ind loop

#write.table(data,paste(outPath,Family,"_logRsMed_WaveCorr_Mean.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

print("GC-corrected file was written")

logRsMed4

}#end function

