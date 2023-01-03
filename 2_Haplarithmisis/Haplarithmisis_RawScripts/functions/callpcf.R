#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					    							Function: callPCF 												   	%%%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): logRs gamma plateau
# 
#
#(->): Segmented logRs
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
callpcf <- function(logRs,gammaSC,gammaMC,plateau,Family,outPath){

#source("fastPCF.R")

print("PCF segmentation is applying...")

SegGenomes <- NULL

for(i in 4:ncol(logRs)){
	
	if(sum(grep("_MC",colnames(logRs)[i]))==1){gamma=gammaMC}else if(sum(grep("E*_Bl*",colnames(logRs)[i]))==1 | sum(grep("S",colnames(logRs)[i]))==1){gamma=gammaSC}else{gamma=gammaMC}
		
        SegGenome <- NULL
	for(chr in unique(logRs[,"Chr"])){
		
		logRChr <- logRs[as.character(logRs$Chr)==chr,i]
		#if(sum(is.na(logRChr))>=1){warning(paste("There are",sum(is.na(logRChr),"missing values in logR-values of Chr.",chr))}
		while(sum(is.na(logRChr))>=1){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}
		sdev <- getMad(logRChr,k=plateau)
		res <- selectFastPcf(logRChr,3,gamma*sdev,T)
		SegChr <- res$yhat
		SegGenome <- c(SegGenome,SegChr)
		
	}#end chr loop
	
	SegGenomes<-cbind(SegGenomes,SegGenome)
	
	print(paste(colnames(logRs)[i],"==> gamma",gamma, "is applied"))
	
}#end file loop

SegLogRs <- cbind(logRs[,c("Name","Chr","Position")],SegGenomes)
colnames(SegLogRs)<- colnames(logRs)
write.table(SegLogRs,paste(outPath,Family,"_gammaSc",gammaSC,"_gammaMc",gammaMC,"_logRseg.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
show(paste("Segmented logR-values (",Family, "_logRseg.txt ) is written on",outPath))

SegLogRs

}#end function
