#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					    							Function: avgint 												   	%%%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): logRs Int
# 
#
#(->): Avg logRs
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avgint <- function(logRs,Int,Family,outPath){

show("Averaging over the aberrant region is applying...")

for(i in 4:ncol(logRs)){
	
#	if(sum(grep("E*_Bl*",colnames(logRs)[i]))==1){gamma=gammaSC}else{gamma=gammaMC}

	
	
	for(chr in c(1:22,"X","Y")){
		
		logRChr <- logRs[as.character(logRs$Chr)==chr,i]
		AnnotChr <- logRs[as.character(logRs$Chr)==chr,c("Name","Chr","Position")]
		#if(sum(is.na(logRChr))>=1){print(paste("There are",sum(is.na(logRChr),"missing values in logR-values of Chr.",chr))}
		while(sum(is.na(logRChr))>=1){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}
		
		AvgLogRChr <- matrix(0,1,length(logRChr))
		 
		if(sum(chr==as.character(Int[,1]))==1){
			
			IntChr <- as.matrix(Int[as.character(Int[,1])==chr,])
			
			AvgLog1 <- mean(na.omit(logRChr[AnnotChr$Position>=as.numeric(IntChr[2]) & AnnotChr$Position<=as.numeric(IntChr[3])])) #Average of abnormal region
			AvgLog2 <- mean(na.omit(c(logRChr[AnnotChr$Position < as.numeric(IntChr[2])], logRChr[AnnotChr$Position > as.numeric(IntChr[3])]))) #Average of normal region		
			AvgLogRChr[] <- AvgLog2
			AvgLogRChr[AnnotChr$Position>=as.numeric(IntChr[2]) & AnnotChr$Position<=as.numeric(IntChr[3])] <- AvgLog1
						
			print(paste("Aberrant region of chromosome",chr, "is incorporated" ))
		}else{
			AvgLog <- mean(na.omit(logRChr))
			AvgLogRChr[] <- AvgLog}
		    #show(paste(chr,"===>",length(AvgLogRChr)))

		if(chr==1){AvgLogR <- AvgLogRChr}else{AvgLogR <- c(AvgLogR,AvgLogRChr)}
		
	}#end chr loop
	
	if(i==4){AvgLogRs <- AvgLogR}else{AvgLogRs<-cbind(AvgLogRs,AvgLogR)}
	
	
}#end i loop

AvgLogRs <- cbind(logRs[,c("Name","Chr","Position")],AvgLogRs)
colnames(AvgLogRs) <-colnames(logRs)
write.table(AvgLogRs,paste(outPath,Family,"_logRavg.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
show(paste("Averaged logR-values (",Family, "_logRavg.txt) is written on",outPath,sep=""))

AvgLogRs

}#end function
