patscore2 <- function(dataPo,QC,Chroms,Gtypes,Family){
	
	
	AvgMDA<-2.2862648
	SdMDA<-0.4721314
	MinCallRate <-16.38
	#MinCallRate <-12.38
	
	print("Computing parental scores...")
	dataCallRate <- QC[["CallRateChrsInd"]]
			
#	colnames(dataCallRate) <- c(Chroms,"Y","Genome")
	CallRates <- dataCallRate[,Chroms]

	Pos <- dataPo[,-c(1:3)]
	
	if(ncol(dataPo)==4){
	Pos <- as.matrix(Pos)
	colnames(Pos)<-colnames(dataPo)[4]
	}
	SamplesPo <- vector("list",ncol(Pos))
	names(SamplesPo)<- colnames(Pos)
		
	for(ind in colnames(Pos)){

		PoScore <- matrix(NA,length(Chroms),6)
		rownames(PoScore) <- Chroms
		for(chr in Chroms){ 
		
		PatGtype <- Gtypes[Gtypes[,"Chr"]==chr,grep(paste("Father_",Family,sep=""),colnames(Gtypes))]
		MatGtype <- Gtypes[Gtypes[,"Chr"]==chr,grep(paste("Mother_",Family,sep=""),colnames(Gtypes))]
		InfPo <- sum((MatGtype== "AA" & PatGtype == "BB") | 
							    (MatGtype== "BB" & PatGtype == "AA") | 
							    (MatGtype== "AB" & PatGtype == "AA") | 
							    (MatGtype== "AA" & PatGtype == "AB") | 
							    (MatGtype== "BB" & PatGtype == "AB") | 
							    (MatGtype== "AB" & PatGtype == "BB"))
		
		dataPoChr <- Pos[dataPo[,"Chr"]==chr ,ind]
		Pat <- (abs(sum(dataPoChr[dataPoChr<0]))/InfPo)*100
		Mat <- (sum(dataPoChr[dataPoChr>0])/InfPo)*100
		if(Pat==0){Pat=0.0001} ; if(Mat==0){Mat=0.0001}
	
		PatMat <- Pat/(Mat+Pat)
		MatPat <- Mat/(Mat+Pat)
		
		if(Pat<=(AvgMDA-(3*SdMDA)) & Mat<=(AvgMDA-(3*SdMDA))){PatMat<-0;MatPat<-0}#else if(Pat<=(AvgMDA-(3*SdMDA))){PatMat<-0};if(Mat<=(AvgMDA+(3*SdMDA))){MatPat<-0}
	#p<- t.test(AvgRsq[AvgRsq[,"Chr"]==chr,"AvgRsq"],logRs[data[,"Chr"]==chr,ind])$p.value
	#to be removed
	ind2 <- gsub("E0_Bl0_","",ind)
	ind2 <- gsub("E0Bl0_","",ind2)

	PoScore[chr,]<-cbind(Pat,Mat,PatMat,MatPat,CallRates[ind2,chr],NA)
	#print(paste(ind,"==>",chr,p,round(dataPoPatChr,digits=2),round(dataPoMatChr,digits=2),round(PatMat,digits=2)))
	
	}
	
		colnames(PoScore) <- c("Pat","Mat","PatMat","MatPat","CallRate","Par")
		SamplesPo[[ind]]<-PoScore
		print(ind)
	
}#end ind loop


	SamplesPo2 <- vector("list",ncol(Pos))
	names(SamplesPo2)<- colnames(Pos)

	for(ind in names(SamplesPo)){

		Sample <- SamplesPo[[ind]]
		#for(c in 1:ncol(Sample)){Sample[,c]<- as.numeric(as.character(Sample[,c]))}
		#Sample[Sample[,"CallRate"]>20 & (Sample[,"PatMat"]<0.9 & Sample[,"MatPat"]<0.9),"Par"]<- 12
		Sample[Sample[,"CallRate"]>MinCallRate & ((Sample[,"Pat"]< (AvgMDA+SdMDA)& Sample[,"Mat"]< (AvgMDA+SdMDA)) | (Sample[,"PatMat"]<0.9 & Sample[,"MatPat"]<0.9)),"Par"]<- 12
		Sample[Sample[,"Pat"]> (AvgMDA+SdMDA) & Sample[,"CallRate"]>MinCallRate & Sample[,"PatMat"]>=0.9 ,"Par"]<- 11
		Sample[Sample[,"Mat"]> (AvgMDA+SdMDA) & Sample[,"CallRate"]>MinCallRate & Sample[,"MatPat"]>=0.9  ,"Par"]<- 22
		Sample[Sample[,"CallRate"]<MinCallRate & (Sample[,"PatMat"]<0.5 & Sample[,"MatPat"]<0.5) ,"Par"]<- 00
 		SamplesPo2[[ind]] <- Sample

	}# end ind loop

SamplesPo2

}#end function


#A <- rbind(ParScore2[["S10018_ebv2_MDA_REF_24h"]][1:22,],ParScore2[["S10018_ebv10_MDA_REF_24h"]][1:22,],ParScore2[["S10016_ebv2_MDA_REF_24h"]][1:22,],ParScore2[["S10016_ebv5_MDA_REF_24h"]][1:22,])

#MeanSd <- rbind(apply(A[,1:5],2,mean),apply(A[,1:5],2,sd))
#c(mean(c(A[,1],A[,2])),sd(c(A[,1],A[,2]))) #==>2.2862648 0.4721314
#c(mean(c(A[,3],A[,4])),sd(c(A[,3],A[,4]))) #==>0.50000000 0.05229585


