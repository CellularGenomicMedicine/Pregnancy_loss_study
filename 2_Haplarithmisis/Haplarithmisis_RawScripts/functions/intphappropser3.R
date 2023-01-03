intphappropser3 <- function(outPath,Family,ParScore,Haps,SibPattern){


	dataHap <- Haps[["dataHap"]]
	dataHapRaw <- Haps[["dataHapRaw"]]

	MinLength = 300
	MinFrac = 0.6
	dataHap1 <- dataHap
	Chroms <- c(1:22,"X")
	Samps <- vector("list",ncol(dataHap)-3)
	names(Samps)<-colnames(dataHap[,-c(1:3)])

	ParScore2 <- ParScore
	Inds <- gsub("_Mat","",colnames(dataHap)[-c(1:3)])
	Inds <- gsub("_Pat","",Inds)
	Inds <- unique(Inds)
	
	for(ind in Inds){
				
		ParScore[[ind]] <- ParScore[[ind]][!is.na(ParScore[[ind]][,"Par"]),]   
		if(sum(nrow(ParScore[[ind]]))>1){
 		PatChr <- rownames(ParScore[[ind]])[ParScore[[ind]][,"Par"] ==11]
		MatChr <- rownames(ParScore[[ind]])[ParScore[[ind]][,"Par"] ==22]
		NullChr <- rownames(ParScore[[ind]])[ParScore[[ind]][,"Par"] ==0 ]
		if(length(PatChr)!=0){for(chr in PatChr){dataHap[dataHap[,"Chr"]==chr,paste(ind,"_Mat",sep="")]<- 0}}
		if(length(MatChr)!=0){for(chr in MatChr){dataHap[dataHap[,"Chr"]==chr,paste(ind,"_Pat",sep="")]<- 0}}
		if(length(NullChr)!=0){for(chr in NullChr){dataHap[dataHap[,"Chr"]==chr,paste(ind,"_Pat",sep="")]<- 0
		dataHap[dataHap[,"Chr"]==chr,paste(ind,"_Mat",sep="")]<- 0}}
		print(ind)
		}else{print(paste(ind,"sample sbould be excluded for diagnosis, as the origin of the chromosomes could not be determined"))}
	
	}#end ind loop

ParScore <- ParScore2
print("--------------------------------------------------------------")
print("Probability computation...")

Blasts <- colnames(dataHap)[grep(SibPattern,colnames(dataHap))]
for(samp in Blasts){

	Indiv <- gsub("_Pat","",samp)
	Indiv <- gsub("_Mat","",Indiv)
		
	for(chr in Chroms){
		
		if(!is.na(ParScore[[Indiv]][chr,"Par"])){
		dataHapChr <- dataHap[dataHap$Chr==chr,]

		IntBlock <- rle(dataHap[dataHap$Chr==chr,samp])
		IntBlock <- cbind(IntBlock$values,IntBlock$lengths)
		colnames(IntBlock)<- c("values","lengths")
		IntBlock <- data.frame(rbind(IntBlock,c(0,0)))

		CumLength <- rbind(cumsum(IntBlock$lengths),IntBlock$values)

		if(nrow(IntBlock)==0){
			CumLength<- cbind(0,0,0,0)
			}else{
			InfSNPs <- rep(NA,ncol(CumLength))
			InfSNPsDisc <- rep(NA,ncol(CumLength))
			Frac <- rep(NA,ncol(CumLength))
			Annot <- matrix(NA,ncol(CumLength),3)
			
			for(i in 1:ncol(CumLength)){
				if(i==1){start=1}else{start=CumLength[1,i-1]+1}
				if(start>CumLength[1,i]){start=CumLength[1,i]}
				InfSNPs[i] <- sum(dataHapRaw[dataHap$Chr==chr,samp][start:CumLength[1,i]]==CumLength[2,i])
				InfSNPsDisc[i] <- sum(dataHapRaw[dataHap$Chr==chr,samp][start:CumLength[1,i]]!=CumLength[2,i] & dataHapRaw[dataHap$Chr==chr,samp][start:CumLength[1,i]]!=0)
				Frac[i] <- round(InfSNPs[i]/(InfSNPs[i]+InfSNPsDisc[i]),digits=3)
				Annot[i,] <- cbind(chr,dataHapChr$Position[start],dataHapChr$Position[CumLength[1,i]])
				if(Frac[i]=="NaN"){Frac[i] <- 0}#when no informative SNP is present for the smoothed haplotype block
				if( Frac[i] < MinFrac){dataHap[dataHap$Chr==chr,samp][start:CumLength[1,i]] <- 0}


			}#end i loop
			
				CumLength<- cbind.data.frame(Annot,IntBlock$lengths,t(CumLength),InfSNPs,InfSNPsDisc,Frac)
		}
		Chroms <- rownames(ParScore[[Indiv]])
		
		if(chr==Chroms[which(ParScore[[Indiv]][,"Par"]!="NA")[1]]){CumLengths<-CumLength}else{CumLengths<-rbind(CumLengths,CumLength)}
		}else{print(paste("The origin of chr.",chr, "of", Indiv,"could not be determined"))}
	
	}#end chr loop

	if(exists("CumLengths")==FALSE){CumLengths <- matrix(0,2,9)}	
	colnames(CumLengths) <- c("Chr","Start","Stop","Length","CumLength","Value","#InfSNPs","#InfSNPsDisc","Frac")
	CumLengths[CumLengths[,"Frac"]<MinFrac,"Value"]<-0
	CumLengths <- CumLengths[CumLengths[,"Value"]!=0,]
	Samps[[samp]]<-CumLengths
	print(samp)
	rm(list="CumLengths")
}#end samp loop


MDAs <- do.call(rbind,Samps[grep(SibPattern,names(Samps))])
 

X <-MDAs[which(MDAs$Chr=="X"),]
MDAs <- rbind(MDAs[MDAs$Chr!="X",],X[grep("Mat",rownames(X)),])
MDAs <- MDAs[MDAs[,"Value"]!=0,]
MDAs <- MDAs[MDAs[,"Length"]>=MinLength,]

#SubBls <- read.table(paste("/uz/data/avalok/symbiosys/raw/Masoud/Htyping/GCcorrectLogRs/SubBls_bin10kb_20130624_184618.txt",sep=""),header=F,sep="\t")

#for(s in as.character(SubBls[,1])){if(sum(grep(s,rownames(MDAs)))!=0){MDAs <- MDAs[-c(grep(s,rownames(MDAs))),];print(paste(s,"is substandard and excluded"))}}

#BadSampsPic <- c("E02_Bl364","E04_Bl357","E06_Bl372","E09_Bl360","E11_Bl352","E15_Bl362")
#for(i in BadSampsPic){if(sum(grep(i,rownames(MDAs)))!=0){MDAs <- MDAs[-c(grep(i,rownames(MDAs))),];print(paste(i,"is amplified with PicoPlex and excluded"))}}


MDAs[,"Chr"] <- as.character(MDAs[,"Chr"])
MDAs[,"Start"] <- as.numeric(as.character(MDAs[,"Start"]))
MDAs[,"Stop"] <- as.numeric(as.character(MDAs[,"Stop"]))


LmFit <- rlm(MDAs[,"#InfSNPs"]~MDAs[,"Length"], data=MDAs)
DataCompl<-cbind(MDAs,LmFit$fitted.values,LmFit$resid,rep(NA,nrow(MDAs)))
colnames(DataCompl)[(ncol(DataCompl)-2):ncol(DataCompl)]<-cbind("Fitted","Resid","Prob")
NegRes <- subset(DataCompl, Resid < 0)


zscores<-scale(DataCompl$Resid,center=TRUE,scale=TRUE)
datawithzscores<-cbind(DataCompl,zscores)
DataCompl[,"Prob"] <- round(pnorm(-abs(datawithzscores$zscores))+0.5,digits=3)
DataCompl[DataCompl[,"Resid"]>0,"Prob"] <- 1


write.table(DataCompl,paste(outPath,Family,"_PerBlockRelScores2.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
write.table(dataHap,paste(outPath,Family,"_Itp2.hap",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

Intp <- vector("list",2)
names(Intp) <- c("dataHap","DataCompl")
Intp[["dataHap"]] <- dataHap 
Intp[["DataCompl"]] <- DataCompl
Intp

}

