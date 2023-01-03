#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					   Function: QC for single-cell genotypes 												   	%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): Gtypes
# 
#
#(->): Chr-sepcific QC
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qcgtype <- function(Gtypes,ChrPos,Family,SibPattern,outPath){


print(paste("Chr-specific QC analysis for",Family,"family"))

#chroms <- c(1:22,"X","Y")
chroms <- unique(Gtypes$Chr)

Sibs <- c(grep(SibPattern,colnames(Gtypes)),grep("Affected",colnames(Gtypes))) 
ColNames<-gsub("E0_Bl0_","",colnames(Gtypes))

QC <- vector("list",4)

print("#1 ==> Call-rate computation...")

for(ind in 4:ncol(Gtypes)){
	
	Gtype <- as.character(Gtypes[,ind])
	CallRate <- callrate(Gtypes[,ind])
	
	for(chr in chroms){

		CallRateChr <- callrate(Gtypes[ChrPos$Chr==chr,ind])
		
		if(chr == chroms[1]){
			CallRateChrs <- CallRateChr[1]
			CallRateChrsHet <- CallRateChr[4]
			CallRateChrsHom <- CallRateChr[5]
		}else{
			CallRateChrs <- cbind(CallRateChrs,CallRateChr[1])
			CallRateChrsHet <- cbind(CallRateChrsHet,CallRateChr[4])
			CallRateChrsHom <- cbind(CallRateChrsHom,CallRateChr[5])
		}
		
	}#end chr loop
	
	CallRateChrs <- cbind(CallRateChrs,CallRate[1])
	CallRateChrsHet <- cbind(CallRateChrsHet,CallRate[4])
	CallRateChrsHom <- cbind(CallRateChrsHom,CallRate[5])
	
	if(ind==4){
		CallRates<-CallRate
		CallRateChrsInd<-CallRateChrs
		CallRateChrsHetInd<-CallRateChrsHet
		CallRateChrsHomInd<-CallRateChrsHom
	}else{
		CallRates<-cbind(CallRates,CallRate)
		CallRateChrsInd<-rbind(CallRateChrsInd,CallRateChrs)
		CallRateChrsHetInd<-rbind(CallRateChrsHetInd,CallRateChrsHet)
		CallRateChrsHomInd<-rbind(CallRateChrsHomInd,CallRateChrsHom)

		}
	print(colnames(Gtypes)[ind])
		
}#end ind loop

rownames(CallRateChrsInd) <- ColNames[4:ncol(Gtypes)]
rownames(CallRateChrsHetInd) <- ColNames[4:ncol(Gtypes)]
rownames(CallRateChrsHomInd) <- ColNames[4:ncol(Gtypes)]

colnames(CallRateChrsInd) <- c(chroms,"Genome")
colnames(CallRateChrsHetInd) <- c(chroms,"Genome")
colnames(CallRateChrsHomInd) <- c(chroms,"Genome")


print("#2 ==> Mendelian inconsistency computation...")

Sibs <- unique(Sibs)
for(sib in Sibs){
	

	Child <- Gtypes[,sib]
	#Genome-wide mendelina inconsistency rate for autosomes
	MendIncRateAut <- mendinc(Gtypes$Father[Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y"],Gtypes$Mother[Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y"],Child[Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y"])
		
	#Chr-specific mendelina inconsistency rate for autosomes
	for(chr in chroms){
		
		MendIncRateChr <- mendinc(Gtypes$Father[Gtypes$Chr==chr],Gtypes$Mother[Gtypes$Chr==chr],Child[Gtypes$Chr==chr])
		
		if(chr==chroms[1]){
			MendIncRateChrs <- MendIncRateChr
		}else{
			MendIncRateChrs <- cbind(MendIncRateChrs,MendIncRateChr)
		}	
	
	}#end chr loop
	
	MendIncRateChrs <- cbind(MendIncRateChrs,MendIncRateAut)
	
	if(sib==Sibs[1]){MendIncRateChrsInd <-  MendIncRateChrs}else{MendIncRateChrsInd <- rbind(MendIncRateChrsInd,MendIncRateChrs)}

	print(colnames(Gtypes)[sib])

}#end sib loop 

if(length(Sibs)==1){MendIncRateChrsInd <- as.matrix(MendIncRateChrsInd)}
rownames(MendIncRateChrsInd) <- ColNames[Sibs]
colnames(MendIncRateChrsInd) <- c(chroms,"GenomeAut")

QC[[1]]<- CallRateChrsInd
QC[[2]]<- CallRateChrsHetInd
QC[[3]]<- CallRateChrsHomInd
QC[[4]]<- MendIncRateChrsInd
names(QC)<-c("CallRateChrsInd","CallRateChrsHetInd","CallRateChrsHomInd","MendIncRateChrsInd")
write.table(MendIncRateChrsInd,paste(outPath,Family,"_ChrSpec_MendInc.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
write.table(CallRates,paste(outPath,Family,"_CallRates.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
write.table(CallRateChrsInd,paste(outPath,Family,"_ChrSpec_CallRates.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
write.table(CallRateChrsHetInd,paste(outPath,Family,"_ChrSpec_Het_CallRates.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
write.table(CallRateChrsHomInd,paste(outPath,Family,"_ChrSpec_Hom_CallRates.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
print(paste("QC files were written on",outPath))

QC
}#end SNPcallQC

