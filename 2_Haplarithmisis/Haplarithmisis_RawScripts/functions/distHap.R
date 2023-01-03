disthap<-function(Family,Int,outPath,SibPattern){

Haps <- read.table(paste(outPath,"/",Family,"_Itp2.hap",sep=""),header=T,sep="\t",stringsAsFactors =F)



Chroms <-Int[Int[,1]!="0",1]
Scs <- colnames(Haps)[grep(SibPattern,colnames(Haps))]


for(ind in Scs){
	
	for(chr in Chroms){
		
		IntChr <- Int[Int[,1]==chr,]	
		#ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",chr,sep=""),2])
		HapChr <- Haps[Haps[,"Chr"]==chr,c("Position",ind)]	
		UP <- HapChr[HapChr$Position<=as.numeric(IntChr[2]),]
		Down <- HapChr[HapChr$Position>=as.numeric(IntChr[3]),]

		UPRle <- rle(UP[,ind])
		DownRle <- rle(Down[,ind])

		#if(UPRle$values[length(UPRle$values)]==DownRle$values[1]){
			Dist <- cbind(ind,chr,UPRle$lengths[length(UPRle$values)],DownRle$lengths[1],UPRle$values[length(UPRle$values)], DownRle$values[1])#}

		
		if(chr == Chroms[1]){Dists <- Dist}else{Dists<- rbind(Dists,Dist)}
	
	}#end chr
	
	if(ind==Scs[1]){DistFam <- Dists }else{DistFam <-rbind(DistFam,Dists)}
	print(ind)
		
}#end ind loop

colnames(DistFam) <- c("Haplotype","Chr","LengthUp","LengthDown","ValueUp","ValueDown")

write.table(DistFam,paste(outPath,"DistanceToHR_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
}
