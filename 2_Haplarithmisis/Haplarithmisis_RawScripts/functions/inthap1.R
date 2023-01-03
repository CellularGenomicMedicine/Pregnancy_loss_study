#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      		        			          Function: inthap            	   											     %%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): raw haplotypes
# 
#
#(->): interpreted haplotypes
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inthap1 <- function(dataHap,Window,VgTh){
	
print("Interpreting the raw haplotypes...")

TipChr = 100#currently is not used 
#Window = 450
#VgTh = 200
Chroms <- c(1:22,"X")
#Paternity <- c("Pat","Mat")

for(h in colnames(dataHap)[-c(1:3)]){
	

	for(chr in Chroms){
		
		HapChr <- dataHap[dataHap$Chr==chr,h]
		
		Par1 <- which(HapChr==1);
		Par2 <- which(HapChr==2);
		
		if(sum(Par1)==0){
			HapSm <- rep(2,1,Par2[1])
			print("Probably the haplotype of the reference individual is reconstructed!!")
			print(paste(h,"Chr.",chr))
		}else if(sum(Par2)==0){
			HapSm <- rep(1,1,Par1[1])
			print("Probably the haplotype of the reference individual is reconstructed!!")
			print(paste(h,"Chr.",chr))
		}else if(Par1[1]<Par2[1]){
			HapSm <- rep(1,1,Par1[1])		
		}else if(Par2[1]<Par1[1]){
			HapSm <- rep(2,1,Par2[1])	
		}else{
			print("No sparse blocks!")}	
		
		HapInt <- rep(0,length(HapChr))
		
		for(k in (length(HapSm)+1):length(HapChr)){
			
			if(k== (length(HapSm)+1) & HapChr[k] ==0){
				HapInt[k] <- HapSm[k-1]
			}else if(HapChr[k]==0){
				HapInt[k] <- HapInt[k-1]
			}else{
				HapInt[k] <- HapChr[k]}
			
		}#end k loop
		
		if(length(unique(HapChr))==1){HapInt=HapChr}
		if(chr==1){HapIntGenome<-HapInt}else{HapIntGenome<-c(HapIntGenome,HapInt)}
		
		}#end chr loop

#===================================
#                        First cycle


for(chr in Chroms){
	
	HapIntChr <- HapIntGenome[dataHap$Chr==chr]
	
	#Defining chr-specific window size 
	WindowChr <- round((Window*length(dataHap$Chr==chr))/length(dataHap$Chr==1))
	if((WindowChr%%2)==0){WindowChr<-WindowChr+1}
	HapBlocks <- rle(HapIntChr)
	if(length(HapBlocks$values)>2){
	for(k in 2:(length(HapBlocks$values)-1)){
		
		if(HapBlocks$values[k]==0 & HapBlocks$values[k+1] == HapBlocks$values[k-1]){
			HapBlocks$value[k] <- HapBlocks$values[k-1]
		}		
		
	}#end k loop
}
#===================================
#                        Second cycle
	
	HapInt <- inverse.rle(HapBlocks)
	HapInt <- medfilt1(HapInt,WindowChr)
	HapBlocks <- rle(as.numeric(HapInt))
	HapBlocks$values[HapBlocks$values==0.5] <- 1
    HapInt <- inverse.rle(HapBlocks)
    HapBlocks <- rle(HapInt)
	
	#Notice: A comment out section is present in MATLAB script

#===================================
#                        Third cycle
	
	#HapIntSm <- inverse.rle(HapBlocks)
	CumLengthBlocks <- cumsum(HapBlocks$lengths)
	VagueBlocks <- which(CumLengthBlocks <= VgTh) 
	
	#Decreasing vagueness
if(vb>1){
	for(vb in  VagueBlocks){
		
		if(vb==1 & HapBlocks$values[vb] != as.numeric(which.max(table(HapChr[1:CumLengthBlocks[vb]])))){
			HapBlocks$Values[vb] <- HapBlocks$Values[vb+1] 
		}else if(vb==length(HapBlocks$lengths) &  HapBlocks$values[vb] != as.numeric(which.max(table(HapChr[CumLengthBlocks[vb-1]:CumLengthBlocks[vb]]))) ){
			HapBlocks$Values[vb] <-  HapBlocks$Values[vb-1]
		}else if( HapBlocks$values[vb] != as.numeric(which.max(table(HapChr[CumLengthBlocks[vb-1]:CumLengthBlocks[vb]]))) & HapBlocks$values[vb-1]==HapBlocks$values[vb+1]){
			HapBlocks$Values[vb] <-  HapBlocks$Values[vb-1]}
		
	}#end vb loop
}
	SmHap <- inverse.rle(HapBlocks)
	HapBlocksSm <- rle(SmHap)
	NotAssBocks <- which(HapBlocksSm$values==1.5) # Blocks which are not assigned and most probably occured in between of 1 and 2 assigned haplotypes (i.e. putative HR-sites)
	
	for(ns  in NotAssBocks){
		
		if(ns!=1 & ns!=length(HapBlocksSm$values) & (HapBlocksSm$values[ns-1]==HapBlocksSm$values[ns+1])){
			HapBlocksSm$values[ns] <- HapBlocksSm$values[ns-1]
		}else if(ns!=length(HapBlocksSm$values)){
			HapBlocksSm$values[ns] <- HapBlocksSm$values[ns-1]
		}else if (ns==1){
			HapBlocksSm$values[ns] <- HapBlocksSm$values[ns+1]}
		
	}#end ns loop

	HapBlocksChrIntFinal <- inverse.rle(HapBlocksSm)
	
	if(chr==1){HapBlocksIntGenome <- HapBlocksChrIntFinal}else{HapBlocksIntGenome <- c(HapBlocksIntGenome,HapBlocksChrIntFinal)}


}#end chr loop

if(h==colnames(dataHap)[-c(1:3)]){ HapIntFamily <- HapBlocksIntGenome}else{HapIntFamily <- cbind(HapIntFamily,HapBlocksIntGenome)}


print(paste("Haplotypes of",h,"is interpreted"))
}#end h loop	

colnames(HapIntFamily)<- colnames(dataHap)[-c(1:3)]
dataHapInt <- cbind(dataHap[,1:3],HapIntFamily)

}#end function