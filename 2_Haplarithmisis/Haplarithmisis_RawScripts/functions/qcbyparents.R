#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					                         Function: qcbyparents            	   									%%%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): Family genotypes
# 
#
#(->): ADO and ADI based on parental genotypes
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qcbyparents <- function(Gtypes,SibPattern){

print("QC by parents analysis...")	

#chroms <- c(1:22,"X","Y")
chroms <- unique(Gtypes$Chr)
Sibs <- c(grep(SibPattern,colnames(Gtypes)),grep("Affected",colnames(Gtypes)))
Sibs <- unique(Sibs)
ColNames<-gsub("E0_Bl0_","",colnames(Gtypes))
 
QC_par <- vector("list",2)
names(QC_par)<-c("ADO","ADI")

for(sib in Sibs){

	#Genome-wide QC by parents
	Father <- Gtypes$Father[Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y"]
	Mother <- Gtypes$Mother[Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y"]
	Sib <- Gtypes[Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y",sib]
	
	ADO_AutGenome <- (sum((Father=="AB" & Mother=="BB" & Sib=="AA") |
						   				   (Father=="AB" & Mother=="AA" & Sib=="BB") | 
						    			   (Father=="BB" & Mother=="AB" & Sib=="AA") |
						    			   (Father=="AA" & Mother=="AB" & Sib=="BB") |
						    			   (Father=="AA" & Mother=="BB" & (Sib=="AA" |  Sib=="BB")) |
						    			   (Father=="BB" & Mother=="AA" & (Sib=="AA" |  Sib=="BB")))/sum((Sib=="AA" | Sib=="BB") & (Father!="NC" | Mother!="NC")))*100
	
	ADI_AutGenome <- (sum((Father=="AA" & Mother=="AA" & Sib=="AB") |
   						 				  (Father=="BB" & Mother=="BB" & Sib=="AB") )/sum(Sib=="AB" & (Father!="NC" | Mother!="NC")))*100

	
	#Chr-specific QC by parents
	for(chr in chroms){
	
		FatherChr <- Gtypes$Father[Gtypes$Chr==chr]
		MotherChr <- Gtypes$Mother[Gtypes$Chr==chr]
		SibChr <- Gtypes[Gtypes$Chr==chr,sib]
		
		ADO <- (sum((FatherChr=="AB" & MotherChr=="BB" & SibChr=="AA") |
						    (FatherChr=="AB" & MotherChr=="AA" & SibChr=="BB") | 
						    (FatherChr=="BB" & MotherChr=="AB" & SibChr=="AA") |
						    (FatherChr=="AA" & MotherChr=="AB" & SibChr=="BB") |
						    (FatherChr=="AA" & MotherChr=="BB" & (SibChr=="AA" |  SibChr=="BB")) |
						    (FatherChr=="BB" & MotherChr=="AA" & (SibChr=="AA" |  SibChr=="BB")))/sum((Sib=="AA" | Sib=="BB") & (Father!="NC" | Mother!="NC")))*100
	
		ADI <- (sum((FatherChr=="AA" & MotherChr=="AA" & SibChr=="AB") |
   						   (FatherChr=="BB" & MotherChr=="BB" & SibChr=="AB") )/sum(Sib=="AB" & (Father!="NC" | Mother!="NC")))*100
						  
	if(chr==chroms[1]){ADOGenome <- ADO; ADIGenome <- ADI }else{ADOGenome <- cbind(ADOGenome,ADO); ADIGenome <- cbind(ADIGenome,ADI) }	

	}#end chr loop
	
	ADOGenome <- cbind(ADOGenome,ADO_AutGenome)
	ADIGenome <- cbind(ADIGenome,ADI_AutGenome)

	if(sib==Sibs[1]){ADOs <- ADOGenome; ADIs <- ADIGenome}else{ADOs <- rbind(ADOs,ADOGenome); ADIs <- rbind(ADIs,ADIGenome)}

	print(ColNames[sib])
}#end sib loop

	
colnames(ADOs) <- c(chroms,"GenomeAut")
colnames(ADIs) <- c(chroms,"GenomeAut")
rownames(ADOs) <- ColNames[Sibs]
rownames(ADIs) <- ColNames[Sibs]
	
QC_par[["ADO"]] <- ADOs
QC_par[["ADI"]] <- ADIs 

QC_par	

}#end function 


