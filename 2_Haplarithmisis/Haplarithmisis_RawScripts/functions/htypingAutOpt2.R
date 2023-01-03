#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					                         Function: htypingAutOpt2            	   									%%%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): Trio ==> Father Mother RefSib Sibs ChrPos
# 
#
#(->): haplotypes of autosomes
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
htypingAutOpt2 <- function(Gtypes,Father,Mother,Ref,Sibs,SibPattern){

print("*******************************************")
print("****     Option2 Htype analysis...     ****")

print("-------------------------------------------")
print("#1Phasing parental genotypes...")
Parents <- phparents(Father,Mother,Ref)
Father <- Parents[,1]
Mother <- Parents[,2]
print("-------------------------------------------")

print("#2 Phasing autosomes of the siblings...")
InfSNPs <- rep(0,length(Father))

#Determination of parental informative SNPs, showing patern of homozygousity in one parent and heterozygousity in the other parent 

InfSNPs[Father=="AB" & Mother=="AA"] <- 1
InfSNPs[Father=="AB" & Mother=="BB"] <- 2
InfSNPs[Father=="BA" & Mother=="AA"] <- 3
InfSNPs[Father=="BA" & Mother=="BB"] <- 4
InfSNPs[Father=="AA" & Mother=="AB"] <- 5
InfSNPs[Father=="BB" & Mother=="AB"] <- 6
InfSNPs[Father=="AA" & Mother=="BA"] <- 7
InfSNPs[Father=="BB" & Mother=="BA"] <- 8



#Determination of paternal and maternal haplotypes

	for(sib in 1:ncol(Sibs)){
	
	Sib <- Sibs[,sib]
	PatHap <- rep(0,length(Father))
	MatHap <- rep(0,length(Mother))

	PatHap[(InfSNPs==2 & Sib=="AB") | (InfSNPs==2 & Sib=="AA") | (InfSNPs==1 & Sib=="AA") |
				 (InfSNPs==3 & Sib=="AB") | (InfSNPs==3 & Sib=="BB") | (InfSNPs==4 & Sib=="BB") ] <- 1
				
	PatHap[(InfSNPs==1 & Sib=="AB") | (InfSNPs==1 & Sib=="BB") | (InfSNPs==2 & Sib=="BB") |
				 (InfSNPs==4 & Sib=="AB") | (InfSNPs==4 & Sib=="AA") | (InfSNPs==3 & Sib=="AA") ] <- 2
	
		
	MatHap[(InfSNPs==6 & Sib=="AB") | (InfSNPs==6 & Sib=="AA") | (InfSNPs==5 & Sib=="AA") |
			   	 (InfSNPs==7 & Sib=="AB") | (InfSNPs==7 & Sib=="BB") | (InfSNPs==8 & Sib=="BB") ] <- 1
				
	MatHap[(InfSNPs==5 & Sib=="AB") | (InfSNPs==5 & Sib=="BB") | (InfSNPs==6 & Sib=="BB") |
				 (InfSNPs==8 & Sib=="AB") | (InfSNPs==8 & Sib=="AA") | (InfSNPs==7 & Sib=="AA") ] <- 2
	
	print(paste("Sibling",gsub(SibPattern,"",colnames(Sibs)[sib]),"is phased"))
	if(sib==1){MatHaps <- MatHap;PatHaps <- PatHap}else{MatHaps <- cbind(MatHaps,MatHap);PatHaps <- cbind(PatHaps,PatHap)}
	
	}#end sib loop

PatHaps <- as.matrix(PatHaps)
MatHaps <- as.matrix(MatHaps)

#colnames(PatHaps)<-paste(gsub(SibPattern,"",colnames(Sibs)),"_Pat",sep="")
#colnames(MatHaps)<-paste(gsub(SibPattern,"",colnames(Sibs)),"_Mat",sep="")
colnames(PatHaps)<-paste(colnames(Sibs),"_Pat",sep="")
colnames(MatHaps)<-paste(colnames(Sibs),"_Mat",sep="")

HapsAut <- cbind(Gtypes[Gtypes$Chr!="X",1:3],PatHaps[Gtypes$Chr!="X",],MatHaps[Gtypes$Chr!="X",])
if(ncol(Sibs)==1){colnames(HapsAut)<- c(colnames(Gtypes)[1:3],colnames(PatHaps),colnames(MatHaps))}

print("-------------------------------------------")

print("*******************************************")

HapsAut

}# end function
