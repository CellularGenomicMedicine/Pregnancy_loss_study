#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      		        			          Function: phparents            	   											     %%%%%%%%%%%%%
#

#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): Trio ==> Father Mother RefSib 
# 
#
#(->): Phased parental genotypes
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phgrandparents <- function(Gtypes,Parent1,Family){

MotherChrX <- Gtypes[Gtypes$Chr=="X",paste("Mother_",Family,sep="")]	
GtypesChrX <- Gtypes[Gtypes$Chr=="X",]

PhasedPar <- Gtypes[,Parent1]

PhasedPar[(Gtypes$Grandfather=="AB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB") |
				   (Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB")]<-"BA"

PhasedPar[(Gtypes$Grandfather=="AA" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB") |
				   (Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="BB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes$Grandfather=="AA" & Gtypes$Grandmother=="BB" & Gtypes[,Parent1]=="BB") | 
				   (Gtypes$Grandfather=="AA" & Gtypes$Grandmother=="BB" & Gtypes[,Parent1]=="AA") | 
				   ((Gtypes$Grandfather=="NC" | Gtypes$Grandmother=="NC") & Gtypes[,Parent1]=="AB") | 
				   (Gtypes$Grandfather=="AB" & Gtypes$Grandmother=="AB")] <- "NC"

Gtypes[,Parent1] <- PhasedPar

#MotherX <- Gtypes[Gtypes$Chr=="X",paste("Mother_",Family,sep="")]

if(Parent1==paste("Mother_",Family,sep="")){

	MotherChrX[(GtypesChrX$Grandfather=="BB" & GtypesChrX$Mother=="AB")]<-"BA"

	MotherChrX[(GtypesChrX$Grandfather=="NC" & GtypesChrX$Mother=="AB") |
					   (GtypesChrX$Grandfather=="BB" & GtypesChrX$Grandmother=="BB" & GtypesChrX$Mother=="AB") | 

					   (GtypesChrX$Grandfather=="AA" & GtypesChrX$Grandmother=="AA" & GtypesChrX$Mother=="AB") ]<-"NC"
	Gtypes[Gtypes$Chr=="X",Parent1] <- MotherChrX
	}	   

if(sum(grep("Father",Parent1))!=0 ){Father=Gtypes[,Parent1];Mother=Gtypes$Mother}else{Father=Gtypes$Father;Mother=Gtypes[,Parent1]}

Parents <- cbind(Father,Mother)

Parents

}#end function
