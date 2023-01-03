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

phparents <- function(Father,Mother,Ref){
	
print("Phasing parental genotypes...")
print("SNP-calls of the reference sibling that are in violation with Mendelian inheritance, are corrected and will be treated as NoCalls...")

	
InfSNPs <- rep(0,length(Father))

#Determination of parental informative SNPs, showing patern of homozygousity in one parent and heterozygousity in the other parent 
InfSNPs[Father=="AA" & Mother=="AB"] <- 1
InfSNPs[Father=="BB" & Mother=="AB"] <- 2

InfSNPs[Father=="AB" & Mother=="AA"] <- 3
InfSNPs[Father=="AB" & Mother=="BB"] <- 4
InfSNPs[Father=="AB" & Mother=="AB"] <- 5


#Assuming that Pat1 and Mat1 homologous chromosomes were transmitted to the reference sibling
Mother[Ref=="AB" & InfSNPs==1] <- "BA" 	
Mother[Ref=="BB" & InfSNPs==1] <-"NC"#Not possible 	
Mother[Ref=="BB" & InfSNPs==2] <- "BA"
Mother[Ref=="AA" & InfSNPs==2] <- "NC"#Not possible
Mother[Ref=="BB" & InfSNPs==5] <- "BA"
Mother[Ref=="AB" & InfSNPs==5] <- "NC"

Father[Ref=="AB" & InfSNPs==3 ] <- "BA"
Father[Ref=="BB" & InfSNPs==3 ] <- "NC"#Not possible
Father[Ref=="BB" & InfSNPs==4 ] <- "BA" 	
Father[Ref=="AA" & InfSNPs==4 ] <- "NC"#Not possible 	
Father[Ref=="BB" & InfSNPs==5 ] <- "BA" 	
Father[Ref=="AB" & InfSNPs==5 ] <- "NC"

Mother[Ref=="NC"] <- "NC" 
Father[Ref=="NC"] <- "NC" 

 	

Parents <- cbind(Father,Mother)
Parents
	
}#end function

