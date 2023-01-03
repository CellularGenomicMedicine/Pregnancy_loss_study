phparentsx <- function(Father,Mother,Ref,RefSex){
#print("Phasing parental genotypes...")
InfSNPs <- rep(0,length(Father))

if(RefSex=="male"){

	Mother[Ref=="BB" & Mother=="AB"] <- "BA" 	
}else if(RefSex=="female"){

	Mother[Ref=="BB" & Father=="BB" & Mother=="AB"] <- "BA" 	
	Mother[Ref=="AB" & Father=="AA" & Mother=="AB"] <- "BA" 	

}else{
	print("Sex of ref. individual is not detemined and chr. X htype is not reliable")
}

Mother[Ref=="NC" & Mother=="AB"] <- "NC" 	

Parents <- cbind(Father,Mother)
Parents
}#end function
