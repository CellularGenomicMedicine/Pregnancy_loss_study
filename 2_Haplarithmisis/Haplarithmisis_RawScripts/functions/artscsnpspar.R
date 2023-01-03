artscsnpspar <- function(Gtypes,SibPattern,outPath){

	Father <- Gtypes$Father
	Mother <- Gtypes$Mother	
	
	Artfs <- matrix(NA,nrow(Gtypes),ncol(Gtypes)-3)

	Artfs <- cbind(Gtypes[,1:3],Artfs)
	colnames(Artfs) <- colnames(Gtypes)

		
	for(sc in colnames(Gtypes)[c(grep(SibPattern,colnames(Gtypes)),grep("Affected",colnames(Gtypes)))]){
			
				print(paste("Sample ",sc, "is processing..." ))
				GtypeSC<- Gtypes[,sc]
				
									Artfs[((GtypeSC=="AA" | GtypeSC=="BB") & (Father =="AB" & Mother =="AB") ) | 
												(GtypeSC=="AA" & (Father =="AA" & Mother =="AA") ) | 
												(GtypeSC=="BB" & (Father =="BB" & Mother =="BB") ) |
												(GtypeSC=="AA" & ((Father =="AB" & Mother =="AA") | (Father =="AA" & Mother =="AB"))) |
												(GtypeSC=="BB" & ((Father =="AB" & Mother =="BB") | (Father =="BB" & Mother =="AB"))) |
												(GtypeSC=="AB" & (Father=="AA" & (Mother =="BB" | Mother == "AB" ))) |
												(GtypeSC=="AB" & (Father=="BB" & (Mother =="AA" | Mother == "AB" ))) |
												(GtypeSC=="AB" & (Mother=="AA" & (Father =="BB" | Father == "AB" ))) |
												(GtypeSC=="AB" & (Mother=="BB" & (Father =="AA" | Father == "AB" ))) |
												(GtypeSC=="AB" & ((Father=="AA" | Father=="BB" | Father=="AB") & Mother =="AB" )) |
										  		(GtypeSC=="AB" & (Father=="AB" & (Mother =="AA" | Mother =="BB" | Mother=="AB"))),sc] <- 0 #Consistent			
							Artfs[((GtypeSC=="AA" | GtypeSC=="BB") & ((Father =="AA" & Mother =="BB") | (Father =="BB" & Mother =="AA"))) | 
						  			(GtypeSC=="AA" & ((Father =="AB" & Mother =="BB") | (Father =="BB" & Mother =="AB"))) |
								    (GtypeSC=="BB" & ((Father =="AB" & Mother =="AA") | (Father =="AA" & Mother =="AB"))),sc] <- 1 #ADO
						  
							Artfs[(GtypeSC=="AB" & (Father =="BB" & Mother =="BB") ) | 
							    	 (GtypeSC=="AB" & (Father =="AA" & Mother =="AA")),sc] <- 2 #ADI
							Artfs[(GtypeSC=="AA" & (Father =="BB" & Mother =="BB")) |
									 (GtypeSC=="BB" & (Father =="AA" & Mother =="AA")),sc] <- 3 #ADOADI
			Artfs[GtypeSC=="NC",sc] <- 4 # NoCall in the single cell

			Artfs[GtypeSC!="NC" & (Father=="NC" | Mother=="NC"),sc]  <- 5 #NoCall in the parents			
	}#end sc loop
	

	#write.table(ADOHtzMC,paste(outPath,PlatAlg,".rts",sep=""),quote=F,col.names=T,row.names=F)
	Artfs <- Artfs[,-c(grep("ther",colnames(Artfs)))]
	write.table(Artfs,paste(outPath,"ArtefactsByParents_illuminaRapid.txt",sep=""),sep="\t",quote=F,col.names=T,row.names=F)
	Artfs

}#end function


