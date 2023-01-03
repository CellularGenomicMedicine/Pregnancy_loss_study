#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%           Function: Parent of Origin             %%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): Trio ==> Father Mother Child
# 
#
#(->): POcalls
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

po2 <- function(Gtypes,Family,SibPattern,outPath){
	
	print("Parent of origin analysis...")
	PatGtype <- Gtypes$Father
	MatGtype <- Gtypes$Mother

	for(sib in unique(c(grep(SibPattern,colnames(Gtypes)),grep("Affect",colnames(Gtypes))))){
	
	ChildGtype <- Gtypes[,sib]
	P <- matrix(0,nrow = nrow(Gtypes), ncol = 1)
	
	P[MatGtype== "AA" & PatGtype == "BB" & ChildGtype == "AA"] <- 1
	P[MatGtype== "AA" & PatGtype == "BB" & ChildGtype == "BB"] <- -1
	P[MatGtype== "BB" & PatGtype == "AA" & ChildGtype == "AA"] <- -1
	P[MatGtype== "BB" & PatGtype == "AA" & ChildGtype == "BB"] <- 1
	P[MatGtype== "AB" & PatGtype == "AA" & ChildGtype == "BB"] <- 0.5
	P[MatGtype== "AA" & PatGtype == "AB" & ChildGtype == "BB"] <- -0.5
	P[MatGtype== "BB" & PatGtype == "AB" & ChildGtype == "AA"] <- -0.5
	P[MatGtype== "AB" & PatGtype == "BB" & ChildGtype == "AA"] <- 0.5

	colnames(P)<-c(colnames(Gtypes)[sib])
	if(sib==grep(SibPattern,colnames(Gtypes))[1]){Ps <- P}else{Ps <- cbind(Ps,P)}

	print(colnames(Gtypes)[sib])	

	}#end sibloop
	
	dataPo <- cbind(Gtypes[,c("Name","Chr","Position")],Ps)
	dataOut <- dataPo
	colnames(dataOut) <- gsub("E0_Bl0_","",colnames(dataPo))
	write.table(dataOut,paste(outPath,Family,".poo",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
	print(paste("Parent of origin origin file (",Family,".poo) is written on: ",outPath,sep=""))
	dataPo
	
			
} #end PO function
