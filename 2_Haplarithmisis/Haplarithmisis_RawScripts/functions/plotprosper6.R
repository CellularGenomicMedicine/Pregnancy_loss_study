#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      					                         Function: plotall            	   									%%%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): outPath chroms logRs matHaps PatHaps ChrsLength AvgLogRs SegLogRs
# 
#
#(->): pdf files and possibly converted tiff files
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotprosper6 <- function(outPath,chroms,logRs,AvgLogRs,SegLogRs,Intp,dataPo,QC,SnpSpecArts,REF_QC,REF_24h_QC,tiffgen){

dataHap <- Intp[["dataHap"]]

logRs_RefAll[,"Chr"] <- as.character(logRs_RefAll[,"Chr"])
logRs_RefAll <- logRs_RefAll[,-c(grep("ther",colnames(logRs_RefAll)),grep("_MC",colnames(logRs_RefAll)))]
PlotFolder <- "Plots_All"
if (file.exists(paste(outPath,PlotFolder,sep=""))){
    outPathPlot <- file.path(outPath, PlotFolder)
	} else {
    dir.create(file.path(outPath, PlotFolder))
    outPathPlot <- file.path(outPath, PlotFolder)}

colnames(dataPo)<-gsub("E0_Bl0_","",colnames(dataPo))		
show(paste("The plot files will be stored at",outPathPlot))

Lab <- 2.4
Main <- 4
Ax <- 2.5
Conf <- 0.95
Alpha <- 0.05

if(length(chroms)>=22){mode="Genome"}else if(length(chroms)==1){mode=paste("Chr",chroms,"Only",sep="")}else{mode="SelChrs"}
show(paste("mode ====>",mode))

logRs <- logRs[,-c(grep("ther",colnames(logRs)))]

for(ind in 4:ncol(logRs)){
	
	print(ind)	
	print(colnames(logRs)[ind])
	dataPCA <- cbind(logRs_RefAll[,-c(1:3)],logRs[,colnames(logRs)[ind]])
	colnames(dataPCA)[ncol(dataPCA)] <- colnames(logRs)[ind]
	if(sum(grep("Affected",colnames(logRs)[ind]))!=1){
	RelScorePat <-  Intp[["DataCompl"]][grep(paste(colnames(logRs)[ind],"_Pat",sep=""),rownames(Intp[["DataCompl"]])),]
	RelScoreMat <-  Intp[["DataCompl"]][grep(paste(colnames(logRs)[ind],"_Mat",sep=""),rownames(Intp[["DataCompl"]])),]
	}
	#print(RelScorePat)
	#print(RelScoreMat)
	#pdf(paste(outPath,"Plots/",Family,"_",colnames(logRs)[ind],"_Genome2.pdf",sep=""),w=20,h=18,bg="transparent",family="Helvetica",title ="MZE plot")
	pdf(paste(outPathPlot,"/",Family,"_",colnames(logRs)[ind],"_",mode,".pdf",sep=""),w=25,h=20,bg="transparent",family="Helvetica",title ="MZE plot")

	
	for(chr in chroms){
		
		print(chr)
		ArtsChr <- SnpSpecArts[SnpSpecArts[,"Chr"]==chr,]
		dataPCAChr <-  dataPCA[as.character(logRs$Chr)==chr,]
		if(sum(grep("Affected",colnames(logRs)[ind]))!=1){
		RelScorePatChr <- RelScorePat[RelScorePat[,"Chr"]==chr,]
		RelScoreMatChr <- RelScoreMat[RelScoreMat[,"Chr"]==chr,]
		}
		
		layout(rbind(matrix(10,1,10),matrix(1,3,10),matrix(3,3,10),matrix(4,3,10),matrix(5,2,10),matrix(6,4,10),matrix(2,8,10),matrix(c(0,7,7,0,8,8,0,9,9,0),7,10,byrow=T)))
	
		ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",chr,sep=""),2])
		LogRChr <- logRs[as.character(logRs$Chr)==chr,]
		while(sum(is.na(LogRChr))>=1){LogRChr[which(is.na(LogRChr))]<-LogRChr[which(is.na(LogRChr))-1]}
		AvgLogRChr <-  AvgLogRs[as.character(AvgLogRs$Chr)==chr,]
		SegmentsChr <-  SegLogRs[as.character(SegLogRs$Chr)==chr,]

		print(nrow(RelScorePatChr) )
		if(chr=="X"){Chr=23}else{Chr <- as.numeric(chr)}
		
		par(mar=c(6,6,1,2))
		plotCyto(chr,BandName=T,ChrsLength,ideogram)
		
		par(mar=c(6,6,3,2))
		
		if(sum(grep("Affected",colnames(logRs)[ind]))==1){AxLim = "c(-2,2)" }else{AxLim = "c(-4,4)"}
		plot(AvgLogRChr$Position,LogRChr[,colnames(logRs)[ind]],col="black",pch=19,cex=0.1, xlab="",ylab=c("logR"),main="",cex.main=1,frame=F,ylim=eval(parse(text=paste(AxLim))),xlim=c(0,ChrLength),cex.lab=Lab,cex.main=Main,cex.axis=Ax)
		abline(h=seq(-2,2,by=1),lty=2,col="grey")
		points(AvgLogRChr$Position,AvgLogRChr[,colnames(logRs[ind])],col="yellow",pch=19,cex=1)
		points(SegmentsChr$Position,SegmentsChr[,colnames(logRs[ind])],col="orange",pch=19,cex=1)
		lines(lowess(LogRChr[,colnames(logRs[ind])] ~ AvgLogRChr$Position,f=1/50),lwd=2,col="red")
		legend((max(logRs[as.character(logRs[,"Chr"])==chr,3])/20)*19,4,c("Avg","PCF","Lowess"),col=c("yellow","orange","red"),border="white",pch=15,bty="n")

	 	if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="red")}
 		
 		par(mar=c(2,6,6,2))	
 		plot(0, xlab = "",ylab="Pat", type = "n", col = "white", xlim = c(0,ChrLength),yaxt="n",xaxt="n",cex.lab=Lab,cex.main=Main,cex.axis=Ax,frame=F,,ylim=c(0,1))
 		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(colnames(logRs)[ind],"_Pat",sep="")]==1], col = "blue")
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(colnames(logRs)[ind],"_Pat",sep="")]==2], col ="cornflowerblue")#blue
 		if(sum(grep("Affected",colnames(logRs)[ind]))!=1 ){
 		if( nrow(RelScorePatChr)!=0){
 			print((RelScorePatChr[i,"Start"]+RelScorePatChr[i,"Stop"])/2)
 		for(b in 1:nrow(RelScorePatChr)){text((RelScorePatChr[b,"Start"]+RelScorePatChr[b,"Stop"])/2,0.5, as.character(round(RelScorePatChr[b,"Prob"],digits=2)),cex=Ax,font=2,pos=3)}
		} 		
 		}
 		
 		par(mar=c(6,6,2,2))
		plot(0, xlab = "",ylab="Mat", type = "n", col = "white", xlim = c(0,ChrLength),yaxt="n",cex.lab=Lab,cex.main=Main,cex.axis=Ax,frame=F,ylim=c(0,1))
 		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(colnames(logRs)[ind],"_Mat",sep="")]==1], col = "red")
		abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(colnames(logRs)[ind],"_Mat",sep="")]==2], col ="pink")
		if(sum(grep("Affected",colnames(logRs)[ind]))!=1 ){
			if(nrow(RelScoreMatChr)!=0){
		 for(b in 1:nrow(RelScoreMatChr)){text((RelScoreMatChr[b,"Start"]+RelScoreMatChr[b,"Stop"])/2,0.5, as.character(round(RelScoreMatChr[b,"Prob"],digits=2)),cex=Ax,font=2,pos=3)}
		 }
		 }
		
		par(mar=c(2,6,4,2))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength),,cex.lab=Lab,cex.main=Main,cex.axis=Ax)
 		abline(v=ArtsChr[ArtsChr[,ind]==4,"Position"], col = "black")
		abline(v=ArtsChr[ArtsChr[,ind]==5,"Position"], col = "grey")
 		abline(v=ArtsChr[ArtsChr[,ind]==0,"Position"], col = "green")	
		abline(v=ArtsChr[ArtsChr[,ind]==3,"Position"], col = "red")			
		abline(v=ArtsChr[ArtsChr[,ind]==1,"Position"], col = "orange")
		abline(v=ArtsChr[ArtsChr[,ind]==2,"Position"], col = "purple")




		par(mar=c(4,6,5,2))
		plot(dataPo$Position[as.character(dataPo$Chr) == chr & dataPo[colnames(logRs)[ind]]>0], dataPo[as.character(dataPo$Chr) == chr& dataPo[colnames(logRs)[ind]]>0,colnames(logRs)[ind]], "h", xlab="", xlim = c(0,ChrLength),ylab="PO", col = "pink",frame=F, ylim = c(-1,1),cex.lab=Lab,cex.main=Main,cex.axis=Ax)
		points(dataPo$Position[as.character(dataPo$Chr) == chr & dataPo[colnames(logRs)[ind]]<0], dataPo[as.character(dataPo$Chr) == chr& dataPo[colnames(logRs)[ind]]<0,colnames(logRs)[ind]], "h", col = "blue")
		points(dataPo$Position[as.character(dataPo$Chr) == chr],rep(0,sum(dataPo[,2] == chr)), "p", col = "grey", pch = 8)
		abline(h=seq(-1,1,by=0.5),lty=2,col="grey")
	#legend((max(logRs[as.character(logRs[,2])==chr,3])/20)*19,4,c("Mat","Pat"),col=c("pink","blue"),border="white",pch=15,bty="n")

			
		Chr <- as.numeric(chr)
	#QC Plots
		par(mar=c(6,6,6,6))
		plot(REF_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_QC[["MendIncRateChrsInd"]])),25] ,
			  REF_QC[["CallRateChrsHetInd"]][grep("Pic",rownames(REF_QC[["CallRateChrsHetInd"]])),25], 
			  xlim=c(0,60),ylim=c(0,40),col="lightpink 1",pch=10,ylab="Call rate (%)",xlab="Mendelian inconsistencies (%)",main="Htz SNP-calls",
			  cex=2,cex.lab=Lab,cex.main=3,cex.axis=Ax,frame=F,yaxt="n",xaxt="n")
		points(REF_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_QC[["MendIncRateChrsInd"]])),25] ,
			 	  REF_QC[["CallRateChrsHetInd"]][grep("MDA",rownames(REF_QC[["CallRateChrsHetInd"]])),25] ,col="lightblue 2",pch=10,cex=2)
		points(REF_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_QC[["MendIncRateChrsInd"]])),25] ,
			      REF_QC[["CallRateChrsHetInd"]][grep("MC",rownames(REF_QC[["CallRateChrsHetInd"]])),25] ,col="darkseagreen 2",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),25] ,
			  REF_24h_QC[["CallRateChrsHetInd"]][grep("Pic",rownames(REF_24h_QC[["CallRateChrsHetInd"]])),25] ,col="lightpink 2",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),25] ,
			 	  REF_24h_QC[["CallRateChrsHetInd"]][grep("MDA",rownames(REF_24h_QC[["CallRateChrsHetInd"]])),25] ,col="lightblue 3",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),25] ,
			      REF_24h_QC[["CallRateChrsHetInd"]][grep("MC",rownames(REF_24h_QC[["CallRateChrsHetInd"]])),25] ,col="darkseagreen 3",pch=10,cex=2)
			      
			      
		points(REF_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_QC[["MendIncRateChrsInd"]])),Chr] ,
			  	  REF_QC[["CallRateChrsHetInd"]][grep("Pic",rownames(REF_QC[["CallRateChrsHetInd"]])),Chr],col="red",pch=10,cex=2)
		points(REF_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_QC[["MendIncRateChrsInd"]])),Chr] ,
			 	  REF_QC[["CallRateChrsHetInd"]][grep("MDA",rownames(REF_QC[["CallRateChrsHetInd"]])),Chr] ,col="blue",pch=10,cex=2)
		points(REF_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_QC[["MendIncRateChrsInd"]])),Chr] ,
			      REF_QC[["CallRateChrsHetInd"]][grep("MC",rownames(REF_QC[["CallRateChrsHetInd"]])),Chr] ,col="green",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),Chr] ,
			  REF_24h_QC[["CallRateChrsHetInd"]][grep("Pic",rownames(REF_24h_QC[["CallRateChrsHetInd"]])),Chr] ,col="red3",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),Chr] ,
			 	  REF_24h_QC[["CallRateChrsHetInd"]][grep("MDA",rownames(REF_24h_QC[["CallRateChrsHetInd"]])),Chr] ,col="blue3",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),Chr] ,
			      REF_24h_QC[["CallRateChrsHetInd"]][grep("MC",rownames(REF_24h_QC[["CallRateChrsHetInd"]])),Chr] ,col="green3",pch=10,cex=2)
		
		points(QC[["MendIncRateChrsInd"]][colnames(logRs)[ind],Chr],
		   		  QC[["CallRateChrsHetInd"]][colnames(logRs)[ind],Chr] , col="black",pch=10,cex=6)
		axis(1, tck = 1, col = "#BABABA50", lty = 2,cex.axis=Ax)
 		axis(2, tck = 1, col = "#BABABA50", lty = 2,cex.axis=Ax)

		par(mar=c(6,6,6,6))
		plot(REF_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_QC[["MendIncRateChrsInd"]])),25] ,
			  REF_QC[["CallRateChrsHomInd"]][grep("Pic",rownames(REF_QC[["CallRateChrsHomInd"]])),25] ,col="lightpink 1",
			  xlim=c(0,60),ylim=c(0,100),pch=10,ylab="Call rate (%)",xlab="Mendelian inconsistencies (%)",main="Hmz SNP-calls",
			  cex=2,cex.lab=Lab,cex.main=3,cex.axis=Ax,frame=F,yaxt="n",xaxt="n")
		points(REF_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_QC[["MendIncRateChrsInd"]])),25] ,
			 	  REF_QC[["CallRateChrsHomInd"]][grep("MDA",rownames(REF_QC[["CallRateChrsHomInd"]])),25] ,col="lightblue 2",pch=10,cex=2)
		points(REF_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_QC[["MendIncRateChrsInd"]])),25] ,
			  REF_QC[["CallRateChrsHomInd"]][grep("MC",rownames(REF_QC[["CallRateChrsHomInd"]])),25] ,col="darkseagreen 2",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),25] ,
			  REF_24h_QC[["CallRateChrsHomInd"]][grep("Pic",rownames(REF_24h_QC[["CallRateChrsHomInd"]])),25] ,col="lightpink 3",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),25] ,
			 	  REF_24h_QC[["CallRateChrsHomInd"]][grep("MDA",rownames(REF_24h_QC[["CallRateChrsHomInd"]])),25] ,col="lightblue 3",pch=10,cex=3)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),25] ,
			  REF_24h_QC[["CallRateChrsHomInd"]][grep("MC",rownames(REF_24h_QC[["CallRateChrsHomInd"]])),25] ,col="darkseagreen 3",pch=10,cex=2)
		
		points(REF_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_QC[["MendIncRateChrsInd"]])),Chr] ,
			  REF_QC[["CallRateChrsHomInd"]][grep("Pic",rownames(REF_QC[["CallRateChrsHomInd"]])),Chr] ,col="red",pch=10,cex=2)
		points(REF_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_QC[["MendIncRateChrsInd"]])),Chr] ,
			 	  REF_QC[["CallRateChrsHomInd"]][grep("MDA",rownames(REF_QC[["CallRateChrsHomInd"]])),Chr] ,col="blue",pch=10,cex=2)
		points(REF_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_QC[["MendIncRateChrsInd"]])),Chr] ,
			  REF_QC[["CallRateChrsHomInd"]][grep("MC",rownames(REF_QC[["CallRateChrsHomInd"]])),Chr] ,col="green",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("Pic",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),Chr] ,
			  REF_24h_QC[["CallRateChrsHomInd"]][grep("Pic",rownames(REF_24h_QC[["CallRateChrsHomInd"]])),Chr] ,col="red3",pch=10,cex=2)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MDA",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),Chr] ,
			 	  REF_24h_QC[["CallRateChrsHomInd"]][grep("MDA",rownames(REF_24h_QC[["CallRateChrsHomInd"]])),Chr] ,col="blue3",pch=10,cex=3)
		points(REF_24h_QC[["MendIncRateChrsInd"]][grep("MC",rownames(REF_24h_QC[["MendIncRateChrsInd"]])),Chr] ,
			  REF_24h_QC[["CallRateChrsHomInd"]][grep("MC",rownames(REF_24h_QC[["CallRateChrsHomInd"]])),Chr] ,col="green3",pch=10,cex=2)
				  
		points(QC[["MendIncRateChrsInd"]][colnames(logRs)[ind],Chr],
		          QC[["CallRateChrsHomInd"]][colnames(logRs)[ind],Chr] ,col="black",pch=10,cex=6)
		          
		axis(1, tck = 1, col = "#BABABA50", lty = 2,cex.axis=Ax)
 		axis(2, tck = 1, col = "#BABABA50", lty = 2,cex.axis=Ax)

		par(mar=c(6,6,6,6))
		if(Family!="REF"){
		PCA_Chr <- pcachr(dataPCAChr,Conf,Alpha,colnames(logRs)[ind])
		}else{plot(0, axes = FALSE, xlab = "",ylab="", type = "n", col = "white")}
		par(mar=c(1,6,1,1))
		plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
		text((ChrLength/2),0, paste("Chromosome",chr,"(",colnames(logRs)[ind],")"),cex=Main,font=2)

		}#end chr loop



show(colnames(logRs)[ind])

dev.off()

	}#end ind loop

show("pdf files were generated")

if(tiffgen==T){
	
	show("Converting pdf files to tiff files...")
	pdfFiles <- list.files(outPathPlot,pattern=paste("\\",mode,".pdf",sep=""))
	for(file in pdfFiles){
		try(system(paste("convert ",outPathPlot,"/",file," ",outPathPlot,"/",file,".tiff",sep="")))
		show(paste("File",file, "is converted" ))
		}#end file loop
	
	show("tiff files were generated")

	}else{show("The user didn't ask for tiff files ==> notice that pdf files are press quality files and sometimes huge")}
	
}#end function
