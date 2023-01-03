chrspecbafplot <- function(dataHap,dataHapRaw,dataPo,BAFs,logRs,logRsSeg,PhBAF,ChrsLength,ideogram,Family,outPath,Chroms,Int){

#dataHap    <- Haps[["dataHap"]]
#dataHapRaw <- Haps[["dataHapRaw"]]
#logRsSeg   <- AvgLogRs

P1    <- PhBAF[["P1"]]
P1Seg <- PhBAF[["P1Seg"]]
P2    <- PhBAF[["P2"]]
P2Seg <- PhBAF[["P2Seg"]]
M1    <- PhBAF[["M1"]]
M1Seg <- PhBAF[["M1Seg"]]
M2    <- PhBAF[["M2"]]
M2Seg <- PhBAF[["M2Seg"]]

require("plotrix")
Lab <- 2.8
Main <- 6
Ax <- 2.8

for(ind in colnames(BAFs)[-c(1:3,grep("ther",colnames(BAFs)))]){

print(ind)

for(Chr in Chroms){

Time <- gsub(" ","_",gsub(":","",gsub("-","",Sys.time())))

#jpeg(paste(outPath,Family,"_",ind,"_",Time,"_Chr",Chr,".jpg",sep=""),width=20000,height=18000,res=600)



jpeg_fn <- paste(outPath,Family,"_",ind,"_",Time,"_Chr",Chr,".jpg",sep="")
jpeg(jpeg_fn,width=20000,height=18000,res=600)




#layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(6,4,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(14,1,10),matrix(15,1.5,10),matrix(6,4,10),matrix(16,1,10),matrix(17,1.5,10),matrix(7,4,10),matrix(8,5,10),matrix(9,4,10),matrix(10,5,10),matrix(11,4,10),matrix(12,5,10),matrix(13,2,10)))

ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",Chr,sep=""),2])


par(mar=c(0,6,2,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
text((ChrLength/2),0, paste("Chromosome",Chr," (",ind,")"),cex=Main,font=2)

par(mar=c(2.5,6,1.5,1))
plotCyto(Chr,BandName=T,ChrsLength,ideogram)
#plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))

head(BAFs)
print(Chr)
chr = Chr
par(mar=c(0,6,0,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))

par(mar=c(0.7,6,0,1))
plot(dataPo$Position[dataPo$Chr == Chr & dataPo[,ind]>0], dataPo[dataPo$Chr == Chr& dataPo[,ind]>0,ind], "h",xaxt="n",yaxt="n",xlab="", xlim = c(0,ChrLength),ylab="PO", col = "pink",frame=F, ylim = c(-1,1),cex.lab=Lab,cex.main=Main,cex.axis=Ax)                
points(dataPo$Position[dataPo$Chr == Chr & dataPo[,ind]<0], dataPo[dataPo$Chr == Chr& dataPo[,ind]<0,ind], "h", col = "blue")
points(dataPo$Position[dataPo$Chr == Chr],rep(0,sum(dataPo[,2] == Chr)), "p", col = "grey", pch = 8,cex=0.1)
abline(h=seq(-1,1,by=0.5),lty=2,col="grey")
axis(side=2,at=c(-1,0,1),labels=c("-1","0","1"),cex.axis=Ax)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}


C=2
		 par(mar=c(0,6,0,1))
    	plot(0,ylab="BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
    	points(BAFs[BAFs[,"Chr"]==Chr,"Position"],BAFs[BAFs[,"Chr"]==Chr,ind],pch=20,col="#00000050",cex=C)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

		 par(mar=c(0,6,0,1))
    	plot(0,ylab="Pat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    			abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
    	points(P1[P1[,"Chr"]==Chr,"Position"],P1[P1[,"Chr"]==Chr,ind],pch=20,col="#0000ff10",cex=C)
		points(P2[P2[,"Chr"]==Chr,"Position"],P2[P2[,"Chr"]==Chr,ind],pch=20,col="#ff000010",cex=C)
		points(P1Seg[P1Seg[,"Chr"]==Chr,"Position"],P1Seg[P1Seg[,"Chr"]==Chr,ind],pch=15,col="#0000ff80",cex=C)
		points(P2Seg[P2Seg[,"Chr"]==Chr,"Position"],P2Seg[P2Seg[,"Chr"]==Chr,ind],pch=15,col="#ff000080",cex=C)
	if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}
 
	par(mar=c(0,6,0,1))
    	plot(0,ylab="Mat-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")
    	points(M1[M1[,"Chr"]==Chr,"Position"],M1[M1[,"Chr"]==Chr,ind],pch=20,col="#0000ff10",cex=C)
		points(M2[M2[,"Chr"]==Chr,"Position"],M2[M2[,"Chr"]==Chr,ind],pch=20,col="#ff000010",cex=C)
		points(M1Seg[M1Seg[,"Chr"]==Chr,"Position"],M1Seg[M1Seg[,"Chr"]==Chr,ind],pch=15,col="#0000ff80",cex=C)
		points(M2Seg[M2Seg[,"Chr"]==Chr,"Position"],M2Seg[M2Seg[,"Chr"]==Chr,ind],pch=15,col="#ff000080",cex=C)
 if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}
	
	par(mar=c(0,6,0,1))
    	plot(0,ylab="logR",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5),labels=c("","-2","","-1","","0","","1","","2",""),cex.axis=Ax)
    		abline(h=seq(-2,2,0.5),lty=2,ylim=c(0,1),col="grey")

    	points(logRs[logRs[,"Chr"]==Chr,"Position"],logRs[logRs[,"Chr"]==Chr,ind],pch=20,col="#00000030",cex=C)
    	points(logRsSeg[logRsSeg[,"Chr"]==Chr,"Position"],logRsSeg[logRsSeg[,"Chr"]==Chr,ind],pch=20,col="#ff000080",cex=C)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

		 par(mar=c(0,6,0,1))
    	plot(0,ylab="Father-BAF",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1),col="grey")

    	points(BAFs[BAFs[,"Chr"]==Chr,"Position"],BAFs[BAFs[,"Chr"]==Chr,paste("Father_",Family,sep="")],pch=20,col="#00000030",cex=C)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}	
		par(mar=c(0,6,0,1))
    	plot(0,ylab="Father-logR",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5),labels=c("","-2","","-1","","0","","1","","2",""),cex.axis=Ax)
    		abline(h=seq(-2,2,0.5),lty=2,ylim=c(0,1),col="grey")

    	points(logRs[logRs[,"Chr"]==Chr,"Position"],logRs[logRs[,"Chr"]==Chr,grep("Father_",colnames(logRs))],pch=20,col="#00000030",cex=C)
    	points(logRsSeg[logRsSeg[,"Chr"]==Chr,"Position"],logRsSeg[logRsSeg[,"Chr"]==Chr,grep("Father_",colnames(logRs))],pch=20,col="#ff000080",cex=C)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}
		 par(mar=c(0,6,0,1))
    	plot(0,ylab="Mother-BAF",col="white",main="",xlab="Position (Gb)",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,yaxt="n",xaxt="n")
		axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
		title(xlab="Position (Gb)")	
    		abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

    	points(BAFs[BAFs[,"Chr"]==Chr,"Position"],BAFs[BAFs[,"Chr"]==Chr,paste("Mother_",Family,sep="")],pch=20,col="#00000030",cex=C)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}
		par(mar=c(0,6,0,1))
    	plot(0,ylab="Mother-logR",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,ChrLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
		axis(side=2,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),labels=c("","-1","","0","","1",""),cex.axis=Ax)
    		abline(h=seq(-1,1,0.5),lty=2,ylim=c(0,1))

    	points(logRs[logRs[,"Chr"]==Chr,"Position"],logRs[logRs[,"Chr"]==Chr,grep("Mother_",colnames(logRs))],pch=20,col="#00000030",cex=C)
    	points(logRsSeg[logRsSeg[,"Chr"]==Chr,"Position"],logRsSeg[logRsSeg[,"Chr"]==Chr,grep("Mother_",colnames(logRs))],pch=20,col="#ff000080",cex=C)

  	axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
	XaxPos = round(seq(0,ChrLength,round(ChrLength/4))/1000000,digits=2)
        axis(side=1,at=XaxPos*1000000,labels=as.character(XaxPos),cex.axis=Ax)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

par(mar=c(0.5,6,2,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
text((ChrLength/2),0, paste("Position (Mb)"),cex=Ax,font=1)

par(mar=c(1,6,1,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==1], col = "blue")
abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Pat",sep="")]==2], col ="cornflowerblue")#blue
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

par(mar=c(0.3,6,0.3,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")]==1], col = "blue")
abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Pat",sep="")]==2], col ="cornflowerblue")#blue
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

par(mar=c(1,6,1,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==1], col = "red")
abline(v=dataHapRaw$Position[as.character(dataHapRaw$Chr)==chr & dataHapRaw[,paste(ind,"_Mat",sep="")]==2], col ="pink")#blue
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

par(mar=c(0.3,6,0.3,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")]==1], col = "red")
abline(v=dataHap$Position[as.character(dataHap$Chr)==chr & dataHap[,paste(ind,"_Mat",sep="")]==2], col = "pink")

if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}

dev.off()

jpeg_fn_resize <- paste(jpeg_fn,".resize.jpg",sep="")
print(jpeg_fn_resize)
system(paste0("/usr/bin/convert -resize 20% ",jpeg_fn, " ", jpeg_fn_resize))
system(paste0("mv ",jpeg_fn_resize, " ", jpeg_fn))




}#end chr loop


}#end ind loop

print(paste("The plots were saved at",outPath))

}#end function
