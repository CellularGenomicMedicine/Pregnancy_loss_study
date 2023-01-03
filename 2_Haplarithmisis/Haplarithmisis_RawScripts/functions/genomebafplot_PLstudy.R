genomebafplot_PLstudy <- function(BAFs,logRs,logRsSeg,PhBAF,ChrsLength,ideogram,Family,outPath){


P1    <- PhBAF[["P1"]]
P1Seg <- PhBAF[["P1Seg"]]
P2    <- PhBAF[["P2"]]
P2Seg <- PhBAF[["P2Seg"]]
M1    <- PhBAF[["M1"]]
M1Seg <- PhBAF[["M1Seg"]]
M2    <- PhBAF[["M2"]]
M2Seg <- PhBAF[["M2Seg"]]

dataHap    <- Haps[["dataHap"]]
dataHapRaw <- Haps[["dataHapRaw"]]
#logRsSeg   <- AvgLogRs

require("plotrix")
Lab <- 1.2
Main <- 2.2
Ax <- 1.2
BAFs <- BAFs[order(BAFs[,"Chr"], BAFs[,"Position"]),]
logRs <- logRs[order(logRs[,"Chr"], logRs[,"Position"]),]
logRsSeg <- logRsSeg[order(logRsSeg[,"Chr"], logRsSeg[,"Position"]),]

CumLengths <- ChrsLength
CumLengths[,"Length"] <- cumsum(as.numeric(CumLengths[,2]))
GenomeLength <- as.numeric(CumLengths[CumLengths[,"Chromosome"]=="chrX","Length"])

P1 <- na.omit(P1)
P1Seg <- na.omit(P1Seg)

P2 <- na.omit(P2)
P2Seg <- na.omit(P2Seg)

M1 <- na.omit(M1)
M1Seg <- na.omit(M1Seg)

M2 <- na.omit(M2)
M2Seg <- na.omit(M2Seg)

for(chr in CumLengths[,"Chromosome"][2:nrow(CumLengths)]){

print(chr)    
    ToAdd <- as.numeric(CumLengths[grep(paste(chr,"$",sep=""),CumLengths[,"Chromosome"])-1,"Length"])
    
    dataPo[dataPo$Chr==gsub("chr","",chr),"Position"] <- dataPo[dataPo$Chr==gsub("chr","",chr),"Position"] +ToAdd
    BAFs[BAFs$Chr==gsub("chr","",chr),"Position"] <- BAFs[BAFs$Chr==gsub("chr","",chr),"Position"] +ToAdd
    logRs[logRs$Chr==gsub("chr","",chr),"Position"] <- logRs[logRs$Chr==gsub("chr","",chr),"Position"] +ToAdd
    logRsSeg[logRsSeg$Chr==gsub("chr","",chr),"Position"] <- logRsSeg[logRsSeg$Chr==gsub("chr","",chr),"Position"] +ToAdd

    P1[P1$Chr==gsub("chr","",chr),"Position"] <- P1[P1$Chr==gsub("chr","",chr),"Position"] +ToAdd
    P1Seg[P1Seg$Chr==gsub("chr","",chr),"Position"] <- P1Seg[P1Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd
    P2[P2$Chr==gsub("chr","",chr),"Position"] <- P2[P2$Chr==gsub("chr","",chr),"Position"] +ToAdd
    P2Seg[P2Seg$Chr==gsub("chr","",chr),"Position"] <- P2Seg[P2Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd

    M1[M1$Chr==gsub("chr","",chr),"Position"] <- M1[M1$Chr==gsub("chr","",chr),"Position"] +ToAdd
    M1Seg[M1Seg$Chr==gsub("chr","",chr),"Position"] <- M1Seg[M1Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd
    M2[M2$Chr==gsub("chr","",chr),"Position"] <- M2[M2$Chr==gsub("chr","",chr),"Position"] +ToAdd
    M2Seg[M2Seg$Chr==gsub("chr","",chr),"Position"] <- M2Seg[M2Seg$Chr==gsub("chr","",chr),"Position"] +ToAdd
#print(chr)
}#end chr loop


Poses <- rep(0,23)
Poses[1] <- (as.numeric(ChrsLength[1,"Length"])/2)

for(i in 2:23){
    
        Poses[i] <- as.numeric(CumLengths[i-1,"Length"])+(as.numeric(ChrsLength[i,"Length"])/2)
}


#for(ind in colnames(BAFs)[-c(1:3,grep("ther",colnames(BAFs)))]){

Time <- gsub(" ","_",gsub(":","",gsub("-","",Sys.time())))

#jpeg(paste(outPath,Family,"_",ind,"_",Time,"_GenomeMultiProfile.jpg",sep=""),width=10000,height=6500,res=600)

jpeg_fn <- paste(outPath,Family,"_",Time,"_GenomeMultiProfile.jpg",sep="")
jpeg(jpeg_fn,width=12000,height=9000,res=600)


#print(ind)

#A= layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,4,10),matrix(5,4,10),matrix(6,4,10),matrix(7,5,10),matrix(8,4,10),matrix(9,5,10),matrix(10,4,10),matrix(11,5,10),matrix(12,2,10)))
#layout.show(A)

layout(rbind(matrix(1,2,10),matrix(2,3.5,10),matrix(3,1,10),matrix(4,2,10),matrix(5,4,10),matrix(6,4,10),matrix(7,4,10),matrix(8,5,10),matrix(9,2,10),matrix(10,4,10),matrix(11,4,10),matrix(12,4,10),matrix(13,5,10),matrix(14,4,10),matrix(15,5,10),matrix(16,4,10),matrix(17,5,10),matrix(18,2,10)))
#layout.show(B)


ChrsLength<- ChrsLength[ChrsLength[,"Chromosome"]!="chrY",]
ideogram<- ideogram[ideogram[,"Chromosome"]!="chrY",]

#1
par(mar=c(0,6,2,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,GenomeLength))
text((GenomeLength/2),0, paste("Whole-genome profile (",Family,")"),cex=Main,font=2)


#2
#source("/plotCytoGenome.R")
par(mar=c(2.5,6,1.5,1))
plotCytoGenome(BandName=T,ChrsLength,ideogram)


#3
par(mar=c(0,6,0,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,GenomeLength))

for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){
        rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430")
        }# end chr loop

text(Poses,0, c(1:22,"X"),cex=Lab,font=2)

#4
#/////
par(mar=c(0.7,6,0,1))
ind = gsub("_Pat","",colnames(dataHap)[grep("E",colnames(dataHap))])[1]
plot(dataPo$Position[dataPo[,ind]>0], dataPo[dataPo[,ind]>0,ind], "h",xaxt="n",yaxt="n",xlab="", xlim = c(0,GenomeLength),ylab="PO", col = "pink",frame=F, ylim = c(-1,1),cex.lab=Lab,cex.main=Main,cex.axis=Ax)                
points(dataPo$Position[dataPo[,ind]<0], dataPo[dataPo[,ind]<0,ind], "h", col = "blue")
points(dataPo$Position,rep(0,nrow(dataPo)), "p", col = "grey", pch = 8,cex=0.1)
abline(h=seq(-1,1,by=0.5),lty=2,col="grey")
axis(side=2,at=c(-1,0,1),labels=c("-1","0","1"),cex.axis=Ax)
if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}
#/////

C=0.4
#5
         par(mar=c(0,6,0,1))
        plot(0,ylab="BAF-EM",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
    }# end chr loop
        abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(BAFs[,"Position"],BAFs[,grep("E",colnames(BAFs))],pch=20,col="#00000050",cex=C)

#6
         par(mar=c(0,6,0,1))
        plot(0,ylab="Pat-EM",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)

    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
    }# end chr loop
                abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(P1[,"Position"],P1[,grep("E",colnames(P1))],pch=20,col="#0000ff10",cex=C)
        points(P2[,"Position"],P2[,grep("E",colnames(P2))],pch=20,col="#ff000010",cex=C)
        points(P1Seg[,"Position"],P1Seg[,grep("E",colnames(P1Seg))],pch=15,col="#0000ff80",cex=C*1.5)
        points(P2Seg[,"Position"],P2Seg[,grep("E",colnames(P1Seg))],pch=15,col="#ff000080",cex=C*1.5)

#7
         par(mar=c(0,6,0,1))
        plot(0,ylab="Mat-EM",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(M1[,"Position"],M1[,grep("E",colnames(M1))],pch=20,col="#0000ff10",cex=C)
        points(M2[,"Position"],M2[,grep("E",colnames(M2))],pch=20,col="#ff000010",cex=C)
        points(M1Seg[,"Position"],M1Seg[,grep("E",colnames(M1Seg))],pch=15,col="#0000ff80",cex=C*1.5)
        points(M2Seg[,"Position"],M2Seg[,grep("E",colnames(M2Seg))],pch=15,col="#ff000080",cex=C*1.5)

         #par(mar=c(0,6,0,1))
        #plot(0,ylab="logR-EM",col="white",main="",xlab="",frame=F,ylim=c(-3,3),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
         #axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""),cex.axis=Ax)
#8
         par(mar=c(0,6,0,1))
         plot(0,ylab="logR-EM",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
         axis(side=2,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),labels=c("","-1","","0","","1",""),cex.axis=Ax)
         abline(h=seq(-1.5,1.5,0.5),lty=2,ylim=c(0,1),col="grey")
         
    
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(-2,2,1),lty=2,ylim=c(0,1))

        points(logRs[,"Position"],logRs[,grep("E",colnames(logRs))],pch=20,col="#00000030",cex=C)
        points(logRsSeg[,"Position"],logRsSeg[,grep("E",colnames(logRsSeg))],pch=20,col="orange",cex=C*1.5)

##################################################### ExEMes - START

 #9       
        #/////
        par(mar=c(0.7,6,0,1))
        ind = gsub("_Pat","",colnames(dataHap)[grep("Affected",colnames(dataHap))])[1]
        plot(dataPo$Position[dataPo[,ind]>0], dataPo[dataPo[,ind]>0,ind], "h",xaxt="n",yaxt="n",xlab="", xlim = c(0,GenomeLength),ylab="PO", col = "pink",frame=F, ylim = c(-1,1),cex.lab=Lab,cex.main=Main,cex.axis=Ax)                
        points(dataPo$Position[dataPo[,ind]<0], dataPo[dataPo[,ind]<0,ind], "h", col = "blue")
        points(dataPo$Position,rep(0,nrow(dataPo)), "p", col = "grey", pch = 8,cex=0.1)
        abline(h=seq(-1,1,by=0.5),lty=2,col="grey")
        axis(side=2,at=c(-1,0,1),labels=c("-1","0","1"),cex.axis=Ax)
        if(sum(Int[,1]==chr)>=1){abline(v=c(Int[Int[,1]==chr,2],Int[Int[,1]==chr,3]),lty=2,lwd=4,col="#ff800040")}
        #/////        
        
        
C=0.4
#10
         par(mar=c(0,6,0,1))
        plot(0,ylab="BAF-CV",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
    }# end chr loop
        abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(BAFs[,"Position"],BAFs[,grep("Affected",colnames(P1))],pch=20,col="#00000050",cex=C)

#11
         par(mar=c(0,6,0,1))
        plot(0,ylab="Pat-CV",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)

    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
    }# end chr loop
                abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(P1[,"Position"],P1[, grep("Affected_",colnames(P1))],pch=20,col="#0000ff10",cex=C)
        points(P2[,"Position"],P2[, grep("Affected_",colnames(P2))],pch=20,col="#ff000010",cex=C)
        points(P1Seg[,"Position"],P1Seg[, grep("Affected_",colnames(P1Seg))],pch=15,col="#0000ff80",cex=C*1.5)
        points(P2Seg[,"Position"],P2Seg[, grep("Affected_",colnames(P2Seg))],pch=15,col="#ff000080",cex=C*1.5)

#12
         par(mar=c(0,6,0,1))
        plot(0,ylab="Mat-CV",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(M1[,"Position"],M1[, grep("Affected_",colnames(M1))],pch=20,col="#0000ff10",cex=C)
        points(M2[,"Position"],M2[, grep("Affected_",colnames(M2))],pch=20,col="#ff000010",cex=C)
        points(M1Seg[,"Position"],M1Seg[, grep("Affected_",colnames(M1Seg))],pch=15,col="#0000ff80",cex=C*1.5)
        points(M2Seg[,"Position"],M2Seg[, grep("Affected_",colnames(M2Seg))],pch=15,col="#ff000080",cex=C*1.5)

#13
         #par(mar=c(0,6,0,1))
         #plot(0,ylab="logR-CV",col="white",main="",xlab="",frame=F,ylim=c(-3,3),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
         #axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""),cex.axis=Ax)
  
         par(mar=c(0,6,0,1))
         plot(0,ylab="logR-CV",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
         axis(side=2,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),labels=c("","-1","","0","","1",""),cex.axis=Ax)
         abline(h=seq(-1.5,1.5,0.5),lty=2,ylim=c(0,1),col="grey")
         
         
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(-2,2,1),lty=2,ylim=c(0,1))

        points(logRs[,"Position"],logRs[, grep("Affected_",colnames(logRs))],pch=20,col="#00000030",cex=C)
        points(logRsSeg[,"Position"],logRsSeg[, grep("Affected_",colnames(logRsSeg))],pch=20,col="orange",cex=C*1.5)


#####################################################.  ExEMes - END

#14
         par(mar=c(0,6,0,1))
        plot(0,ylab="BAF-Father",col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
    
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-1,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],2,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(BAFs[,"Position"],BAFs[,paste("Father_",Family,sep="")],pch=20,col="#00000030",cex=C)
    
#15        
        #par(mar=c(0,6,0,1))
        #plot(0,ylab="Father-logR",col="white",main="",xlab="",frame=F,ylim=c(-3,3),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        #axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""),cex.axis=Ax)

        par(mar=c(0,6,0,1))
        plot(0,ylab="logR-Father",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
        axis(side=2,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),labels=c("","-1","","0","","1",""),cex.axis=Ax)
        abline(h=seq(-1.5,1.5,0.5),lty=2,ylim=c(0,1),col="grey")
        
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(-2,2,1),lty=2,ylim=c(0,1))

        points(logRs[,"Position"],logRs[,grep("Father_",colnames(logRs))],pch=20,col="#00000030",cex=C)
        points(logRsSeg[,"Position"],logRsSeg[,grep("Father_",colnames(logRs))],pch=20,col="orange",cex=C*1.5)

#16
         par(mar=c(0,6,0,1))
        plot(0,ylab="BAF-Mother",col="white",main="",xlab="Position (Gb)",frame=F,ylim=c(-0.1,1.1),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,yaxt="n",xaxt="n")
        axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
        #axis(side=1,at=seq(0,3000000000,500000000),labels=c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),cex.axis=Ax,xlab="Position (Gb)")    
        title(xlab="Position (Gb)")    
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))

        points(BAFs[,"Position"],BAFs[,paste("Mother_",Family,sep="")],pch=20,col="#00000030",cex=C)

#17
        #par(mar=c(0,6,0,1))
        #plot(0,ylab="logR-Mother",col="white",main="",xlab="",frame=F,ylim=c(-3,3),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
    #axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""),cex.axis=Ax)

    par(mar=c(0,6,0,1))
    plot(0,ylab="logR-Mother",col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim=c(0,GenomeLength),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
    axis(side=2,at=c(-1.5,-1,-0.5,0,0.5,1,1.5),labels=c("","-1","","0","","1",""),cex.axis=Ax)
    abline(h=seq(-1.5,1.5,0.5),lty=2,ylim=c(0,1),col="grey")
    
    for(chr in CumLengths[,"Chromosome"][seq(2,22,2)]){                 
    rect((as.numeric(CumLengths[CumLengths[,"Chromosome"]==chr,"Length"])-as.numeric(ChrsLength[ChrsLength[,"Chromosome"]==chr,"Length"])),-5,CumLengths[CumLengths[,"Chromosome"]==chr,"Length"],5,col="#84848430",border="#84848430") 
    }# end chr loop
            abline(h=seq(-2,2,1),lty=2,ylim=c(0,1))

        points(logRs[,"Position"],logRs[,grep("Mother_",colnames(logRs))],pch=20,col="#00000030",cex=C)
        points(logRsSeg[,"Position"],logRsSeg[,grep("Mother_",colnames(logRs))],pch=20,col="orange",cex=C*1.5)

        axis(side=1,at=seq(0,3000000000,500000000),labels=c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),cex.axis=Ax,xlab="Position (Gb)")

#18
par(mar=c(0.5,6,2,1))
plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,GenomeLength))
text((GenomeLength/2),0, paste("Position (Gb)"),cex=Ax,font=1)

dev.off()


jpeg_fn_resize <- paste(jpeg_fn,".resize.jpg",sep="")
system(paste0("/usr/bin/convert -resize 20% ",jpeg_fn, " ", jpeg_fn_resize))
system(paste0("mv ",jpeg_fn_resize, " ", jpeg_fn))

#}#end ind loop

print(paste("The plots were saved at",outPath))

}#end function

