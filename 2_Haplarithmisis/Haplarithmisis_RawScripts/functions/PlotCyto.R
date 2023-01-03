plotCyto <-function(chr,BandName=T,ChrsLength,ideogram){
 	
#load("/Users/u0066087/Data/refData/Ideogram_hg19.rda") 	


chromName<-paste("chr",chr,sep="")
xRange<-c(0,as.numeric(as.character(ChrsLength[ChrsLength[,1]==chromName,2])))


plot(0, axes = F, xlab = "",ylab="", type = "n", col = "gray", xlim = xRange)
	

up <-0.5
bottom <- -0.5

   cylindrect(xRange[1],(bottom/5),xRange[2],(up/5),col="#84848450",gradient="y",nslices=100)

#segments(0,0,as.numeric(ChrsLength[ChrsLength[,1]==chromName,2]),0,col="#84848440",cex=200,lwd=5) 
bandes<-ideogram[ideogram$Chromosome==chromName,]
bandes$posGrafic <- bandes$Start+(bandes$Stop-bandes$Start)/2
numBandes<-dim(bandes)[1]
gieStain<-levels(bandes$Color) 

colorsBandes<-c("red","gray100","black","gray25","gray50","gray75","gray50","gray50")
#japintat<-0

for(b in 1:numBandes){
   chromStart<-bandes$Start[b]
   chromEnd<-bandes$Stop[b]
   gBand<-bandes$Color[b]
   colorBanda<-colorsBandes[gBand]
   if(gBand == "acen"){
     #segments(chromStart,bottom,chromStart,up,col="black")
    } else{ if(gBand == "gvar"){
      rect(chromStart,bottom,chromEnd,up,col=colorBanda,angle=90,density=0.9)
      cylindrect(chromStart,bottom,chromEnd,up,col=colorBanda,gradient="y",nslices=50)
     }else{
      rect(chromStart,bottom,chromEnd,up,col=colorBanda,border="gray50")
      cylindrect(chromStart,bottom,chromEnd,up,col=colorBanda,gradient="y",nslices=50)
     } 
    }
#   if(japintat==0 & gBand == "acen"){
  #   japintat <- 1
     #points(chromEnd, (up+bottom)/2, pch=19, col = "#84848440", xpd = TRUE, cex=2)
   #}
 }#end b loop

  
# Noms de les bandes
if (BandName){
  	text(bandes$posGrafic,rep(bottom-0.2,numBandes),bandes[,4],srt=-90,adj=0,cex=1.2,family="Helvetica",xpd=TRUE)
}


}#end Function

