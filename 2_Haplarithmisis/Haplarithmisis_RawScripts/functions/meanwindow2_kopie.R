meanwindow2 <- function(logRsRaw,GC,window,Func,Family,ParScore,outPath,SibPattern){
  
Rs<- logRsRaw[,-c(1:3)]

MedRs3 <- NULL
for(ind in colnames(Rs)){  
 MedR3 <- NULL
 print(paste(ind,"is being processed"))
 for(chr in unique(logRsRaw$Chr)){
  print(paste("chr",chr,"is being processed"))
  RsChr3    <- Rs[logRsRaw$Chr==chr,ind]
  MedRChr3  <- slideFunct(RsChr3, window, 1)
  MedRChr3  <- c(rep(mean(as.numeric(na.omit(RsChr3[1:window]))),window),MedRChr3[,1],rep(mean(as.numeric(na.omit(RsChr3[(length(RsChr3)-window):length(RsChr3)]))),window))
  MedR3     <- c(MedR3,MedRChr3)
 }#end chr loop
 if(ind==colnames(Rs)[1]){ MedRs3<- data.frame(MedR3,stringsAsFactors=FALSE) } else{MedRs3<- data.frame(MedRs3,MedR3,stringsAsFactors=FALSE) }
 print(paste(ind,"is normalized"))
}#end ind loop


logRsMed3 <- data.frame(logRsRaw[,c(1:3)],MedRs3,stringsAsFactors=FALSE)
colnames(logRsMed3)<- colnames(logRsRaw)


#ParScore<-load(paste(outPath,"ParScore_",Family,".rda",sep=""))
#ParScore<-eval(parse(text=ParScore))


#--------------------------------------
#---------- window on GC

for(chr in unique(logRsRaw$Chr)){
  print(paste("Chromosome",chr))
  MedGcChr3 <- NULL
  GcChr3 <- GC[logRsRaw$Chr==chr,6]
  MedGcChr3  <- slideFunct(GcChr3, window, 1)
  MedGcChr3  <- c(rep(mean(as.numeric(na.omit(GcChr3[1:window]))),window),MedGcChr3[,1],rep(mean(as.numeric(na.omit(GcChr3[(length(GcChr3)-window):length(GcChr3)]))),window))

  if(chr==unique(logRsRaw$Chr)[1]){MedGc3 <- MedGcChr3}else{MedGc3 <- c(MedGc3,MedGcChr3)}
  

}#end chr loop


GC3 <- data.frame(logRsRaw[,c(1:3)],MedGc3,stringsAsFactors=FALSE)
colnames(GC3)[4]<-"GC"

rowsTot3 <- intersect(rownames(logRsMed3),rownames(GC3))

logRsMed4 <- logRsMed3[rowsTot3,]

GC3        <- data.frame(GC3[rowsTot3,],stringsAsFactors=FALSE)
GC3[,"GC"] <- as.numeric(GC3[,"GC"])

logRs <- wavecorrtmean2.1(logRsMed4,ParScore,GC3,Family,outPath,SibPattern)
logRs

}#end function

