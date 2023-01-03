phasebaf2 <- function(Gtypes,dataraw,Family,Parent1,gamma,outPath,SibPattern){ 

print("Computing phased BAF...")

if(sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB")>=1){
warning(paste("On the ChrX  of the father",sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB"),"(",sum(Gtypes[Gtypes[,"Chr"]=="X", paste("Father_",Family,sep="")]=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!","These are treated as NoCalls"))

Gtypes[Gtypes[,"Chr"]=="X" & Gtypes[,paste("Father_",Family,sep="")]=="AB",paste("Father_",Family,sep="")] <- "NC"#There could not be htz SNPs on chromosome X of the father
}
  
plateau <-100
PhBAF <- vector("list",8)
names(PhBAF) <- c("P1","P2","M1","M2","P1Seg","P2Seg","M1Seg","M2Seg")
  
Parents <- phgrandparents(Gtypes,Parent1,Family)
Father <- Parents[,1]
Mother <- Parents[,2]

BAFs <- dataraw[,c(1:3,grep(".B.Allele.",colnames(dataraw)))]
BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))] <- 1 - BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))]
colnames(BAFs) <- gsub(".B.Allele.Freq","",colnames(BAFs))
colnames(BAFs)[4:ncol(BAFs)] <- paste(colnames(BAFs)[4:ncol(BAFs)],"_",Family,sep="")

ToRep <- BAFs 
ToRep[,4:ncol(BAFs)] <- matrix(-1,ncol(BAFs)-3,nrow(BAFs)) 

BAFs <- BAFs[Gtypes$Grandfather!="NC" | Gtypes$Grandmother!="NC",]
Father <- Father[Gtypes$Grandfather!="NC" | Gtypes$Grandmother!="NC"]
Mother <- Mother[Gtypes$Grandfather!="NC" | Gtypes$Grandmother!="NC"]


PhBAF[["P1"]] <-  BAFs[((Father=="AB"& Mother=="AA") | (Father=="BA"& Mother=="BB")),]
PhBAF[["P2"]] <-  BAFs[((Father=="AB"& Mother=="BB") | (Father=="BA"& Mother=="AA")),]

PhBAF[["M1"]] <- BAFs[((Mother=="AB"& Father=="AA") | (Mother=="BA"& Father=="BB")),]
PhBAF[["M2"]] <- BAFs[((Mother=="AB"& Father=="BB") | (Mother=="BA"& Father=="AA")),]

if(sum(grep("Father",Parent1))!=0){PhBAF[["M1"]]<-ToRep;PhBAF[["M2"]]<-ToRep;print("Parent1 is the Father")}else{PhBAF[["P1"]]<-ToRep;PhBAF[["P2"]]<-ToRep;print("Parent1 is the Mother")}

gammaSC <- gamma
gammaMC <- gamma

PhBAF[["P1Seg"]] <- callpcfBAF(PhBAF[["P1"]],gammaSC,gammaMC,plateau,SibPattern)
PhBAF[["P2Seg"]] <- callpcfBAF(PhBAF[["P2"]],gammaSC,gammaMC,plateau,SibPattern)
PhBAF[["M1Seg"]] <- callpcfBAF(PhBAF[["M1"]],gammaSC,gammaMC,plateau,SibPattern)
PhBAF[["M2Seg"]]<- callpcfBAF(PhBAF[["M2"]],gammaSC,gammaMC,plateau,SibPattern)

write.table(PhBAF[["P1"]],paste(outPath,"P1_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
write.table(PhBAF[["P1Seg"]],paste(outPath,"P1Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
write.table(PhBAF[["P2"]],paste(outPath,"P2_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
write.table(PhBAF[["P2Seg"]],paste(outPath,"P2Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

write.table(PhBAF[["M1"]],paste(outPath,"M1_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
write.table(PhBAF[["M1Seg"]],paste(outPath,"M1Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
write.table(PhBAF[["M2"]],paste(outPath,"M2_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
write.table(PhBAF[["M2Seg"]],paste(outPath,"M2Seg_",Family,".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)

PhBAF

}#end function
