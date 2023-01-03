phasebaf4_mod <- function(Gtypes,dataraw,Family,gamma,ParScore,outPath,SibPattern){ 

print("Computing phased BAF...")
 
plateau <-100
PhBAF <- vector("list",8)
names(PhBAF) <- c("P1","P2","M1","M2","P1Seg","P2Seg","M1Seg","M2Seg")
  
Father <- Gtypes[,paste("Father_",Family,sep="")]
Mother <- Gtypes[,paste("Mother_",Family,sep="")]
Ref <- Gtypes[,grep("Affected",colnames(Gtypes))]

FatherX <- Father[Gtypes$Chr=="X"]
MotherX <- Mother[Gtypes$Chr=="X"] 
RefX <- Ref[Gtypes$Chr=="X"]

if(sum(FatherX=="AB")>=1){
warning(paste("On the chr. x  of the father",sum(FatherX=="AB"),"(",sum(FatherX=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!"))
print("These will be treated as NoCalls")
FatherX[FatherX=="AB"] <- "NC"#There could not be htz SNPs on chromosome X of the father
}


sc <- colnames(Gtypes)[grep("Affected",colnames(Gtypes))]
if(is.na(ParScore[[sc]]["X","Par"])){print(paste("Sex of ",sc ,"could not be determined!!"))
		}else if(ParScore[[sc]]["X","Par"]==12){RefSex <- "female" ;print(paste(sc,"is female"))
		} else if(ParScore[[sc]]["X","Par"]==22){RefSex<- "male";print(paste(sc,"is male"))
		}else if(ParScore[[sc]]["X","Par"]==11){RefSex <- "female";print(paste(sc,"is a female but without maternal X chromosome"))
		}else if(ParScore[[sc]]["X","Par"]==0){print(paste("Chr.X nullisomy of ",sc))
		}else{print(paste("Sex of ",sc ,"could not be determined!!"))}

ParentsX <- phparentsx (FatherX,MotherX,RefX,RefSex)

Parents <- phparents(Father,Mother,Ref)
Father <- Parents[,1]
Mother <- Parents[,2]
Father[Gtypes$Chr=="X"] <- ParentsX[,"Father"]
Mother[Gtypes$Chr=="X"] <- ParentsX[,"Mother"]

BAFs <- dataraw[,c(1:3,grep(".B.Allele.",colnames(dataraw)))]
BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))] <- 1 - BAFs[Father=="BA" | Mother== "BA",grep(".B.Allele.Freq",colnames(BAFs))]
colnames(BAFs) <- gsub(".B.Allele.Freq","",colnames(BAFs))
colnames(BAFs)[4:ncol(BAFs)] <- paste(colnames(BAFs)[4:ncol(BAFs)],"_",Family,sep="")

BAFs <- BAFs[Ref!="NC",]
Father <- Father[Ref!="NC"]
Mother <- Mother[Ref!="NC"]

PhBAF[["P1"]] <-  BAFs[((Father=="AB"& Mother=="AA") | (Father=="BA"& Mother=="BB")),-c(grep("ther",colnames(BAFs)))]
PhBAF[["P2"]] <-  BAFs[((Father=="AB"& Mother=="BB") | (Father=="BA"& Mother=="AA")),-c(grep("ther",colnames(BAFs)))]

PhBAF[["M1"]] <- BAFs[((Mother=="AB"& Father=="AA") | (Mother=="BA"& Father=="BB")),-c(grep("ther",colnames(BAFs)))]
PhBAF[["M2"]] <- BAFs[((Mother=="AB"& Father=="BB") | (Mother=="BA"& Father=="AA")),-c(grep("ther",colnames(BAFs)))]

gammaSC <- gamma
gammaMC <- gamma

PhBAF[["P1Seg"]] <- callpcfBAF_mod(PhBAF[["P1"]],gammaSC,gammaMC,plateau,SibPattern)
PhBAF[["P2Seg"]] <- callpcfBAF_mod(PhBAF[["P2"]],gammaSC,gammaMC,plateau,SibPattern)
PhBAF[["M1Seg"]] <- callpcfBAF_mod(PhBAF[["M1"]],gammaSC,gammaMC,plateau,SibPattern)
PhBAF[["M2Seg"]]<- callpcfBAF_mod(PhBAF[["M2"]],gammaSC,gammaMC,plateau,SibPattern)

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

