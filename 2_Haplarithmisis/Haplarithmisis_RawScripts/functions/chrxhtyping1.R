#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%   	      nction: chrxHtyping            	   									%%%%%%%%%%%%%%
#
#Author: MZE (RGLab)
#cDate: 
#mDtae: 
#
#(<-): Gtypes dataPo
# 
#
#(->): X-chromosome haplotype
#
#Description: 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chrxhtyping1 <- function(Gtypes,ParScore,SibPattern){


print("Haplotype reconstrucion of the sex chromosome...")
print("Ammended chrxhtyping")

Father <- Gtypes$Father[Gtypes$Chr=="X"]
Mother <- Gtypes$Mother[Gtypes$Chr=="X"]

if(sum(Father=="AB")>=1){
warning(paste("On the chr. x  of the father",sum(Father=="AB"),"(",sum(Father=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!"))
print("These will be treated as NoCalls")
Father[Father=="AB"] <- "NC"#There could not be htz SNPs on chromosome X of the father
}
	
RefSib <- Gtypes[Gtypes$Chr=="X",grep("Affected",colnames(Gtypes))]
ScSexes <- matrix(NA,length(ParScore),2)
rownames(ScSexes)<-names(ParScore)


for(sc in names(ParScore)){
	
	ScSexes[sc,1] <- sc
	if(is.na(ParScore[[sc]]["X","Par"])){print(paste("Sex of ",sc ,"could not be determined!!"))
		}else if(ParScore[[sc]]["X","Par"]==12){ScSexes[sc,2] <- "female" ;print(paste(sc,"is female"))
		} else if(ParScore[[sc]]["X","Par"]==22){ScSexes[sc,2] <- "male";print(paste(sc,"is male"))
		}else if(ParScore[[sc]]["X","Par"]==11){ScSexes[sc,2] <- "female";print(paste(sc,"is a female but without maternal X chromosome"))
		}else if(ParScore[[sc]]["X","Par"]==0){print(paste("Chr.X nullisomy of ",sc))
		}else{print(paste("Sex of ",sc ,"could not be determined!!"))}
}



	
	if(ScSexes[grep("Affected",ScSexes[,1]),2]=="female"){
	RefSibSex="female";print(paste("Referece sibling is female"))
	}else{
	RefSibSex="male";print(paste("Referece sibling is male"))}
	
#<---------------------------------------------------------------------------------
#Phasing genotype of the mother based on the affected sibling genotype

if(RefSibSex=="male"){
	Mother[RefSib=="BB" & Mother=="AB"] = "BA" #M1 is the carrier chromosome
	RefSib[RefSib=="AB"]="NC"#There could not be htz SNPs on chromosome X of male affected
	}else{
	#RefSib[RefSib=="AB" & Father=="BB"]="BA"#There could not be htz SNPs on chromosome X of male affected
	Mother[RefSib=="AB" & Father== "AA" & Mother=="AB"] = "BA" #M1 is the carrier chromosome
	Mother[RefSib=="BB" & Mother=="AB"] = "BA" #M1 is the carrier chromosome
	Mother[RefSib=="BB" & Father== "AA" & Mother=="AB"] = "NC"
	Mother[RefSib=="AA" & Father== "BB" & Mother=="AB"] = "NC"

	}
#>-------------------------------------------------------------------------------------
	
MatHaps <- do.call("rbind",strsplit(Mother,split=""))
PatHap <-	as.matrix(do.call("rbind",strsplit(Father,split=""))[,1])

M1<-as.matrix(MatHaps[,1])
M2<-as.matrix(MatHaps[,2])

#if(sum(grep("REF",Family))!=1){
Bls<- Gtypes[Gtypes$Chr=="X",ScSexes[unique(c(grep(SibPattern,ScSexes[,1]),grep("Affected",ScSexes[,1]))),1]]#Blastomeres
#}else{
#Bls<- Gtypes[Gtypes$Chr=="X",ScSexes[c(grep(SibPattern,ScSexes[,1])),1]]#Blastomeres
#}
Bls <- as.matrix(Bls)

for(bl in 1:ncol(Bls)){
	
	
	Bl<-as.matrix(Bls[,bl])
	UnInf<-which(Bl=="NC" | Mother=="NC"| Father=="NC"| RefSib=="NC" | Mother=="AA" | Mother=="BB" )

	B1<-matrix(0,nrow(Bl),1)
	B2<-matrix(0,nrow(Bl),1)

	#InfBl <- which(Bl!="00")
	if(is.na(ScSexes[bl,2])){
	 	print(paste("The Chr.X haplotype of Sample",colnames(Bls)[bl],"could not be reconstructed as the origin of this chromosome could not be determined..."))
	}else if(ScSexes[bl,2]=="male"){
		#<----------------------------------------------------
		#if the embryo is male
		Bl[Bl=="AB"]<-"NC"
		BlHaps<- do.call("rbind",strsplit(Bl,split=""))

		B1[M1==BlHaps[,1]]<-1
		B1[M2==BlHaps[,1]]<-2
		B1[UnInf]<-0				
		print(paste("The",colnames(Bls)[bl],"blastomere is from a male Embryo"))
		#>----------------------------------------------------
	}else if (ScSexes[bl,2]=="female"){
		#<----------------------------------------------------
		#if the embryo is female
		Bl[PatHap=="B" & Bl=="AB"]="BA"#Phasing of blastomere based on the father's genotype
		BlHaps<- do.call("rbind",strsplit(Bl,split=""))	   	
	  	B1[M1==BlHaps[,2]]<-1
		B1[M2==BlHaps[,2]]<-2
		B2[,1]<-1
		B1[UnInf]<-0				

		print(paste("The",colnames(Bls)[bl],"blastomere is from a female Embryo"))
		#>----------------------------------------------------
	}
	
	#B1[MatUnInf]<-0
	PhBl<-cbind(B1,B2)
	if(bl==1){
		PhBls<-PhBl
	}else{
		PhBls<-cbind(PhBls,PhBl)
		}

}#end bl loop



matHaps <-  PhBls[,seq(1,ncol(PhBls),2)]
patHaps <-  PhBls[,seq(2,ncol(PhBls),2)]

matHaps <- as.matrix(matHaps)
patHaps <- as.matrix(patHaps)


#if(sum(grep("REF",Family))!=1){
colnames(matHaps)<-paste(colnames(Gtypes)[unique(c(grep(SibPattern,colnames(Gtypes)),grep("Affected",colnames(Gtypes))))],"_Mat",sep="")
colnames(patHaps)<-paste(colnames(Gtypes)[unique(c(grep(SibPattern,colnames(Gtypes)),grep("Affected",colnames(Gtypes))))],"_Pat",sep="")
#}else{
#colnames(matHaps)<-paste(colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))],"_Mat",sep="")
#colnames(patHaps)<-paste(colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))],"_Pat",sep="")

#}

HapsChrX <- cbind(Gtypes[Gtypes$Chr=="X",1:3],patHaps,matHaps)
#colnames(HapsChrX)<-gsub("E0_Bl0_","",colnames(HapsChrX))

#write.table(Gtypes[,1:3],patHaps,matHaps),paste(outPath,"Famil".hap",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
HapsChrX

}#end function

