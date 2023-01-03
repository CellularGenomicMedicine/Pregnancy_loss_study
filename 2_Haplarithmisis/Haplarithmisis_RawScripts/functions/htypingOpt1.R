htypingOpt1 <- function(Gtypes,dataPo,ParScore,Parent1,SibPattern,Family,outPath){

print("*******************************************")
print("****     Option1Htype analysis...     ****")
	
ScSexes <- matrix(NA,length(ParScore),2)
rownames(ScSexes)<-names(ParScore)

for(sc in names(ParScore)){

        ScSexes[sc,1] <- sc
        if(is.na(ParScore[[sc]]["X","Par"])){print(paste("Sex of ",sc ,"could not be determined!!"))
                }else if(ParScore[[sc]]["X","Par"]==12){ScSexes[sc,2] <- "female" ;print(paste(sc,"is female"
))
                } else if(ParScore[[sc]]["X","Par"]==22){ScSexes[sc,2] <- "male";print(paste(sc,"is male"))
                }else if(ParScore[[sc]]["X","Par"]==11){ScSexes[sc,2] <- "female";print(paste(sc,"is a female
 but without maternal X chromosome"))
                }else if(ParScore[[sc]]["X","Par"]==0){ScSexes[sc,2] <- "ND";print(paste("Chr.X nullisomy of ",sc))
                }else{print(paste("Sex of ",sc ,"could not be determined!!"))}
}

Father <- Gtypes[Gtypes[,"Chr"]=="X",paste("Father_",Family,sep="")]

if(sum(Father=="AB")>=1){
	warning(paste("On the chr. X  of the father",sum(Father=="AB"),"(",sum(Father=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!"))
	print("These will be treated as NoCalls")
	Father[Father=="AB"] <- "NC"#There could not be htz SNPs on chromosome X of the father
	Gtypes[Gtypes[,"Chr"]=="X",paste("Father_",Family,sep="")] <- Father
}


#------------------------------------------------------------
#------------			Phasing chr X			-------------

chrxhtypingOpt1 <- function(Gtypes,Parent1,ScSexes,Family){

GtypesChrX <- Gtypes[Gtypes$Chr=="X",]	
MotherChrX <- GtypesChrX[,paste("Mother_",Family,sep="")]
FatherChrX <- GtypesChrX[,paste("Father_",Family,sep="")]
GrandfatherChrX <- GtypesChrX[,paste("Grandfather_",Family,sep="")]
GrandmotherChrX <- GtypesChrX[,paste("Grandmother_",Family,sep="")]

HapXMat <- matrix(0,nrow(GtypesChrX),ncol(as.matrix(GtypesChrX[,grep(SibPattern,colnames(GtypesChrX))])))
colnames(HapXMat)<-paste(colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))])

if(Parent1==paste("Mother_",Family,sep="")){

	MotherChrX[(GtypesChrX$Grandfather=="BB" & GtypesChrX$Mother=="AB")]<-"BA"
	
	MotherChrX[(GtypesChrX$Grandfather=="NC" & GtypesChrX$Mother=="AB") |
					   (GtypesChrX$Grandfather=="BB" & GtypesChrX$Grandmother=="BB" & GtypesChrX$Mother=="AB") | 
					   (GtypesChrX$Grandfather=="AA" & GtypesChrX$Grandmother=="AA" & GtypesChrX$Mother=="AB") ]<-"NC"
	}	   

na.omit(ScSexes[ScSexes[,2]=="male",1])


for(s in as.character(na.omit(ScSexes[ScSexes[,2]=="male",1]))){
	
	print(paste(s,"is a male sample but has",sum(GtypesChrX[,s]=="AB"),"heterozygous SNP-calls; these are treated as NoCalls"))			   
	GtypesChrX[GtypesChrX[,s]=="AB",s]<-"NC"
	HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "BB") |
				   (MotherChrX=="BA" & GtypesChrX[,s] == "AA"),s] <- 2
	HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "AA") |
					(MotherChrX=="BA" & GtypesChrX[,s] == "BB"),s] <- 1	
	HapXMat[GrandfatherChrX=="NC" | GrandmotherChrX=="NC" ,s] <- 0						
						
	print(s)							
 }#end s loop



for(s in as.character(na.omit(ScSexes[ScSexes[,2]=="female",1]))){
	
	HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "BB") |
				   (MotherChrX=="BA" & GtypesChrX[,s] == "AA") |
				   (MotherChrX=="AB" & FatherChrX == "AA" & GtypesChrX[,s] == "AB") |
				   (MotherChrX=="BA" & FatherChrX == "BB" & GtypesChrX[,s] == "AB") ,s] <- 2
	HapXMat[(MotherChrX=="AB" & GtypesChrX[,s] == "AA") |
				    (MotherChrX=="BA" & GtypesChrX[,s] == "BB") |
					(MotherChrX=="BA" & FatherChrX == "AA" & GtypesChrX[,s] == "AB")|
					(MotherChrX=="AB" & FatherChrX == "BB" & GtypesChrX[,s] == "AB"),s] <- 1
	HapXMat[GrandfatherChrX=="NC" | GrandmotherChrX=="NC" | FatherChrX=="NC" ,s] <- 0						
	
	print(s)							
 }#end s loop
	
	HapXMat

}#end chrxhtypingOpt1 function

 HapXMat <- chrxhtypingOpt1(Gtypes,Parent1,ScSexes,Family)


PhasedPar <- Gtypes[,Parent1]
if(Parent1==paste("Mother_",Family,sep="")){Parent2 = paste("Father_",Family,sep="")}else{Parent2 = paste("Mother_",Family,sep="")}



PhasedPar[ (Gtypes$Grandfather=="AB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB")   |
				   (Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB") ] <-"BA"

PhasedPar[(Gtypes$Grandfather=="AA" & Gtypes$Grandmother=="AA" & Gtypes[,Parent1]=="AB") |
				   (Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="BB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes$Grandfather=="AA" & Gtypes$Grandmother=="BB" & Gtypes[,Parent1]=="BB") | 
				   (Gtypes$Grandfather=="AA" & Gtypes$Grandmother=="BB" & Gtypes[,Parent1]=="AA") | 
				   ((Gtypes$Grandfather=="AA" | Gtypes$Grandmother=="AA") & Gtypes[,Parent1]=="BB") |
				   ((Gtypes$Grandfather=="AA" | Gtypes$Grandmother=="AB") & Gtypes[,Parent1]=="BB") |
                                   ((Gtypes$Grandfather=="BB" | Gtypes$Grandmother=="AA") & Gtypes[,Parent1]=="BB") |
                                   ((Gtypes$Grandfather=="BB" | Gtypes$Grandmother=="AA") & Gtypes[,Parent1]=="AA") |
                                   ((Gtypes$Grandfather=="BB" | Gtypes$Grandmother=="AB") & Gtypes[,Parent1]=="AA") |
				   ((Gtypes$Grandfather=="BB" | Gtypes$Grandmother=="BB") & Gtypes[,Parent1]=="AA") |
				   ((Gtypes$Grandfather=="NC" | Gtypes$Grandmother=="NC") & Gtypes[,Parent1]=="AB") 	] <- "NC"

 
Gtypes[,Parent1] <- PhasedPar

Sibs <-  colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))]

for(sib in Sibs){

Hap1 <- rep(0,nrow(Gtypes))
Hap2 <- rep(0,nrow(Gtypes))


print(sib)

GtypeChild <- Gtypes[,sib]

#-------------------------------------------------------------------------------------------------------
Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") | 
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA")] <- 1



Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AA") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AB") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="AA") |  
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,sib]=="BB") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="AB") | 
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,sib]=="BB")] <- 2
			 
Hap1[Gtypes$Grandfather=="AB" & Gtypes$Grandmother=="AB"] <-0
 			 
#----------------------------------------------------------------------------------------------------------

if(sib==Sibs[1]){ParHaps <-Hap1 }else{ParHaps <-cbind(ParHaps,Hap1)}			 
			 
}#end sib loop			 

colnames(ParHaps)<- Sibs


	#for(s in as.character(ScSexes[ScSexes[,2]=="male",1])){
	
		#print(paste(s,"is a male sample but has",sum(GtypesChrX[,s]=="AB"),"heterozygous SNP-calls; these are treated as NoCalls"))			   
		#GtypesChrX[GtypesChrX[,s]=="AB",s]<-"NC"

	#}#end s loop

#------------------------------------------------------------
Haps <- vector("list",2)
names(Haps) <- c("dataHapRaw","dataHap")

ParHaps[Gtypes$Chr=="X",] <- HapXMat

ParHapsAll <- cbind(Gtypes[,c("Name","Chr","Position")],ParHaps)
dataHapRaw <- ParHapsAll

dataHap1<-inthapnew1(dataHapRaw,Window,Int)

if(Parent1==paste("Father_",Family,sep="")){ par1="_Pat"; par2="_Mat" }else{ par1="_Mat"; par2="_Pat" }

Haps[["dataHapRaw"]] <- cbind(ParHapsAll,matrix(0,nrow(ParHapsAll),(ncol(ParHapsAll)-3)))
colnames(Haps[["dataHapRaw"]])<- c("Name","Chr","Position",paste(colnames(dataHap1)[-c(1:3)],par1,sep=""),paste(colnames(dataHap1)[-c(1:3)],par2,sep=""))

Haps[["dataHap"]] <- cbind(dataHap1,matrix(0,nrow(dataHap1),(ncol(dataHap1)-3)))
colnames(Haps[["dataHap"]])<- colnames(Haps[["dataHapRaw"]])


Haps

}#end function


