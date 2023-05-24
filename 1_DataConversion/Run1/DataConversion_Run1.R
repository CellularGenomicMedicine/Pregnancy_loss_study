###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To convert SNP genotyping raw data into data that is compatible with haplarithmisis. To extract individual family data files from the raw SNP genotyping data file.

# input: Raw SNP genotyping data from the 1st run (out of 3 total), available on GEO (GSE228151).

# output: Data files for individual pregnancy loss families that are compatible with haplarithmisis (format: .adj). 

###########################################################################################################################

#Production of intermediate file that contains the column "Fams1" stating which sample ID belongs to which family ID and 
# individual (Mother or Father) or fetal tissue type (Extramebryonic mesoderm or Chorionic villi). 
#This file is required for the proper extraction of individual family data from the raw input files where all families are combined in a single file.  

rm(list=ls(all=T))

data <- read.table("/Projects/PregnancyLoss/DataConversion/Run1/Callrates_2020.txt",sep="\t",header=T,stringsAsFactors=F)
Peds <- read.table("/Projects/PregnancyLoss/DataConversion/Run1/FamilyNames_2020.txt",sep="\t",header=T,stringsAsFactors=F)

Fams <-rep(NA,nrow(data))
Fams[grep("-1",data[,2])] <- paste("Father",gsub("-1","",paste0("PL",data[grep("-1",data[,2]),2])),sep="_")
Fams[grep("-2",data[,2])] <- paste("Mother",gsub("-2","",paste0("PL",data[grep("-2",data[,2]),2])),sep="_")

for(f in Peds[,1]){ Fams[data[,"Sample.ID"]==f] <- paste(Peds[Peds[,"Sample"]==f,"Sample"],paste0("PL",Peds[Peds[,"Sample"]==f,"Family"]),sep="_")}

A <- do.call("rbind",strsplit(Fams,"_"))

Fams1 <- Fams

Fams1[duplicated(Fams1)] <- paste("ExtMes",Fams[duplicated(Fams1)],sep="")

FamsU <- unique(A[,2])
data2 <- data.frame(data,Fams1)
data2[,"Fams1"] <- as.character(data2[,"Fams1"])

for(r in 1:nrow(data2)){if((sum(grep("Father",data2[r,"Fams1"]))==0 & sum(grep("Mother",data2[r,"Fams1"]))==0 & sum(grep("ExtMes",data2[r,"Fams1"]))==0)==T){data2[r,"Fams1"]<-paste0("ChrVil",data2[r,"Fams1"])}}
data2[duplicated(data2[,"Sample.ID"]),"Sample.ID"] <- paste0(data2[duplicated(data2[,"Sample.ID"]),"Sample.ID"],".1")


###Extraction of individual families from the raw SNP genotyping data file. 
#Individual family data files are given the appropriate column names required for hapalrithmisis. 
#Column names consist of "Name", "Chr", "Position", "E01_Bl01" (fetal extraembryonic mesoderm DNA), "Affected" (fetal chorionic villi DNA), "Mother", "Father". 

dataPath<- "/Projects/PregnancyLoss/Data/Raw/"
outPath <- "/Projects/PregnancyLoss/Data/Processed/"


data <- read.table(paste0(dataPath,"Families_2020.txt"),sep="\t",header=T,stringsAsFactors=F)

A <- gsub("X","S",colnames(data))

A2 <- gsub(".GType.1","_b.GType",A)
A2 <- gsub(".B.Allele.Freq.1","_b.B.Allele.Freq",A2)
A2 <- gsub(".Log.R.Ratio.1","_b.Log.R.Ratio",A2)

A2 <- gsub(".1.GType","_F.GType",A2,fixed=T)
A2 <- gsub(".1.B.Allele.Freq","_F.B.Allele.Freq",A2,fixed=T)
A2 <- gsub(".1.Log.R.Ratio","_F.Log.R.Ratio",A2,fixed=T)
A2 <- gsub(".2.GType","_M.GType",A2,fixed=T)
A2 <- gsub(".2.B.Allele.Freq","_M.B.Allele.Freq",A2,fixed=T)
A2 <- gsub(".2.Log.R.Ratio","_M.Log.R.Ratio",A2,fixed=T)

FamData <- data2
FamData[,2] <- gsub("-1$","_F",FamData[,2])
FamData[,2] <- gsub("-2$","_M",FamData[,2])
FamData[,2] <- paste0("S",FamData[,2])
FamData[,2] <- gsub(".1$","_b",FamData[,2])

FamData[1:10,c("Sample.ID","Fams1")]

head(cbind(colnames(data),A2),10)

data2 <- data
FamData_b <- FamData[grep("_b",FamData[,"Sample.ID"]),]
FamData_a <- FamData[-c(grep("_b",FamData[,"Sample.ID"])),]

for(c in FamData_b[,"Sample.ID"]){A2 <- gsub(FamData_b[FamData_b[,"Sample.ID"]==c,"Sample.ID"], FamData_b[FamData_b[,2]==c,"Fams1"],A2)}
for(c in FamData_a[,"Sample.ID"]){A2 <- gsub(FamData_a[FamData_a[,"Sample.ID"]==c,"Sample.ID"], FamData_a[FamData_a[,2]==c,"Fams1"],A2)}

colnames(data2) =A2

colnames(data2) <- gsub("ExtMes","E01_Bl",colnames(data2))
colnames(data2) <- gsub("ChrVil","Affected",colnames(data2))

for(f in gsub("X","PL",unique(do.call("rbind",strsplit(FamData[,"Fams1"],"_"))[,2]))){ 
	dataOut<- data2[,c(1:3,grep(f,colnames(data2)))]
	colnames(dataOut) <- gsub(paste0("_",f),"",colnames(dataOut))
	colnames(dataOut) <- gsub("Affected.*.GType","Affected.Gtype", colnames(dataOut)) 
	colnames(dataOut) <- gsub("Affected.*.B.Allele.Freq","Affected.B.Allele.Freq", colnames(dataOut))
	colnames(dataOut) <- gsub("Affected.*.Log.R.Ratio","Affected.Log.R.Ratio", colnames(dataOut))
	write.table(dataOut, paste0(outPath,f,".adj"),col.names=T,row.names=F,quote=F,sep="\t") 
	}

