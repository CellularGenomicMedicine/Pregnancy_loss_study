###File conversion and fambuilder for round1 of pregnancy loss study families.

rm(list=ls(all=T))
dataPath<- "/ifs/data/research/CGM/Projects/Miscarr/MiscarriageStudy/Data/Run1/"
outPath <- "/ifs/data/research/CGM/Projects/Miscarr/MiscarriageStudy/Data/Run1/FamiliesRun1/"

data <- read.table(paste0(dataPath,"Families_2020.txt"),sep="\t",header=T,stringsAsFactors=F)


A <- gsub("X","S",colnames(data))

A2 <- gsub(".GType.1","_b.GType",A)
A2 <- gsub(".B.Allele.Freq.1","_b.B.Allele.Freq",A2)
A2 <- gsub(".Log.R.Ratio.1","_b.Log.R.Ratio",A2)

#A2 <- gsub(".1$","_1",A)
A2 <- gsub(".1.GType","_F.GType",A2,fixed=T)
A2 <- gsub(".1.B.Allele.Freq","_F.B.Allele.Freq",A2,fixed=T)
A2 <- gsub(".1.Log.R.Ratio","_F.Log.R.Ratio",A2,fixed=T)
A2 <- gsub(".2.GType","_M.GType",A2,fixed=T)
A2 <- gsub(".2.B.Allele.Freq","_M.B.Allele.Freq",A2,fixed=T)
A2 <- gsub(".2.Log.R.Ratio","_M.Log.R.Ratio",A2,fixed=T)

FamData <- read.table("/ifs/data/research/CGM/Projects/Miscarr/MiscarriageStudy/Data/Run1/FamilyNames_2020_2.txt",sep="\t",header=T,stringsAsFactors=F)
FamData[,2] <- gsub("-1$","_F",FamData[,2])
FamData[,2] <- gsub("-2$","_M",FamData[,2])
FamData[,2] <- paste0("S",FamData[,2])
FamData[,2] <- gsub(".1$","_b",FamData[,2])

FamData[1:10,c("Sample.ID","Fams1")]

head(cbind(colnames(data),A2),10)

#FamsData2 <- FamsData2[]
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
	write.table(dataOut, paste0(outPath,f,".adj"),col.names=T,row.names=F,quote=F,sep="\t") }#end f loop


