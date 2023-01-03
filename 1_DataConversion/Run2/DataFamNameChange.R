args <- commandArgs(TRUE)
Family <- args[1]

dataPath <- "/ifs/data/research/CGM/Projects/Miscarr/2021/Families/"

dataOut <- read.table(paste0(dataPath,Family),sep="\t",header=T)

f <- gsub(".adj","",Family)

colnames(dataOut) <- gsub(paste0("_",f),"",colnames(dataOut))
colnames(dataOut) <- gsub("EM", "E01_Bl01-",colnames(dataOut))
colnames(dataOut) <- gsub("CV", "Affected-",colnames(dataOut))
colnames(dataOut) <- gsub("Father.", "Father",colnames(dataOut))
colnames(dataOut) <- gsub("Mother.", "Mother",colnames(dataOut))

colnames(dataOut) <- gsub("-.*?\\.","",colnames(dataOut))
colnames(dataOut) <- gsub("GType",".GType",colnames(dataOut))
colnames(dataOut) <- gsub("B.Allele.Freq",".B.Allele.Freq",colnames(dataOut))
colnames(dataOut) <- gsub("Log.R.Ratio",".Log.R.Ratio",colnames(dataOut))

write.table(dataOut,paste0(dataPath,Family),col.names=T,row.names=F,quote=F,sep="\t")
~
