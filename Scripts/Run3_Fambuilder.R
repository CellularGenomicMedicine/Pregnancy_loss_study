args <- commandArgs(TRUE)
Family <- args[1]

dataPath <- "/ifs/data/research/CGM/Projects/Miscarr/2022/Families/FamsRound3/Test/"

dataPL <- read.table(paste0(dataPath,Family),sep="\t",header=T)

colnames(dataPL) <- gsub("father","Father",colnames(dataPL))
colnames(dataPL) <- gsub("mother","Mother",colnames(dataPL))
colnames(dataPL) <- gsub("chorionicvilli","Affected",colnames(dataPL))
colnames(dataPL) <- gsub("mesoderma","E01_Bl01",colnames(dataPL))
colnames(dataPL) <- gsub(Family,"",colnames(dataPL))

write.table(dataPL,paste0(dataPath,Family),col.names=T,row.names=F,quote=F,sep="\t")
