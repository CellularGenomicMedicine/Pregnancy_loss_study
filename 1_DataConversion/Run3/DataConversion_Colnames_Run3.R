###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To convert SNP genotyping raw data into data that is compatible with haplarithmisis.

# input: .adj file containing SNP genotyping data for each individual pregnancy loss family. 

# output: Data files for individual pregnancy loss families that are compatible with haplarithmisis (format: .adj). 

###########################################################################################################################

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

