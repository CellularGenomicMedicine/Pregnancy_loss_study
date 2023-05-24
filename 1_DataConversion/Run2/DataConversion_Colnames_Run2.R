###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To convert SNP genotyping raw data into data that is compatible with haplarithmisis. 

# input: .adj file containing SNP genotyping data for each individual pregnancy loss family. 

# output: Data files for individual pregnancy loss families that are compatible with haplarithmisis (format: .adj). 

###########################################################################################################################

#Changes consist of appropriate colnames.

args <- commandArgs(TRUE)
Family <- args[1]

dataPath <- "/Projects/PregnancyLoss/Data/Processed/"

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

