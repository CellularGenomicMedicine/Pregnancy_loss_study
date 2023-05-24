###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To convert SNP genotyping raw data into data that is compatible with haplarithmisis. To extract individual family data files from the raw SNP genotyping data file.

# input: Raw SNP genotyping data from the 3rd run (out of 3 total), available on GEO (GSE228151).

# output: Data files for individual pregnancy loss families that are compatible with haplarithmisis (format: .adj). 

###########################################################################################################################

dataPath<- "/Projects/PregnancyLoss/Data/Raw/"
outPath <- "/Projects/PregnancyLoss/Data/Processed/"

data <- read.table(paste0(dataPath,"MiscarrSamples_2022_02_14.txt"),sep="\t",header=T,nrow=1)
dataAll <- read.table(paste0(dataPath,"MiscarrSamples_2022_02_14.txt"),sep="\t",header=T)

Cols <- gsub(".Log.R.Ratio","",colnames(data))
Cols <- gsub(".GType","",Cols)
Cols <- gsub(".B.Allele.Freq","",Cols)

Cols <- gsub("father","",Cols)
Cols <- gsub("mother","",Cols)
Cols <- gsub("chorionicvilli","",Cols)
Cols <- gsub("mesoderma","",Cols)
Cols <- gsub("chorionic.villi","",Cols)

Fams <- unique(do.call("rbind",strsplit(Cols,"_"))[,1])[-c(1:3)]

for(f in Fams){ 
  dataOut<- dataAll[,c(1:3,grep(f,colnames(data)))]
  
  print(f)
  write.table(dataOut, paste0(outPath,f,".adj"),col.names=T,row.names=F,quote=F,sep="\t") 
}
