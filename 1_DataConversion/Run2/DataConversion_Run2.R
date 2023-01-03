
###Samples_2021_02_19
dataPath<- "/ifs/data/research/CGM/Projects/Miscarr/2021/"
outPath <- "/ifs/data/research/CGM/Projects/Miscarr/2021/Families/"

data <- read.table(paste0(dataPath,"MiscarrSamples_2021_02_19.txt"),sep="\t",header=T,nrow=1)
dataAll <- read.table(paste0(dataPath,"MiscarrSamples_2021_02_19.txt"),sep="\t",header=T)

Cols <- gsub(".Log.R.Ratio","",colnames(data))
Cols <- gsub(".GType","",Cols)
Cols <- gsub(".B.Allele.Freq","",Cols)

Fams <- unique(do.call("rbind",strsplit(Cols,"_"))[,2])[-c(1:3)]

for(f in Fams){ 
  dataOut<- dataAll[,c(1:3,grep(f,colnames(data)))]

  print(f)
  write.table(dataOut, paste0(outPath,f,".adj"),col.names=T,row.names=F,quote=F,sep="\t") 
}#end f loop

