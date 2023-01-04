
rm(list=ls(all=T))

data <- read.table("/ifs/data/research/CGM/Projects/Miscarr/MiscarriageStudy/Data/Run1/Callrates_2020.txt",sep="\t",header=T,stringsAsFactors=F)
Peds <- read.table("/ifs/data/research/CGM/Projects/Miscarr/MiscarriageStudy/Data/Run1/FamilyNames_2020.txt",sep="\t",header=T,stringsAsFactors=F)

Fams <-rep(NA,nrow(data))
Fams[grep("-1",data[,2])] <- paste("Father",gsub("-1","",paste0("PL",data[grep("-1",data[,2]),2])),sep="_")
Fams[grep("-2",data[,2])] <- paste("Mother",gsub("-2","",paste0("PL",data[grep("-2",data[,2]),2])),sep="_")

for(f in Peds[,1]){ Fams[data[,"Sample.ID"]==f] <- paste(Peds[Peds[,"Sample"]==f,"Sample"],paste0("PL",Peds[Peds[,"Sample"]==f,"Family"]),sep="_")}#end f loop

A <- do.call("rbind",strsplit(Fams,"_"))

Fams1 <- Fams

Fams1[duplicated(Fams1)] <- paste("ExtMes",Fams[duplicated(Fams1)],sep="")

FamsU <- unique(A[,2])
data2 <- data.frame(data,Fams1)
data2[,"Fams1"] <- as.character(data2[,"Fams1"])

for(r in 1:nrow(data2)){if((sum(grep("Father",data2[r,"Fams1"]))==0 & sum(grep("Mother",data2[r,"Fams1"]))==0 & sum(grep("ExtMes",data2[r,"Fams1"]))==0)==T){data2[r,"Fams1"]<-paste0("ChrVil",data2[r,"Fams1"])}}
data2[duplicated(data2[,"Sample.ID"]),"Sample.ID"] <- paste0(data2[duplicated(data2[,"Sample.ID"]),"Sample.ID"],".1")

write.table(data2,"/ifs/data/research/CGM/Projects/Miscarr/MiscarriageStudy/Data/Run1/FamilyNames_2020_2.txt",quote=F,sep="\t",col.names=T,row.names=F)


