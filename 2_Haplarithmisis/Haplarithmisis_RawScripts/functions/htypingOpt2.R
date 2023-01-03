htypingOpt2 <- function(Gtypes,Window,Int,dataPo,ParScore,SibPattern,Family,outPath){

Haps <- vector("list",2)
names(Haps) <- c("dataHapRaw","dataHap")

Children <- unique(c(grep(SibPattern,colnames(Gtypes)),grep("Affected",colnames(Gtypes))))

Sibs <- Gtypes[,Children]

if(length(Children)==1){Sibs<-as.matrix(Sibs);colnames(Sibs)<-colnames(Gtypes)[grep(SibPattern,colnames(Gtypes))]} 

HapsAut <- htypingAutOpt2(Gtypes,Gtypes[,paste("Father_",Family,sep="")],Gtypes[,paste("Mother_",Family,sep="")],Gtypes[,grep("Affected",colnames(Gtypes))],Sibs,SibPattern)
print("Autosmes are haplotyped")
HapsChrX <- chrxhtyping1(Gtypes,ParScore,SibPattern)

Haps[["dataHapRaw"]] <- rbind(HapsAut,HapsChrX)
Haps[["dataHap"]]<-inthapnew1(Haps[["dataHapRaw"]],Window,Int)
Haps

}#end function

