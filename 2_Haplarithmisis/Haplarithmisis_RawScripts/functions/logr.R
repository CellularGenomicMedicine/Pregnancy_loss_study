logr <- function(BAF,Name){
	
	BAF <- BAF[BAF[,"count"]>=2,]
	logR <- log2(BAF[,"count"]/median(BAF[,"count"]))
	Out <- data.frame(BAF[,c("chr","cord","BAF")],logR,stringsAsFactors=F)
	colnames(Out) <- c("Chr","Position",paste0(Name,".B.Allele.Freq"),paste0(Name,".Log.R.Ratio"))
	Out
}#end logr function
