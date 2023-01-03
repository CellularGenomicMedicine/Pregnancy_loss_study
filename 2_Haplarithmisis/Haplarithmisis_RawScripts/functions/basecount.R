basecount <- function(df,sample_ID){
 t(pbsapply(1:nrow(df), function(w) {
   tbl <- data.frame(table(unlist(strsplit(as.character(df[w,paste(sample_ID,"base_calls",sep=".")]),""))),stringsAsFactors=FALSE)
   ref_df <- data.frame(df[w,],stringsAsFactors=FALSE)
   ref_df$ref  <- sum(as.integer(tbl[tbl[,]==".","Freq"]),as.integer(tbl[tbl[,1]==",","Freq"]))
   ref_df$a    <- if (length(tbl[toupper(tbl[,1])=="A",1])==0) 0 else if (length(tbl[toupper(tbl[,1])=="A",1])==1) as.integer(tbl[toupper(tbl[,1])=="A","Freq"]) else sum(as.integer(tbl[tbl[,1]=="a","Freq"]),as.integer(tbl[tbl[,1]=="A","Freq"]))
   ref_df$c    <- if (length(tbl[toupper(tbl[,1])=="C",1])==0) 0 else if (length(tbl[toupper(tbl[,1])=="C",1])==1) as.integer(tbl[toupper(tbl[,1])=="C","Freq"]) else sum(as.integer(tbl[tbl[,1]=="c","Freq"]),as.integer(tbl[tbl[,1]=="C","Freq"]))
   ref_df$g    <- if (length(tbl[toupper(tbl[,1])=="G",1])==0) 0 else if (length(tbl[toupper(tbl[,1])=="G",1])==1) as.integer(tbl[toupper(tbl[,1])=="G","Freq"]) else sum(as.integer(tbl[tbl[,1]=="g","Freq"]),as.integer(tbl[tbl[,1]=="G","Freq"]))
   ref_df$t    <- if (length(tbl[toupper(tbl[,1])=="T",1])==0) 0 else if (length(tbl[toupper(tbl[,1])=="T",1])==1) as.integer(tbl[toupper(tbl[,1])=="T","Freq"]) else sum(as.integer(tbl[tbl[,1]=="t","Freq"]),as.integer(tbl[tbl[,1]=="T","Freq"]))
   if (toupper(df[w,"REFallele"])=="A") ref_df$a <- sum(ref_df$a,ref_df$ref)
   if (toupper(df[w,"REFallele"])=="C") ref_df$c <- sum(ref_df$c,ref_df$ref)
   if (toupper(df[w,"REFallele"])=="G") ref_df$g <- sum(ref_df$g,ref_df$ref)
   if (toupper(df[w,"REFallele"])=="T") ref_df$t <- sum(ref_df$t,ref_df$ref)
   ref_df$Deletions  <- if (length(tbl[tbl[,1]=="+",1])==0) 0 else as.integer(tbl[tbl[,1]=="+","Freq"])
   ref_df$Insertions <- if (length(tbl[tbl[,1]=="-",1])==0) 0 else as.integer(tbl[tbl[,1]=="-","Freq"])
   ref_df$star       <- if (length(tbl[tbl[,1]=="*",1])==0) 0 else as.integer(tbl[tbl[,1]=="*","Freq"])
   ref_df$alels      <- paste(colnames(ref_df[,c("a","c","g","t")])[ref_df[,c("a","c","g","t")]>0],collapse='')
   ref_df$Status     <- if(nchar(ref_df$alels)<2) "ADO" else {if(nchar(ref_df$alels)==2) "NORMAL" else  "ADI" }
   ref_df_allels     <- ref_df[,c("a","c","g","t")]
   ref_df$BAF        <- BAF_func(ref_df$Status,ref_df_allels,ref_df)
   rownames(ref_df) <- paste(ref_df[,"chr"],ref_df[,"cord"],sep="_")
   ref_df
   }
 ))
}
