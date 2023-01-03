BAF_func <- function(alel_id,ref_df_allels,ref_df){
if(alel_id=="ADO")    { if(toupper(ref_df$alels)==toupper(ref_df$REFallele)) {BAF_call <- 0} else {BAF_call <- 1}} 
if(alel_id=="NORMAL") {if(toupper(unlist(strsplit(ref_df$alels,""))[1])==toupper(ref_df$REFallel)) {BAF_call <- ref_df_allels[,colnames(ref_df_allels)==(unlist(strsplit(ref_df$alels,""))[2])]/sum(ref_df_allels[ref_df_allels>0]) } else {BAF_call <- ref_df_allels[,colnames(ref_df_allels)==(unlist(strsplit(ref_df$alels,""))[1])]/sum(ref_df_allels[ref_df_allels>0])} } 
if(alel_id=="ADI") {if(toupper(unlist(strsplit(ref_df$alels,""))[1])==toupper(ref_df$REFallel)) {BAF_call <- sum(ref_df_allels[,colnames(ref_df_allels)!=(unlist(strsplit(ref_df$alels,""))[1])])/ sum(ref_df_allels[ref_df_allels>0]) } else {BAF_call <- ref_df_allels[,colnames(ref_df_allels)==(unlist(strsplit(ref_df$alels,""))[1])]/sum(ref_df_allels[ref_df_allels>0])} }

BAF_call
}

