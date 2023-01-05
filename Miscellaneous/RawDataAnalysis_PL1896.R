###PL1896 case where embryo tissue shows only paternal contirbution to chr. X 
###Could be Turner syndrome with paternal X retained
###Could be Klinefalter syndrome with a rescue mechanism leading to paternal UPD chr XY
###To deduce which it is; check chr. Y SNP values in the SNP genotyping files 
###PL1063 used as comparison, normal plot with male embryo
###PL1246 used as comparison, normal plot with female embryo

PL1896 <- read.delim("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/PL1896_Y_chrom/PL1896.adj")
PL1063 <- read.delim("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/PL1896_Y_chrom/PL1063.adj")
PL1246 <- read.delim("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/PL1896_Y_chrom/PL1246.adj")

EmbryoY_PL1896 <- PL1896[PL1896$Chr == "Y", c("Chr","E01_Bl01.GType","E01_Bl01.B.Allele.Freq","E01_Bl01.Log.R.Ratio")]
EmbryoY_PL1063 <- PL1063[PL1063$Chr == "Y", c("Chr","E01_Bl01.GType","E01_Bl01.B.Allele.Freq","E01_Bl01.Log.R.Ratio")]
EmbryoY_PL1246 <- PL1246[PL1246$Chr == "Y", c("Chr","E01_Bl01.GType","E01_Bl01.B.Allele.Freq","E01_Bl01.Log.R.Ratio")]

View(EmbryoY_PL1896)
View(EmbryoY_PL1063)
View(EmbryoY_PL1246)

###Conclusion: no Y SNPs visible, therefore cannot be Klinefilter pat UPD, but is instead confirmed Turner syndrome. 
###Confirmed Turner syndrome