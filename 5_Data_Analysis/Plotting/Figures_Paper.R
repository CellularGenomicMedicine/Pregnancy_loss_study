library(ggplot2)
library(ggpattern)

###All chroms for MII, MI, Mitotic 
{
Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv", sep=";")

Ordering3 <- c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic I","Meiotic II")

Affected_Chromosomes$Aff_chr <- factor(Affected_Chromosomes$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Chrom_vs_SegOrigin <- ggplot(data=Affected_Chromosomes, aes(x=Aff_chr,fill = factor(Seg_origin, levels = Ordering3, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7, aes(y = (..count..)/sum(..count..))) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.151)) + 
  scale_fill_discrete(breaks = Ordering3) +
  ylab("Proportion") 
Chrom_vs_SegOrigin
}

###All chroms for Meiotic vs mitotic 
{
Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv", sep=";")

Ordering4 <- c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic")

Affected_Chromosomes$Aff_chr <- factor(Affected_Chromosomes$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Chrom_vs_SegOrigin2 <- ggplot(data=Affected_Chromosomes, aes(x=Aff_chr,fill = factor(Seg_origin2, levels = Ordering4, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7, aes(y = (..count..)/sum(..count..))) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.151)) + 
  scale_fill_discrete(breaks = Ordering4) +
  scale_fill_manual(values = c("Meiotic" = "#048F45", "Mitotic" = "#0CB702", "Mitotic & Meiotic " = "#92D895", "NA" = "#8C8C8C")) +
  ylab("Proportion") +
  xlab("Chromosome")
Chrom_vs_SegOrigin2
}

###All chroms for Parental_error_origin
{
Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv", sep=";")

Ordering5 <- c(NA,"Paternal","Maternal")

Affected_Chromosomes$Aff_chr <- factor(Affected_Chromosomes$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Chrom_vs_ParOrigin <- ggplot(data=Affected_Chromosomes, aes(x=Aff_chr,fill = factor(Parental_error_origin, levels = Ordering5, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7, aes(y = (..count..)/sum(..count..))) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.151)) + 
  scale_fill_manual(values = c("Maternal" = "#C77CFF", "Paternal" = "#F37735", "NA" = "#8C8C8C")) +
  ylab("Proportion") +
  xlab("Chromosome")
Chrom_vs_ParOrigin
}

###SPL & RPL comparison Parental error origin
{
Parental_origin_plot <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv", sep=";")

Ordering2 <- c(NA,"Paternal","Maternal")

###Position = fill 
ParOrigin_Plot3 <- ggplot(data = Parental_origin_plot, aes(x=Group, fill = factor(Parental_error_origin, levels = Ordering2, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.5) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank(), legend.position = "none") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(breaks = Ordering2) +
  scale_fill_manual(values = c("Maternal" = "#C77CFF", "Paternal" = "#F37735", "NA" = "#8C8C8C")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("")
ParOrigin_Plot3  
}

###Statistics###
#PatMat origin all aberrations
{
data <- data.frame("Maternal" = c(16,9), "Paternal" = c(4,8), row.names = c("RPL","SPL"), stringsAsFactors = F)
data

#Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
chisq.test(data)$expected
#Maternal Paternal
#RPL 13.51351 6.486486
#SPL 11.48649 5.513514
###Conclusion = chisquared is appropriate

#With Yates' continuity correction
chisq.test(data)
#Pearson's Chi-squared test with Yates' continuity correction
#data:  data
#X-squared = 1.9596, df = 1, p-value = 0.1616

#Without Yates' continuity correction
chisq.test(data, correct = FALSE)
#Pearson's Chi-squared test
#data:  data
#X-squared = 3.0703, df = 1, p-value = 0.07974
}

###SPL & RPL comparison Segregational error origin
{
Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv", sep=";")

Affected_Chromosomes$Seg_origin2 <- factor(Affected_Chromosomes$Seg_origin2, levels = c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic"))

Ordering3 <- c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic")

SegOrigin3_Plot <- ggplot(data=Affected_Chromosomes, aes(x=Group, fill = factor(Seg_origin2, levels =Ordering3, exclude = NULL))) + 
  geom_bar(stat="count",position = "fill", color =  "black",width=0.5) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.position = "none") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Meiotic" = "#048F45", "Mitotic" = "#0CB702", "Mitotic & Meiotic " = "#92D895", "NA" = "#8C8C8C")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("")
SegOrigin3_Plot
}

###Statistics###
#Mitotic/meiotic only for all aberrations 
{
data <- data.frame("Mitotic" = c(6,8), "Meiotic" = c(20,10), row.names = c("RPL","SPL"), stringsAsFactors = F)
data

###Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
chisq.test(data)$expected
#Conclusion = chisquared is appropriate

###With Yates' continuity correction
chisq.test(data)
#X-squared = 1.3619, df = 1, p-value = 0.2432

###Without Yates' continuity correction
chisq.test(data, correct = FALSE)
#X-squared = 2.2385, df = 1, p-value = 0.1346


#Mitotic/meiotic only for trisomy cases 

data <- data.frame("Mitotic" = c(5,4), "Meiotic" = c(11,5), row.names = c("RPL","SPL"), stringsAsFactors = F)
data

data <- data.frame("Mitotic" = c(5,4), "Meiotic" = c(11,5), row.names = c("RPL","SPL"), stringsAsFactors = F)
data

###Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
chisq.test(data)$expected
#Mitotic Meiotic
#RPL    5.76   10.24
#SPL    3.24    5.76
###Conclusion = chisquared approximation may be incorrect

###Fisher's exact test 
fisher.test(data)
#Fisher's Exact Test for Count Data

#data:  data
#p-value = 0.6707
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.07847419 4.30133893
#sample estimates:
#odds ratio 
#  0.581578 

###Chi-squared With Yates' continuity correction
chisq.test(data)
#Pearson's Chi-squared test with Yates' continuity correction
#data:  data
#X-squared = 0.050938, df = 1, p-value = 0.8214

###Chi-squared Without Yates' continuity correction
chisq.test(data, correct = FALSE)
#Pearson's Chi-squared test

#data:  data
#X-squared = 0.43523, df = 1, p-value = 0.5094



#Mitotic/meiotic only for only trisomy cases

MitMei_TrisOnly <- subset(Affected_Chromosomes, Indication_w_translocation == "Trisomy")
MitMei_TrisOnly_RPL <- subset(MitMei_TrisOnly, Group == "RPL")
MitMei_TrisOnly_SPL <- subset(MitMei_TrisOnly, Group == "SPL")
table(MitMei_TrisOnly_RPL$Seg_origin2)
table(MitMei_TrisOnly_SPL$Seg_origin2)

data <- data.frame("Mitotic" = c(4,4), "Meiotic" = c(10,5), row.names = c("RPL","SPL"), stringsAsFactors = F)
data

###Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
chisq.test(data)$expected
#Mitotic  Meiotic
#RPL 4.869565 9.130435
#SPL 3.130435 5.869565
###Conclusion = chisquared approximation may be incorrect

###Fisher's exact test 
fisher.test(data)
#Fisher's Exact Test for Count Data

#data:  data
#p-value = 0.657
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 0.06251808 4.06861482
#sample estimates:
#odds ratio 
# 0.5158756 

###With Yates' continuity correction
chisq.test(data)
#Pearson's Chi-squared test with Yates' continuity correction

#data:  data
#X-squared = 0.1099, df = 1, p-value = 0.7403

#Without Yates' continuity correction
chisq.test(data, correct = FALSE)
#Pearson's Chi-squared test

#data:  data
#X-squared = 0.60847, df = 1, p-value = 0.4354
}