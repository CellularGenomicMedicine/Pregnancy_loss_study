rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

###All chroms for segregational and parental origin combined. 

Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")

Newfile <- rbind(Affected_Chromosomes, Affected_Chromosomes)

OriginCol <- data.frame(Origin = c(as.character(Affected_Chromosomes$Parental_error_origin), as.character(Affected_Chromosomes$Seg_origin2)))

Newfile <- cbind(Newfile, OriginCol) 

Newfile$Origin_Type <- c(rep("Parental origin", times = 39), rep("Segregational origin", times = 39))

Newfile$Aff_chr <- factor(Newfile$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Newfile$Origin <- factor(Newfile$Origin, levels = c(NA, "Mitotic & Meiotic ", "Mitotic", "Meiotic", "Paternal","Maternal"))

Order <- c(NA, "Mitotic & Meiotic ", "Mitotic", "Meiotic", "Paternal","Maternal")

Chroms_Seg_Par_Combined <- ggplot(data=Newfile, aes(x=Origin_Type,fill = factor(Origin, levels = Order, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7) + 
  scale_x_discrete(drop=F, expand = expansion(add = 1)) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), strip.background = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.spacing.x = unit(0.2, "lines"), strip.placement = "outside",
        legend.position = "none") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(breaks = Order) +
  scale_fill_manual(values = c("Maternal" = color_Maternal, "Paternal" = color_Paternal, "Meiotic" = color_Meiotic, 
                               "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  facet_grid(~ Aff_chr, switch = "x", drop = F, scales = "free_x") +  
  ylab("Genomic aberrations (#)") +
  xlab("Chromosome") 
Chroms_Seg_Par_Combined

ggsave(plot = Chroms_Seg_Par_Combined, width = 11, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Chroms_Seg_Par4.pdf")
ggsave(plot = Chroms_Seg_Par_Combined, width = 11, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Chroms_Seg_Par4.png")


###Combine ParOr, SegOr, GainNeutralLoss
Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")

Newfile <- rbind(Affected_Chromosomes, Affected_Chromosomes, Affected_Chromosomes)

OriginCol <- data.frame(Origin = c(as.character(Affected_Chromosomes$Parental_error_origin), as.character(Affected_Chromosomes$Seg_origin2)
                                   , as.character(Affected_Chromosomes$Loss_gain_or_neutral)))

Newfile <- cbind(Newfile, OriginCol) 

Newfile$Origin_Type <- c(rep("1_Parental origin", times = 39), rep("2_Segregational origin", times = 39), rep("3_GainNeutralLoss", times = 39))

Newfile$Aff_chr <- factor(Newfile$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Newfile$Origin <- factor(Newfile$Origin, levels = c("Loss","Neutral","Gain",NA, "Mitotic & Meiotic ", "Mitotic", 
                                                    "Meiotic", "Paternal","Maternal"))

Order <- c("Loss","Neutral","Gain", NA, "Mitotic & Meiotic ", "Mitotic", "Meiotic", "Paternal","Maternal")

Chroms_Seg_Par_Combined <- ggplot(data=Newfile, aes(x=Origin_Type,fill = factor(Origin, levels = Order, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7) + 
  scale_x_discrete(drop=F, expand = expansion(add = 1)) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), strip.background = element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), panel.spacing.x = unit(0, "null"), strip.placement = "outside") + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(breaks = Order) +
  scale_fill_manual(values = c("Maternal" = color_Maternal, "Paternal" = color_Paternal, "Meiotic" = color_Meiotic, 
                               "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA,
                               "Gain" = color_Gain, "Neutral" = color_Neutral, "Loss" = color_Loss)) +
  facet_grid(~ Aff_chr, switch = "x", drop = F, scales = "free_x") +  
  ylab("Genomic aberrations (#)") +
  xlab("Chromosome") 
Chroms_Seg_Par_Combined



###All chroms for MII, MI, Mitotic 
{
Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")

Ordering3 <- c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic I","Meiotic II")

Affected_Chromosomes$Aff_chr <- factor(Affected_Chromosomes$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Chrom_vs_SegOrigin <- ggplot(data=Affected_Chromosomes, aes(x=Aff_chr,fill = factor(Seg_origin, levels = Ordering3, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7, aes(y = (..count..)/sum(..count..))) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.151)) + 
  scale_fill_discrete(breaks = Ordering3) +
  scale_fill_manual(values = c("Meiotic II" = color_Meiotic_II, "Meiotic I" = color_Meiotic, "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  ylab("Proportion") 
Chrom_vs_SegOrigin
}

###All chroms for Meiotic vs mitotic 
{
  Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
Ordering4 <- c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic")

Affected_Chromosomes$Aff_chr <- factor(Affected_Chromosomes$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Chrom_vs_SegOrigin2 <- ggplot(data=Affected_Chromosomes, aes(x=Aff_chr,fill = factor(Seg_origin2, levels = Ordering4, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7) + #, aes(y = (..count..)/sum(..count..))) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6)) + 
  scale_fill_discrete(breaks = Ordering4) +
  scale_fill_manual(values = c("Meiotic" = color_Meiotic, "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  ylab("Proportion") #+
  #xlab("Chromosome")
Chrom_vs_SegOrigin2

ggsave(plot = Chrom_vs_SegOrigin2, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/SegOr_Chroms.pdf")
ggsave(plot = Chrom_vs_SegOrigin2, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/SegOr_Chroms.png")


}

###All chroms for Parental_error_origin
{
  Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
Ordering5 <- c(NA,"Paternal","Maternal")

Affected_Chromosomes$Aff_chr <- factor(Affected_Chromosomes$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Chrom_vs_ParOrigin <- ggplot(data=Affected_Chromosomes, aes(x=Aff_chr,fill = factor(Parental_error_origin, levels = Ordering5, exclude = NULL))) + 
  geom_bar(stat="count",color =  "black",width=0.7) + #, aes(y = (..count..)/sum(..count..))) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), 
        legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6)) + 
  scale_fill_manual(values = c("Maternal" = color_Maternal, "Paternal" = color_Paternal, "NA" = color_NA)) +
  ylab("Proportion") #+
  #xlab("Chromosome")
Chrom_vs_ParOrigin

ggsave(plot = Chrom_vs_ParOrigin, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/ParOr_Chroms.pdf")
ggsave(plot = Chrom_vs_ParOrigin, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/ParOr_Chroms.png")


}

###SPL & RPL comparison Parental error origin
{
  
Parental_origin_plot <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
Ordering2 <- c(NA,"Paternal","Maternal")

###Position = fill 
ParOrigin_Plot3 <- ggplot(data = Parental_origin_plot, aes(x=Group, fill = factor(Parental_error_origin, levels = Ordering2, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.7) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank(), legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(breaks = Ordering2) +
  scale_fill_manual(values = c("Maternal" = color_Maternal, "Paternal" = color_Paternal, "NA" = color_NA)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") #+
  #xlab("")
ParOrigin_Plot3  

ggsave(plot = ParOrigin_Plot3, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/ParOr_SPL_v_RPL.pdf")
ggsave(plot = ParOrigin_Plot3, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/ParOr_SPL_v_RPL.png")

}

###Statistics###
#PatMat origin all aberrations
{
rm(list=ls(all=T))
  
  Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
Affected_Chromosomes_RPL <- subset(Affected_Chromosomes, Group == "RPL")
Affected_Chromosomes_SPL <- subset(Affected_Chromosomes, Group == "SPL")
table(Affected_Chromosomes_RPL$Parental_error_origin)
table(Affected_Chromosomes_SPL$Parental_error_origin)
  
data <- data.frame("Maternal" = c(16,11), "Paternal" = c(4,7), row.names = c("RPL","SPL"), stringsAsFactors = F)
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

}

###SPL & RPL comparison Segregational error origin
{
  Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
  source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")
  
Affected_Chromosomes$Seg_origin2 <- factor(Affected_Chromosomes$Seg_origin2, levels = c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic"))

Ordering3 <- c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic")

SegOrigin3_Plot <- ggplot(data=Affected_Chromosomes, aes(x=Group, fill = factor(Seg_origin2, levels =Ordering3, exclude = NULL))) + 
  geom_bar(stat="count",position = "fill", color =  "black",width=0.7) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Meiotic" = color_Meiotic, "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") #+
  #xlab("")
SegOrigin3_Plot

ggsave(plot = SegOrigin3_Plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/SegOr_SPL_v_RPL.pdf")
ggsave(plot = SegOrigin3_Plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/SegOr_SPL_v_RPL.png")

}

###Statistics###
#Mitotic/meiotic only for all aberrations 
{
rm(list=ls(all=T))
  
  Affected_Chromosomes <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
Affected_Chromosomes_RPL <- subset(Affected_Chromosomes, Group == "RPL")
Affected_Chromosomes_SPL <- subset(Affected_Chromosomes, Group == "SPL")
table(Affected_Chromosomes_RPL$Seg_origin2)
table(Affected_Chromosomes_SPL$Seg_origin2)

data <- data.frame("Mitotic" = c(7,11), "Meiotic" = c(11,5), row.names = c("RPL","SPL"), stringsAsFactors = F)
data

###Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
chisq.test(data)$expected
#Mitotic Meiotic
#RPL    5.76   10.24
#SPL    3.24    5.76
###Conclusion = chisquared approximation is appropriate

###Chi-squared With Yates' continuity correction
chisq.test(data)
#Pearson's Chi-squared test with Yates' continuity correction

#data:  data
#X-squared = 2.378, df = 1, p-value = 0.1231


#Mitotic/meiotic only for only trisomy cases

MitMei_TrisOnly <- subset(Affected_Chromosomes, Indication == "Trisomy")
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

}



#### Aberration type plotting 
rm(list=ls(all=T))

library(ggplot2)
library(tidyr)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

###Gain and loss, RPL vs SPL 
#PerAberration <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv")
PerAberration <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")

##
PerAberration_RPL <- subset(PerAberration, Group == "RPL")
PerAberration_SPL <- subset(PerAberration, Group == "SPL")
##


Order <- c("Loss","Neutral","Gain")

PerAberration$Loss_gain_or_neutral <- factor(PerAberration$Loss_gain_or_neutral, levels = c("Loss","Neutral","Gain"))

GainLossNeutral_plot <- ggplot(PerAberration, aes(x=Group,fill = factor(Loss_gain_or_neutral, levels = Order, exclude = NULL))) + 
  geom_bar(stat="count", position = "fill", color =  "black",width=0.7) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.00)) + 
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  scale_fill_manual(values = c("Loss" = color_Loss, "Neutral" = color_Neutral, "Gain" = color_Gain)) +
  ylab("Proportion") 
GainLossNeutral_plot

ggsave(plot = GainLossNeutral_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GainLossNeutral_SPL_v_RPL.pdf")
ggsave(plot = GainLossNeutral_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GainLossNeutral_SPL_v_RPL.png")


###Gain and loss, per chromosome 
PerAberration <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv")

Order <- c("Loss","Neutral","Gain")

PerAberration$Loss_gain_or_neutral <- factor(PerAberration$Loss_gain_or_neutral, levels = c("Loss","Neutral","Gain"))

PerAberration$Aff_chr <- factor(PerAberration$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

GainLossNeutral_Chrom_plot <- ggplot(PerAberration, aes(x=Aff_chr,fill = factor(Loss_gain_or_neutral, levels = Order, exclude = NULL))) + 
  geom_bar(stat="count", color =  "black",width=0.7) +
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  scale_x_discrete(drop=F) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6)) + 
  scale_fill_manual(values = c("Loss" = color_Loss, "Neutral" = color_Neutral, "Gain" = color_Gain)) +
  ylab("Proportion") 
GainLossNeutral_Chrom_plot

ggsave(plot = GainLossNeutral_Chrom_plot, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GainLossNeutral_Chroms.pdf")
ggsave(plot = GainLossNeutral_Chrom_plot, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GainLossNeutral_Chroms.png")


###Gain and loss vs par origin
PerAberration <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv")

Order <- c("Loss","Neutral","Gain")

PerAberration$Loss_gain_or_neutral <- factor(PerAberration$Loss_gain_or_neutral, levels = c("Loss","Neutral","Gain"))

GainLossNeutral_ParOr_plot <- ggplot(data = subset(PerAberration, !is.na(Parental_error_origin)), aes(x=Parental_error_origin,fill = factor(Loss_gain_or_neutral, levels = Order))) + 
  geom_bar(stat="count", position = "fill", color =  "black",width=0.7) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.00)) + 
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  scale_fill_manual(values = c("Loss" = color_Loss, "Neutral" = color_Neutral, "Gain" = color_Gain)) +
  ylab("Proportion") 
GainLossNeutral_ParOr_plot

ggsave(plot = GainLossNeutral_ParOr_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GainLossNeutral_ParOr.pdf")
ggsave(plot = GainLossNeutral_ParOr_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GainLossNeutral_ParOr.png")


###Affected region plotting
rm(list=ls(all=T))

library(ggplot2)
library(tidyr)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

###Aff region, RPL vs SPL 

#PerAberration <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv")
PerAberration <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")

Order <- c("Segmental","Chromosomal","Genome wide")

PerAberration$Aff_region <- factor(PerAberration$Aff_region, levels = c("Segmental","Chromosomal","Genome wide"))

Aff_Region_plot <- ggplot(PerAberration, aes(x=Group,fill = factor(Aff_region, levels = Order, exclude = NULL))) + 
  geom_bar(stat="count", position = "fill", color =  "black",width=0.7) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.00)) + 
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  scale_fill_manual(values = c("Segmental" = color_Segmental, "Chromosomal" = color_Chromosomal, "Genome wide" = color_GenomeWide)) +
  ylab("Proportion") 
Aff_Region_plot

ggsave(plot = Aff_Region_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Aff_Region_SPL_v_RPL.pdf")
ggsave(plot = Aff_Region_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Aff_Region_SPL_v_RPL.png")


###Aff region, per chromosome 
PerAberration <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerAberration.csv")

Order <- c("Segmental","Chromosomal","Genome wide")

PerAberration$Aff_region <- factor(PerAberration$Aff_region, levels = c("Segmental","Chromosomal","Genome wide"))

PerAberration$Aff_chr <- factor(PerAberration$Aff_chr, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","GW"))

Aff_Region_Chrom_plot <- ggplot(PerAberration, aes(x=Aff_chr,fill = factor(Aff_region, levels = Order, exclude = NULL))) + 
  geom_bar(stat="count", color =  "black",width=0.7) +
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
        legend.position = "none") + 
  scale_x_discrete(drop=F) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6)) + 
  scale_fill_manual(values = c("Segmental" = color_Segmental, "Chromosomal" = color_Chromosomal, "Genome wide" = color_GenomeWide)) +
  ylab("Proportion") 
Aff_Region_Chrom_plot

ggsave(plot = Aff_Region_Chrom_plot, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Aff_Region_Chroms.pdf")
ggsave(plot = Aff_Region_Chrom_plot, width = 10, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Aff_Region_Chroms.png")


###Aff region vs par origin
PerAberration <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerAberration_Source_data.csv", sep=";")

Order <- c("Segmental","Chromosomal","Genome wide")

PerAberration$Aff_region <- factor(PerAberration$Aff_region, levels = c("Segmental","Chromosomal","Genome wide"))

Aff_Region_ParOr_plot <- ggplot(data = subset(PerAberration, !is.na(Parental_error_origin)), aes(x=Parental_error_origin,fill = factor(Aff_region, levels = Order))) + 
  geom_bar(stat="count", position = "fill", color =  "black",width=0.7) + 
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.00)) + 
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  scale_fill_manual(values = c("Segmental" = color_Segmental, "Chromosomal" = color_Chromosomal, "Genome wide" = color_GenomeWide)) +
  ylab("Proportion") 
Aff_Region_ParOr_plot

ggsave(plot = Aff_Region_ParOr_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Aff_Region_ParOr.pdf")
ggsave(plot = Aff_Region_ParOr_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Aff_Region_ParOr.png")




