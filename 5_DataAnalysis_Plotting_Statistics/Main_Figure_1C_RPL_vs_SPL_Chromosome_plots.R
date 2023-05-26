###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: Firstly, to make barplots that compare multiple factors (parental origin, segregational origin, aberration size, 
#                 and copy number change) between sporadic pregnancy loss (SPL) and recurrent pregnancy loss (RPL)
#                 subcohorts, analyzed by genome haplarithmisis. Secondly, to make a chromosome specific barplot showing 
#                 parental and segregational origin of aberrations detected by genome haplarithmisis. 

# input: PL_PerAberration_Source_data.csv, this file contains information on each individual aberration detected by genome 
#        haplarithmisis in miscarried POCs included in the pregnancy loss study. 

# output: Main figure 1 C; parental origin (maternal or paternal), segregational origin (meiotic, mitotic, or combination), 
#         aberration size (segmental, chromosomal, genome wide), and copy number change (loss, gain, neutral) between 
#         sporadic pregnancy loss (SPL) and recurrent pregnancy loss (RPL) subcohorts. 
#         Main figure 1 D; chromosome specific barplots showing parental (maternal or paternal) and segregational 
#         origin (meiotic, mitotic, or combination)

###########################################################################################################################

rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

###All chroms for segregational and parental origin combined. 

Affected_Chromosomes <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")

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
       "/Projects/PregnancyLoss/Figures/Chroms_Seg_Par4.pdf")
ggsave(plot = Chroms_Seg_Par_Combined, width = 11, height = 5,
       "/Projects/PregnancyLoss/Figures/Chroms_Seg_Par4.png")


###All chroms for Meiotic vs mitotic 

  Affected_Chromosomes <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
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
       "/Projects/PregnancyLoss/Figures/SegOr_Chroms.pdf")
ggsave(plot = Chrom_vs_SegOrigin2, width = 10, height = 5,
       "/Projects/PregnancyLoss/Figures/SegOr_Chroms.png")




###All chroms for Parental_error_origin
Affected_Chromosomes <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
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
       "/Projects/PregnancyLoss/Figures/ParOr_Chroms.pdf")
ggsave(plot = Chrom_vs_ParOrigin, width = 10, height = 5,
       "/Projects/PregnancyLoss/Figures/ParOr_Chroms.png")




###SPL & RPL comparison Parental error origin
Parental_origin_plot <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
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
  ylab("Proportion") 
ParOrigin_Plot3  

ggsave(plot = ParOrigin_Plot3, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/ParOr_SPL_v_RPL.pdf")
ggsave(plot = ParOrigin_Plot3, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/ParOr_SPL_v_RPL.png")



###SPL & RPL comparison Segregational error origin
Affected_Chromosomes <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

Affected_Chromosomes$Seg_origin2 <- factor(Affected_Chromosomes$Seg_origin2, levels = c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic"))

Ordering3 <- c(NA,"Mitotic & Meiotic ","Mitotic","Meiotic")

SegOrigin3_Plot <- ggplot(data=Affected_Chromosomes, aes(x=Group, fill = factor(Seg_origin2, levels =Ordering3, exclude = NULL))) + 
  geom_bar(stat="count",position = "fill", color =  "black",width=0.7) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank(), legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Meiotic" = color_Meiotic, "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") 
  #xlab("")
SegOrigin3_Plot

ggsave(plot = SegOrigin3_Plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/SegOr_SPL_v_RPL.pdf")
ggsave(plot = SegOrigin3_Plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/SegOr_SPL_v_RPL.png")




#### Aberration type plotting 
rm(list=ls(all=T))

library(ggplot2)
library(tidyr)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

###Gain and loss, RPL vs SPL 
PerAberration <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")

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
       "/Projects/PregnancyLoss/Figures/GainLossNeutral_SPL_v_RPL.pdf")
ggsave(plot = GainLossNeutral_plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/GainLossNeutral_SPL_v_RPL.png")


###Affected region plotting
rm(list=ls(all=T))

library(ggplot2)
library(tidyr)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

###Aff region, RPL vs SPL 

PerAberration <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")

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
       "/Projects/PregnancyLoss/Figures/Aff_Region_SPL_v_RPL.pdf")
ggsave(plot = Aff_Region_plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Aff_Region_SPL_v_RPL.png")










