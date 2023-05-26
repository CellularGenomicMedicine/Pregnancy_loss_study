###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To produce barplots showcasing the ratio of miscarried POCs that were found genetically normal and abnormal 
#                 by genome haplarithmisis, divided into sporadic pregnancy loss (SPL) and recurrent pregnancy loss (RPL) 
#                 subcohorts. 

# input: PL_PerPOC_Source_data.csv, this file contains information on each miscarried POC included in the pregnancy loss study.

# output: Extended data figure 2 A; normal vs. abnormal ratio for SPL, 2 B; normal vs. abnormal ratio for RPL, 2 C; the 
#         number of pregnancy losses (1-5) and parental origin ()

###########################################################################################################################

rm(list=ls(all=T))

library(ggplot2)
library(scales)

###SPL
###Plot colors for indications
source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
Demographic_data_SNPHapla     <- subset(Demographic_data_TOTAL, Group == "SNPHapla" & Indication_SNPHapla != "Excluded")
Demographic_data_SNPHapla_SPL     <- subset(Demographic_data_SNPHapla, Group2 == "SPL")


Demographic_data_SNPHapla_SPL$Indication_SNPHapla <- factor(Demographic_data_SNPHapla_SPL$Indication_SNPHapla, levels = c("Normal","Abnormal"))
Order <- c("Abnormal","Normal")
Demographic_data_SNPHapla_SPL$Group <- factor(Demographic_data_SNPHapla_SPL$Group, levels = "SNPHapla")

Barplot_SNPHapla_Normal_Abnormal_SPL <- ggplot(data = Demographic_data_SNPHapla_SPL, aes(x=Group, fill = factor(Indication_SNPHapla, levels = Order, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.3) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label = paste0("\n", (..count..*2)," ","(",percent((..count..*2)/84, accuracy = 0.1),")"), x = 1.35), 
            position = position_fill(vjust = 0.5), stat = "count", size = 3) +
  scale_fill_manual(values = c("Normal" = color_Normal, "Abnormal" = color_Abnormal)) +
  ylab("Proportion") +
  xlab("SPL")
Barplot_SNPHapla_Normal_Abnormal_SPL

ggsave(plot = Barplot_SNPHapla_Normal_Abnormal_SPL, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Normal_Abnormal_SPL.png")
ggsave(plot = Barplot_SNPHapla_Normal_Abnormal_SPL, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Normal_Abnormal_SPL.pdf")


###RPL
rm(list=ls(all=T))

library(ggplot2)
library(scales)

###Plot colors for indications
source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
Demographic_data_SNPHapla     <- subset(Demographic_data_TOTAL, Group == "SNPHapla" & Indication_SNPHapla != "Excluded")
Demographic_data_SNPHapla_RPL     <- subset(Demographic_data_SNPHapla, Group2 == "RPL")

Demographic_data_SNPHapla_RPL <- rbind(Demographic_data_SNPHapla_RPL, Demographic_data_SNPHapla_RPL)
Demographic_data_SNPHapla_RPL$Tissue <- c(rep("EM", times = 52), rep("CV", times = 52))

Demographic_data_SNPHapla_RPL[9,14] <- "Normal"

Demographic_data_SNPHapla_RPL$Indication_SNPHapla <- factor(Demographic_data_SNPHapla_RPL$Indication_SNPHapla, levels = c("Normal","Abnormal"))
Order <- c("Abnormal","Normal")
Demographic_data_SNPHapla_RPL$Group <- factor(Demographic_data_SNPHapla_RPL$Group, levels = "SNPHapla")

Barplot_SNPHapla_Normal_Abnormal_RPL <- ggplot(data = Demographic_data_SNPHapla_RPL, aes(x=Group, fill = factor(Indication_SNPHapla, levels = Order, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.3) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label = paste0("\n", (..count..)," ","(",percent(((..count..))/104, accuracy = 0.1),")"), x = 1.35), 
            position = position_fill(vjust = 0.5), stat = "count", size = 3) +
  scale_fill_manual(values = c("Normal" = color_Normal, "Abnormal" = color_Abnormal)) +
  ylab("Proportion") +
  xlab("RPL")
Barplot_SNPHapla_Normal_Abnormal_RPL

ggsave(plot = Barplot_SNPHapla_Normal_Abnormal_RPL, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Normal_Abnormal_RPL.png")
ggsave(plot = Barplot_SNPHapla_Normal_Abnormal_RPL, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Normal_Abnormal_RPL.pdf")


rm(list=ls(all=T))

library(ggplot2)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

PerFamily <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")

NumPLs_SNPHapla <- subset(PerFamily, Group == "SNPHapla" & Indication_SNPHapla == "Abnormal")

Seg_Or_Order <- c("Mitotic", "Mitotic & Meiotic ", "Meiotic")
Par_Or_Order <- c("Paternal","Maternal & Paternal ", "Maternal")

###Position = "fill"
NumPL_Seg_plot <- ggplot(data = subset(NumPLs_SNPHapla, !is.na(Seg_origin2)), aes(x = Total_PLs, fill = factor(Seg_origin2, levels = Seg_Or_Order))) +
  geom_bar(stat="count", color =  "black", width=0.7, position = "fill") +
  scale_fill_manual(values = c("Meiotic" = color_Meiotic, "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("Number of pregnancy losses")
NumPL_Seg_plot

ggsave(plot = NumPL_Seg_plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/NumPL_SegOr_fill.pdf")
ggsave(plot = NumPL_Seg_plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/NumPL_SegOr_fill.png")


###Parental origin
###Position = "fill"
NumPL_PatMat_plot <- ggplot(data = subset(NumPLs_SNPHapla, !is.na(Parental_origin)), aes(x = Total_PLs, fill = factor(Parental_origin, levels = Par_Or_Order))) +
  geom_bar(stat="count", color =  "black", width=0.7, position = "fill") +
  scale_fill_manual(values = c("Maternal" = color_Maternal, "Paternal" = color_Paternal, "Maternal & Paternal " = color_Maternal_and_Paternal)) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("Number of pregnancy losses")
NumPL_PatMat_plot

ggsave(plot = NumPL_PatMat_plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/NumPL_ParOr_fill.pdf")
ggsave(plot = NumPL_PatMat_plot, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/NumPL_ParOr_fill.png")
