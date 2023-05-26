###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To produce barplots for the pregnancy loss study

# input: PL_PerPOC_Source_data.csv, this file contains information on each miscarried POC included in the pregnancy loss study. 

# output: Barplots that were used in the pregnancy loss paper. 
#         Main figure 1A; barplots used to indicate normal vs. abnormal number (or percentage) of PL cases assessed by 
#         conventional karyotyping or genome haplarithmisis. As well as the distribution of specific aberration types for the
#         PL cases found to be genetically abnormal. 

###########################################################################################################################


###Barplots ConvKaryo and SNPHapla indications
rm(list=ls(all=T))

library(ggplot2)
library(scales)

###Plot colors for indications
source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

###Conv Karyo Indications normal and abnormal
#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")

#Demographic_data_TOTAL              <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Source_data\ /20221027_Miscarr_PerPOC.csv", sep=";")
Demographic_data_ConvKaryo          <- subset(Demographic_data_TOTAL, Group == "ConvKaryo" | Group == "SNPHapla")

Demographic_data_ConvKaryo_Normal <- subset(Demographic_data_ConvKaryo, Indication_ConvKaryo == "Normal")
Demographic_data_ConvKaryo_Abnormal <- subset(Demographic_data_ConvKaryo, Indication_ConvKaryo == "Abnormal")


Demographic_data_ConvKaryo$Indication_ConvKaryo <- factor(Demographic_data_ConvKaryo$Indication_ConvKaryo, levels = c("Normal","Abnormal"))
Order <- c("Abnormal","Normal")

Demographic_data_ConvKaryo$Group <- gsub("SNPHapla", "ConvKaryo",Demographic_data_ConvKaryo$Group)

Barplot_ConvKaryo_Normal_Abnormal <- ggplot(data = Demographic_data_ConvKaryo, aes(x=Group, fill = factor(Indication_ConvKaryo, levels = Order, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.3) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label = paste0("\n", ..count..," ","(",percent(..count../1745, accuracy = 0.1),")"), x = 1.35), 
            position = position_fill(vjust = 0.5), stat = "count", size = 3) +
  scale_fill_manual(values = c("Normal" = color_Normal, "Abnormal" = color_Abnormal)) +
  ylab("Proportion") +
  xlab("")
Barplot_ConvKaryo_Normal_Abnormal

ggsave(plot = Barplot_ConvKaryo_Normal_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_ConvKaryo_Normal_Abnormal.png")
ggsave(plot = Barplot_ConvKaryo_Normal_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_ConvKaryo_Normal_Abnormal.pdf")
ggsave(plot = Barplot_ConvKaryo_Normal_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_ConvKaryo_Normal_Abnormal.svg")


###Conv Karyo Indications only abnormal individual classification
#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
Demographic_data_ConvKaryo_Abnormal     <- subset(Demographic_data_TOTAL, Indication_ConvKaryo == "Abnormal")

Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly <- gsub("[?]", "",Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly)
Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly <- gsub("Structural rearrangements", "Structural rearrangement",Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly)

Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly <- factor(Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly, levels = c("Autosomal monosomy","Structural rearrangement","Others","Combined abnormalities","Gonosomal aneuploidy","Polyploidy","Autosomal trisomy"))
Order <- c("Autosomal monosomy","Structural rearrangement","Others","Combined abnormalities","Gonosomal aneuploidy","Polyploidy","Autosomal trisomy")

Demographic_data_ConvKaryo_Abnormal$Group <- gsub("SNPHapla", "ConvKaryo",Demographic_data_ConvKaryo_Abnormal$Group)

Barplot_ConvKaryo_Abnormal <- ggplot(data = Demographic_data_ConvKaryo_Abnormal, aes(x=Group, fill = factor(Classification_ConvKaryo_Poly, levels = Order, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.5) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label = paste0("\n", ..count..," ","(",percent(..count../879, accuracy = 0.1),")"), x = 1.35), 
            position = position_fill(vjust = 0.5), stat = "count", size = 3) +
  scale_fill_manual(values = c("Autosomal monosomy" = color_Aut_mon,"Structural rearrangement" = color_Str,"Others" = color_Oth,
                               "Combined abnormalities" = color_Com_ab,"Gonosomal aneuploidy" = color_Gon,"Polyploidy" = color_Tetra,
                               "Autosomal trisomy" = color_Aut_tri)) +
  ylab("Proportion") +
  xlab("")
Barplot_ConvKaryo_Abnormal

ggsave(plot = Barplot_ConvKaryo_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_ConvKaryo_Abnormal.png")
ggsave(plot = Barplot_ConvKaryo_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_ConvKaryo_Abnormal.pdf")
ggsave(plot = Barplot_ConvKaryo_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_ConvKaryo_Abnormal.svg")


###SNPHapla Indications normal and abnormal
#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
Demographic_data_SNPHapla     <- subset(Demographic_data_TOTAL, Group == "SNPHapla" & Indication_SNPHapla != "Excluded")

Demographic_data_SNPHapla$Indication_SNPHapla <- factor(Demographic_data_SNPHapla$Indication_SNPHapla, levels = c("Normal","Abnormal"))
Order <- c("Abnormal","Normal")
Demographic_data_SNPHapla$Group <- factor(Demographic_data_SNPHapla$Group, levels = "SNPHapla")

Barplot_SNPHapla_Normal_Abnormal <- ggplot(data = Demographic_data_SNPHapla, aes(x=Group, fill = factor(Indication_SNPHapla, levels = Order, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.3) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label = paste0("\n", ..count..," ","(",percent(..count../94, accuracy = 0.1),")"), x = 1.35), 
            position = position_fill(vjust = 0.5), stat = "count", size = 3) +
  scale_fill_manual(values = c("Normal" = color_Normal, "Abnormal" = color_Abnormal)) +
  ylab("Proportion") +
  xlab("")
Barplot_SNPHapla_Normal_Abnormal

ggsave(plot = Barplot_SNPHapla_Normal_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Normal_Abnormal.png")
ggsave(plot = Barplot_SNPHapla_Normal_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Normal_Abnormal.pdf")
ggsave(plot = Barplot_SNPHapla_Normal_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Normal_Abnormal.svg")

###SNP Hapla Indications only abnormal individual classification
#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
Demographic_data_SNPHapla_Abnormal     <- subset(Demographic_data_TOTAL, Indication_SNPHapla == "Abnormal")

Demographic_data_SNPHapla_Abnormal$Classification_SNPHapla <- factor(Demographic_data_SNPHapla_Abnormal$Classification_SNPHapla, levels = c("Autosomal monosomy","Polyploidy","Structural rearrangement","Gonosomal aneuploidy","Combined abnormalities","Autosomal trisomy"))
Order <- c("Autosomal monosomy","Polyploidy","Structural rearrangement","Gonosomal aneuploidy","Combined abnormalities","Autosomal trisomy")
Demographic_data_SNPHapla_Abnormal$Group <- factor(Demographic_data_SNPHapla_Abnormal$Group, levels = "SNPHapla")

Barplot_SNPHapla_Abnormal <- ggplot(data = Demographic_data_SNPHapla_Abnormal, aes(x=Group, fill = factor(Classification_SNPHapla, levels = Order, exclude = NULL))) +
  geom_bar(stat="count", position = 'fill', color =  "black",width=0.5) + 
  scale_x_discrete(drop=FALSE) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label = paste0("\n", ..count..," ","(",percent(..count../33, accuracy = 0.1),")"), x = 1.4), 
            position = position_fill(vjust = 0.5), stat = "count", size = 3) +
  scale_fill_manual(values = c("Autosomal monosomy" = color_Aut_mon,"Structural rearrangement" = color_Str,"Others" = color_Oth,
                               "Combined abnormalities" = color_Com_ab,"Gonosomal aneuploidy" = color_Gon,"Polyploidy" = color_Tetra,
                               "Autosomal trisomy" = color_Aut_tri)) +
  ylab("Proportion") +
  xlab("")
Barplot_SNPHapla_Abnormal

ggsave(plot = Barplot_SNPHapla_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Abnormal.png")
ggsave(plot = Barplot_SNPHapla_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Abnormal.pdf")
ggsave(plot = Barplot_SNPHapla_Abnormal, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/Barplot_SNPHapla_Abnormal.svg")



