###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To produce a combination of half-violin and scatter plots comparing maternal, paternal and gestational age
#                 for genetically normal and abnormal PLs.  

# input: PL_PerPOC_Source_data.csv, this file contains information on each miscarried POC included in the pregnancy loss study. 

# output: Main figure 1B; Maternal, paternal, and gestational age comparison between genetically normal and abnormal cases for 
#         conventional karyotyping as well as haplarithmisis. 

###########################################################################################################################

rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)
library(ggsignif)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")
source("/Projects/PregnancyLoss/Figures/functions/Flat_Violin.R")


###Preprocessing of data files
  PerFamily <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
  PerFamily <- subset(PerFamily, Group == "ConvKaryo" | Group == "SNPHapla")
  
  ConvKaryoConversion <- subset(PerFamily, Group == "SNPHapla")
  
  ConvKaryoConversion$Group <- c(rep("ConvKaryo", times = 114))
  
  All_Ages <- data.frame(All_Ages = c(ConvKaryoConversion$Maternal.age, ConvKaryoConversion$Paternal.age, ConvKaryoConversion$Gestational.age))
  ConvKaryoConversion <- rbind(ConvKaryoConversion, ConvKaryoConversion, ConvKaryoConversion)
  ConvKaryoConversion$Age_Type <- c(rep("Maternal age", times = 114), rep("Paternal age", times = 114), rep("Gestational age", times = 114))
  
  IndicationCol <- subset(ConvKaryoConversion, Age_Type == "Maternal age")
  IndicationCol2 <- subset(ConvKaryoConversion, Age_Type == "Paternal age")
  IndicationCol3 <- subset(ConvKaryoConversion, Age_Type == "Gestational age")
  Indication <- data.frame(Indication = c(as.character(IndicationCol$Indication_ConvKaryo), 
                                          as.character(IndicationCol2$Indication_ConvKaryo), 
                                          as.character(IndicationCol3$Indication_ConvKaryo)))
  
  IndicationConvKaryoConversion <- cbind(ConvKaryoConversion, All_Ages, Indication)
  
  
  All_Ages <- data.frame(All_Ages = c(PerFamily$Maternal.age, PerFamily$Paternal.age, PerFamily$Gestational.age))
  PerFamily <- rbind(PerFamily, PerFamily, PerFamily)
  PerFamily$Age_Type <- c(rep("Maternal age", times = 1745), rep("Paternal age", times = 1745), rep("Gestational age", times = 1745))
  
  IndicationCol4 <- subset(PerFamily, Group == "SNPHapla" & Age_Type == "Maternal age")
  IndicationCol5 <- subset(PerFamily, Group == "ConvKaryo" & Age_Type == "Paternal age")
  IndicationCol6 <- subset(PerFamily, Group == "SNPHapla" & Age_Type == "Gestational age")
  IndicationCol7 <- subset(PerFamily, Group == "ConvKaryo" & Age_Type == "Gestational age")
  IndicationCol8 <- subset(PerFamily, Group == "SNPHapla" & Age_Type == "Gestational age")
  IndicationCol9 <- subset(PerFamily, Group == "ConvKaryo" & Age_Type == "Gestational age")
  
  Indication2 <- data.frame(Indication = c(as.character(IndicationCol4$Indication_SNPHapla), 
                                           as.character(IndicationCol5$Indication_ConvKaryo), 
                                           as.character(IndicationCol6$IndicationSNPHapla),
                                           as.character(IndicationCol7$IndicationConvKaryo),
                                           as.character(IndicationCol8$IndicationSNPHapla),
                                           as.character(IndicationCol9$IndicationConvKaryo)))
  
  PerFamily <- cbind(PerFamily, All_Ages, Indication2)
  
  PerFamily <- rbind(PerFamily, IndicationConvKaryoConversion)
  PerFamily <- subset(PerFamily, select = -c(Maternal.age, Paternal.age, Gestational.age))
  
  
  
  ###Total without gestational age 
  PerFamilyPatMat <- PerFamily
  
  PerFamilyPatMat$Group <- gsub('SNPHapla', 'SNP haplotyping \n (n = 91)',PerFamilyPatMat$Group)
  PerFamilyPatMat$Group <- gsub('ConvKaryo', 'Conventional karyotyping \n (n = 1745)',PerFamilyPatMat$Group)
  
  PerFamilyPatMat <- subset(PerFamilyPatMat, Indication != "Excluded")
  PerFamilyPatMat <- subset(PerFamilyPatMat, Age_Type != "Gestational age")
  
  PerFamilyPatMat$Age_Type               <- factor(PerFamilyPatMat$Age_Type, levels = c("Maternal age","Paternal age"))
  PerFamilyPatMat$Indication  <- factor(PerFamilyPatMat$Indication, levels = c("Normal","Abnormal"))


###ConvKaryo individual plot 
###Violin
  PerFamilyPatMatConvKaryo <- subset(PerFamilyPatMat, Group == "Conventional karyotyping \n (n = 1745)")
  
  TOTAL_AgesPlot_PatMat_ConvKaryo <- 
    ggplot(data=subset(PerFamilyPatMatConvKaryo, !is.na(All_Ages)), aes(x=Age_Type, y=All_Ages, fill=Indication, color = Indication)) + 
    scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
    scale_color_manual(values = c("black","black")) +
    geom_flat_violin(width = 1, position = position_dodge(width = 0.7), alpha = 0.5) +
    geom_dotplot(stackdir = "up", binaxis = "y", binwidth = 0.08, alpha = 0.5, position = position_jitterdodge(0.2), ) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) +
    theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), legend.title = element_blank(), 
          axis.title.x = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
          legend.key = element_rect(fill = "transparent"), axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13)) +
    ggtitle("Conventional karyotyping \n (n = 1745)") +
    ylim(10,63) +
    ylab("Age (years)") +
    geom_signif(y_position = c(62, 62), xmin = c(0.82, 1.82), xmax = c(1.18,2.18), annotation = c("p = 6.7e-5", "p = 4.6e-3"), 
                tip_length = 0.01, textsize = 4)
  TOTAL_AgesPlot_PatMat_ConvKaryo
  

###SNPHapla individual plot 
###Violin
  PerFamilyPatMatSNPHapla <- subset(PerFamilyPatMat, Group == "SNP haplotyping \n (n = 91)")
  
  TOTAL_AgesPlot_PatMat_SNPHapla <- 
    ggplot(data=subset(PerFamilyPatMatSNPHapla, !is.na(All_Ages)), aes(x=Age_Type, y=All_Ages, fill=Indication, color = Indication)) + 
    scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
    scale_color_manual(values = c("black","black")) +
    geom_flat_violin(width = 1, position = position_dodge(width = 0.7), alpha = 0.5) +
    geom_dotplot(stackdir = "up", binaxis = "y", binwidth = 0.08, alpha = 0.5, position = position_jitterdodge(0.2), ) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) + 
    theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), legend.title = element_blank(), 
          axis.title.x = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
          legend.key = element_rect(fill = "transparent"), axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13)) +
    ggtitle("Genome haplarithmisis \n (n = 91)") +
    ylim(10,63) +
    ylab("Age (years)") +
    geom_signif(y_position = c(62, 62), xmin = c(0.82, 1.82), xmax = c(1.18,2.18), annotation = c("P = 1.81e-7", "P = 7.64e-6"), 
                tip_length = 0.01, textsize = 4)
  TOTAL_AgesPlot_PatMat_SNPHapla

  
  
###Gestational age ConvKaryo
###Violin
  PerFamilyGest <- PerFamily
  
  PerFamilyGest$Group <- gsub('SNPHapla', 'SNP Haplotyping \n (n = 91)',PerFamilyGest$Group)
  PerFamilyGest$Group <- gsub('ConvKaryo', 'Conventional karyotyping \n (n = 1745)',PerFamilyGest$Group)
  
  PerFamilyGest <- subset(PerFamilyGest, Indication != "Excluded")
  PerFamilyGest <- subset(PerFamilyGest, Age_Type != "Paternal age" & Age_Type !="Maternal age")
  
  PerFamilyGest$Age_Type               <- factor(PerFamilyGest$Age_Type, levels = c("Gestational age"))
  PerFamilyGest$Indication  <- factor(PerFamilyGest$Indication, levels = c("Normal","Abnormal"))
  
  PerFamilyPatMatConvKaryoGest <- subset(PerFamilyGest, Group == "Conventional karyotyping \n (n = 1745)")
  
  TOTAL_AgesPlot_PatMat_ConvKaryo_Gest <- 
    ggplot(data=subset(PerFamilyPatMatConvKaryoGest, !is.na(All_Ages)), aes(x=Age_Type, y=All_Ages, fill=Indication, color = Indication)) + 
    scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
    scale_color_manual(values = c("black","black")) +
    geom_flat_violin(width = 1, position = position_dodge(width = 0.7), alpha = 0.5) +
    geom_dotplot(stackdir = "up", binaxis = "y", binwidth = 0.08, alpha = 0.5, position = position_jitterdodge(0.2), ) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) +  
    theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), legend.title = element_blank(), 
          axis.title.x = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = "transparent"), axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13)) +
    ggtitle("Conventional karyotyping \n (n = 1745)") +
    ylim(2,30) +
    ylab("Age (weeks)") +
    geom_signif(y_position = 30, xmin = 0.82, xmax = 1.18, annotation = "NS", tip_length = 0.01, textsize = 4) 
  TOTAL_AgesPlot_PatMat_ConvKaryo_Gest


###Gestational age SNPHapla
###Violin
  PerFamilyPatMatSNPHaplaGest <- subset(PerFamilyGest, Group == "SNP Haplotyping \n (n = 91)")
  
  TOTAL_AgesPlot_PatMat_SNPHapla_Gest <- 
    ggplot(data=subset(PerFamilyPatMatSNPHaplaGest, !is.na(All_Ages)), aes(x=Age_Type, y=All_Ages, fill=Indication, color = Indication)) + 
    scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
    scale_color_manual(values = c("black","black")) +
    geom_flat_violin(width = 1, position = position_dodge(width = 0.7), alpha = 0.5) +
    geom_dotplot(stackdir = "up", binaxis = "y", binwidth = 0.08, alpha = 0.5, position = position_jitterdodge(0.2), ) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) +
    theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), legend.title = element_blank(), 
          axis.title.x = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5),
          legend.key = element_rect(fill = "transparent"), axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13)) +
    ggtitle("Genome Haplarithmisis \n (n = 91)") +
    ylim(2,30) +
    ylab("Age (weeks)") +
    geom_signif(y_position = 30, xmin = 0.82, xmax = 1.18, annotation = "NS", tip_length = 0.01, textsize = 4) 
  TOTAL_AgesPlot_PatMat_SNPHapla_Gest


###Save all plots same height and width
  ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/PatMatAge_ConvKaryo_flat_Violin.pdf")
  ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/PatMatAge_ConvKaryo_flat_Violin.png")
  
  ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/PatMatAge_SNPHapla_flat_Violin.pdf")
  ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/PatMatAge_SNPHapla_flat_Violin.png")
  
  ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo_Gest, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/GestAge_ConvKaryo_flat_Violin.pdf")
  ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo_Gest, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/GestAge_ConvKaryo_flat_Violin.png")
  
  ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla_Gest, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/GestAge_SNPHapla_flat_Violin.pdf")
  ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla_Gest, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/GestAge_SNPHapla_flat_Violin.png")
