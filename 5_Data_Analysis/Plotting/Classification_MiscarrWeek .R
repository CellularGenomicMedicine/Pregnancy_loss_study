rm(list=ls(all=T))

library(ggplot2)
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")

MiscarrWeekBinning <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv", stringsAsFactors = F)
MiscarrWeekBinning_SNPHapla <- subset(MiscarrWeekBinning, Group == "SNPHapla" & Indication_SNPHapla != "Excluded")

MiscarrWeekBinning_SNPHapla$Week_bin <- as.numeric(cut_width(MiscarrWeekBinning_SNPHapla$Gestational.age, 1))

MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 1] <- "4-5"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 2] <- "6-7"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 3] <- "6-7"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 4] <- "8-9"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 5] <- "8-9"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 6] <- "10-13"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 7] <- "10-13"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 8] <- "10-13"
MiscarrWeekBinning_SNPHapla$Week_bin[MiscarrWeekBinning_SNPHapla$Week_bin == 9] <- "10-13"


###RPL and SPL 
MiscarrWeekBinning_SNPHapla$Week_bin <- factor(MiscarrWeekBinning_SNPHapla$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))
MiscarrWeekBinning_SNPHapla$Indication_SNPHapla <- factor(MiscarrWeekBinning_SNPHapla$Indication_SNPHapla, levels = c("Normal","Abnormal"))


MiscarrWeekBin <- ggplot(data= subset(MiscarrWeekBinning_SNPHapla, !is.na(Week_bin)), aes(fill =  Indication_SNPHapla, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Gestational age (in weeks)") +
  ggtitle("Genome haplarithmisis")
MiscarrWeekBin

ggsave(plot = MiscarrWeekBin, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MiscarrWeekBin.pdf")
ggsave(plot = MiscarrWeekBin, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MiscarrWeekBin.png")


###RPL
MiscarrWeekBinning_RPL <- MiscarrWeekBinning_SNPHapla[which(MiscarrWeekBinning_SNPHapla$Group2 == "RPL"),]

MiscarrWeekBinning_RPL$Week_bin <- factor(MiscarrWeekBinning_RPL$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))

MiscarrWeekBin_RPL <- ggplot(data= subset(MiscarrWeekBinning_RPL, !is.na(Week_bin)), aes(fill =  Indication_SNPHapla, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
  ggtitle("RPL") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Gestational age (in weeks)")
MiscarrWeekBin_RPL

ggsave(plot = MiscarrWeekBin_RPL, width = 6, height = 5,
        "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_RPL.pdf")
ggsave(plot = MiscarrWeekBin_RPL, width = 6, height = 5,
        "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_RPL.png")


###SPL
MiscarrWeekBinning_SPL <- MiscarrWeekBinning_SNPHapla[which(MiscarrWeekBinning_SNPHapla$Group2 == "SPL"),]

MiscarrWeekBinning_SPL$Week_bin <- factor(MiscarrWeekBinning_SPL$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))

MiscarrWeekBin_SPL <- ggplot(data = subset(MiscarrWeekBinning_SPL, !is.na(Week_bin)), aes(fill =  Indication_SNPHapla, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
  ggtitle("SPL") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Gestational age (in weeks)") 
MiscarrWeekBin_SPL

ggsave(plot = MiscarrWeekBin_SPL, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_SPL.pdf")
ggsave(plot = MiscarrWeekBin_SPL, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_SPL.png")

