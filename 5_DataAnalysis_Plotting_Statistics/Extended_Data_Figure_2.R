rm(list=ls(all=T))

library(ggplot2)
library(scales)

###SPL
###Plot colors for indications
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerPOC_Source_data.csv")
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/Barplot_SNPHapla_Normal_Abnormal_SPL.png")
ggsave(plot = Barplot_SNPHapla_Normal_Abnormal_SPL, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/Barplot_SNPHapla_Normal_Abnormal_SPL.pdf")


###RPL
rm(list=ls(all=T))

library(ggplot2)
library(scales)

###Plot colors for indications
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

#Load and alter input data
Demographic_data_TOTAL <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerPOC_Source_data.csv")
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/Barplot_SNPHapla_Normal_Abnormal_RPL.png")
ggsave(plot = Barplot_SNPHapla_Normal_Abnormal_RPL, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/Barplot_SNPHapla_Normal_Abnormal_RPL.pdf")


rm(list=ls(all=T))

library(ggplot2)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerPOC_Source_data.csv")

#PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")
NumPLs_SNPHapla <- subset(PerFamily, Group == "SNPHapla" & Indication_SNPHapla == "Abnormal")

Seg_Or_Order <- c("Mitotic", "Mitotic & Meiotic ", "Meiotic")
Par_Or_Order <- c("Paternal","Maternal & Paternal ", "Maternal")

###Position = "fill"
{
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_SegOr_fill.pdf")
ggsave(plot = NumPL_Seg_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_SegOr_fill.png")
}

###Binning of number of pregnancy losses
###Input file conversion 

###Bin 1 and 2 into 1-2 
{
NumPLs_SNPHapla$Total_PLs_Bin <- with(NumPLs_SNPHapla, ifelse(Total_PLs == 1 | Total_PLs == 2, "1-2","3-5"))

NumPL_Seg_bin_plot <- ggplot(data = subset(NumPLs_SNPHapla, !is.na(Seg_origin2)), aes(x = Total_PLs_Bin, fill = factor(Seg_origin2, levels = Seg_Or_Order))) +
  geom_bar(stat="count", color =  "black", width=0.7, position = "fill") +
  scale_fill_manual(values = c("Meiotic" = color_Meiotic, "Mitotic" = color_Mitotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("Number of pregnancy losses")
NumPL_Seg_bin_plot

ggsave(plot = NumPL_Seg_bin_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_SegOr_Bin_fill.pdf")
ggsave(plot = NumPL_Seg_bin_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_SegOr_Bin_fill.png")
}

###Parental origin
###Position = "fill"
{
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_ParOr_fill.pdf")
ggsave(plot = NumPL_PatMat_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_ParOr_fill.png")
}

###Binning of number of pregnancy losses
###Input file conversion 

###Bin 1 and 2 into 1-2 
{
NumPLs_SNPHapla$Total_PLs_Bin <- with(NumPLs_SNPHapla, ifelse(Total_PLs == 1 | Total_PLs == 2, "1-2","3-5"))

NumPL_PatMat_bin_plot <- ggplot(data = subset(NumPLs_SNPHapla, !is.na(Parental_origin)), aes(x = Total_PLs_Bin, fill = factor(Parental_origin, levels = Par_Or_Order))) +
  geom_bar(stat="count", color =  "black", width=0.7, position = "fill") +
  scale_fill_manual(values = c("Maternal" = color_Maternal, "Paternal" = color_Paternal, "Maternal & Paternal " = color_Maternal_and_Paternal)) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("Number of pregnancy losses")
NumPL_PatMat_bin_plot

ggsave(plot = NumPL_PatMat_bin_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_ParOr_Bin_fill.pdf")
ggsave(plot = NumPL_PatMat_bin_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/NEW/Extended_Figures/NumPL_ParOr_Bin_fill.png")
}
