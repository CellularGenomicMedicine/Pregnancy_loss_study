rm(list=ls(all=T))

library(ggplot2)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")

PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_SegOr_fill.pdf")
ggsave(plot = NumPL_Seg_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_SegOr_fill.png")
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_SegOr_Bin_fill.pdf")
ggsave(plot = NumPL_Seg_bin_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_SegOr_Bin_fill.png")
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_ParOr_fill.pdf")
ggsave(plot = NumPL_PatMat_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_ParOr_fill.png")
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
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_ParOr_Bin_fill.pdf")
ggsave(plot = NumPL_PatMat_bin_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/NumPL_ParOr_Bin_fill.png")
}
