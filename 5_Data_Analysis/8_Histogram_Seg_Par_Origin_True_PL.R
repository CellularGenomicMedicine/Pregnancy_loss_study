###True RPL (no previous live born) vs assumed RPL (previous live born, but still 2 consecutive PLs)

rm(list=ls(all=T))

library(ggplot2)
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")

PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")
RPL_types_SNPHapla <- subset(PerFamily, Group == "SNPHapla" & Group2 == "RPL" & Indication_SNPHapla == "Abnormal")

###Segregational origin
Seg_Or_Order <- c("Mitotic", "Meiotic")

###Position = "fill"
RPL_types_Seg_Or_plot <- ggplot(data = subset(RPL_types_SNPHapla, !is.na(Seg_origin2)), aes(x = RPL_status_based_on_Livebornchild, fill = factor(Seg_origin2, levels = Seg_Or_Order))) +
  geom_bar(stat="count", color =  "black", width=0.7, position = "fill") +
  scale_fill_manual(values = c("Meiotic" = color_Meiotic, "Mitotic" = color_Mitotic)) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("")
RPL_types_Seg_Or_plot

ggsave(plot = RPL_types_Seg_Or_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/RPL_types_SegOr_fill.pdf")
ggsave(plot = RPL_types_Seg_Or_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/RPL_types_SegOr_fill.png")


###Parental origin
Par_Or_Order <- c("Paternal", "Maternal & Paternal ", "Maternal")

RPL_types_Par_Or_plot <- ggplot(data = subset(RPL_types_SNPHapla, !is.na(Parental_origin)), aes(x = RPL_status_based_on_Livebornchild, fill = factor(Parental_origin, levels = Par_Or_Order))) +
  geom_bar(stat="count", color =  "black", width=0.7, position = "fill") +
  scale_fill_manual(values = c("Maternal" = color_Maternal, "Paternal" = color_Paternal, "Maternal & Paternal " = color_Maternal_and_Paternal)) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank()) + 
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  ylab("Proportion") +
  xlab("")
RPL_types_Par_Or_plot

ggsave(plot = RPL_types_Par_Or_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/RPL_types_ParOr_fill.pdf")
ggsave(plot = RPL_types_Par_Or_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/RPL_types_ParOr_fill.png")
