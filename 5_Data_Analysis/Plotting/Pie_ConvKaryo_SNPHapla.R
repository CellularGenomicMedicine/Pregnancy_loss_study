###Circle plots
rm(list=ls(all=T))

library(ggplot2)
library(ggrepel)
library(scales)
library(ggpp)


###ConvKaryo plot
{
#Plot colors for indications
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")

#Load and alter input data
Demographic_data_TOTAL              <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv", sep=";")
Demographic_data_ConvKaryo_Abnormal     <- subset(Demographic_data_TOTAL, Indication_ConvKaryo == "Abnormal")

Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly <- gsub("[?]", "",Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly)
Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly <- gsub("Structural rearrangements", "Structural rearrangement",Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly)

Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly <- factor(Demographic_data_ConvKaryo_Abnormal$Classification_ConvKaryo_Poly, levels = c("Autosomal monosomy","Structural rearrangement","Others","Combined abnormalities","Gonosomal aneuploidy","Polyploidy","Autosomal trisomy"))
Order <- c("Autosomal \n monosomy","Structural \n rearrangement","Others","Combined \n abnormalities","Gonosomal \n aneuploidy","Polyploidy","Autosomal \n trisomy")

#Load function position_stack_and_nudge
sapply("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/PositionStackandNudge.R", FUN=source)

#Plotting function
Demographic_data_ConvKaryo_PieChart <- ggplot(Demographic_data_ConvKaryo_Abnormal, aes(x="", fill = Classification_ConvKaryo_Poly)) +
  geom_bar(stat = "count") +
  coord_polar("y", start = 0) +
  ggtitle("Coventional karyotyping \n (n = 879)") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16)) +
  geom_text(aes(label = paste0("\n", percent(..count../879)), x = 1.7), position = position_stack(vjust = 0.5), stat = "count", size = 6) +
  scale_fill_manual(values = c(color_Aut_mon,color_Str,color_Oth,color_Com_ab,color_Gon,color_Tetra,color_Aut_tri)) 
Demographic_data_ConvKaryo_PieChart

ggsave(plot = Demographic_data_ConvKaryo_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/ConvKaryo_Pie_Abnormal.pdf")
ggsave(plot = Demographic_data_ConvKaryo_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/ConvKaryo_Pie_Abnormal.png")
ggsave(plot = Demographic_data_ConvKaryo_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/ConvKaryo_Pie_Abnormal.svg")
}


###SNPHapla plot 
{
rm(list=ls(all=T))

#Plot colors for indications
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")
  
#Load and alter input data
Demographic_data_TOTAL              <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv", sep=";")
Demographic_data_SNPHapla_Abnormal     <- subset(Demographic_data_TOTAL, Indication_SNPHapla == "Abnormal")

Demographic_data_SNPHapla_Abnormal$Classification_SNPHapla <- factor(Demographic_data_SNPHapla_Abnormal$Classification_SNPHapla, levels = c("Autosomal monosomy","Polyploidy","Structural rearrangement","Gonosomal aneuploidy","Combined abnormalities","Autosomal trisomy"))
Order <- c("Autosomal \n monosomy","Polyploidy","Structural \n rearrangement","Gonosomal \n aneuploidy","Combined \n abnormalities","Autosomal \n trisomy")

#Load function position_stack_and_nudge
sapply("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/PositionStackandNudge.R", FUN=source)

#Plotting function
Demographic_data_SNPHapla_PieChart <- ggplot(Demographic_data_SNPHapla_Abnormal, aes(x="", fill = Classification_SNPHapla)) +
  geom_bar(stat = "count") +
  coord_polar("y", start = 0) +
  ggtitle("Genome haplarithmisis \n (n = 33)") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16)) +
  geom_text(aes(label = paste0("\n", percent(..count../33)), x = 1.7), stat = "count", 
            position = position_stack(vjust = 0.5), size = 6) +
  scale_fill_manual(values = c(color_Aut_mon,color_Tetra,color_Str,color_Gon,color_Com_ab,color_Aut_tri)) 
Demographic_data_SNPHapla_PieChart

ggsave(plot = Demographic_data_SNPHapla_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_Pie_Abnormal.pdf")
ggsave(plot = Demographic_data_SNPHapla_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_Pie_Abnormal.png")
ggsave(plot = Demographic_data_SNPHapla_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_Pie_Abnormal.svg")

}

###SNPHapla SPL plot
{
rm(list=ls(all=T))

#Plot colors for indications
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")
  
#Load and alter input data
Demographic_data_TOTAL              <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv", sep=";")
Demographic_data_SNPHapla_Abnormal      <- subset(Demographic_data_TOTAL, Indication_SNPHapla == "Abnormal")
Demographic_data_SNPHapla_Abnormal_SPL  <- subset(Demographic_data_SNPHapla_Abnormal, Group2 == "SPL")

Demographic_data_SNPHapla_Abnormal_SPL$Classification_SNPHapla <- factor(Demographic_data_SNPHapla_Abnormal_SPL$Classification_SNPHapla, levels = c("Structural rearrangement","Combined abnormalities","Autosomal trisomy"))
Order <- c("Structural \n rearrangement","Combined \n abnormalities","Autosomal \n trisomy")

#Load function position_stack_and_nudge
sapply("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/PositionStackandNudge.R", FUN=source)

#Plotting function
Demographic_data_SNPHapla_SPL_PieChart <- ggplot(Demographic_data_SNPHapla_Abnormal_SPL, aes(x="", fill = Classification_SNPHapla)) +
  geom_bar(stat = "count") +
  coord_polar("y", start = 0) +
  ggtitle("SNP Haplotyping SPL \n (n = 14)") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16)) +
  geom_text(aes(label = paste0("\n", percent(..count../14)), x = 1.7), stat = "count", 
            position = position_stack(vjust = 0.5), size = 6) +
  scale_fill_manual(values = c(color_Str,color_Com_ab,color_Aut_tri)) 
Demographic_data_SNPHapla_SPL_PieChart

ggsave(plot = Demographic_data_SNPHapla_SPL_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_SPL_Pie_Abnormal.pdf")
ggsave(plot = Demographic_data_SNPHapla_SPL_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_SPL_Pie_Abnormal.png")
ggsave(plot = Demographic_data_SNPHapla_SPL_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_SPL_Pie_Abnormal.svg")

}

###SNPHapla RPL plot
{
rm(list=ls(all=T))

#Plot colors for indications
source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")
  

#Load and alter input data
Demographic_data_TOTAL              <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv", sep=";")
Demographic_data_SNPHapla_Abnormal      <- subset(Demographic_data_TOTAL, Indication_SNPHapla == "Abnormal")
Demographic_data_SNPHapla_Abnormal_RPL  <- subset(Demographic_data_SNPHapla_Abnormal, Group2 == "RPL")

Demographic_data_SNPHapla_Abnormal_RPL$Classification_SNPHapla <- factor(Demographic_data_SNPHapla_Abnormal_RPL$Classification_SNPHapla, levels = c("Autosomal monosomy","Combined abnormalities","Polyploidy","Gonosomal aneuploidy","Autosomal trisomy"))
Order <- c("Autosomal \n monosomy","Combined \n abnormalities","Polyploidy","Gonosomal \n aneuploidy","Autosomal \n trisomy")

#Load function position_stack_and_nudge
sapply("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/PositionStackandNudge.R", FUN=source)

#Plotting function
Demographic_data_SNPHapla_RPL_PieChart <- ggplot(Demographic_data_SNPHapla_Abnormal_RPL, aes(x="", fill = Classification_SNPHapla)) +
  geom_bar(stat = "count") +
  coord_polar("y", start = 0) +
  ggtitle("SNP Haplotyping RPL \n (n = 19)") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16)) +
  geom_text(aes(label = paste0("\n", percent(..count../19)), x = 1.7), stat = "count", 
            position = position_stack(vjust = 0.5), size = 6) +
  scale_fill_manual(values = c(color_Aut_mon,color_Com_ab,color_Tetra,color_Gon,color_Aut_tri)) 
Demographic_data_SNPHapla_RPL_PieChart


ggsave(plot = Demographic_data_SNPHapla_RPL_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_RPL_Pie_Abnormal.pdf")
ggsave(plot = Demographic_data_SNPHapla_RPL_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_RPL_Pie_Abnormal.png")
ggsave(plot = Demographic_data_SNPHapla_RPL_PieChart, width = 12.1, height = 8.18,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/2023/SNPHapla_RPL_Pie_Abnormal.svg")

}
