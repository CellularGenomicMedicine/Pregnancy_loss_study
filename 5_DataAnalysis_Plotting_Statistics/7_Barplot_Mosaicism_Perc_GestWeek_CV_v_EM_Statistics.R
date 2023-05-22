rm(list=ls(all=T))

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerTissue_Source_data.csv")

PerTissue_BioRel <- subset(PerTissue, Biologically_Relevant_to_Include_MosDiff == "y")


###Statistics mosaicism percentage and CV vs EM 
PerTissue_CV <- subset(PerTissue_BioRel, Group2 == "CV")
PerTissue_EM <- subset(PerTissue_BioRel, Group2 == "EM")

summary(PerTissue_EM$Mosaicism_Perc_Av, na.rm = T)
sd(PerTissue_EM$Mosaicism_Perc_Av, na.rm = T)

summary(PerTissue_CV$Mosaicism_Perc_Av, na.rm = T)
sd(PerTissue_CV$Mosaicism_Perc_Av, na.rm = T)

wilcox.test(PerTissue_CV$Mosaicism_Perc_Av, PerTissue_EM$Mosaicism_Perc_Av, paired = T)


###Total Barchart
rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

PerTissue$Seg_Origin <- factor(PerTissue$Seg_Origin, levels = c(NA, "Mitotic & Meiotic ", "Mitotic", "Meiotic"))
Order <- c(NA, "Mitotic & Meiotic ", "Mitotic", "Meiotic")

PerTissue$DummyVar <- c(rep("x", times = 72))

AllAbb_SegOrigin <- ggplot(data=PerTissue, aes(x=DummyVar,fill = factor(Seg_Origin, levels = Order, exclude = NULL))) + 
  geom_bar(stat="count", position = 'stack', color =  "black",width=0.3) +
  scale_x_discrete(drop=T) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank(),
        legend.position = "none") + 
  geom_text(aes(label = ..count..), 
            position = position_stack(vjust = 0.5), stat = "count", size = 3) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Mitotic" = color_Mitotic, "Meiotic" = color_Meiotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) +
  ylab("Genomic aberrations (#)") +
  xlab("") +
  coord_flip()
AllAbb_SegOrigin

ggsave(plot = AllAbb_SegOrigin, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/AllAbb_SegOrigin.png")
ggsave(plot = AllAbb_SegOrigin, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/AllAbb_SegOrigin.pdf")


###BioRel Barchart 
rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

PerTissue_BioRel <- subset(PerTissue, Biologically_Relevant_to_Include_MosDiff == "y")

PerTissue_BioRel$Seg_Origin <- factor(PerTissue_BioRel$Seg_Origin , levels = c("Mitotic", "Meiotic"))

Ordering3 <- c("Mitotic","Meiotic")



BioRel_SegOrigin <- ggplot(data=PerTissue_BioRel, aes(x=Biologically_Relevant_to_Include_MosDiff,fill = factor(Seg_Origin, levels = Ordering3))) + 
  geom_bar(stat="count", position = 'stack', color =  "black",width=0.3) +
  scale_x_discrete(drop=T) + 
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank(),
        legend.position = "none") + 
  geom_text(aes(label = ..count..), 
            position = position_stack(vjust = 0.5), stat = "count", size = 3) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("Mitotic" = color_Mitotic, "Meiotic" = color_Meiotic)) +
  ylab("Genomic aberrations (#)") +
  xlab("") +
  coord_flip()
BioRel_SegOrigin

ggsave(plot = BioRel_SegOrigin, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/BioRel_SegOrigin2.png")
ggsave(plot = BioRel_SegOrigin, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/BioRel_SegOrigin2.pdf")


###Only mitotic excluding the 2 meiotic cases. --> result still significant but less so. 
rm(list=ls(all=T))

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

PerTissue_BioRel <- subset(PerTissue, Biologically_Relevant_to_Include_MosDiff == "y")
PerTissue_BioRel <- subset(PerTissue_BioRel, Seg_Origin == "Mitotic")

PerTissue_CV <- subset(PerTissue_BioRel, Group2 == "CV")
PerTissue_EM <- subset(PerTissue_BioRel, Group2 == "EM")

summary(PerTissue_EM$Mosaicism_Perc_Av, na.rm = T)
sd(PerTissue_EM$Mosaicism_Perc_Av, na.rm = T)

summary(PerTissue_CV$Mosaicism_Perc_Av, na.rm = T)
sd(PerTissue_CV$Mosaicism_Perc_Av, na.rm = T)

wilcox.test(PerTissue_CV$Mosaicism_Perc_Av, PerTissue_EM$Mosaicism_Perc_Av, paired = T)



##Coord_flip
{
  rm(list=ls(all=T))
  
  library(ggplot2)
  library(ggsignif)
  library(ggdist)
  
  source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")
  source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Flat_Violin.R")
  
  
  PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")
  
  PerTissue_BioRel <- subset(PerTissue, Biologically_Relevant_to_Include_MosDiff == "y")
  
  PerTissue_BioRel$Mosaicism_Perc_Av[PerTissue_BioRel$Mosaicism_Perc_Av == 0] <- 0.5 
  
  PerTissue_BioRel$Group2  <- factor(PerTissue_BioRel$Group2, levels = c("EM","CV"))
  
  Order <- c("EM", "CV")
  
  MosaicismPercPlot <- 
    ggplot(data=subset(PerTissue_BioRel, !is.na(Mosaicism_Perc_Av)), aes(x=factor(Group2, levels = Order), y=Mosaicism_Perc_Av, fill = Group2)) + 
    scale_fill_manual(values = c("CV" = color_CV, "EM"  = color_EM)) +
    #option1
    geom_flat_violin(width = -0.9, alpha = 0.5) +
    stat_dots(side = "left", binwidth = 1, justification = 1.1, shape = 16, alpha = 1, aes( fill = "black")) + #, aes(color = "black", fill = "black") +
    #geom_dotplot(stackdir = "up", binaxis = "y", binwidth = 1, alpha = 0.3, position = position_jitterdodge(1.2), stackratio = 1.5) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) +
    #option2
    #geom_violin(width = 0.5, position = position_dodge(width = 0.7), alpha = 0.5) +
    #geom_jitter(shape = 16, size = 0.4, alpha = 0.3, position = position_jitterdodge(0.2)) + 
    #geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) +  
    #option3
    #geom_boxplot(width = 0.5, position = position_dodge(width = 0.7), outlier.shape = NA) +  
    #geom_jitter(size = 0.5, alpha = 0.3, position = position_jitterdodge(0.2)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), legend.title = element_blank(), 
          axis.title.y.left = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
          axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13), 
          legend.position = "none") +
    geom_signif(y_position = 99, xmin = 1, xmax = 2, annotation = "", 
                tip_length = 0.01, textsize = 4) +
    #ggtitle("Tissue comparison, abormal cases \n (n = 32)") +
    ylab("Mosaicism (%)") +
    coord_flip()
  #geom_signif(y_position = 99, xmin = 0.82, xmax = 1.18, annotation = "NS", tip_length = 0.01, textsize = 4)
  #to add figure legend: remove legend.position = none, add: legend.key = element_rect(fill = "transparent"), legend.key.size = unit(0.8,'cm'),
  MosaicismPercPlot
  
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM3.pdf")
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM3.png")
}

###Mosaicism % plot CV vs EM
###Alter input file so value 0% is 0.5, only for plotting purposes not statistics! 
{
  rm(list=ls(all=T))
  
  library(ggplot2)
  library(ggsignif)
  library(ggdist)
  
  source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")
  source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Flat_Violin.R")
  
  
  PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")
  
  PerTissue_BioRel <- subset(PerTissue, Biologically_Relevant_to_Include_MosDiff == "y")
  
  PerTissue_BioRel$Mosaicism_Perc_Av[PerTissue_BioRel$Mosaicism_Perc_Av == 0] <- 0.5 
  
  PerTissue_BioRel$Group2  <- factor(PerTissue_BioRel$Group2, levels = c("CV","EM"))
  
  MosaicismPercPlot <- 
    ggplot(data=subset(PerTissue_BioRel, !is.na(Mosaicism_Perc_Av)), aes(x=Group2, y=Mosaicism_Perc_Av, fill = Group2)) + 
    scale_fill_manual(values = c("CV" = color_CV, "EM"  = color_EM)) +
    #option1
    geom_flat_violin(width = 0.9, position = position_dodge(width = 0.7), alpha = 0.5) +
    stat_dots(side = "right", binwidth = 1, justification = -0.1, shape = 16, alpha = 1, aes( fill = "black")) + #, aes(color = "black", fill = "black") +
    #geom_dotplot(stackdir = "up", binaxis = "y", binwidth = 1, alpha = 0.3, position = position_jitterdodge(1.2), stackratio = 1.5) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) +
    #option2
    #geom_violin(width = 0.5, position = position_dodge(width = 0.7), alpha = 0.5) +
    #geom_jitter(shape = 16, size = 0.4, alpha = 0.3, position = position_jitterdodge(0.2)) + 
    #geom_boxplot(width = 0.1, position = position_dodge(width = 0.7), outlier.shape = NA) +  
    #option3
    #geom_boxplot(width = 0.5, position = position_dodge(width = 0.7), outlier.shape = NA) +  
    #geom_jitter(size = 0.5, alpha = 0.3, position = position_jitterdodge(0.2)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
    theme(panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), legend.title = element_blank(), 
          axis.title.x = element_blank(), strip.background = element_blank(),
          strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
          axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13), 
          legend.position = "none") +
    geom_signif(y_position = 99, xmin = 1, xmax = 2, annotation = "", 
                tip_length = 0.01, textsize = 4) +
    #ggtitle("Tissue comparison, abormal cases \n (n = 32)") +
    ylab("Mosaicism (%)") 
  #geom_signif(y_position = 99, xmin = 0.82, xmax = 1.18, annotation = "NS", tip_length = 0.01, textsize = 4)
  #to add figure legend: remove legend.position = none, add: legend.key = element_rect(fill = "transparent"), legend.key.size = unit(0.8,'cm'),
  MosaicismPercPlot
  
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM2.pdf")
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM2.png")
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM2.svg")
  
}

###Mosdif vs tissue type plot 10% bins
{
  rm(list=ls(all=T))
  
  library(ggplot2)
  
  PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")
  
  All_MossDiff_Tissue_plot <- ggplot(data= subset(PerTissue, !is.na(Mosaicism_bin_range_10)), aes(fill =  Mosaicism_bin_range_10, x = Group2)) +
    geom_bar(position = 'fill', stat = 'count', width = 0.5) +
    scale_fill_brewer(palette = "Blues") +
    scale_y_continuous(expand = c(0,0)) +
    theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
    scale_x_discrete(labels=c("CV","EM")) +
    geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
    ylab("Proportion") +
    xlab("Tissue type") 
  All_MossDiff_Tissue_plot
  
  ggsave(plot = All_MossDiff_Tissue_plot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_Tissue_10.pdf")
  ggsave(plot = All_MossDiff_Tissue_plot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_Tissue_10.png")
}

###Mosdif vs tissue type plot 20% bins
{
  rm(list=ls(all=T))
  
  library(ggplot2)
  
  PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")
  
  All_MossDiff_Tissue_plot <- ggplot(data= subset(PerTissue, !is.na(Mosaicism_bin_range_20)), aes(fill =  Mosaicism_bin_range_20, x = Group2)) +
    geom_bar(position = 'fill', stat = 'count', width = 0.5) +
    scale_fill_brewer(palette = "Blues") +
    scale_y_continuous(expand = c(0,0)) +
    theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
    scale_x_discrete(labels=c("CV","EM")) +
    geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
    ylab("Proportion") +
    xlab("Tissue type")
  All_MossDiff_Tissue_plot
  
  ggsave(plot = All_MossDiff_Tissue_plot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_Tissue_20.pdf")
  ggsave(plot = All_MossDiff_Tissue_plot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_Tissue_20.png")
}

###MosDiff vs gestational age plot
{
  rm(list=ls(all=T))
  
  library(ggplot2)
  
  PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")
  
  PerTissue$Week_bin <- factor(PerTissue$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))
  
  All_MossDiffplot <- ggplot(data= subset(PerTissue, !is.na(Week_bin) & !is.na(Mosaicism_bin_range_10)), aes(fill =  Mosaicism_bin_range_10, x = Week_bin)) +
    geom_bar(position = 'fill', stat = 'count', width = 0.5) +
    scale_fill_brewer(palette = "Blues") +
    scale_y_continuous(expand = c(0,0)) +
    theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
    scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
    geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
    ylab("Proportion") +
    xlab("Gestational age (in weeks)")
  All_MossDiffplot
  
  ggsave(plot = All_MossDiffplot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_GestWeek.pdf")
  ggsave(plot = All_MossDiffplot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_GestWeek.png")
}



### Distribution mosaicism percentage (10%) between EM and CV tissues, All cases including meiotic origiin. 
rm(list=ls(all=T))

library(ggplot2)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

MosDiff_CV_v_EM_comp <- ggplot(data= subset(PerTissue, !is.na(Mosaicism_bin_range_10)), aes(x = Mosaicism_bin_range_10, fill = Group2)) +
  geom_bar(stat = 'count', position = position_dodge(width = 0.8, preserve = 'single'), color =  "black",width=0.7) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title = element_blank(), axis.text.x = element_text(angle = 330, vjust = 0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,12.5)) +
  scale_fill_manual(values = c("CV" = color_CV, "EM" = color_EM)) +
  ylab("Genomic aberrations (#)") + 
  xlab("Mosaicism")
MosDiff_CV_v_EM_comp

#geom_text(aes(label= ..count..), stat = 'count', position = position_dodge(width = 0.8)) +

ggsave(plot = MosDiff_CV_v_EM_comp, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM_comp.pdf")
ggsave(plot = MosDiff_CV_v_EM_comp, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM_comp.png")



### Distribution mosaicism percentage (10%) between EM and CV tissues, ONLY biologically relevant (e.g. mitotic)
rm(list=ls(all=T))

library(ggplot2)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

PerTissue_BioRel <- subset(PerTissue, Biologically_Relevant_to_Include_MosDiff == "y")

PerTissue_BioRel$Mosaicism_bin_range_10 <- factor(PerTissue_BioRel$Mosaicism_bin_range_10, levels = c("0-10%", "11-20%", "21-30%", "31-40%", "41-50%", "51-60%", "61-70%", "71-80%", "81-90%", "91-100%"))

MosDiff_CV_v_EM_comp_BioRel <- ggplot(data= subset(PerTissue_BioRel, !is.na(Mosaicism_bin_range_10)), aes(x = Mosaicism_bin_range_10, fill = Group2)) +
  geom_bar(stat = 'count', position = position_dodge(width = 0.8, preserve = 'single'), color =  "black",width=0.7) +
  theme(legend.key.size = unit(0.5,'cm'), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.title = element_blank(), axis.text.x = element_text(angle = 330, vjust = 0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5)) +
  scale_x_discrete(drop = F) +
  scale_fill_manual(values = c("CV" = color_CV, "EM" = color_EM)) +
  ylab("Genomic aberrations (#)") + 
  xlab("Mosaicism")
MosDiff_CV_v_EM_comp_BioRel

ggsave(plot = MosDiff_CV_v_EM_comp_BioRel, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM_comp_BioRel2.pdf")
ggsave(plot = MosDiff_CV_v_EM_comp_BioRel, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_CV_v_EM_comp_BioRel2.png")



### Pie chart MI, MII, Mitotic 
rm(list=ls(all=T))

library(ggplot2)
library(scales)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/Indication_colors.R")

#Load function position_stack_and_nudge
sapply("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/functions/PositionStackandNudge.R", FUN=source)

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Source_data/PL_PerTissue_Source_data.csv")

PerTissue$Seg_Origin <- factor(PerTissue$Seg_Origin, levels = c(NA, "Mitotic & Meiotic ", "Meiotic", "Mitotic"))
Order <- c(NA, "Mitotic & Meiotic ", "Meiotic", "Mitotic")

#Plotting function
PieChart_SegOr_Hapla <- ggplot(PerTissue, aes(x="", fill = factor(Seg_Origin, levels = Order, exclude = NULL))) +
  geom_bar(stat = "count") +
  coord_polar("y", start = 0) +
  #ggtitle("Genome haplarithmisis \n (n = 72 tissues)") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 16)) +
  geom_text(aes(label = paste0("\n", percent(..count../74)), x = 1.7), stat = "count", 
            position = position_stack(vjust = 0.5), size = 6) +
  scale_fill_manual(values = c("Mitotic" = color_Mitotic, "Meiotic" = color_Meiotic, "Mitotic & Meiotic " = color_Mitotic_and_meiotic, "NA" = color_NA)) 
PieChart_SegOr_Hapla

ggsave(plot = PieChart_SegOr_Hapla, width = 12.1, height = 8.18, 
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_SegOr_PieChartBlue.pdf")
ggsave(plot = PieChart_SegOr_Hapla, width = 12.1, height = 8.18, 
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_SegOr_PieChartBlue.png") 


