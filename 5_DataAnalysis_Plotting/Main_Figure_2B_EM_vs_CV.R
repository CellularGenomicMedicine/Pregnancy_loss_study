###########################################################################################################################
# Author: Rick Essers
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To produce barplots showcasing the segregational origin (meiotic, mitotic, or combination) and the 
#                 difference in mosaicism degrees (in percentage) between fetal tissue DNAs (extraembryonic mesoderm vs. 
#                 chorionic villi)

# input: PL_PerTissue_Source_data.csv, this file contains information on each individual tissue type in miscarried 
#        POCs included in the pregnancy loss study. 

# output: Main figure 2 B; 2 barplots indicating the segregational origin of aberrations in miscarried POCs and the subset 
#         (mitotic cases + meiotic with >10% difference in mosaicism between the two tissue types) of cases utilized for 
#         the mosaicism comparison between EM and CV. A combination of half violin and scatterplot to indicate the difference 
#         in mosaicism degree of aberration in the EM anc CV tissues as well as a barplot with the number of genomic 
#         aberrations and their mosaicism degrees. 

###########################################################################################################################



###Total Barchart
rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

PerTissue <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerTissue_Source_data.csv")

PerTissue$Seg_Origin <- factor(PerTissue$Seg_Origin, levels = c(NA, "Mitotic & Meiotic ", "Mitotic", "Meiotic"))
Order <- c(NA, "Mitotic & Meiotic ", "Mitotic", "Meiotic")

PerTissue$DummyVar <- c(rep("x", times = 74))

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
       "/Projects/PregnancyLoss/Figures/AllAbb_SegOrigin.png")
ggsave(plot = AllAbb_SegOrigin, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/AllAbb_SegOrigin.pdf")


###BioRel Barchart 
rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

PerTissue <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerTissue_Source_data.csv")

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
       "/Projects/PregnancyLoss/Figures/BioRel_SegOrigin2.png")
ggsave(plot = BioRel_SegOrigin, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/BioRel_SegOrigin2.pdf")


##Coord_flip
  rm(list=ls(all=T))
  
  library(ggplot2)
  library(ggsignif)
  library(ggdist)
  
  source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")
  source("/Projects/PregnancyLoss/Figures/functions/Flat_Violin.R")
  
  
  PerTissue <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerTissue_Source_data.csv")
  
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
         "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM3.pdf")
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM3.png")


###Mosaicism % plot CV vs EM
###Alter input file so value 0% is 0.5, only for plotting purposes not statistics! 

  rm(list=ls(all=T))
  
  library(ggplot2)
  library(ggsignif)
  library(ggdist)
  
  source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")
  source("/Projects/PregnancyLoss/Figures/functions/Flat_Violin.R")
  
  
  PerTissue <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerTissue_Source_data.csv")
  
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
         "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM2.pdf")
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM2.png")
  ggsave(plot = MosaicismPercPlot, width = 6, height = 5,
         "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM2.svg")
  


### Distribution mosaicism percentage (10%) between EM and CV tissues, All cases including meiotic origiin. 
rm(list=ls(all=T))

library(ggplot2)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

PerTissue <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerTissue_Source_data.csv")

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
       "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM_comp.pdf")
ggsave(plot = MosDiff_CV_v_EM_comp, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM_comp.png")



### Distribution mosaicism percentage (10%) between EM and CV tissues, ONLY biologically relevant (e.g. mitotic)
rm(list=ls(all=T))

library(ggplot2)

source("/Projects/PregnancyLoss/Figures/functions/Indication_colors.R")

PerTissue <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerTissue_Source_data.csv")

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
       "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM_comp_BioRel2.pdf")
ggsave(plot = MosDiff_CV_v_EM_comp_BioRel, width = 6, height = 5,
       "/Projects/PregnancyLoss/Figures/MosDiff_CV_v_EM_comp_BioRel2.png")

