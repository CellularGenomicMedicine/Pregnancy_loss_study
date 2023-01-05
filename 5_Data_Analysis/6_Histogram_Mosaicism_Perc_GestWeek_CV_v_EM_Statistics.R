rm(list=ls(all=T))

library(ggplot2)

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

PerTissue_CV <- subset(PerTissue, Group2 == "CV")
PerTissue_EM <- subset(PerTissue, Group2 == "EM")

summary(PerTissue_EM$Mosaicism_Perc2, na.rm = T)
sd(PerTissue_EM$Mosaicism_Perc2, na.rm = T)

summary(PerTissue_CV$Mosaicism_Perc2, na.rm = T)
sd(PerTissue_CV$Mosaicism_Perc2, na.rm = T)

t.test(PerTissue_CV$Mosaicism_Perc2, PerTissue_EM$Mosaicism_Perc2, paired = T)

#Paired t-test

#data:  PerTissue_CV$Mosaicism_Perc2 and PerTissue_EM$Mosaicism_Perc2
#t = -2.2782, df = 33, p-value = 0.02932
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -16.1463453  -0.9124783
#sample estimates:
#  mean of the differences 
#-8.529412 

###Mosaicism % plot EM vs CV (Unfinished)
{
MosaicismPercPlot <- 
  ggplot(data=subset(PerTissue, !is.na(Mosaicism_Perc2)), aes(x=Group2, y=Mosaicism_Perc2, fill = Group2)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank(), 
        axis.title.x = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = "transparent"), 
        axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13)) +
  ggtitle("Conventional karyotyping \n (n = 1745)") +
  #ylim(0,100) +
  ylab("Mosaicism (%)") 
#+
 # geom_signif(y_position = 99, xmin = 0.82, xmax = 1.18, annotation = "NS", tip_length = 0.01, textsize = 4)
MosaicismPercPlot
}

###Mosdif vs tissue type plot
rm(list=ls(all=T))

library(ggplot2)

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

All_MossDiff_Tissue_plot <- ggplot(data= subset(PerTissue, !is.na(Mosaicism_bin)), aes(fill =  Mosaicism_bin, x = Group2)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("CV","EM")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Tissue type")
All_MossDiff_Tissue_plot

ggsave(plot = All_MossDiff_Tissue_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_Tissue.pdf")
ggsave(plot = All_MossDiff_Tissue_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_Tissue.png")

###MosDiff vs gestational age plot
rm(list=ls(all=T))

library(ggplot2)

PerTissue <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerTissue.csv")

PerTissue$Week_bin <- factor(PerTissue$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))

All_MossDiffplot <- ggplot(data= subset(PerTissue, !is.na(Week_bin) & !is.na(Mosaicism_bin)), aes(fill =  Mosaicism_bin, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Gestational age (in weeks)")
All_MossDiffplot

ggsave(plot = All_MossDiffplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_GestWeek.pdf")
ggsave(plot = All_MossDiffplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/MosDiff_GestWeek.png")

