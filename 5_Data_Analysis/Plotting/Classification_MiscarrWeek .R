library(ggplot2)

MiscarrWeekBinning <- read.csv("~/Surfdrive/ClinicalGenetics/Miscarriage/Paper/MiscarrWeek_AllBin.csv", sep=";")

###RPL and SPL 
MiscarrWeekBinning$Week_bin <- factor(MiscarrWeekBinning$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))
MiscarrWeekBinning$Classification <- factor(MiscarrWeekBinning$Classification, levels = c("Normal","Abnormal"))

MiscarrWeekBin <- ggplot(data= subset(MiscarrWeekBinning, !is.na(Week_bin)), aes(fill =  Classification, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
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
MiscarrWeekBinning_RPL <- MiscarrWeekBinning[which(MiscarrWeekBinning$Group == "RPL"),]

MiscarrWeekBinning_RPL$Week_bin <- factor(MiscarrWeekBinning_RPL$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))

MiscarrWeekBin_RPL <- ggplot(data=MiscarrWeekBinning_RPL, aes(fill =  Classification, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  ggtitle("RPL") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Gestational age (in weeks)")
MiscarrWeekBin_RPL

ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_RPL.pdf")
ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_RPL.png")

###SPL
MiscarrWeekBinning_SPL <- MiscarrWeekBinning[which(MiscarrWeekBinning$Group == "SPL"),]

MiscarrWeekBinning_SPL$Week_bin <- factor(MiscarrWeekBinning_SPL$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))

MiscarrWeekBin_SPL <- ggplot(data=MiscarrWeekBinning_SPL, aes(fill =  Classification, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  ggtitle("SPL") +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Gestational age (in weeks)") 
MiscarrWeekBin_SPL

ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_SPL.pdf")
ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/WeekOfMiscarrPlot/MiscarrWeekBin_SPL.png")



##3MosDiff 
All_MossDef <- read.csv("~/Surfdrive/ClinicalGenetics/Miscarriage/Paper/All_MossDef.csv", sep=";")

All_MossDef$Week_bin <- factor(All_MossDef$Week_bin, levels = c("4-5","6-7","8-9","10-13","NA"))

All_MossDefplot <- ggplot(data= subset(All_MossDef, !is.na(Week_bin) & !is.na(Mosaicism_bin)), aes(fill =  Mosaicism_bin, x = Week_bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c("4-5","6-7","8-9","10-13", "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Gestational age (in weeks) \n ")
All_MossDefplot

ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/Mosaicism_differencePlot/MiscarrWeek_Mosaicism.pdf")
ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/Mosaicism_differencePlot/MiscarrWeek_Mosaicism.png")

