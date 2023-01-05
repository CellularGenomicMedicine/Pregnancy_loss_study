rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)
library(ggsignif)

source("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Scripts/Github/5_DataAnalysis_Plotting_Statistics/Plotting/functions/Indication_colors.R")

###Preprocessing of data files
{
PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")

ConvKaryoConversion <- subset(PerFamily, Group == "SNPHapla")

ConvKaryoConversion$Group <- c(rep("ConvKaryo", times = 111))

All_Ages <- data.frame(All_Ages = c(ConvKaryoConversion$Maternal.age, ConvKaryoConversion$Paternal.age, ConvKaryoConversion$Gestational.age))
ConvKaryoConversion <- rbind(ConvKaryoConversion, ConvKaryoConversion, ConvKaryoConversion)
ConvKaryoConversion$Age_Type <- c(rep("Maternal age", times = 111), rep("Paternal age", times = 111), rep("Gestational age", times = 111))

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
}

###ConvKaryo individual plot 
{
PerFamilyPatMatConvKaryo <- subset(PerFamilyPatMat, Group == "Conventional karyotyping \n (n = 1745)")

TOTAL_AgesPlot_PatMat_ConvKaryo <- 
  ggplot(data=subset(PerFamilyPatMatConvKaryo, !is.na(All_Ages)), aes(x=Age_Type, y=All_Ages, fill=Indication, color = Indication)) + 
  geom_jitter(size = 0.5, alpha = 0.3, position = position_jitterdodge(0.2)) + 
  scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
  scale_color_manual(values = c("black","black")) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.7), outlier.shape = NA) +
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank(), 
        axis.title.x = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), 
        legend.key = element_rect(fill = "transparent"), axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13)) +
  ggtitle("Conventional karyotyping \n (n = 1745)") +
  ylim(10,63) +
  ylab("Age (years)") +
  geom_signif(y_position = c(62, 62), xmin = c(0.82, 1.82), xmax = c(1.18,2.18), annotation = c("p = 6.72e-5", "p = 4.57e-3"), 
              tip_length = 0.01, textsize = 4)
TOTAL_AgesPlot_PatMat_ConvKaryo
}

###SNPHapla individual plot 
{
PerFamilyPatMatSNPHapla <- subset(PerFamilyPatMat, Group == "SNP haplotyping \n (n = 91)")

TOTAL_AgesPlot_PatMat_SNPHapla <- 
  ggplot(data=subset(PerFamilyPatMatSNPHapla, !is.na(All_Ages)), aes(x=Age_Type, y=All_Ages, fill=Indication, color = Indication)) + 
  geom_jitter(size = 0.5, alpha = 0.3, position = position_jitterdodge(0.2)) + 
  scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
  scale_color_manual(values = c("black","black")) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.7), outlier.shape = NA) +
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
}

###Gestational age ConvKaryo
{
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
  geom_jitter(size = 0.5, alpha = 0.3, position = position_jitterdodge(0.2)) + 
  scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
  scale_color_manual(values = c("black","black")) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.7), outlier.shape = NA) +
  theme(legend.key.size = unit(0.8,'cm'), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.title = element_blank(), 
        axis.title.x = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 12), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = "transparent"), axis.title.y = element_text(size = 13), axis.text.x = element_text(size = 13)) +
  ggtitle("Conventional karyotyping \n (n = 1745)") +
  ylim(2,30) +
  ylab("Age (weeks)") +
  geom_signif(y_position = 30, xmin = 0.82, xmax = 1.18, annotation = "NS", tip_length = 0.01, textsize = 4) 
TOTAL_AgesPlot_PatMat_ConvKaryo_Gest
}

###Gestational age SNPHapla
{
PerFamilyPatMatSNPHaplaGest <- subset(PerFamilyGest, Group == "SNP Haplotyping \n (n = 91)")

TOTAL_AgesPlot_PatMat_SNPHapla_Gest <- 
  ggplot(data=subset(PerFamilyPatMatSNPHaplaGest, !is.na(All_Ages)), aes(x=Age_Type, y=All_Ages, fill=Indication, color = Indication)) + 
  geom_jitter(size = 0.5, alpha = 0.3, position = position_jitterdodge(0.2)) + 
  scale_fill_manual(values = c("Normal" = color_Normal,"Abnormal" = color_Abnormal)) +
  scale_color_manual(values = c("black","black")) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.7), outlier.shape = NA) +
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
}

###Save all plots same height and width
{
ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/PatMatAge_ConvKaryo.pdf")
ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/PatMatAge_ConvKaryo.png")

ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/PatMatAge_SNPHapla.pdf")
ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/PatMatAge_SNPHapla.png")

ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo_Gest, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GestAge_ConvKaryo.pdf")
ggsave(plot = TOTAL_AgesPlot_PatMat_ConvKaryo_Gest, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GestAge_ConvKaryo.png")

ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla_Gest, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GestAge_SNPHapla.pdf")
ggsave(plot = TOTAL_AgesPlot_PatMat_SNPHapla_Gest, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/GestAge_SNPHapla.png")
}

###Statistics
{
###Preparing input files for statistics, ConvKaryo vs SNPHapla, and Normal vs Abnormal within the groups. 
Demographic_data_TOTAL <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")
Demographic_data_SNPHapla           <- subset(Demographic_data_TOTAL, Group == "SNPHapla")

Demographic_data_TOTAL_Normal       <- subset(Demographic_data_TOTAL, Indication_ConvKaryo == "Normal")
Demographic_data_TOTAL_Abnormal     <- subset(Demographic_data_TOTAL, Indication_ConvKaryo == "Abnormal")

Demographic_data_SNPHapla_Normal    <- subset(Demographic_data_SNPHapla, Indication_SNPHapla == "Normal")
Demographic_data_SNPHapla_Abnormal  <- subset(Demographic_data_SNPHapla, Indication_SNPHapla == "Abnormal")

#Demographic_data_SNPHapla_SPL      <- subset(Demographic_data_SNPHapla, Group2 == "SPL")
#Demographic_data_SNPHapla_RPL      <- subset(Demographic_data_SNPHapla, Group2 == "RPL")

###Calculate mean + sd for gestational age in SNPHapla group (n= 91)
Demographic_data_SNPHapla_Exc       <- rbind(Demographic_data_SNPHapla_Normal, Demographic_data_SNPHapla_Abnormal)
summary(Demographic_data_SNPHapla_Exc$Gestational.age)
sd(Demographic_data_SNPHapla_Exc$Gestational.age, na.rm = T)

Demographic_data_SNPHapla_SPL      <- subset(Demographic_data_SNPHapla_Exc, Group2 == "SPL")
Demographic_data_SNPHapla_RPL      <- subset(Demographic_data_SNPHapla_Exc, Group2 == "RPL")

###Calculate mean + sd for all ages in SPL SNPHapla group (n= 42)
summary(Demographic_data_SNPHapla_SPL$Maternal.age)
sd(Demographic_data_SNPHapla_SPL$Maternal.age, na.rm = T)

summary(Demographic_data_SNPHapla_SPL$Paternal.age)
sd(Demographic_data_SNPHapla_SPL$Paternal.age, na.rm = T)

summary(Demographic_data_SNPHapla_SPL$Gestational.age)
sd(Demographic_data_SNPHapla_SPL$Gestational.age, na.rm = T)

###Calculate mean + sd for all ages in RPL SNPHapla group (n= 49)
summary(Demographic_data_SNPHapla_RPL$Maternal.age)
sd(Demographic_data_SNPHapla_RPL$Maternal.age, na.rm = T)

summary(Demographic_data_SNPHapla_RPL$Paternal.age)
sd(Demographic_data_SNPHapla_RPL$Paternal.age, na.rm = T)

summary(Demographic_data_SNPHapla_RPL$Gestational.age)
sd(Demographic_data_SNPHapla_RPL$Gestational.age, na.rm = T)

summary(Demographic_data_SNPHapla_RPL$Total_PLs)
sd(Demographic_data_SNPHapla_RPL$Total_PLs, na.rm = T)

###Maternal age between ConvKaryo vs SNPHapla
t.test(Demographic_data_TOTAL$Maternal.age, Demographic_data_SNPHapla$Maternal.age, paired = F)
#Welch Two Sample t-test

#data:  Demographic_data_TOTAL$Maternal.age and Demographic_data_SNPHapla$Maternal.age
#t = -2.0356, df = 122.81, p-value = 0.04394
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -2.41208431 -0.03374638
#sample estimates:
#  mean of x mean of y 
#28.40094  29.62385 

###Paternal age
t.test(Demographic_data_TOTAL$Paternal.age, Demographic_data_SNPHapla$Paternal.age)
#Welch Two Sample t-test

#data:  Demographic_data_TOTAL$Paternal.age and Demographic_data_SNPHapla$Paternal.age
#t = -0.36377, df = 118.97, p-value = 0.7167
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.4487222  0.9990432
#sample estimates:
#  mean of x mean of y 
#30.89397  31.11881 

###Gestational age
t.test(Demographic_data_TOTAL$Gestational.age, Demographic_data_SNPHapla$Gestational.age)
#Welch Two Sample t-test

#data:  Demographic_data_TOTAL$Gestational.age and Demographic_data_SNPHapla$Gestational.age
#t = -0.64862, df = 100.95, p-value = 0.5181
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.5376068  0.2726721
#sample estimates:
#  mean of x mean of y 
#7.478376  7.610843 

###ConvKaryo normal vs abnormal, maternal age
t.test(Demographic_data_TOTAL_Normal$Maternal.age, Demographic_data_TOTAL_Abnormal$Maternal.age)

#Welch Two Sample t-test

#data:  Demographic_data_TOTAL_Normal$Maternal.age and Demographic_data_TOTAL_Abnormal$Maternal.age
#t = -4.0083, df = 1702.1, p-value = 6.379e-05
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.7809359 -0.6106713
#sample estimates:
#  mean of x mean of y 
#27.79042  28.98622 

###ConvKaryo normal vs abnormal, paternal age
t.test(Demographic_data_TOTAL_Normal$Paternal.age, Demographic_data_TOTAL_Abnormal$Paternal.age)
#Welch Two Sample t-test

#data:  Demographic_data_TOTAL_Normal$Paternal.age and Demographic_data_TOTAL_Abnormal$Paternal.age
#t = -2.837, df = 1329.5, p-value = 0.004623
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.7082326 -0.3115724
#sample estimates:
#  mean of x mean of y 
#30.32558  31.33548 

###ConvKaryo normal vs abnormal, gestational age
t.test(as.numeric(Demographic_data_TOTAL_Normal$Gestational.age), as.numeric(Demographic_data_TOTAL_Abnormal$Gestational.age))
#Welch Two Sample t-test

#data:  as.numeric(Demographic_data_TOTAL_Normal$Gestational.age) and as.numeric(Demographic_data_TOTAL_Abnormal$Gestational.age)
#t = -0.37164, df = 953.53, p-value = 0.7102
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.3098213  0.2111601
#sample estimates:
#  mean of x mean of y 
#7.450206  7.499536 

###SNPHapla normal vs abnormal, maternal age
t.test(Demographic_data_SNPHapla_Normal$Maternal.age, Demographic_data_SNPHapla_Abnormal$Maternal.age)
#Welch Two Sample t-test

#data:  Demographic_data_SNPHapla_Normal$Maternal.age and Demographic_data_SNPHapla_Abnormal$Maternal.age
#t = -5.7912, df = 70.251, p-value = 1.806e-07
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -8.879161 -4.330269
#sample estimates:
#  mean of x mean of y 
#27.61404  34.21875 

###SNPHapla normal vs abnormal, paternal age
t.test(Demographic_data_SNPHapla_Normal$Paternal.age, Demographic_data_SNPHapla_Abnormal$Paternal.age)
#Welch Two Sample t-test

#data:  Demographic_data_SNPHapla_Normal$Paternal.age and Demographic_data_SNPHapla_Abnormal$Paternal.age
#t = -4.8423, df = 68.594, p-value = 7.638e-06
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -7.999698 -3.331106
#sample estimates:
#  mean of x mean of y 
#29.43137  35.09677 

###SNPHapla normal vs abnormal, Gestational age
t.test(as.numeric(Demographic_data_SNPHapla_Normal$Gestational.age), as.numeric(Demographic_data_SNPHapla_Abnormal$Gestational.age))
#Welch Two Sample t-test

#data:  as.numeric(Demographic_data_SNPHapla_Normal$Gestational.age) and as.numeric(Demographic_data_SNPHapla_Abnormal$Gestational.age)
#t = 0.29286, df = 54.657, p-value = 0.7707
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.7134178  0.9575737
#sample estimates:
#  mean of x mean of y 
#7.611364  7.489286 

###SNPHapla SPL vs RPL 
t.test(Demographic_data_SNPHapla_SPL$Maternal.age, Demographic_data_SNPHapla_RPL$Maternal.age)

t.test(Demographic_data_SNPHapla_SPL$Paternal.age, Demographic_data_SNPHapla_RPL$Paternal.age)

t.test(Demographic_data_SNPHapla_SPL$Gestational.age, Demographic_data_SNPHapla_RPL$Gestational.age)

}
