rm(list=ls(all=T))

library(ggplot2)
library(ggpattern)
library(ggsignif)

PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")

Ages_PerFamily <- subset(PerFamily, Group == "SNPHapla")

Ages_PerFamily$Group <- gsub('SNPHapla', 'SNP Haplotyping \n (n = 91)',Ages_PerFamily$Group)

Ages_PerFamily <- subset(Ages_PerFamily, Indication_SNPHapla != "Excluded")


###Make ages bin columns 
{
UniqueMatAge <- unique(Ages_PerFamily$Maternal.age)
UniqueMatAge[order(UniqueMatAge)]

Ages_PerFamily$Mat_Age_Bin <- ifelse(Ages_PerFamily$Maternal.age == 17, "LowerThan24", 
                              ifelse(Ages_PerFamily$Maternal.age == 18, "LowerThan24",       
                              ifelse(Ages_PerFamily$Maternal.age == 20, "LowerThan24",       
                              ifelse(Ages_PerFamily$Maternal.age == 21, "LowerThan24",       
                              ifelse(Ages_PerFamily$Maternal.age == 22, "LowerThan24",       
                              ifelse(Ages_PerFamily$Maternal.age == 23, "LowerThan24",       
                              ifelse(Ages_PerFamily$Maternal.age == 24, "LowerThan24",       
                              ifelse(Ages_PerFamily$Maternal.age == 25, "25-29",       
                              ifelse(Ages_PerFamily$Maternal.age == 26, "25-29",       
                              ifelse(Ages_PerFamily$Maternal.age == 27, "25-29",       
                              ifelse(Ages_PerFamily$Maternal.age == 28, "25-29",       
                              ifelse(Ages_PerFamily$Maternal.age == 29, "25-29",       
                              ifelse(Ages_PerFamily$Maternal.age == 30, "30-34",       
                              ifelse(Ages_PerFamily$Maternal.age == 31, "30-34",       
                              ifelse(Ages_PerFamily$Maternal.age == 32, "30-34",       
                              ifelse(Ages_PerFamily$Maternal.age == 33, "30-34",       
                              ifelse(Ages_PerFamily$Maternal.age == 34, "30-34",       
                              ifelse(Ages_PerFamily$Maternal.age == 35, "35-39",       
                              ifelse(Ages_PerFamily$Maternal.age == 36, "35-39",       
                              ifelse(Ages_PerFamily$Maternal.age == 37, "35-39",       
                              ifelse(Ages_PerFamily$Maternal.age == 38, "35-39",       
                              ifelse(Ages_PerFamily$Maternal.age == 39, "35-39",       
                              ifelse(Ages_PerFamily$Maternal.age == 40, "HigherThan40",       
                              ifelse(Ages_PerFamily$Maternal.age == 41, "HigherThan40",       
                              ifelse(Ages_PerFamily$Maternal.age == 42, "HigherThan40",       
                              ifelse(Ages_PerFamily$Maternal.age == 44, "HigherThan40",       
                                NA))))))))))))))))))))))))))
  
UniquePatAge <- unique(Ages_PerFamily$Paternal.age)
UniquePatAge[order(UniquePatAge)]

Ages_PerFamily$Pat_Age_Bin <- ifelse(Ages_PerFamily$Paternal.age == 18, "LowerThan24", 
                              ifelse(Ages_PerFamily$Paternal.age == 20, "LowerThan24",       
                              ifelse(Ages_PerFamily$Paternal.age == 22, "LowerThan24",       
                              ifelse(Ages_PerFamily$Paternal.age == 23, "LowerThan24",       
                              ifelse(Ages_PerFamily$Paternal.age == 24, "LowerThan24",       
                              ifelse(Ages_PerFamily$Paternal.age == 25, "25-29",       
                              ifelse(Ages_PerFamily$Paternal.age == 26, "25-29",       
                              ifelse(Ages_PerFamily$Paternal.age == 27, "25-29",       
                              ifelse(Ages_PerFamily$Paternal.age == 28, "25-29",       
                              ifelse(Ages_PerFamily$Paternal.age == 29, "25-29",       
                              ifelse(Ages_PerFamily$Paternal.age == 30, "30-34",       
                              ifelse(Ages_PerFamily$Paternal.age == 31, "30-34",       
                              ifelse(Ages_PerFamily$Paternal.age == 32, "30-34",       
                              ifelse(Ages_PerFamily$Paternal.age == 33, "30-34",       
                              ifelse(Ages_PerFamily$Paternal.age == 34, "30-34",       
                              ifelse(Ages_PerFamily$Paternal.age == 35, "35-39",       
                              ifelse(Ages_PerFamily$Paternal.age == 36, "35-39",       
                              ifelse(Ages_PerFamily$Paternal.age == 37, "35-39",       
                              ifelse(Ages_PerFamily$Paternal.age == 38, "35-39",       
                              ifelse(Ages_PerFamily$Paternal.age == 39, "35-39",       
                              ifelse(Ages_PerFamily$Paternal.age == 40, "HigherThan40",       
                              ifelse(Ages_PerFamily$Paternal.age == 41, "HigherThan40",       
                              ifelse(Ages_PerFamily$Paternal.age == 44, "HigherThan40",       
                              ifelse(Ages_PerFamily$Paternal.age == 45, "HigherThan40",
                                NA))))))))))))))))))))))))

UniqueGestAge <- unique(Ages_PerFamily$Gestational.age)
UniqueGestAge[order(UniqueGestAge)]

Ages_PerFamily$Gest_Age_Bin <- ifelse(Ages_PerFamily$Gestational.age == 4.5, "4-5", 
                               ifelse(Ages_PerFamily$Gestational.age == 5.0, "4-5",       
                               ifelse(Ages_PerFamily$Gestational.age == 5.5, "4-5",       
                               ifelse(Ages_PerFamily$Gestational.age == 6.0, "6-7", 
                               ifelse(Ages_PerFamily$Gestational.age == 6.1, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 6.2, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 6.4, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 6.5, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 6.6, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 6.7, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 7.0, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 7.1, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 7.4, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 7.5, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 7.6, "6-7",       
                               ifelse(Ages_PerFamily$Gestational.age == 8.0, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 8.1, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 8.2, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 8.4, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 8.5, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 8.6, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 9.0, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 9.1, "8-9",       
                               ifelse(Ages_PerFamily$Gestational.age == 9.2, "8-9",
                               ifelse(Ages_PerFamily$Gestational.age == 9.3, "8-9",
                               ifelse(Ages_PerFamily$Gestational.age == 9.4, "8-9",
                               ifelse(Ages_PerFamily$Gestational.age == 9.5, "8-9",
                               ifelse(Ages_PerFamily$Gestational.age == 10.3, "10-13",
                               ifelse(Ages_PerFamily$Gestational.age == 11.0, "10-13",
                               ifelse(Ages_PerFamily$Gestational.age == 12.0, "10-13",
                               ifelse(Ages_PerFamily$Gestational.age == 13.0, "10-13",
                                 NA)))))))))))))))))))))))))))))))
}


###RPL and SPL mat age
{
Ages_PerFamily$Mat_Age_Bin <- factor(Ages_PerFamily$Mat_Age_Bin, levels = c("LowerThan24","25-29","30-34","35-39","HigherThan40"))
Ages_PerFamily$Indication_SNPHapla <- factor(Ages_PerFamily$Indication_SNPHapla, levels = c("Normal","Abnormal"))
             
All_age_Mother_plot <- ggplot(data= subset(Ages_PerFamily, !is.na(Mat_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Mat_Age_Bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c(expression(""<=24),"25-29","30-34","35-39",expression("">=40))) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Maternal age (in years)") +
  ggtitle("Genome haplarithmisis") 
All_age_Mother_plot

ggsave(plot = All_age_Mother_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Mother2.pdf")
ggsave(plot = All_age_Mother_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Mother2.png")


###RPL and SPL split
{
All_age_Mother$Age_bin <- factor(All_age_Mother$Age_bin, levels = c("LowerThan24","25-29","30-34","35-39","HigherThan40"))
All_age_Mother$Indication <- factor(All_age_Mother$Indication, levels = c("Normal","Abnormal"))

plottest <- ggplot(data= subset(All_age_Mother, !is.na(Age_bin)), aes(fill =  Indication, x = Age_bin)) +
              geom_bar(position = 'fill', stat = 'count', width = 0.5) +
              facet_wrap( ~Group) +
              scale_fill_manual(values = c("#00BFC4","#F8766D")) +
              scale_y_continuous(expand = c(0,0), labels = scales::percent) +
              theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
              scale_x_discrete(labels=c(expression(""<=24),"25-29","30-34","35-39",expression("">=40))) +
              geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
              xlab("Maternal age (in years)")
plottest

ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/AgePlot/Age_Split_Mother.pdf")
ggsave("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/Figures/AgePlot/Age_Split_Mother.png")
}

###RPL mat
Ages_PerFamily_RPL <- Ages_PerFamily[which(Ages_PerFamily$Group2 == "RPL"),]

Ages_PerFamily_RPL$Mat_Age_Bin <- factor(Ages_PerFamily_RPL$Mat_Age_Bin, levels = c("LowerThan24","25-29","30-34","35-39","HigherThan40"))

All_age_Mother_RPLplot <- ggplot(data=subset(Ages_PerFamily_RPL,!is.na(Mat_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Mat_Age_Bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  ggtitle("RPL") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c(expression(""<=24),"25-29","30-34","35-39",expression("">=40))) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Maternal age (in years)")
All_age_Mother_RPLplot

ggsave(plot = All_age_Mother_RPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Mother_RPL.pdf")
ggsave(plot = All_age_Mother_RPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Mother_RPL.png")


###SPL mat 
Ages_PerFamily_SPL <- Ages_PerFamily[which(Ages_PerFamily$Group2 == "SPL"),]

Ages_PerFamily_SPL$Mat_Age_Bin <- factor(Ages_PerFamily_SPL$Mat_Age_Bin, levels = c("LowerThan24","25-29","30-34","35-39","HigherThan40"))

All_age_Mother_SPLplot <- ggplot(data= subset(Ages_PerFamily_SPL, !is.na(Mat_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Mat_Age_Bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  ggtitle("SPL") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c(expression(""<=24),"25-29","30-34","35-39",expression("">=40), "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Maternal age (in years)")
All_age_Mother_SPLplot

ggsave(plot = All_age_Mother_SPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Mother_SPL.pdf")
ggsave(plot = All_age_Mother_SPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Mother_SPL.png")
}


###RPL and SPL pat age
{
Ages_PerFamily$Pat_Age_Bin <- factor(Ages_PerFamily$Pat_Age_Bin, levels = c("LowerThan24","25-29","30-34","35-39","HigherThan40"))
Ages_PerFamily$Indication_SNPHapla <- factor(Ages_PerFamily$Indication_SNPHapla, levels = c("Normal","Abnormal"))

All_age_Father_plot <- ggplot(data= subset(Ages_PerFamily, !is.na(Pat_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Pat_Age_Bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c(expression(""<=24),"25-29","30-34","35-39",expression("">=40))) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Paternal age (in years)") +
  ggtitle("Genome haplarithmisis") 
All_age_Father_plot

ggsave(plot = All_age_Father_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Father2.pdf")
ggsave(plot = All_age_Father_plot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Father2.png")


###RPL Pat
Ages_PerFamily_RPL <- Ages_PerFamily[which(Ages_PerFamily$Group2 == "RPL"),]

Ages_PerFamily_RPL$Pat_Age_Bin <- factor(Ages_PerFamily_RPL$Pat_Age_Bin, levels = c("LowerThan24","25-29","30-34","35-39","HigherThan40"))

All_age_Father_RPLplot <- ggplot(data=subset(Ages_PerFamily_RPL,!is.na(Pat_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Pat_Age_Bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  ggtitle("RPL") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c(expression(""<=24),"25-29","30-34","35-39",expression("">=40))) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Paternal age (in years)")
All_age_Father_RPLplot

ggsave(plot = All_age_Father_RPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Father_RPL.pdf")
ggsave(plot = All_age_Father_RPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Father_RPL.png")


###SPL Pat 
Ages_PerFamily_SPL <- Ages_PerFamily[which(Ages_PerFamily$Group2 == "SPL"),]

Ages_PerFamily_SPL$Pat_Age_Bin <- factor(Ages_PerFamily_SPL$Pat_Age_Bin, levels = c("LowerThan24","25-29","30-34","35-39","HigherThan40"))

All_age_Father_SPLplot <- ggplot(data= subset(Ages_PerFamily_SPL, !is.na(Pat_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Pat_Age_Bin)) +
  geom_bar(position = 'fill', stat = 'count', width = 0.5) +
  scale_y_continuous(expand = c(0,0), labels = scales::percent) +
  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
  ggtitle("SPL") +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
  scale_x_discrete(labels=c(expression(""<=24),"25-29","30-34","35-39",expression("">=40), "NA")) +
  geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
  xlab("Paternal age (in years)")
All_age_Father_SPLplot

ggsave(plot = All_age_Father_SPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Father_SPL.pdf")
ggsave(plot = All_age_Father_SPLplot, width = 6, height = 5,
       "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Father_SPL.png")
}


###RPL and SPL Gest age
{
  Ages_PerFamily$Gest_Age_Bin <- factor(Ages_PerFamily$Gest_Age_Bin, levels = c("4-5","6-7","8-9","10-13"))
  Ages_PerFamily$Indication_SNPHapla <- factor(Ages_PerFamily$Indication_SNPHapla, levels = c("Normal","Abnormal"))
  
  All_age_Gest_plot <- ggplot(data= subset(Ages_PerFamily, !is.na(Gest_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Gest_Age_Bin)) +
    geom_bar(position = 'fill', stat = 'count', width = 0.5) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) +
    scale_fill_manual(values = c("#00BFC4","#F8766D")) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
    scale_x_discrete(labels=c("4-5","6-7","8-9","10-13")) +
    geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
    xlab("Gestational age (in weeks)") +
    ggtitle("Genome haplarithmisis") 
  All_age_Gest_plot
  
  ggsave(plot = All_age_Gest_plot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Gest2.pdf")
  ggsave(plot = All_age_Gest_plot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Gest2.png")
  
  
  ###RPL Gest
  Ages_PerFamily_RPL <- Ages_PerFamily[which(Ages_PerFamily$Group2 == "RPL"),]
  
  Ages_PerFamily_RPL$Gest_Age_Bin <- factor(Ages_PerFamily_RPL$Gest_Age_Bin, levels = c("4-5","6-7","8-9","10-13"))
  
  All_age_Gest_RPLplot <- ggplot(data=subset(Ages_PerFamily_RPL,!is.na(Gest_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Gest_Age_Bin)) +
    geom_bar(position = 'fill', stat = 'count', width = 0.5) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) +
    scale_fill_manual(values = c("#00BFC4","#F8766D")) +
    ggtitle("RPL") +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
    scale_x_discrete(labels=c("4-5","6-7","8-9","10-13")) +
    geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
    xlab("Gestational age (in weeks)")
  All_age_Gest_RPLplot
  
  ggsave(plot = All_age_Gest_RPLplot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Gest_RPL.pdf")
  ggsave(plot = All_age_Gest_RPLplot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Gest_RPL.png")
  
  
  ###SPL Gest 
  Ages_PerFamily_SPL <- Ages_PerFamily[which(Ages_PerFamily$Group2 == "SPL"),]
  
  Ages_PerFamily_SPL$Gest_Age_Bin <- factor(Ages_PerFamily_SPL$Gest_Age_Bin, levels = c("4-5","6-7","8-9","10-13"))
  
  All_age_Gest_SPLplot <- ggplot(data= subset(Ages_PerFamily_SPL, !is.na(Gest_Age_Bin)), aes(fill =  Indication_SNPHapla, x = Gest_Age_Bin)) +
    geom_bar(position = 'fill', stat = 'count', width = 0.5) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent) +
    scale_fill_manual(values = c("#00BFC4","#F8766D")) +
    ggtitle("SPL") +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black")) + 
    scale_x_discrete(labels=c("4-5","6-7","8-9","10-13")) +
    geom_text(aes(label= ..count..), stat = 'count', position = position_fill(vjust=0.5)) +
    xlab("Gestational age (in weeks)")
  All_age_Gest_SPLplot
  
  ggsave(plot = All_age_Gest_SPLplot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Gest_SPL.pdf")
  ggsave(plot = All_age_Gest_SPLplot, width = 6, height = 5,
         "/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Figures/Supplementary/Age_All_Gest_SPL.png")
}

#1
((879 + (866/(91/33))) * 100) / 1745

#2
((879 + (866*33/91))*100)/1745

#3
NewAbnormal <- (879 + (866*33/91))
NewAbnormal *100 / 1745
