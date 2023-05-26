
###Statistics
{
  ###Preparing input files for statistics, ConvKaryo vs SNPHapla, and Normal vs Abnormal within the groups. 
  Demographic_data_TOTAL <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
  
  
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


### logistic regression analysis paternal and maternal age 
rm(list=ls(all=T))

library(ggplot2)

#Load datafile
PerFamily <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
#Select cases included in the haplarithmisis cohort (n = 91)
PerFamily <- subset(PerFamily, Group == "SNPHapla" & Indication_SNPHapla!= "Excluded")
#Split based on RPL or SPL 
PerFamily_RPL <- subset(PerFamily, Group2 == "RPL")
PerFamily_SPL <- subset(PerFamily, Group2 == "SPL")

#Summarize the number of PLs per group 
summary(PerFamily_RPL$Total_PLs)
summary(PerFamily_SPL$Total_PLs)

##################

rm(list=ls(all=T))

PerFamily <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerPOC_Source_data.csv")
PerFamily <- subset(PerFamily, Group == "ConvKaryo" | Group == "SNPHapla")

###Multivariate logistic regression
fit <- glm(Indication_ConvKaryo_.Normal_0_Abnormal_1.~Paternal.age+Maternal.age+Gestational.age, data=PerFamily, family = binomial(), na.action = na.omit)
summary(fit)
pchisq(4.1,3, lower.tail = F)
confint(fit)
exp(coef(fit))
exp(confint(fit))

fit <- glm(Indication_ConvKaryo_.Normal_0_Abnormal_1.~Paternal.age+Maternal.age, data=PerFamily, family = binomial(), na.action = na.omit)
summary(fit)
pchisq(9.1,2, lower.tail = F)
confint(fit)
exp(coef(fit))
exp(confint(fit))

fit <- glm(Indication_ConvKaryo_.Normal_0_Abnormal_1.~Paternal.age, data=PerFamily, family = binomial(), na.action = na.omit)
summary(fit)
pchisq(7.9,1, lower.tail = F)
confint(fit)
exp(coef(fit))
exp(confint(fit))

fit <- glm(Indication_ConvKaryo_.Normal_0_Abnormal_1.~Maternal.age, data=PerFamily, family = binomial(), na.action = na.omit)
summary(fit)
pchisq(15.9,1, lower.tail = F)
confint(fit)
exp(coef(fit))
exp(confint(fit))

fit <- glm(Indication_ConvKaryo_.Normal_0_Abnormal_1.~Gestational.age, data=PerFamily, family = binomial(), na.action = na.omit)
summary(fit)
pchisq(0.2,1, lower.tail = F)
confint(fit)
exp(coef(fit))
exp(confint(fit))


###t.test
PerFamily_ConvKaryo_Normal    <- subset(PerFamily, Indication_ConvKaryo == "Normal")
PerFamily_ConvKaryo_Abnormal  <- subset(PerFamily, Indication_ConvKaryo == "Abnormal")

t.test(PerFamily_ConvKaryo_Normal$Maternal.age, PerFamily_ConvKaryo_Abnormal$Maternal.age)



#RPL vs SPL
###Statistics###
#PatMat origin all aberrations
  rm(list=ls(all=T))
  
  Affected_Chromosomes <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
  Affected_Chromosomes_RPL <- subset(Affected_Chromosomes, Group == "RPL")
  Affected_Chromosomes_SPL <- subset(Affected_Chromosomes, Group == "SPL")
  table(Affected_Chromosomes_RPL$Parental_error_origin)
  table(Affected_Chromosomes_SPL$Parental_error_origin)
  
  data <- data.frame("Maternal" = c(16,11), "Paternal" = c(4,7), row.names = c("RPL","SPL"), stringsAsFactors = F)
  data
  
  #Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
  chisq.test(data)$expected
  #Maternal Paternal
  #RPL 13.51351 6.486486
  #SPL 11.48649 5.513514
  ###Conclusion = chisquared is appropriate
  
  #With Yates' continuity correction
  chisq.test(data)
  #Pearson's Chi-squared test with Yates' continuity correction
  #data:  data
  #X-squared = 1.9596, df = 1, p-value = 0.1616
  



###Statistics###
#Mitotic/meiotic only for all aberrations 
  rm(list=ls(all=T))
  
  Affected_Chromosomes <- read.csv("/Projects/PregnancyLoss/Source_data/PL_PerAberration_Source_data.csv", sep=";")
  
  Affected_Chromosomes_RPL <- subset(Affected_Chromosomes, Group == "RPL")
  Affected_Chromosomes_SPL <- subset(Affected_Chromosomes, Group == "SPL")
  table(Affected_Chromosomes_RPL$Seg_origin2)
  table(Affected_Chromosomes_SPL$Seg_origin2)
  
  data <- data.frame("Mitotic" = c(7,11), "Meiotic" = c(11,5), row.names = c("RPL","SPL"), stringsAsFactors = F)
  data
  
  ###Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
  chisq.test(data)$expected
  #Mitotic Meiotic
  #RPL    5.76   10.24
  #SPL    3.24    5.76
  ###Conclusion = chisquared approximation is appropriate
  
  ###Chi-squared With Yates' continuity correction
  chisq.test(data)
  #Pearson's Chi-squared test with Yates' continuity correction
  
  #data:  data
  #X-squared = 2.378, df = 1, p-value = 0.1231
  
  
  #Mitotic/meiotic only for only trisomy cases
  
  MitMei_TrisOnly <- subset(Affected_Chromosomes, Indication == "Trisomy")
  MitMei_TrisOnly_RPL <- subset(MitMei_TrisOnly, Group == "RPL")
  MitMei_TrisOnly_SPL <- subset(MitMei_TrisOnly, Group == "SPL")
  table(MitMei_TrisOnly_RPL$Seg_origin2)
  table(MitMei_TrisOnly_SPL$Seg_origin2)
  
  data <- data.frame("Mitotic" = c(4,4), "Meiotic" = c(10,5), row.names = c("RPL","SPL"), stringsAsFactors = F)
  data
  
  ###Check if all expected values are >5 (If not Fisher's exact test is appropriate instead of chisquared)
  chisq.test(data)$expected

  fisher.test(data)

  


### Statistics for mosaicism difference between extraembryonic mesoderm and chorionic villi DNAs of miscarried POCs 
rm(list=ls(all=T))

PerTissue <- read.csv2("/Projects/PregnancyLoss/Source_data/PL_PerTissue_Source_data.csv")

PerTissue_BioRel <- subset(PerTissue, Biologically_Relevant_to_Include_MosDiff == "y")


###Statistics mosaicism percentage and CV vs EM 
PerTissue_CV <- subset(PerTissue_BioRel, Group2 == "CV")
PerTissue_EM <- subset(PerTissue_BioRel, Group2 == "EM")

summary(PerTissue_EM$Mosaicism_Perc_Av, na.rm = T)
sd(PerTissue_EM$Mosaicism_Perc_Av, na.rm = T)

summary(PerTissue_CV$Mosaicism_Perc_Av, na.rm = T)
sd(PerTissue_CV$Mosaicism_Perc_Av, na.rm = T)

wilcox.test(PerTissue_CV$Mosaicism_Perc_Av, PerTissue_EM$Mosaicism_Perc_Av, paired = T)

