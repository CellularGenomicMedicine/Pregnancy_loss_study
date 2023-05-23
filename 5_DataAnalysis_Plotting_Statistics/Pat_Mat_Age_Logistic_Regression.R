
rm(list=ls(all=T))

library(ggplot2)

#Load datafile
PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")
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

PerFamily <- read.csv2("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Data/20221027_Miscarr_PerFamily.csv")
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

