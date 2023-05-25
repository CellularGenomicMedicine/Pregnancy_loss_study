###########################################################################################################################
# Author: Rick Essers 
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To calculate mosaicism degree (in percentage) of an aberration. 
#                 P1, P2, M1, M2 haplotyping values (from haplarithmisis and "1_Automated_Segmentation_PL.R") are compared
#                 to tables from Conlin, et al. 2010 that contain mosaicism degrees corresponding to specific paternal and 
#                 maternal haplotyping values. 
#                 This script is specific for aberrations affected by maternal contamination. This script can compensate 
#                 for maternal contamination by utilizing input haplotyping values. 

# input: P1, P2, M1, M2 haplotyping values (from haplarithmisis and "1_Automated_Segmentation_PL.R").
#        Input concerning the type of aberration ("Gain", "Loss"), and parental contribution to the aberration 
#        ("Father" or "Mother) is required. Indicate the parent with the highest chromosomal contribution to the genomic 
#        location of interest, e.g. trisomy; 2 chromosomes from the mother, 1 from the father --> "Mother", 
#        or for monosomy; 1 chromosome from the father, 0 from the mother --> "Father"

# output: A percentage of mosaicism corresponding to parental haplotyping values of the genomic location of interest. 

###########################################################################################################################

rm(list=ls(all=T))
options(digits = 16)

###Insert the P1, P2, M1, M2 values for the genomic location of interest, these values are extracted from "1_Automated_Segmentation_PL.R" script.
P1value <- "0.00735407103794966"
P2value <- "0.443787887368421"
M1value <- "0.0192494415618391"
M2value <- "0.5643018808"

### Aberration type can be either "Gain" or "Loss". 
AberrationType <- "Gain"

### More contibution to genome from mother or father?
Parent <- "Father"

###Compensate for maternal contamination (ONLY use when mat cont is present, otherwise you unnecessarily compensate)
P2value <- as.numeric(P2value) - as.numeric(M1value) 


if (AberrationType == "Loss"){

###When aberration is a loss

##Load mosaicism percentage reference files from Conlin et al. 2010.
MosDet_Loss_PatMat <- read.csv("/Projects/PregnancyLoss/Scripts/Haplarithmisis/Mosaicism_Detection/Refs/MosDet_Loss.csv", row.names=1)
    
  if (Parent == "Mother"){
  
    zeroAB_ABB <- MosDet_Loss_PatMat[3,1] - as.numeric(P2value)
    zeroAB_ABA <- MosDet_Loss_PatMat[2,1] - as.numeric(M2value)
    zeroTotal <- abs(zeroAB_ABA) + abs(zeroAB_ABB)
  
    fiveAB_ABB <- MosDet_Loss_PatMat[3,2] - as.numeric(P2value)
    fiveAB_ABA <- MosDet_Loss_PatMat[2,2] - as.numeric(M2value)
    fiveTotal <- abs(fiveAB_ABA) + abs(fiveAB_ABB)
  
    tenAB_ABB <- MosDet_Loss_PatMat[3,3] - as.numeric(P2value)
    tenAB_ABA <- MosDet_Loss_PatMat[2,3] - as.numeric(M2value)
    tenTotal <- abs(tenAB_ABA) + abs(tenAB_ABB)
  
    fifteenAB_ABB <- MosDet_Loss_PatMat[3,4] - as.numeric(P2value)
    fifteenAB_ABA <- MosDet_Loss_PatMat[2,4] - as.numeric(M2value)
    fifteenTotal <- abs(fifteenAB_ABA) + abs(fifteenAB_ABB)
  
    twentyAB_ABB <- MosDet_Loss_PatMat[3,5] - as.numeric(P2value)
    twentyAB_ABA <- MosDet_Loss_PatMat[2,5] - as.numeric(M2value)
    twentyTotal <- abs(twentyAB_ABA) + abs(twentyAB_ABB)
  
    twentyfiveAB_ABB <- MosDet_Loss_PatMat[3,6] - as.numeric(P2value)
    twentyfiveAB_ABA <- MosDet_Loss_PatMat[2,6] - as.numeric(M2value)
    twentyfiveTotal <- abs(twentyfiveAB_ABA) + abs(twentyfiveAB_ABB)
  
    thirtyAB_ABB <- MosDet_Loss_PatMat[3,7] - as.numeric(P2value)
    thirtyAB_ABA <- MosDet_Loss_PatMat[2,7] - as.numeric(M2value)
    thirtyTotal <- abs(thirtyAB_ABA) + abs(thirtyAB_ABB)
  
    thirtyfiveAB_ABB <- MosDet_Loss_PatMat[3,8] - as.numeric(P2value)
    thirtyfiveAB_ABA <- MosDet_Loss_PatMat[2,8] - as.numeric(M2value)
    thirtyfiveTotal <- abs(thirtyfiveAB_ABA) + abs(thirtyfiveAB_ABB)
  
    fortyAB_ABB <- MosDet_Loss_PatMat[3,9] - as.numeric(P2value)
    fortyAB_ABA <- MosDet_Loss_PatMat[2,9] - as.numeric(M2value)
    fortyTotal <- abs(fortyAB_ABA) + abs(fortyAB_ABB)
  
    fortyfiveAB_ABB <- MosDet_Loss_PatMat[3,10] - as.numeric(P2value)
    fortyfiveAB_ABA <- MosDet_Loss_PatMat[2,10] - as.numeric(M2value)
    fortyfiveTotal <- abs(fortyfiveAB_ABA) + abs(fortyfiveAB_ABB)
  
    fiftyAB_ABB <- MosDet_Loss_PatMat[3,11] - as.numeric(P2value)
    fiftyAB_ABA <- MosDet_Loss_PatMat[2,11] - as.numeric(M2value)
    fiftyTotal <- abs(fiftyAB_ABA) + abs(fiftyAB_ABB)
  
    fiftyfiveAB_ABB <- MosDet_Loss_PatMat[3,12] - as.numeric(P2value)
    fiftyfiveAB_ABA <- MosDet_Loss_PatMat[2,12] - as.numeric(M2value)
    fiftyfiveTotal <- abs(fiftyfiveAB_ABA) + abs(fiftyfiveAB_ABB)
  
    sixtyAB_ABB <- MosDet_Loss_PatMat[3,13] - as.numeric(P2value)
    sixtyAB_ABA <- MosDet_Loss_PatMat[2,13] - as.numeric(M2value)
    sixtyTotal <- abs(sixtyAB_ABA) + abs(sixtyAB_ABB)
  
    sixtyfiveAB_ABB <- MosDet_Loss_PatMat[3,14] - as.numeric(P2value)
    sixtyfiveAB_ABA <- MosDet_Loss_PatMat[2,14] - as.numeric(M2value)
    sixtyfiveTotal <- abs(sixtyfiveAB_ABA) + abs(sixtyfiveAB_ABB)
  
    seventyAB_ABB <- MosDet_Loss_PatMat[3,15] - as.numeric(P2value)
    seventyAB_ABA <- MosDet_Loss_PatMat[2,15] - as.numeric(M2value)
    seventyTotal <- abs(seventyAB_ABA) + abs(seventyAB_ABB)
  
    seventyfiveAB_ABB <- MosDet_Loss_PatMat[3,16] - as.numeric(P2value)
    seventyfiveAB_ABA <- MosDet_Loss_PatMat[2,16] - as.numeric(M2value)
    seventyfiveTotal <- abs(seventyfiveAB_ABA) + abs(seventyfiveAB_ABB)
  
    eightyAB_ABB <- MosDet_Loss_PatMat[3,17] - as.numeric(P2value)
    eightyAB_ABA <- MosDet_Loss_PatMat[2,17] - as.numeric(M2value)
    eightyTotal <- abs(eightyAB_ABA) + abs(eightyAB_ABB)
  
    eightyfiveAB_ABB <- MosDet_Loss_PatMat[3,18] - as.numeric(P2value)
    eightyfiveAB_ABA <- MosDet_Loss_PatMat[2,18] - as.numeric(M2value)
    eightyfiveTotal <- abs(eightyfiveAB_ABA) + abs(eightyfiveAB_ABB)
  
    ninetyAB_ABB <- MosDet_Loss_PatMat[3,19] - as.numeric(P2value)
    ninetyAB_ABA <- MosDet_Loss_PatMat[2,19] - as.numeric(M2value)
    ninetyTotal <- abs(ninetyAB_ABA) + abs(ninetyAB_ABB)
  
    ninetyfiveAB_ABB <- MosDet_Loss_PatMat[3,20] - as.numeric(P2value)
    ninetyfiveAB_ABA <- MosDet_Loss_PatMat[2,20] - as.numeric(M2value)
    ninetyfiveTotal <- abs(ninetyfiveAB_ABA) + abs(ninetyfiveAB_ABB)
  
    onehundredAB_ABB <- MosDet_Loss_PatMat[3,21] - as.numeric(P2value)
    onehundredAB_ABA <- MosDet_Loss_PatMat[2,21] - as.numeric(M2value)
    onehundredTotal <- abs(onehundredAB_ABA) + abs(onehundredAB_ABB)
  
    ######################
  
    CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
    Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
    Percentage
  
  } else if (Parent == "Father"){
  
  zeroAB_ABB <- MosDet_Loss_PatMat[3,1] - as.numeric(M2value)
  zeroAB_ABA <- MosDet_Loss_PatMat[2,1] - as.numeric(P2value)
  zeroTotal <- abs(zeroAB_ABA) + abs(zeroAB_ABB)
  
  fiveAB_ABB <- MosDet_Loss_PatMat[3,2] - as.numeric(M2value)
  fiveAB_ABA <- MosDet_Loss_PatMat[2,2] - as.numeric(P2value)
  fiveTotal <- abs(fiveAB_ABA) + abs(fiveAB_ABB)
  
  tenAB_ABB <- MosDet_Loss_PatMat[3,3] - as.numeric(M2value)
  tenAB_ABA <- MosDet_Loss_PatMat[2,3] - as.numeric(P2value)
  tenTotal <- abs(tenAB_ABA) + abs(tenAB_ABB)
  
  fifteenAB_ABB <- MosDet_Loss_PatMat[3,4] - as.numeric(M2value)
  fifteenAB_ABA <- MosDet_Loss_PatMat[2,4] - as.numeric(P2value)
  fifteenTotal <- abs(fifteenAB_ABA) + abs(fifteenAB_ABB)
  
  twentyAB_ABB <- MosDet_Loss_PatMat[3,5] - as.numeric(M2value)
  twentyAB_ABA <- MosDet_Loss_PatMat[2,5] - as.numeric(P2value)
  twentyTotal <- abs(twentyAB_ABA) + abs(twentyAB_ABB)
  
  twentyfiveAB_ABB <- MosDet_Loss_PatMat[3,6] - as.numeric(M2value)
  twentyfiveAB_ABA <- MosDet_Loss_PatMat[2,6] - as.numeric(P2value)
  twentyfiveTotal <- abs(twentyfiveAB_ABA) + abs(twentyfiveAB_ABB)
  
  thirtyAB_ABB <- MosDet_Loss_PatMat[3,7] - as.numeric(M2value)
  thirtyAB_ABA <- MosDet_Loss_PatMat[2,7] - as.numeric(P2value)
  thirtyTotal <- abs(thirtyAB_ABA) + abs(thirtyAB_ABB)
  
  thirtyfiveAB_ABB <- MosDet_Loss_PatMat[3,8] - as.numeric(M2value)
  thirtyfiveAB_ABA <- MosDet_Loss_PatMat[2,8] - as.numeric(P2value)
  thirtyfiveTotal <- abs(thirtyfiveAB_ABA) + abs(thirtyfiveAB_ABB)
  
  fortyAB_ABB <- MosDet_Loss_PatMat[3,9] - as.numeric(M2value)
  fortyAB_ABA <- MosDet_Loss_PatMat[2,9] - as.numeric(P2value)
  fortyTotal <- abs(fortyAB_ABA) + abs(fortyAB_ABB)
  
  fortyfiveAB_ABB <- MosDet_Loss_PatMat[3,10] - as.numeric(M2value)
  fortyfiveAB_ABA <- MosDet_Loss_PatMat[2,10] - as.numeric(P2value)
  fortyfiveTotal <- abs(fortyfiveAB_ABA) + abs(fortyfiveAB_ABB)
  
  fiftyAB_ABB <- MosDet_Loss_PatMat[3,11] - as.numeric(M2value)
  fiftyAB_ABA <- MosDet_Loss_PatMat[2,11] - as.numeric(P2value)
  fiftyTotal <- abs(fiftyAB_ABA) + abs(fiftyAB_ABB)
  
  fiftyfiveAB_ABB <- MosDet_Loss_PatMat[3,12] - as.numeric(M2value)
  fiftyfiveAB_ABA <- MosDet_Loss_PatMat[2,12] - as.numeric(P2value)
  fiftyfiveTotal <- abs(fiftyfiveAB_ABA) + abs(fiftyfiveAB_ABB)
  
  sixtyAB_ABB <- MosDet_Loss_PatMat[3,13] - as.numeric(M2value)
  sixtyAB_ABA <- MosDet_Loss_PatMat[2,13] - as.numeric(P2value)
  sixtyTotal <- abs(sixtyAB_ABA) + abs(sixtyAB_ABB)
  
  sixtyfiveAB_ABB <- MosDet_Loss_PatMat[3,14] - as.numeric(M2value)
  sixtyfiveAB_ABA <- MosDet_Loss_PatMat[2,14] - as.numeric(P2value)
  sixtyfiveTotal <- abs(sixtyfiveAB_ABA) + abs(sixtyfiveAB_ABB)
  
  seventyAB_ABB <- MosDet_Loss_PatMat[3,15] - as.numeric(M2value)
  seventyAB_ABA <- MosDet_Loss_PatMat[2,15] - as.numeric(P2value)
  seventyTotal <- abs(seventyAB_ABA) + abs(seventyAB_ABB)
  
  seventyfiveAB_ABB <- MosDet_Loss_PatMat[3,16] - as.numeric(M2value)
  seventyfiveAB_ABA <- MosDet_Loss_PatMat[2,16] - as.numeric(P2value)
  seventyfiveTotal <- abs(seventyfiveAB_ABA) + abs(seventyfiveAB_ABB)
  
  eightyAB_ABB <- MosDet_Loss_PatMat[3,17] - as.numeric(M2value)
  eightyAB_ABA <- MosDet_Loss_PatMat[2,17] - as.numeric(P2value)
  eightyTotal <- abs(eightyAB_ABA) + abs(eightyAB_ABB)
  
  eightyfiveAB_ABB <- MosDet_Loss_PatMat[3,18] - as.numeric(M2value)
  eightyfiveAB_ABA <- MosDet_Loss_PatMat[2,18] - as.numeric(P2value)
  eightyfiveTotal <- abs(eightyfiveAB_ABA) + abs(eightyfiveAB_ABB)
  
  ninetyAB_ABB <- MosDet_Loss_PatMat[3,19] - as.numeric(M2value)
  ninetyAB_ABA <- MosDet_Loss_PatMat[2,19] - as.numeric(P2value)
  ninetyTotal <- abs(ninetyAB_ABA) + abs(ninetyAB_ABB)
  
  ninetyfiveAB_ABB <- MosDet_Loss_PatMat[3,20] - as.numeric(M2value)
  ninetyfiveAB_ABA <- MosDet_Loss_PatMat[2,20] - as.numeric(P2value)
  ninetyfiveTotal <- abs(ninetyfiveAB_ABA) + abs(ninetyfiveAB_ABB)
  
  onehundredAB_ABB <- MosDet_Loss_PatMat[3,21] - as.numeric(M2value)
  onehundredAB_ABA <- MosDet_Loss_PatMat[2,21] - as.numeric(P2value)
  onehundredTotal <- abs(onehundredAB_ABA) + abs(onehundredAB_ABB)
  
  ######################
  
  CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
  Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
  Percentage
  
  }#end of mother or father 
} else if (AberrationType == "Gain"){

###When aberration is a gain 

##Load mosaicism percentage reference files from Conlin et al. 2010.
MosDet_Gain_PatMat <- read.delim("/Projects/PregnancyLoss/Scripts/Haplarithmisis/Mosaicism_Detection/Refs/MosDet_Gain_4_bands.csv", row.names=1)
  
    if (Parent == "Mother"){
    
    zeroAB_ABB <- MosDet_Gain_PatMat[3,1] - as.numeric(P2value)
    zeroAB_ABA <- MosDet_Gain_PatMat[2,1] - as.numeric(M2value)
    zeroTotal <- abs(zeroAB_ABA) + abs(zeroAB_ABB)
    
    fiveAB_ABB <- MosDet_Gain_PatMat[3,2] - as.numeric(P2value)
    fiveAB_ABA <- MosDet_Gain_PatMat[2,2] - as.numeric(M2value)
    fiveTotal <- abs(fiveAB_ABA) + abs(fiveAB_ABB)
    
    tenAB_ABB <- MosDet_Gain_PatMat[3,3] - as.numeric(P2value)
    tenAB_ABA <- MosDet_Gain_PatMat[2,3] - as.numeric(M2value)
    tenTotal <- abs(tenAB_ABA) + abs(tenAB_ABB)
    
    fifteenAB_ABB <- MosDet_Gain_PatMat[3,4] - as.numeric(P2value)
    fifteenAB_ABA <- MosDet_Gain_PatMat[2,4] - as.numeric(M2value)
    fifteenTotal <- abs(fifteenAB_ABA) + abs(fifteenAB_ABB)
    
    twentyAB_ABB <- MosDet_Gain_PatMat[3,5] - as.numeric(P2value)
    twentyAB_ABA <- MosDet_Gain_PatMat[2,5] - as.numeric(M2value)
    twentyTotal <- abs(twentyAB_ABA) + abs(twentyAB_ABB)
    
    twentyfiveAB_ABB <- MosDet_Gain_PatMat[3,6] - as.numeric(P2value)
    twentyfiveAB_ABA <- MosDet_Gain_PatMat[2,6] - as.numeric(M2value)
    twentyfiveTotal <- abs(twentyfiveAB_ABA) + abs(twentyfiveAB_ABB)
    
    thirtyAB_ABB <- MosDet_Gain_PatMat[3,7] - as.numeric(P2value)
    thirtyAB_ABA <- MosDet_Gain_PatMat[2,7] - as.numeric(M2value)
    thirtyTotal <- abs(thirtyAB_ABA) + abs(thirtyAB_ABB)
    
    thirtyfiveAB_ABB <- MosDet_Gain_PatMat[3,8] - as.numeric(P2value)
    thirtyfiveAB_ABA <- MosDet_Gain_PatMat[2,8] - as.numeric(M2value)
    thirtyfiveTotal <- abs(thirtyfiveAB_ABA) + abs(thirtyfiveAB_ABB)
    
    fortyAB_ABB <- MosDet_Gain_PatMat[3,9] - as.numeric(P2value)
    fortyAB_ABA <- MosDet_Gain_PatMat[2,9] - as.numeric(M2value)
    fortyTotal <- abs(fortyAB_ABA) + abs(fortyAB_ABB)
    
    fortyfiveAB_ABB <- MosDet_Gain_PatMat[3,10] - as.numeric(P2value)
    fortyfiveAB_ABA <- MosDet_Gain_PatMat[2,10] - as.numeric(M2value)
    fortyfiveTotal <- abs(fortyfiveAB_ABA) + abs(fortyfiveAB_ABB)
    
    fiftyAB_ABB <- MosDet_Gain_PatMat[3,11] - as.numeric(P2value)
    fiftyAB_ABA <- MosDet_Gain_PatMat[2,11] - as.numeric(M2value)
    fiftyTotal <- abs(fiftyAB_ABA) + abs(fiftyAB_ABB)
    
    fiftyfiveAB_ABB <- MosDet_Gain_PatMat[3,12] - as.numeric(P2value)
    fiftyfiveAB_ABA <- MosDet_Gain_PatMat[2,12] - as.numeric(M2value)
    fiftyfiveTotal <- abs(fiftyfiveAB_ABA) + abs(fiftyfiveAB_ABB)
    
    sixtyAB_ABB <- MosDet_Gain_PatMat[3,13] - as.numeric(P2value)
    sixtyAB_ABA <- MosDet_Gain_PatMat[2,13] - as.numeric(M2value)
    sixtyTotal <- abs(sixtyAB_ABA) + abs(sixtyAB_ABB)
    
    sixtyfiveAB_ABB <- MosDet_Gain_PatMat[3,14] - as.numeric(P2value)
    sixtyfiveAB_ABA <- MosDet_Gain_PatMat[2,14] - as.numeric(M2value)
    sixtyfiveTotal <- abs(sixtyfiveAB_ABA) + abs(sixtyfiveAB_ABB)
    
    seventyAB_ABB <- MosDet_Gain_PatMat[3,15] - as.numeric(P2value)
    seventyAB_ABA <- MosDet_Gain_PatMat[2,15] - as.numeric(M2value)
    seventyTotal <- abs(seventyAB_ABA) + abs(seventyAB_ABB)
    
    seventyfiveAB_ABB <- MosDet_Gain_PatMat[3,16] - as.numeric(P2value)
    seventyfiveAB_ABA <- MosDet_Gain_PatMat[2,16] - as.numeric(M2value)
    seventyfiveTotal <- abs(seventyfiveAB_ABA) + abs(seventyfiveAB_ABB)
    
    eightyAB_ABB <- MosDet_Gain_PatMat[3,17] - as.numeric(P2value)
    eightyAB_ABA <- MosDet_Gain_PatMat[2,17] - as.numeric(M2value)
    eightyTotal <- abs(eightyAB_ABA) + abs(eightyAB_ABB)
    
    eightyfiveAB_ABB <- MosDet_Gain_PatMat[3,18] - as.numeric(P2value)
    eightyfiveAB_ABA <- MosDet_Gain_PatMat[2,18] - as.numeric(M2value)
    eightyfiveTotal <- abs(eightyfiveAB_ABA) + abs(eightyfiveAB_ABB)
    
    ninetyAB_ABB <- MosDet_Gain_PatMat[3,19] - as.numeric(P2value)
    ninetyAB_ABA <- MosDet_Gain_PatMat[2,19] - as.numeric(M2value)
    ninetyTotal <- abs(ninetyAB_ABA) + abs(ninetyAB_ABB)
    
    ninetyfiveAB_ABB <- MosDet_Gain_PatMat[3,20] - as.numeric(P2value)
    ninetyfiveAB_ABA <- MosDet_Gain_PatMat[2,20] - as.numeric(M2value)
    ninetyfiveTotal <- abs(ninetyfiveAB_ABA) + abs(ninetyfiveAB_ABB)
    
    onehundredAB_ABB <- MosDet_Gain_PatMat[3,21] - as.numeric(P2value)
    onehundredAB_ABA <- MosDet_Gain_PatMat[2,21] - as.numeric(M2value)
    onehundredTotal <- abs(onehundredAB_ABA) + abs(onehundredAB_ABB)
    
    ######################
    
    CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
    Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
    Percentage
    
  } else if (Parent == "Father"){
    
    zeroAB_ABB <- MosDet_Gain_PatMat[3,1] - as.numeric(M2value)
    zeroAB_ABA <- MosDet_Gain_PatMat[2,1] - as.numeric(P2value)
    zeroTotal <- abs(zeroAB_ABA) + abs(zeroAB_ABB)
    
    fiveAB_ABB <- MosDet_Gain_PatMat[3,2] - as.numeric(M2value)
    fiveAB_ABA <- MosDet_Gain_PatMat[2,2] - as.numeric(P2value)
    fiveTotal <- abs(fiveAB_ABA) + abs(fiveAB_ABB)
    
    tenAB_ABB <- MosDet_Gain_PatMat[3,3] - as.numeric(M2value)
    tenAB_ABA <- MosDet_Gain_PatMat[2,3] - as.numeric(P2value)
    tenTotal <- abs(tenAB_ABA) + abs(tenAB_ABB)
    
    fifteenAB_ABB <- MosDet_Gain_PatMat[3,4] - as.numeric(M2value)
    fifteenAB_ABA <- MosDet_Gain_PatMat[2,4] - as.numeric(P2value)
    fifteenTotal <- abs(fifteenAB_ABA) + abs(fifteenAB_ABB)
    
    twentyAB_ABB <- MosDet_Gain_PatMat[3,5] - as.numeric(M2value)
    twentyAB_ABA <- MosDet_Gain_PatMat[2,5] - as.numeric(P2value)
    twentyTotal <- abs(twentyAB_ABA) + abs(twentyAB_ABB)
    
    twentyfiveAB_ABB <- MosDet_Gain_PatMat[3,6] - as.numeric(M2value)
    twentyfiveAB_ABA <- MosDet_Gain_PatMat[2,6] - as.numeric(P2value)
    twentyfiveTotal <- abs(twentyfiveAB_ABA) + abs(twentyfiveAB_ABB)
    
    thirtyAB_ABB <- MosDet_Gain_PatMat[3,7] - as.numeric(M2value)
    thirtyAB_ABA <- MosDet_Gain_PatMat[2,7] - as.numeric(P2value)
    thirtyTotal <- abs(thirtyAB_ABA) + abs(thirtyAB_ABB)
    
    thirtyfiveAB_ABB <- MosDet_Gain_PatMat[3,8] - as.numeric(M2value)
    thirtyfiveAB_ABA <- MosDet_Gain_PatMat[2,8] - as.numeric(P2value)
    thirtyfiveTotal <- abs(thirtyfiveAB_ABA) + abs(thirtyfiveAB_ABB)
    
    fortyAB_ABB <- MosDet_Gain_PatMat[3,9] - as.numeric(M2value)
    fortyAB_ABA <- MosDet_Gain_PatMat[2,9] - as.numeric(P2value)
    fortyTotal <- abs(fortyAB_ABA) + abs(fortyAB_ABB)
    
    fortyfiveAB_ABB <- MosDet_Gain_PatMat[3,10] - as.numeric(M2value)
    fortyfiveAB_ABA <- MosDet_Gain_PatMat[2,10] - as.numeric(P2value)
    fortyfiveTotal <- abs(fortyfiveAB_ABA) + abs(fortyfiveAB_ABB)
    
    fiftyAB_ABB <- MosDet_Gain_PatMat[3,11] - as.numeric(M2value)
    fiftyAB_ABA <- MosDet_Gain_PatMat[2,11] - as.numeric(P2value)
    fiftyTotal <- abs(fiftyAB_ABA) + abs(fiftyAB_ABB)
    
    fiftyfiveAB_ABB <- MosDet_Gain_PatMat[3,12] - as.numeric(M2value)
    fiftyfiveAB_ABA <- MosDet_Gain_PatMat[2,12] - as.numeric(P2value)
    fiftyfiveTotal <- abs(fiftyfiveAB_ABA) + abs(fiftyfiveAB_ABB)
    
    sixtyAB_ABB <- MosDet_Gain_PatMat[3,13] - as.numeric(M2value)
    sixtyAB_ABA <- MosDet_Gain_PatMat[2,13] - as.numeric(P2value)
    sixtyTotal <- abs(sixtyAB_ABA) + abs(sixtyAB_ABB)
    
    sixtyfiveAB_ABB <- MosDet_Gain_PatMat[3,14] - as.numeric(M2value)
    sixtyfiveAB_ABA <- MosDet_Gain_PatMat[2,14] - as.numeric(P2value)
    sixtyfiveTotal <- abs(sixtyfiveAB_ABA) + abs(sixtyfiveAB_ABB)
    
    seventyAB_ABB <- MosDet_Gain_PatMat[3,15] - as.numeric(M2value)
    seventyAB_ABA <- MosDet_Gain_PatMat[2,15] - as.numeric(P2value)
    seventyTotal <- abs(seventyAB_ABA) + abs(seventyAB_ABB)
    
    seventyfiveAB_ABB <- MosDet_Gain_PatMat[3,16] - as.numeric(M2value)
    seventyfiveAB_ABA <- MosDet_Gain_PatMat[2,16] - as.numeric(P2value)
    seventyfiveTotal <- abs(seventyfiveAB_ABA) + abs(seventyfiveAB_ABB)
    
    eightyAB_ABB <- MosDet_Gain_PatMat[3,17] - as.numeric(M2value)
    eightyAB_ABA <- MosDet_Gain_PatMat[2,17] - as.numeric(P2value)
    eightyTotal <- abs(eightyAB_ABA) + abs(eightyAB_ABB)
    
    eightyfiveAB_ABB <- MosDet_Gain_PatMat[3,18] - as.numeric(M2value)
    eightyfiveAB_ABA <- MosDet_Gain_PatMat[2,18] - as.numeric(P2value)
    eightyfiveTotal <- abs(eightyfiveAB_ABA) + abs(eightyfiveAB_ABB)
    
    ninetyAB_ABB <- MosDet_Gain_PatMat[3,19] - as.numeric(M2value)
    ninetyAB_ABA <- MosDet_Gain_PatMat[2,19] - as.numeric(P2value)
    ninetyTotal <- abs(ninetyAB_ABA) + abs(ninetyAB_ABB)
    
    ninetyfiveAB_ABB <- MosDet_Gain_PatMat[3,20] - as.numeric(M2value)
    ninetyfiveAB_ABA <- MosDet_Gain_PatMat[2,20] - as.numeric(P2value)
    ninetyfiveTotal <- abs(ninetyfiveAB_ABA) + abs(ninetyfiveAB_ABB)
    
    onehundredAB_ABB <- MosDet_Gain_PatMat[3,21] - as.numeric(M2value)
    onehundredAB_ABA <- MosDet_Gain_PatMat[2,21] - as.numeric(P2value)
    onehundredTotal <- abs(onehundredAB_ABA) + abs(onehundredAB_ABB)
    
    ######################
    
    CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
    Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
    Percentage
    
  }#End of parental split. 
}#End of aberration type split
  
