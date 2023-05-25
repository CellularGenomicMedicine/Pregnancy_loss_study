###########################################################################################################################
# Author: Rick Essers 
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University Medical Center (MUMC+)

# script purpose: To calculate mosaicism degree (in percentage) of an aberration. 
#                 Mean BAF values are compared to tables from Conlin, et al. 2010 that contain mosaicism degrees  
#                 corresponding to specific BAF values. 
#                 This script is specific for aberrations that show 4 bands in the BAF values and can be used for 
#                 gain, loss and neutral copy number aberrations.  

# input: Mean BAF values for each band (4 total). 
#        Input concerning the type of aberration ("Gain", "Loss", "Neutral") is required 

# output: A percentage of mosaicism corresponding to mean BAF values of the genomic location of interest. 

###########################################################################################################################

rm(list=ls(all=T))
options(digits = 16)

###BAF - 4 bands
Band1   <- "0.981111262849162"
Band2   <- "0.668406166666667"
Band3   <- "0.371275215189873"
Band4   <- "0.0177244540798319"

### Aberration type can be either "Gain", "Loss", or "Neutral. 
AberrationType <- "Neutral"


if (AberrationType == "Gain"){

  ##Load mosaicism percentage reference files from Conlin et al. 2010.
  
  MosDet_Gain_4_bands_BAF <- read.delim("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Github/4_MosaicismDetection/Refs/MosDet_Gain_4_bands.csv", row.names=1)
  
  zeroAA_AAA <- MosDet_Gain_4_bands_BAF[1,1] - as.numeric(Band4)
  zeroAB_ABA <- MosDet_Gain_4_bands_BAF[2,1] - as.numeric(Band3)
  zeroAB_ABB <- MosDet_Gain_4_bands_BAF[3,1] - as.numeric(Band2)
  zeroBB_BBB <- MosDet_Gain_4_bands_BAF[4,1] - as.numeric(Band1)
  zeroTotal <- abs(zeroBB_BBB) + abs(zeroAB_ABA) + abs(zeroAB_ABB) + abs(zeroAA_AAA)

  fiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,2] - as.numeric(Band4)
  fiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,2] - as.numeric(Band3)
  fiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,2] - as.numeric(Band2)
  fiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,2] - as.numeric(Band1)
  fiveTotal <- abs(fiveBB_BBB) + abs(fiveAB_ABA) + abs(fiveAB_ABB) + abs(fiveAA_AAA)

  tenAA_AAA <- MosDet_Gain_4_bands_BAF[1,3] - as.numeric(Band4)
  tenAB_ABA <- MosDet_Gain_4_bands_BAF[2,3] - as.numeric(Band3)
  tenAB_ABB <- MosDet_Gain_4_bands_BAF[3,3] - as.numeric(Band2)
  tenBB_BBB <- MosDet_Gain_4_bands_BAF[4,3] - as.numeric(Band1)
  tenTotal <- abs(tenBB_BBB) + abs(tenAB_ABA) + abs(tenAB_ABB) + abs(tenAA_AAA)

  fifteenAA_AAA <- MosDet_Gain_4_bands_BAF[1,4] - as.numeric(Band4)
  fifteenAB_ABA <- MosDet_Gain_4_bands_BAF[2,4] - as.numeric(Band3)
  fifteenAB_ABB <- MosDet_Gain_4_bands_BAF[3,4] - as.numeric(Band2)
  fifteenBB_BBB <- MosDet_Gain_4_bands_BAF[4,4] - as.numeric(Band1)
  fifteenTotal <- abs(fifteenBB_BBB) + abs(fifteenAB_ABA) + abs(fifteenAB_ABB) + abs(fifteenAA_AAA)

  twentyAA_AAA <- MosDet_Gain_4_bands_BAF[1,5] - as.numeric(Band4)
  twentyAB_ABA <- MosDet_Gain_4_bands_BAF[2,5] - as.numeric(Band3)
  twentyAB_ABB <- MosDet_Gain_4_bands_BAF[3,5] - as.numeric(Band2)
  twentyBB_BBB <- MosDet_Gain_4_bands_BAF[4,5] - as.numeric(Band1)
  twentyTotal <- abs(twentyBB_BBB) + abs(twentyAB_ABA) + abs(twentyAB_ABB) + abs(twentyAA_AAA)

  twentyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,6] - as.numeric(Band4)
  twentyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,6] - as.numeric(Band3)
  twentyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,6] - as.numeric(Band2)
  twentyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,6] - as.numeric(Band1)
  twentyfiveTotal <- abs(twentyfiveBB_BBB) + abs(twentyfiveAB_ABA) + abs(twentyfiveAB_ABB) + abs(twentyfiveAA_AAA)

  thirtyAA_AAA <- MosDet_Gain_4_bands_BAF[1,7] - as.numeric(Band4)
  thirtyAB_ABA <- MosDet_Gain_4_bands_BAF[2,7] - as.numeric(Band3)
  thirtyAB_ABB <- MosDet_Gain_4_bands_BAF[3,7] - as.numeric(Band2)
  thirtyBB_BBB <- MosDet_Gain_4_bands_BAF[4,7] - as.numeric(Band1)
  thirtyTotal <- abs(thirtyBB_BBB) + abs(thirtyAB_ABA) + abs(thirtyAB_ABB) + abs(thirtyAA_AAA)

  thirtyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,8] - as.numeric(Band4)
  thirtyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,8] - as.numeric(Band3)
  thirtyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,8] - as.numeric(Band2)
  thirtyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,8] - as.numeric(Band1)
  thirtyfiveTotal <- abs(thirtyfiveBB_BBB) + abs(thirtyfiveAB_ABA) + abs(thirtyfiveAB_ABB) + abs(thirtyfiveAA_AAA)

  fortyAA_AAA <- MosDet_Gain_4_bands_BAF[1,9] - as.numeric(Band4)
  fortyAB_ABA <- MosDet_Gain_4_bands_BAF[2,9] - as.numeric(Band3)
  fortyAB_ABB <- MosDet_Gain_4_bands_BAF[3,9] - as.numeric(Band2)
  fortyBB_BBB <- MosDet_Gain_4_bands_BAF[4,9] - as.numeric(Band1)
  fortyTotal <- abs(fortyBB_BBB) + abs(fortyAB_ABA) + abs(fortyAB_ABB) + abs(fortyAA_AAA)

  fortyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,10] - as.numeric(Band4)
  fortyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,10] - as.numeric(Band3)
  fortyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,10] - as.numeric(Band2)
  fortyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,10] - as.numeric(Band1)
  fortyfiveTotal <- abs(fortyfiveBB_BBB) + abs(fortyfiveAB_ABA) + abs(fortyfiveAB_ABB) + abs(fortyfiveAA_AAA)

  fiftyAA_AAA <- MosDet_Gain_4_bands_BAF[1,11] - as.numeric(Band4)
  fiftyAB_ABA <- MosDet_Gain_4_bands_BAF[2,11] - as.numeric(Band3)
  fiftyAB_ABB <- MosDet_Gain_4_bands_BAF[3,11] - as.numeric(Band2)
  fiftyBB_BBB <- MosDet_Gain_4_bands_BAF[4,11] - as.numeric(Band1)
  fiftyTotal <- abs(fiftyBB_BBB) + abs(fiftyAB_ABA) + abs(fiftyAB_ABB) + abs(fiftyAA_AAA)

  fiftyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,12] - as.numeric(Band4)
  fiftyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,12] - as.numeric(Band3)
  fiftyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,12] - as.numeric(Band2)
  fiftyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,12] - as.numeric(Band1)
  fiftyfiveTotal <- abs(fiftyfiveBB_BBB) + abs(fiftyfiveAB_ABA) + abs(fiftyfiveAB_ABB) + abs(fiftyfiveAA_AAA)

  sixtyAA_AAA <- MosDet_Gain_4_bands_BAF[1,13] - as.numeric(Band4)
  sixtyAB_ABA <- MosDet_Gain_4_bands_BAF[2,13] - as.numeric(Band3)
  sixtyAB_ABB <- MosDet_Gain_4_bands_BAF[3,13] - as.numeric(Band2)
  sixtyBB_BBB <- MosDet_Gain_4_bands_BAF[4,13] - as.numeric(Band1)
  sixtyTotal <- abs(sixtyBB_BBB) + abs(sixtyAB_ABA) + abs(sixtyAB_ABB) + abs(sixtyAA_AAA)

  sixtyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,14] - as.numeric(Band4)
  sixtyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,14] - as.numeric(Band3)
  sixtyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,14] - as.numeric(Band2)
  sixtyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,14] - as.numeric(Band1)
  sixtyfiveTotal <- abs(sixtyfiveBB_BBB) + abs(sixtyfiveAB_ABA) + abs(sixtyfiveAB_ABB) + abs(sixtyfiveAA_AAA)

  seventyAA_AAA <- MosDet_Gain_4_bands_BAF[1,15] - as.numeric(Band4)
  seventyAB_ABA <- MosDet_Gain_4_bands_BAF[2,15] - as.numeric(Band3)
  seventyAB_ABB <- MosDet_Gain_4_bands_BAF[3,15] - as.numeric(Band2)
  seventyBB_BBB <- MosDet_Gain_4_bands_BAF[4,15] - as.numeric(Band1)
  seventyTotal <- abs(seventyBB_BBB) + abs(seventyAB_ABA) + abs(seventyAB_ABB) + abs(seventyAA_AAA)

  seventyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,16] - as.numeric(Band4)
  seventyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,16] - as.numeric(Band3)
  seventyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,16] - as.numeric(Band2)
  seventyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,16] - as.numeric(Band1)
  seventyfiveTotal <- abs(seventyfiveBB_BBB) + abs(seventyfiveAB_ABA) + abs(seventyfiveAB_ABB) + abs(seventyfiveAA_AAA)

  eightyAA_AAA <- MosDet_Gain_4_bands_BAF[1,17] - as.numeric(Band4)
  eightyAB_ABA <- MosDet_Gain_4_bands_BAF[2,17] - as.numeric(Band3)
  eightyAB_ABB <- MosDet_Gain_4_bands_BAF[3,17] - as.numeric(Band2)
  eightyBB_BBB <- MosDet_Gain_4_bands_BAF[4,17] - as.numeric(Band1)
  eightyTotal <- abs(eightyBB_BBB) + abs(eightyAB_ABA) + abs(eightyAB_ABB) + abs(eightyAA_AAA)

  eightyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,18] - as.numeric(Band4)
  eightyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,18] - as.numeric(Band3)
  eightyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,18] - as.numeric(Band2)
  eightyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,18] - as.numeric(Band1)
  eightyfiveTotal <- abs(eightyfiveBB_BBB) + abs(eightyfiveAB_ABA) + abs(eightyfiveAB_ABB) + abs(eightyfiveAA_AAA)

  ninetyAA_AAA <- MosDet_Gain_4_bands_BAF[1,19] - as.numeric(Band4)
  ninetyAB_ABA <- MosDet_Gain_4_bands_BAF[2,19] - as.numeric(Band3)
  ninetyAB_ABB <- MosDet_Gain_4_bands_BAF[3,19] - as.numeric(Band2)
  ninetyBB_BBB <- MosDet_Gain_4_bands_BAF[4,19] - as.numeric(Band1)
  ninetyTotal <- abs(ninetyBB_BBB) + abs(ninetyAB_ABA) + abs(ninetyAB_ABB) + abs(ninetyAA_AAA)

  ninetyfiveAA_AAA <- MosDet_Gain_4_bands_BAF[1,20] - as.numeric(Band4)
  ninetyfiveAB_ABA <- MosDet_Gain_4_bands_BAF[2,20] - as.numeric(Band3)
  ninetyfiveAB_ABB <- MosDet_Gain_4_bands_BAF[3,20] - as.numeric(Band2)
  ninetyfiveBB_BBB <- MosDet_Gain_4_bands_BAF[4,20] - as.numeric(Band1)
  ninetyfiveTotal <- abs(ninetyfiveBB_BBB) + abs(ninetyfiveAB_ABA) + abs(ninetyfiveAB_ABB) + abs(ninetyfiveAA_AAA)

  onehundredAA_AAA <- MosDet_Gain_4_bands_BAF[1,21] - as.numeric(Band4)
  onehundredAB_ABA <- MosDet_Gain_4_bands_BAF[2,21] - as.numeric(Band3)
  onehundredAB_ABB <- MosDet_Gain_4_bands_BAF[3,21] - as.numeric(Band2)
  onehundredBB_BBB <- MosDet_Gain_4_bands_BAF[4,21] - as.numeric(Band1)
  onehundredTotal <- abs(onehundredBB_BBB) + abs(onehundredAB_ABA) + abs(onehundredAB_ABB) + abs(onehundredAA_AAA)

  ######################

  CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))

  Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
  Percentage
  #End copy gain aberration

} else if (AberrationType == "Loss"){
  
  ##Load mosaicism percentage reference files from Conlin et al. 2010.
  
  MosDet_Loss_BAF <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Github/4_MosaicismDetection/Refs/MosDet_Loss.csv", row.names=1)

  zeroAA_A <- MosDet_Loss_BAF[1,1] - as.numeric(Band4)
  zeroAB_A <- MosDet_Loss_BAF[2,1] - as.numeric(Band3)
  zeroAB_B <- MosDet_Loss_BAF[3,1] - as.numeric(Band2)
  zeroBB_B <- MosDet_Loss_BAF[4,1] - as.numeric(Band1)
  zeroTotal <- abs(zeroBB_B) + abs(zeroAB_A) + abs(zeroAB_B) + abs(zeroAA_A)
  
  fiveAA_A <- MosDet_Loss_BAF[1,2] - as.numeric(Band4)
  fiveAB_A <- MosDet_Loss_BAF[2,2] - as.numeric(Band3)
  fiveAB_B <- MosDet_Loss_BAF[3,2] - as.numeric(Band2)
  fiveBB_B <- MosDet_Loss_BAF[4,2] - as.numeric(Band1)
  fiveTotal <- abs(fiveBB_B) + abs(fiveAB_A) + abs(fiveAB_B) + abs(fiveAA_A)
  
  tenAA_A <- MosDet_Loss_BAF[1,3] - as.numeric(Band4)
  tenAB_A <- MosDet_Loss_BAF[2,3] - as.numeric(Band3)
  tenAB_B <- MosDet_Loss_BAF[3,3] - as.numeric(Band2)
  tenBB_B <- MosDet_Loss_BAF[4,3] - as.numeric(Band1)
  tenTotal <- abs(tenBB_B) + abs(tenAB_A) + abs(tenAB_B) + abs(tenAA_A)
  
  fifteenAA_A <- MosDet_Loss_BAF[1,4] - as.numeric(Band4)
  fifteenAB_A <- MosDet_Loss_BAF[2,4] - as.numeric(Band3)
  fifteenAB_B <- MosDet_Loss_BAF[3,4] - as.numeric(Band2)
  fifteenBB_B <- MosDet_Loss_BAF[4,4] - as.numeric(Band1)
  fifteenTotal <- abs(fifteenBB_B) + abs(fifteenAB_A) + abs(fifteenAB_B) + abs(fifteenAA_A)
  
  twentyAA_A <- MosDet_Loss_BAF[1,5] - as.numeric(Band4)
  twentyAB_A <- MosDet_Loss_BAF[2,5] - as.numeric(Band3)
  twentyAB_B <- MosDet_Loss_BAF[3,5] - as.numeric(Band2)
  twentyBB_B <- MosDet_Loss_BAF[4,5] - as.numeric(Band1)
  twentyTotal <- abs(twentyBB_B) + abs(twentyAB_A) + abs(twentyAB_B) + abs(twentyAA_A)
  
  twentyfiveAA_A <- MosDet_Loss_BAF[1,6] - as.numeric(Band4)
  twentyfiveAB_A <- MosDet_Loss_BAF[2,6] - as.numeric(Band3)
  twentyfiveAB_B <- MosDet_Loss_BAF[3,6] - as.numeric(Band2)
  twentyfiveBB_B <- MosDet_Loss_BAF[4,6] - as.numeric(Band1)
  twentyfiveTotal <- abs(twentyfiveBB_B) + abs(twentyfiveAB_A) + abs(twentyfiveAB_B) + abs(twentyfiveAA_A)
  
  thirtyAA_A <- MosDet_Loss_BAF[1,7] - as.numeric(Band4)
  thirtyAB_A <- MosDet_Loss_BAF[2,7] - as.numeric(Band3)
  thirtyAB_B <- MosDet_Loss_BAF[3,7] - as.numeric(Band2)
  thirtyBB_B <- MosDet_Loss_BAF[4,7] - as.numeric(Band1)
  thirtyTotal <- abs(thirtyBB_B) + abs(thirtyAB_A) + abs(thirtyAB_B) + abs(thirtyAA_A)
  
  thirtyfiveAA_A <- MosDet_Loss_BAF[1,8] - as.numeric(Band4)
  thirtyfiveAB_A <- MosDet_Loss_BAF[2,8] - as.numeric(Band3)
  thirtyfiveAB_B <- MosDet_Loss_BAF[3,8] - as.numeric(Band2)
  thirtyfiveBB_B <- MosDet_Loss_BAF[4,8] - as.numeric(Band1)
  thirtyfiveTotal <- abs(thirtyfiveBB_B) + abs(thirtyfiveAB_A) + abs(thirtyfiveAB_B) + abs(thirtyfiveAA_A)
  
  fortyAA_A <- MosDet_Loss_BAF[1,9] - as.numeric(Band4)
  fortyAB_A <- MosDet_Loss_BAF[2,9] - as.numeric(Band3)
  fortyAB_B <- MosDet_Loss_BAF[3,9] - as.numeric(Band2)
  fortyBB_B <- MosDet_Loss_BAF[4,9] - as.numeric(Band1)
  fortyTotal <- abs(fortyBB_B) + abs(fortyAB_A) + abs(fortyAB_B) + abs(fortyAA_A)
  
  fortyfiveAA_A <- MosDet_Loss_BAF[1,10] - as.numeric(Band4)
  fortyfiveAB_A <- MosDet_Loss_BAF[2,10] - as.numeric(Band3)
  fortyfiveAB_B <- MosDet_Loss_BAF[3,10] - as.numeric(Band2)
  fortyfiveBB_B <- MosDet_Loss_BAF[4,10] - as.numeric(Band1)
  fortyfiveTotal <- abs(fortyfiveBB_B) + abs(fortyfiveAB_A) + abs(fortyfiveAB_B) + abs(fortyfiveAA_A)
  
  fiftyAA_A <- MosDet_Loss_BAF[1,11] - as.numeric(Band4)
  fiftyAB_A <- MosDet_Loss_BAF[2,11] - as.numeric(Band3)
  fiftyAB_B <- MosDet_Loss_BAF[3,11] - as.numeric(Band2)
  fiftyBB_B <- MosDet_Loss_BAF[4,11] - as.numeric(Band1)
  fiftyTotal <- abs(fiftyBB_B) + abs(fiftyAB_A) + abs(fiftyAB_B) + abs(fiftyAA_A)
  
  fiftyfiveAA_A <- MosDet_Loss_BAF[1,12] - as.numeric(Band4)
  fiftyfiveAB_A <- MosDet_Loss_BAF[2,12] - as.numeric(Band3)
  fiftyfiveAB_B <- MosDet_Loss_BAF[3,12] - as.numeric(Band2)
  fiftyfiveBB_B <- MosDet_Loss_BAF[4,12] - as.numeric(Band1)
  fiftyfiveTotal <- abs(fiftyfiveBB_B) + abs(fiftyfiveAB_A) + abs(fiftyfiveAB_B) + abs(fiftyfiveAA_A)
  
  sixtyAA_A <- MosDet_Loss_BAF[1,13] - as.numeric(Band4)
  sixtyAB_A <- MosDet_Loss_BAF[2,13] - as.numeric(Band3)
  sixtyAB_B <- MosDet_Loss_BAF[3,13] - as.numeric(Band2)
  sixtyBB_B <- MosDet_Loss_BAF[4,13] - as.numeric(Band1)
  sixtyTotal <- abs(sixtyBB_B) + abs(sixtyAB_A) + abs(sixtyAB_B) + abs(sixtyAA_A)
  
  sixtyfiveAA_A <- MosDet_Loss_BAF[1,14] - as.numeric(Band4)
  sixtyfiveAB_A <- MosDet_Loss_BAF[2,14] - as.numeric(Band3)
  sixtyfiveAB_B <- MosDet_Loss_BAF[3,14] - as.numeric(Band2)
  sixtyfiveBB_B <- MosDet_Loss_BAF[4,14] - as.numeric(Band1)
  sixtyfiveTotal <- abs(sixtyfiveBB_B) + abs(sixtyfiveAB_A) + abs(sixtyfiveAB_B) + abs(sixtyfiveAA_A)
  
  seventyAA_A <- MosDet_Loss_BAF[1,15] - as.numeric(Band4)
  seventyAB_A <- MosDet_Loss_BAF[2,15] - as.numeric(Band3)
  seventyAB_B <- MosDet_Loss_BAF[3,15] - as.numeric(Band2)
  seventyBB_B <- MosDet_Loss_BAF[4,15] - as.numeric(Band1)
  seventyTotal <- abs(seventyBB_B) + abs(seventyAB_A) + abs(seventyAB_B) + abs(seventyAA_A)
  
  seventyfiveAA_A <- MosDet_Loss_BAF[1,16] - as.numeric(Band4)
  seventyfiveAB_A <- MosDet_Loss_BAF[2,16] - as.numeric(Band3)
  seventyfiveAB_B <- MosDet_Loss_BAF[3,16] - as.numeric(Band2)
  seventyfiveBB_B <- MosDet_Loss_BAF[4,16] - as.numeric(Band1)
  seventyfiveTotal <- abs(seventyfiveBB_B) + abs(seventyfiveAB_A) + abs(seventyfiveAB_B) + abs(seventyfiveAA_A)
  
  eightyAA_A <- MosDet_Loss_BAF[1,17] - as.numeric(Band4)
  eightyAB_A <- MosDet_Loss_BAF[2,17] - as.numeric(Band3)
  eightyAB_B <- MosDet_Loss_BAF[3,17] - as.numeric(Band2)
  eightyBB_B <- MosDet_Loss_BAF[4,17] - as.numeric(Band1)
  eightyTotal <- abs(eightyBB_B) + abs(eightyAB_A) + abs(eightyAB_B) + abs(eightyAA_A)
  
  eightyfiveAA_A <- MosDet_Loss_BAF[1,18] - as.numeric(Band4)
  eightyfiveAB_A <- MosDet_Loss_BAF[2,18] - as.numeric(Band3)
  eightyfiveAB_B <- MosDet_Loss_BAF[3,18] - as.numeric(Band2)
  eightyfiveBB_B <- MosDet_Loss_BAF[4,18] - as.numeric(Band1)
  eightyfiveTotal <- abs(eightyfiveBB_B) + abs(eightyfiveAB_A) + abs(eightyfiveAB_B) + abs(eightyfiveAA_A)
  
  ninetyAA_A <- MosDet_Loss_BAF[1,19] - as.numeric(Band4)
  ninetyAB_A <- MosDet_Loss_BAF[2,19] - as.numeric(Band3)
  ninetyAB_B <- MosDet_Loss_BAF[3,19] - as.numeric(Band2)
  ninetyBB_B <- MosDet_Loss_BAF[4,19] - as.numeric(Band1)
  ninetyTotal <- abs(ninetyBB_B) + abs(ninetyAB_A) + abs(ninetyAB_B) + abs(ninetyAA_A)
  
  ninetyfiveAA_A <- MosDet_Loss_BAF[1,20] - as.numeric(Band4)
  ninetyfiveAB_A <- MosDet_Loss_BAF[2,20] - as.numeric(Band3)
  ninetyfiveAB_B <- MosDet_Loss_BAF[3,20] - as.numeric(Band2)
  ninetyfiveBB_B <- MosDet_Loss_BAF[4,20] - as.numeric(Band1)
  ninetyfiveTotal <- abs(ninetyfiveBB_B) + abs(ninetyfiveAB_A) + abs(ninetyfiveAB_B) + abs(ninetyfiveAA_A)
  
  onehundredAA_A <- MosDet_Loss_BAF[1,21] - as.numeric(Band4)
  onehundredAB_A <- MosDet_Loss_BAF[2,21] - as.numeric(Band3)
  onehundredAB_B <- MosDet_Loss_BAF[3,21] - as.numeric(Band2)
  onehundredBB_B <- MosDet_Loss_BAF[4,21] - as.numeric(Band1)
  onehundredTotal <- abs(onehundredBB_B) + abs(onehundredAB_A) + abs(onehundredAB_B) + abs(onehundredAA_A)

  
  ###
  
  CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
  
  Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
  Percentage
  ##End copy loss aberration
  
} else if (AberrationType == "Neutral"){

  ##Load mosaicism percentage reference files from Conlin et al. 2010.
  
  MosDet_Neutral <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive/Submission/Github/4_MosaicismDetection/Refs/MosDet_Neutral.csv", row.names=1)
  
  zeroAA_AA <- MosDet_Neutral[1,1] - as.numeric(Band4)
  zeroAB_AA <- MosDet_Neutral[2,1] - as.numeric(Band3)
  zeroAB_BB <- MosDet_Neutral[3,1] - as.numeric(Band2)
  zeroBB_BB <- MosDet_Neutral[4,1] - as.numeric(Band1)
  zeroTotal <- abs(zeroBB_BB) + abs(zeroAB_AA) + abs(zeroAB_BB) + abs(zeroAA_AA)
  
  fiveAA_AA <- MosDet_Neutral[1,2] - as.numeric(Band4)
  fiveAB_AA <- MosDet_Neutral[2,2] - as.numeric(Band3)
  fiveAB_BB <- MosDet_Neutral[3,2] - as.numeric(Band2)
  fiveBB_BB <- MosDet_Neutral[4,2] - as.numeric(Band1)
  fiveTotal <- abs(fiveBB_BB) + abs(fiveAB_AA) + abs(fiveAB_BB) + abs(fiveAA_AA)
  
  tenAA_AA <- MosDet_Neutral[1,3] - as.numeric(Band4)
  tenAB_AA <- MosDet_Neutral[2,3] - as.numeric(Band3)
  tenAB_BB <- MosDet_Neutral[3,3] - as.numeric(Band2)
  tenBB_BB <- MosDet_Neutral[4,3] - as.numeric(Band1)
  tenTotal <- abs(tenBB_BB) + abs(tenAB_AA) + abs(tenAB_BB) + abs(tenAA_AA)
  
  fifteenAA_AA <- MosDet_Neutral[1,4] - as.numeric(Band4)
  fifteenAB_AA <- MosDet_Neutral[2,4] - as.numeric(Band3)
  fifteenAB_BB <- MosDet_Neutral[3,4] - as.numeric(Band2)
  fifteenBB_BB <- MosDet_Neutral[4,4] - as.numeric(Band1)
  fifteenTotal <- abs(fifteenBB_BB) + abs(fifteenAB_AA) + abs(fifteenAB_BB) + abs(fifteenAA_AA)
  
  twentyAA_AA <- MosDet_Neutral[1,5] - as.numeric(Band4)
  twentyAB_AA <- MosDet_Neutral[2,5] - as.numeric(Band3)
  twentyAB_BB <- MosDet_Neutral[3,5] - as.numeric(Band2)
  twentyBB_BB <- MosDet_Neutral[4,5] - as.numeric(Band1)
  twentyTotal <- abs(twentyBB_BB) + abs(twentyAB_AA) + abs(twentyAB_BB) + abs(twentyAA_AA)
  
  twentyfiveAA_AA <- MosDet_Neutral[1,6] - as.numeric(Band4)
  twentyfiveAB_AA <- MosDet_Neutral[2,6] - as.numeric(Band3)
  twentyfiveAB_BB <- MosDet_Neutral[3,6] - as.numeric(Band2)
  twentyfiveBB_BB <- MosDet_Neutral[4,6] - as.numeric(Band1)
  twentyfiveTotal <- abs(twentyfiveBB_BB) + abs(twentyfiveAB_AA) + abs(twentyfiveAB_BB) + abs(twentyfiveAA_AA)
  
  
  thirtyAA_AA <- MosDet_Neutral[1,7] - as.numeric(Band4)
  thirtyAB_AA <- MosDet_Neutral[2,7] - as.numeric(Band3)
  thirtyAB_BB <- MosDet_Neutral[3,7] - as.numeric(Band2)
  thirtyBB_BB <- MosDet_Neutral[4,7] - as.numeric(Band1)
  thirtyTotal <- abs(thirtyBB_BB) + abs(thirtyAB_AA) + abs(thirtyAB_BB) + abs(thirtyAA_AA)
  
  thirtyfiveAA_AA <- MosDet_Neutral[1,8] - as.numeric(Band4)
  thirtyfiveAB_AA <- MosDet_Neutral[2,8] - as.numeric(Band3)
  thirtyfiveAB_BB <- MosDet_Neutral[3,8] - as.numeric(Band2)
  thirtyfiveBB_BB <- MosDet_Neutral[4,8] - as.numeric(Band1)
  thirtyfiveTotal <- abs(thirtyfiveBB_BB) + abs(thirtyfiveAB_AA) + abs(thirtyfiveAB_BB) + abs(thirtyfiveAA_AA)
  
  fortyAA_AA <- MosDet_Neutral[1,9] - as.numeric(Band4)
  fortyAB_AA <- MosDet_Neutral[2,9] - as.numeric(Band3)
  fortyAB_BB <- MosDet_Neutral[3,9] - as.numeric(Band2)
  fortyBB_BB <- MosDet_Neutral[4,9] - as.numeric(Band1)
  fortyTotal <- abs(fortyBB_BB) + abs(fortyAB_AA) + abs(fortyAB_BB) + abs(fortyAA_AA)
  
  fortyfiveAA_AA <- MosDet_Neutral[1,10] - as.numeric(Band4)
  fortyfiveAB_AA <- MosDet_Neutral[2,10] - as.numeric(Band3)
  fortyfiveAB_BB <- MosDet_Neutral[3,10] - as.numeric(Band2)
  fortyfiveBB_BB <- MosDet_Neutral[4,10] - as.numeric(Band1)
  fortyfiveTotal <- abs(fortyfiveBB_BB) + abs(fortyfiveAB_AA) + abs(fortyfiveAB_BB) + abs(fortyfiveAA_AA)
  
  fiftyAA_AA <- MosDet_Neutral[1,11] - as.numeric(Band4)
  fiftyAB_AA <- MosDet_Neutral[2,11] - as.numeric(Band3)
  fiftyAB_BB <- MosDet_Neutral[3,11] - as.numeric(Band2)
  fiftyBB_BB <- MosDet_Neutral[4,11] - as.numeric(Band1)
  fiftyTotal <- abs(fiftyBB_BB) + abs(fiftyAB_AA) + abs(fiftyAB_BB) + abs(fiftyAA_AA)
  
  fiftyfiveAA_AA <- MosDet_Neutral[1,12] - as.numeric(Band4)
  fiftyfiveAB_AA <- MosDet_Neutral[2,12] - as.numeric(Band3)
  fiftyfiveAB_BB <- MosDet_Neutral[3,12] - as.numeric(Band2)
  fiftyfiveBB_BB <- MosDet_Neutral[4,12] - as.numeric(Band1)
  fiftyfiveTotal <- abs(fiftyfiveBB_BB) + abs(fiftyfiveAB_AA) + abs(fiftyfiveAB_BB) + abs(fiftyfiveAA_AA)
  
  sixtyAA_AA <- MosDet_Neutral[1,13] - as.numeric(Band4)
  sixtyAB_AA <- MosDet_Neutral[2,13] - as.numeric(Band3)
  sixtyAB_BB <- MosDet_Neutral[3,13] - as.numeric(Band2)
  sixtyBB_BB <- MosDet_Neutral[4,13] - as.numeric(Band1)
  sixtyTotal <- abs(sixtyBB_BB) + abs(sixtyAB_AA) + abs(sixtyAB_BB) + abs(sixtyAA_AA)
  
  sixtyfiveAA_AA <- MosDet_Neutral[1,14] - as.numeric(Band4)
  sixtyfiveAB_AA <- MosDet_Neutral[2,14] - as.numeric(Band3)
  sixtyfiveAB_BB <- MosDet_Neutral[3,14] - as.numeric(Band2)
  sixtyfiveBB_BB <- MosDet_Neutral[4,14] - as.numeric(Band1)
  sixtyfiveTotal <- abs(sixtyfiveBB_BB) + abs(sixtyfiveAB_AA) + abs(sixtyfiveAB_BB) + abs(sixtyfiveAA_AA)
  
  seventyAA_AA <- MosDet_Neutral[1,15] - as.numeric(Band4)
  seventyAB_AA <- MosDet_Neutral[2,15] - as.numeric(Band3)
  seventyAB_BB <- MosDet_Neutral[3,15] - as.numeric(Band2)
  seventyBB_BB <- MosDet_Neutral[4,15] - as.numeric(Band1)
  seventyTotal <- abs(seventyBB_BB) + abs(seventyAB_AA) + abs(seventyAB_BB) + abs(seventyAA_AA)
  
  seventyfiveAA_AA <- MosDet_Neutral[1,16] - as.numeric(Band4)
  seventyfiveAB_AA <- MosDet_Neutral[2,16] - as.numeric(Band3)
  seventyfiveAB_BB <- MosDet_Neutral[3,16] - as.numeric(Band2)
  seventyfiveBB_BB <- MosDet_Neutral[4,16] - as.numeric(Band1)
  seventyfiveTotal <- abs(seventyfiveBB_BB) + abs(seventyfiveAB_AA) + abs(seventyfiveAB_BB) + abs(seventyfiveAA_AA)
  
  eightyAA_AA <- MosDet_Neutral[1,17] - as.numeric(Band4)
  eightyAB_AA <- MosDet_Neutral[2,17] - as.numeric(Band3)
  eightyAB_BB <- MosDet_Neutral[3,17] - as.numeric(Band2)
  eightyBB_BB <- MosDet_Neutral[4,17] - as.numeric(Band1)
  eightyTotal <- abs(eightyBB_BB) + abs(eightyAB_AA) + abs(eightyAB_BB) + abs(eightyAA_AA)
  
  eightyfiveAA_AA <- MosDet_Neutral[1,18] - as.numeric(Band4)
  eightyfiveAB_AA <- MosDet_Neutral[2,18] - as.numeric(Band3)
  eightyfiveAB_BB <- MosDet_Neutral[3,18] - as.numeric(Band2)
  eightyfiveBB_BB <- MosDet_Neutral[4,18] - as.numeric(Band1)
  eightyfiveTotal <- abs(eightyfiveBB_BB) + abs(eightyfiveAB_AA) + abs(eightyfiveAB_BB) + abs(eightyfiveAA_AA)
  
  ninetyAA_AA <- MosDet_Neutral[1,19] - as.numeric(Band4)
  ninetyAB_AA <- MosDet_Neutral[2,19] - as.numeric(Band3)
  ninetyAB_BB <- MosDet_Neutral[3,19] - as.numeric(Band2)
  ninetyBB_BB <- MosDet_Neutral[4,19] - as.numeric(Band1)
  ninetyTotal <- abs(ninetyBB_BB) + abs(ninetyAB_AA) + abs(ninetyAB_BB) + abs(ninetyAA_AA)
  
  ninetyfiveAA_AA <- MosDet_Neutral[1,20] - as.numeric(Band4)
  ninetyfiveAB_AA <- MosDet_Neutral[2,20] - as.numeric(Band3)
  ninetyfiveAB_BB <- MosDet_Neutral[3,20] - as.numeric(Band2)
  ninetyfiveBB_BB <- MosDet_Neutral[4,20] - as.numeric(Band1)
  ninetyfiveTotal <- abs(ninetyfiveBB_BB) + abs(ninetyfiveAB_AA) + abs(ninetyfiveAB_BB) + abs(ninetyfiveAA_AA)
  
  onehundredAA_AA <- MosDet_Neutral[1,21] - as.numeric(Band4)
  onehundredAB_AA <- MosDet_Neutral[2,21] - as.numeric(Band3)
  onehundredAB_BB <- MosDet_Neutral[3,21] - as.numeric(Band2)
  onehundredBB_BB <- MosDet_Neutral[4,21] - as.numeric(Band1)
  onehundredTotal <- abs(onehundredBB_BB) + abs(onehundredAB_AA) + abs(onehundredAB_BB) + abs(onehundredAA_AA)
  
  ######################
  
  CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
  
  Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
  Percentage

}  ###End copy neutral aberration