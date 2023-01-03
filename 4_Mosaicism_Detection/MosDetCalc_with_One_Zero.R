###Trisomy MII_Mit (4 BAF bands)###

rm(list=ls(all=T))

##Load mosaicism percentage reference files from Conlin et al. 2010.
MosDet_TRISOMY_MII_Mit <- read.table("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/MosaicismDetection/MosDet_TRISOMY_MII_Mit.csv",sep="\t",header=T,stringsAsFactors=F,row.names = 1)

options(digits = 16)

###BAFTest_EM_PL2682 -- 4 bands
Band1   <- "0.985042076093294"
Band2   <- "0.654957995274585"
Band3   <- "0.317782212368024"
Band4   <- "0.0136834723909699"

###Trisomy MII_Mit (4 BAF bands)###

  zeroAB_ABB <- MosDet_TRISOMY_MII_Mit[3,1] - as.numeric(P2value)
  zeroAB_ABA <- MosDet_TRISOMY_MII_Mit[2,1] - as.numeric(M2value)
  zeroBB_BBB <- MosDet_TRISOMY_MII_Mit[4,1] - as.numeric(P1value)
  zeroAA_AAA <- MosDet_TRISOMY_MII_Mit[1,1] - as.numeric(M1value)
  zeroTotal <- abs(zeroBB_BBB) + abs(zeroAB_ABA) + abs(zeroAB_ABB) + abs(zeroAA_AAA)
  
  fiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,2] - as.numeric(P2value)
  fiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,2] - as.numeric(M2value)
  fiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,2] - as.numeric(P1value)
  fiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,2] - as.numeric(M1value)
  fiveTotal <- abs(fiveBB_BBB) + abs(fiveAB_ABA) + abs(fiveAB_ABB) + abs(fiveAA_AAA)
  
  tenAB_ABB <- MosDet_TRISOMY_MII_Mit[3,3] - as.numeric(P2value)
  tenAB_ABA <- MosDet_TRISOMY_MII_Mit[2,3] - as.numeric(M2value)
  tenBB_BBB <- MosDet_TRISOMY_MII_Mit[4,3] - as.numeric(P1value)
  tenAA_AAA <- MosDet_TRISOMY_MII_Mit[1,3] - as.numeric(M1value)
  tenTotal <- abs(tenBB_BBB) + abs(tenAB_ABA) + abs(tenAB_ABB) + abs(tenAA_AAA)
  
  fifteenAB_ABB <- MosDet_TRISOMY_MII_Mit[3,4] - as.numeric(P2value)
  fifteenAB_ABA <- MosDet_TRISOMY_MII_Mit[2,4] - as.numeric(M2value)
  fifteenBB_BBB <- MosDet_TRISOMY_MII_Mit[4,4] - as.numeric(P1value)
  fifteenAA_AAA <- MosDet_TRISOMY_MII_Mit[1,4] - as.numeric(M1value)
  fifteenTotal <- abs(fifteenBB_BBB) + abs(fifteenAB_ABA) + abs(fifteenAB_ABB) + abs(fifteenAA_AAA)
  
  twentyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,5] - as.numeric(P2value)
  twentyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,5] - as.numeric(M2value)
  twentyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,5] - as.numeric(P1value)
  twentyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,5] - as.numeric(M1value)
  twentyTotal <- abs(twentyBB_BBB) + abs(twentyAB_ABA) + abs(twentyAB_ABB) + abs(twentyAA_AAA)
  
  twentyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,6] - as.numeric(P2value)
  twentyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,6] - as.numeric(M2value)
  twentyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,6] - as.numeric(P1value)
  twentyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,6] - as.numeric(M1value)
  twentyfiveTotal <- abs(twentyfiveBB_BBB) + abs(twentyfiveAB_ABA) + abs(twentyfiveAB_ABB) + abs(twentyfiveAA_AAA)
  
  thirtyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,7] - as.numeric(P2value)
  thirtyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,7] - as.numeric(M2value)
  thirtyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,7] - as.numeric(P1value)
  thirtyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,7] - as.numeric(M1value)
  thirtyTotal <- abs(thirtyBB_BBB) + abs(thirtyAB_ABA) + abs(thirtyAB_ABB) + abs(thirtyAA_AAA)
  
  thirtyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,8] - as.numeric(P2value)
  thirtyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,8] - as.numeric(M2value)
  thirtyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,8] - as.numeric(P1value)
  thirtyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,8] - as.numeric(M1value)
  thirtyfiveTotal <- abs(thirtyfiveBB_BBB) + abs(thirtyfiveAB_ABA) + abs(thirtyfiveAB_ABB) + abs(thirtyfiveAA_AAA)
  
  fortyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,9] - as.numeric(P2value)
  fortyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,9] - as.numeric(M2value)
  fortyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,9] - as.numeric(P1value)
  fortyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,9] - as.numeric(M1value)
  fortyTotal <- abs(fortyBB_BBB) + abs(fortyAB_ABA) + abs(fortyAB_ABB) + abs(fortyAA_AAA)
  
  fortyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,10] - as.numeric(P2value)
  fortyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,10] - as.numeric(M2value)
  fortyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,10] - as.numeric(P1value)
  fortyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,10] - as.numeric(M1value)
  fortyfiveTotal <- abs(fortyfiveBB_BBB) + abs(fortyfiveAB_ABA) + abs(fortyfiveAB_ABB) + abs(fortyfiveAA_AAA)
  
  fiftyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,11] - as.numeric(P2value)
  fiftyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,11] - as.numeric(M2value)
  fiftyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,11] - as.numeric(P1value)
  fiftyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,11] - as.numeric(M1value)
  fiftyTotal <- abs(fiftyBB_BBB) + abs(fiftyAB_ABA) + abs(fiftyAB_ABB) + abs(fiftyAA_AAA)
  
  fiftyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,12] - as.numeric(P2value)
  fiftyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,12] - as.numeric(M2value)
  fiftyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,12] - as.numeric(P1value)
  fiftyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,12] - as.numeric(M1value)
  fiftyfiveTotal <- abs(fiftyfiveBB_BBB) + abs(fiftyfiveAB_ABA) + abs(fiftyfiveAB_ABB) + abs(fiftyfiveAA_AAA)
  
  sixtyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,13] - as.numeric(P2value)
  sixtyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,13] - as.numeric(M2value)
  sixtyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,13] - as.numeric(P1value)
  sixtyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,13] - as.numeric(M1value)
  sixtyTotal <- abs(sixtyBB_BBB) + abs(sixtyAB_ABA) + abs(sixtyAB_ABB) + abs(sixtyAA_AAA)
  
  sixtyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,14] - as.numeric(P2value)
  sixtyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,14] - as.numeric(M2value)
  sixtyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,14] - as.numeric(P1value)
  sixtyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,14] - as.numeric(M1value)
  sixtyfiveTotal <- abs(sixtyfiveBB_BBB) + abs(sixtyfiveAB_ABA) + abs(sixtyfiveAB_ABB) + abs(sixtyfiveAA_AAA)
  
  seventyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,15] - as.numeric(P2value)
  seventyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,15] - as.numeric(M2value)
  seventyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,15] - as.numeric(P1value)
  seventyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,15] - as.numeric(M1value)
  seventyTotal <- abs(seventyBB_BBB) + abs(seventyAB_ABA) + abs(seventyAB_ABB) + abs(seventyAA_AAA)
  
  seventyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,16] - as.numeric(P2value)
  seventyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,16] - as.numeric(M2value)
  seventyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,16] - as.numeric(P1value)
  seventyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,16] - as.numeric(M1value)
  seventyfiveTotal <- abs(seventyfiveBB_BBB) + abs(seventyfiveAB_ABA) + abs(seventyfiveAB_ABB) + abs(seventyfiveAA_AAA)
  
  eightyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,17] - as.numeric(P2value)
  eightyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,17] - as.numeric(M2value)
  eightyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,17] - as.numeric(P1value)
  eightyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,17] - as.numeric(M1value)
  eightyTotal <- abs(eightyBB_BBB) + abs(eightyAB_ABA) + abs(eightyAB_ABB) + abs(eightyAA_AAA)
  
  eightyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,18] - as.numeric(P2value)
  eightyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,18] - as.numeric(M2value)
  eightyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,18] - as.numeric(P1value)
  eightyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,18] - as.numeric(M1value)
  eightyfiveTotal <- abs(eightyfiveBB_BBB) + abs(eightyfiveAB_ABA) + abs(eightyfiveAB_ABB) + abs(eightyfiveAA_AAA)
  
  ninetyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,19] - as.numeric(P2value)
  ninetyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,19] - as.numeric(M2value)
  ninetyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,19] - as.numeric(P1value)
  ninetyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,19] - as.numeric(M1value)
  ninetyTotal <- abs(ninetyBB_BBB) + abs(ninetyAB_ABA) + abs(ninetyAB_ABB) + abs(ninetyAA_AAA)
  
  ninetyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,20] - as.numeric(P2value)
  ninetyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,20] - as.numeric(M2value)
  ninetyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,20] - as.numeric(P1value)
  ninetyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,20] - as.numeric(M1value)
  ninetyfiveTotal <- abs(ninetyfiveBB_BBB) + abs(ninetyfiveAB_ABA) + abs(ninetyfiveAB_ABB) + abs(ninetyfiveAA_AAA)
  
  onehundredAB_ABB <- MosDet_TRISOMY_MII_Mit[3,21] - as.numeric(P2value)
  onehundredAB_ABA <- MosDet_TRISOMY_MII_Mit[2,21] - as.numeric(M2value)
  onehundredBB_BBB <- MosDet_TRISOMY_MII_Mit[4,21] - as.numeric(P1value)
  onehundredAA_AAA <- MosDet_TRISOMY_MII_Mit[1,21] - as.numeric(M1value)
  onehundredTotal <- abs(onehundredBB_BBB) + abs(onehundredAB_ABA) + abs(onehundredAB_ABB) + abs(onehundredAA_AAA)
  
  ######################3
  
  CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
  
  Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
  Percentage
 