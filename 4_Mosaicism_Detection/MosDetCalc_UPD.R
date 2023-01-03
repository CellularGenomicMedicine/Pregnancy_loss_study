rm(list=ls(all=T))

##Load mosaicism percentage reference files from Conlin et al. 2010.
MosDet_UPD <- read.table("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/MosaicismDetection/MosDet_UPD.csv",sep=",",header=T,stringsAsFactors=F,row.names = 1)

options(digits = 16)

###Sample
P1value <- "0.00441444576968892"
P2value <- "0.535782276179428"
M1value <- "0.00637766369634855"
M2value <- "0.466122348587571"

### More contibution to genome from mother or father?
Parent <- "Mother"

###Compensate for maternal contamination (ONLY use when mat cont is present, otherwise you unnecessarily compensate)
P2value <- as.numeric(P2value) - as.numeric(M1value) 


if (Parent == "Mother"){
  
  zeroAB_BB <- MosDet_UPD[3,1] - as.numeric(P2value)
  zeroAB_AA <- MosDet_UPD[2,1] - as.numeric(M2value)
  zeroTotal <- abs(zeroAB_AA) + abs(zeroAB_BB)

  fiveAB_BB <- MosDet_UPD[3,2] - as.numeric(P2value)
  fiveAB_AA <- MosDet_UPD[2,2] - as.numeric(M2value)
  fiveTotal <- abs(fiveAB_AA) + abs(fiveAB_BB)

  tenAB_BB <- MosDet_UPD[3,3] - as.numeric(P2value)
  tenAB_AA <- MosDet_UPD[2,3] - as.numeric(M2value)
  tenTotal <- abs(tenAB_AA) + abs(tenAB_BB)

  fifteenAB_BB <- MosDet_UPD[3,4] - as.numeric(P2value)
  fifteenAB_AA <- MosDet_UPD[2,4] - as.numeric(M2value)
  fifteenTotal <- abs(fifteenAB_AA) + abs(fifteenAB_BB)

  twentyAB_BB <- MosDet_UPD[3,5] - as.numeric(P2value)
  twentyAB_AA <- MosDet_UPD[2,5] - as.numeric(M2value)
  twentyTotal <- abs(twentyAB_AA) + abs(twentyAB_BB)

  twentyfiveAB_BB <- MosDet_UPD[3,6] - as.numeric(P2value)
  twentyfiveAB_AA <- MosDet_UPD[2,6] - as.numeric(M2value)
  twentyfiveTotal <- abs(twentyfiveAB_AA) + abs(twentyfiveAB_BB)

  thirtyAB_BB <- MosDet_UPD[3,7] - as.numeric(P2value)
  thirtyAB_AA <- MosDet_UPD[2,7] - as.numeric(M2value)
  thirtyTotal <- abs(thirtyAB_AA) + abs(thirtyAB_BB)

  thirtyfiveAB_BB <- MosDet_UPD[3,8] - as.numeric(P2value)
  thirtyfiveAB_AA <- MosDet_UPD[2,8] - as.numeric(M2value)
  thirtyfiveTotal <- abs(thirtyfiveAB_AA) + abs(thirtyfiveAB_BB)

  fortyAB_BB <- MosDet_UPD[3,9] - as.numeric(P2value)
  fortyAB_AA <- MosDet_UPD[2,9] - as.numeric(M2value)
  fortyTotal <- abs(fortyAB_AA) + abs(fortyAB_BB)

  fortyfiveAB_BB <- MosDet_UPD[3,10] - as.numeric(P2value)
  fortyfiveAB_AA <- MosDet_UPD[2,10] - as.numeric(M2value)
  fortyfiveTotal <- abs(fortyfiveAB_AA) + abs(fortyfiveAB_BB)

  fiftyAB_BB <- MosDet_UPD[3,11] - as.numeric(P2value)
  fiftyAB_AA <- MosDet_UPD[2,11] - as.numeric(M2value)
  fiftyTotal <- abs(fiftyAB_AA) + abs(fiftyAB_BB)

  fiftyfiveAB_BB <- MosDet_UPD[3,12] - as.numeric(P2value)
  fiftyfiveAB_AA <- MosDet_UPD[2,12] - as.numeric(M2value)
  fiftyfiveTotal <- abs(fiftyfiveAB_AA) + abs(fiftyfiveAB_BB)

  sixtyAB_BB <- MosDet_UPD[3,13] - as.numeric(P2value)
  sixtyAB_AA <- MosDet_UPD[2,13] - as.numeric(M2value)
  sixtyTotal <- abs(sixtyAB_AA) + abs(sixtyAB_BB)

  sixtyfiveAB_BB <- MosDet_UPD[3,14] - as.numeric(P2value)
  sixtyfiveAB_AA <- MosDet_UPD[2,14] - as.numeric(M2value)
  sixtyfiveTotal <- abs(sixtyfiveAB_AA) + abs(sixtyfiveAB_BB)

  seventyAB_BB <- MosDet_UPD[3,15] - as.numeric(P2value)
  seventyAB_AA <- MosDet_UPD[2,15] - as.numeric(M2value)
  seventyTotal <- abs(seventyAB_AA) + abs(seventyAB_BB)

  seventyfiveAB_BB <- MosDet_UPD[3,16] - as.numeric(P2value)
  seventyfiveAB_AA <- MosDet_UPD[2,16] - as.numeric(M2value)
  seventyfiveTotal <- abs(seventyfiveAB_AA) + abs(seventyfiveAB_BB)

  eightyAB_BB <- MosDet_UPD[3,17] - as.numeric(P2value)
  eightyAB_AA <- MosDet_UPD[2,17] - as.numeric(M2value)
  eightyTotal <- abs(eightyAB_AA) + abs(eightyAB_BB)

  eightyfiveAB_BB <- MosDet_UPD[3,18] - as.numeric(P2value)
  eightyfiveAB_AA <- MosDet_UPD[2,18] - as.numeric(M2value)
  eightyfiveTotal <- abs(eightyfiveAB_AA) + abs(eightyfiveAB_BB)

  ninetyAB_BB <- MosDet_UPD[3,19] - as.numeric(P2value)
  ninetyAB_AA <- MosDet_UPD[2,19] - as.numeric(M2value)
  ninetyTotal <- abs(ninetyAB_AA) + abs(ninetyAB_BB)

  ninetyfiveAB_BB <- MosDet_UPD[3,20] - as.numeric(P2value)
  ninetyfiveAB_AA <- MosDet_UPD[2,20] - as.numeric(M2value)
  ninetyfiveTotal <- abs(ninetyfiveAB_AA) + abs(ninetyfiveAB_BB)

  onehundredAB_BB <- MosDet_UPD[3,21] - as.numeric(P2value)
  onehundredAB_AA <- MosDet_UPD[2,21] - as.numeric(M2value)
  onehundredTotal <- abs(onehundredAB_AA) + abs(onehundredAB_BB)

######################

  CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
  Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
  Percentage

  } else if (Parent == "Father"){
  
  zeroAB_BB <- MosDet_UPD[3,1] - as.numeric(M2value)
  zeroAB_AA <- MosDet_UPD[2,1] - as.numeric(P2value)
  zeroTotal <- abs(zeroAB_AA) + abs(zeroAB_BB)
  
  fiveAB_BB <- MosDet_UPD[3,2] - as.numeric(M2value)
  fiveAB_AA <- MosDet_UPD[2,2] - as.numeric(P2value)
  fiveTotal <- abs(fiveAB_AA) + abs(fiveAB_BB)
  
  tenAB_BB <- MosDet_UPD[3,3] - as.numeric(M2value)
  tenAB_AA <- MosDet_UPD[2,3] - as.numeric(P2value)
  tenTotal <- abs(tenAB_AA) + abs(tenAB_BB)
  
  fifteenAB_BB <- MosDet_UPD[3,4] - as.numeric(M2value)
  fifteenAB_AA <- MosDet_UPD[2,4] - as.numeric(P2value)
  fifteenTotal <- abs(fifteenAB_AA) + abs(fifteenAB_BB)
  
  twentyAB_BB <- MosDet_UPD[3,5] - as.numeric(M2value)
  twentyAB_AA <- MosDet_UPD[2,5] - as.numeric(P2value)
  twentyTotal <- abs(twentyAB_AA) + abs(twentyAB_BB)
  
  twentyfiveAB_BB <- MosDet_UPD[3,6] - as.numeric(M2value)
  twentyfiveAB_AA <- MosDet_UPD[2,6] - as.numeric(P2value)
  twentyfiveTotal <- abs(twentyfiveAB_AA) + abs(twentyfiveAB_BB)
  
  thirtyAB_BB <- MosDet_UPD[3,7] - as.numeric(M2value)
  thirtyAB_AA <- MosDet_UPD[2,7] - as.numeric(P2value)
  thirtyTotal <- abs(thirtyAB_AA) + abs(thirtyAB_BB)
  
  thirtyfiveAB_BB <- MosDet_UPD[3,8] - as.numeric(M2value)
  thirtyfiveAB_AA <- MosDet_UPD[2,8] - as.numeric(P2value)
  thirtyfiveTotal <- abs(thirtyfiveAB_AA) + abs(thirtyfiveAB_BB)
  
  fortyAB_BB <- MosDet_UPD[3,9] - as.numeric(M2value)
  fortyAB_AA <- MosDet_UPD[2,9] - as.numeric(P2value)
  fortyTotal <- abs(fortyAB_AA) + abs(fortyAB_BB)
  
  fortyfiveAB_BB <- MosDet_UPD[3,10] - as.numeric(M2value)
  fortyfiveAB_AA <- MosDet_UPD[2,10] - as.numeric(P2value)
  fortyfiveTotal <- abs(fortyfiveAB_AA) + abs(fortyfiveAB_BB)
  
  fiftyAB_BB <- MosDet_UPD[3,11] - as.numeric(M2value)
  fiftyAB_AA <- MosDet_UPD[2,11] - as.numeric(P2value)
  fiftyTotal <- abs(fiftyAB_AA) + abs(fiftyAB_BB)
  
  fiftyfiveAB_BB <- MosDet_UPD[3,12] - as.numeric(M2value)
  fiftyfiveAB_AA <- MosDet_UPD[2,12] - as.numeric(P2value)
  fiftyfiveTotal <- abs(fiftyfiveAB_AA) + abs(fiftyfiveAB_BB)
  
  sixtyAB_BB <- MosDet_UPD[3,13] - as.numeric(M2value)
  sixtyAB_AA <- MosDet_UPD[2,13] - as.numeric(P2value)
  sixtyTotal <- abs(sixtyAB_AA) + abs(sixtyAB_BB)
  
  sixtyfiveAB_BB <- MosDet_UPD[3,14] - as.numeric(M2value)
  sixtyfiveAB_AA <- MosDet_UPD[2,14] - as.numeric(P2value)
  sixtyfiveTotal <- abs(sixtyfiveAB_AA) + abs(sixtyfiveAB_BB)
  
  seventyAB_BB <- MosDet_UPD[3,15] - as.numeric(M2value)
  seventyAB_AA <- MosDet_UPD[2,15] - as.numeric(P2value)
  seventyTotal <- abs(seventyAB_AA) + abs(seventyAB_BB)
  
  seventyfiveAB_BB <- MosDet_UPD[3,16] - as.numeric(M2value)
  seventyfiveAB_AA <- MosDet_UPD[2,16] - as.numeric(P2value)
  seventyfiveTotal <- abs(seventyfiveAB_AA) + abs(seventyfiveAB_BB)
  
  eightyAB_BB <- MosDet_UPD[3,17] - as.numeric(M2value)
  eightyAB_AA <- MosDet_UPD[2,17] - as.numeric(P2value)
  eightyTotal <- abs(eightyAB_AA) + abs(eightyAB_BB)
  
  eightyfiveAB_BB <- MosDet_UPD[3,18] - as.numeric(M2value)
  eightyfiveAB_AA <- MosDet_UPD[2,18] - as.numeric(P2value)
  eightyfiveTotal <- abs(eightyfiveAB_AA) + abs(eightyfiveAB_BB)
  
  ninetyAB_BB <- MosDet_UPD[3,19] - as.numeric(M2value)
  ninetyAB_AA <- MosDet_UPD[2,19] - as.numeric(P2value)
  ninetyTotal <- abs(ninetyAB_AA) + abs(ninetyAB_BB)
  
  ninetyfiveAB_BB <- MosDet_UPD[3,20] - as.numeric(M2value)
  ninetyfiveAB_AA <- MosDet_UPD[2,20] - as.numeric(P2value)
  ninetyfiveTotal <- abs(ninetyfiveAB_AA) + abs(ninetyfiveAB_BB)
  
  onehundredAB_BB <- MosDet_UPD[3,21] - as.numeric(M2value)
  onehundredAB_AA <- MosDet_UPD[2,21] - as.numeric(P2value)
  onehundredTotal <- abs(onehundredAB_AA) + abs(onehundredAB_BB)
  
  ######################
  
  CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))
  Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
  Percentage
  
}



