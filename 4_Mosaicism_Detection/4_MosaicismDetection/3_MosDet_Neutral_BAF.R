###Trisomy MII_Mit (4 BAF bands)###

rm(list=ls(all=T))

##Load mosaicism percentage reference files from Conlin et al. 2010.

MosDet_Neutral <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive\ kopie/Submission/Github/4_MosaicismDetection/Refs/MosDet_Neutral.csv", row.names=1)

options(digits = 16)

###BAF - 4 bands
Band1   <- "0.981111262849162"
Band2   <- "0.668406166666667"
Band3   <- "0.371275215189873"
Band4   <- "0.0177244540798319"

###Trisomy MII_Mit (4 BAF bands)###

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

######################3

CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))

Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
Percentage
