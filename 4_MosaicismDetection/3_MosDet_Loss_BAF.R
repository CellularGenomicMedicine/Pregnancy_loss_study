###Monosomy (4 BAF bands)###

rm(list=ls(all=T))

##Load mosaicism percentage reference files from Conlin et al. 2010.

MosDet_Loss_BAF <- read.csv("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/Paper/1.Paper/Miscarriage_Surfdrive\ kopie/Submission/Github/4_MosaicismDetection/Refs/MosDet_Loss.csv", row.names=1)

#Set level of digits after the comma 
options(digits = 16)

###BAF - 4 bands, only change band 1 and 4 if there is maternal contamination and these are moved away from 1 and 0. 
Band1   <- "0.979718698882682"
Band2   <- "0.669384310126582"
Band3   <- "0.369220026470588"
Band4   <- "0.0203403704546853"

###Monosomy (4 BAF bands)###

{
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
}

###

CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))

Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
Percentage
