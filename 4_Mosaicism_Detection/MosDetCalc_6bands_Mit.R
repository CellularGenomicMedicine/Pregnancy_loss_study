###Trisomy Mitosis (6 BAF bands)###

rm(list=ls(all=T))

##Load mosaicism percentage reference files from Conlin et al. 2010.

MosDet_TRISOMY_Mit <- read.table("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/MosaicismDetection/MosDet_TRISOMY_Mit.csv",sep="\t",header=T,stringsAsFactors=F,row.names = 1)

options(digits = 20)

###BAFTest_CV_PL2682 -- 6 bands
Band1   <- "0.979329175706215"
Band2   <- "0.903565606896552"
Band3   <- "0.584196066666667"
Band4   <- "0.349702738596491"
Band5   <- "0.0968285926666667"
Band6   <- "0.0177600185727622"

###Trisomy MI (6 BAF bands)###

zeroAA_AAA <- MosDet_TRISOMY_Mit[1,1] - as.numeric(Band6)
zeroAA_AAB <- MosDet_TRISOMY_Mit[2,1] - as.numeric(Band5)
zeroAB_ABA <- MosDet_TRISOMY_Mit[3,1] - as.numeric(Band4)
zeroAB_ABB <- MosDet_TRISOMY_Mit[4,1] - as.numeric(Band3)
zeroBB_BBA <- MosDet_TRISOMY_Mit[5,1] - as.numeric(Band2)
zeroBB_BBB <- MosDet_TRISOMY_Mit[6,1] - as.numeric(Band1)
zeroTotal <- abs(zeroBB_BBB) + abs(zeroBB_BBA) + abs(zeroAB_ABB) + abs(zeroAB_ABA) + abs(zeroAA_AAB) + abs(zeroAA_AAA)

fiveAA_AAA <- MosDet_TRISOMY_Mit[1,2] - as.numeric(Band6)
fiveAA_AAB <- MosDet_TRISOMY_Mit[2,2] - as.numeric(Band5)
fiveAB_ABA <- MosDet_TRISOMY_Mit[3,2] - as.numeric(Band4)
fiveAB_ABB <- MosDet_TRISOMY_Mit[4,2] - as.numeric(Band3)
fiveBB_BBA <- MosDet_TRISOMY_Mit[5,2] - as.numeric(Band2)
fiveBB_BBB <- MosDet_TRISOMY_Mit[6,2] - as.numeric(Band1)
fiveTotal <- abs(fiveBB_BBB) + abs(fiveBB_BBA) + abs(fiveAB_ABB) + abs(fiveAB_ABA) + abs(fiveAA_AAB) + abs(fiveAA_AAA)

tenAA_AAA <- MosDet_TRISOMY_Mit[1,3] - as.numeric(Band6)
tenAA_AAB <- MosDet_TRISOMY_Mit[2,3] - as.numeric(Band5)
tenAB_ABA <- MosDet_TRISOMY_Mit[3,3] - as.numeric(Band4)
tenAB_ABB <- MosDet_TRISOMY_Mit[4,3] - as.numeric(Band3)
tenBB_BBA <- MosDet_TRISOMY_Mit[5,3] - as.numeric(Band2)
tenBB_BBB <- MosDet_TRISOMY_Mit[6,3] - as.numeric(Band1)
tenTotal <- abs(tenBB_BBB) + abs(tenBB_BBA) + abs(tenAB_ABB) + abs(tenAB_ABA) + abs(tenAA_AAB) + abs(tenAA_AAA)

fifteenAA_AAA <- MosDet_TRISOMY_Mit[1,4] - as.numeric(Band6)
fifteenAA_AAB <- MosDet_TRISOMY_Mit[2,4] - as.numeric(Band5)
fifteenAB_ABA <- MosDet_TRISOMY_Mit[3,4] - as.numeric(Band4)
fifteenAB_ABB <- MosDet_TRISOMY_Mit[4,4] - as.numeric(Band3)
fifteenBB_BBA <- MosDet_TRISOMY_Mit[5,4] - as.numeric(Band2)
fifteenBB_BBB <- MosDet_TRISOMY_Mit[6,4] - as.numeric(Band1)
fifteenTotal <- abs(fifteenBB_BBB) + abs(fifteenBB_BBA) + abs(fifteenAB_ABB) + abs(fifteenAB_ABA) + abs(fifteenAA_AAB) + abs(fifteenAA_AAA)

twentyAA_AAA <- MosDet_TRISOMY_Mit[1,5] - as.numeric(Band6)
twentyAA_AAB <- MosDet_TRISOMY_Mit[2,5] - as.numeric(Band5)
twentyAB_ABA <- MosDet_TRISOMY_Mit[3,5] - as.numeric(Band4)
twentyAB_ABB <- MosDet_TRISOMY_Mit[4,5] - as.numeric(Band3)
twentyBB_BBA <- MosDet_TRISOMY_Mit[5,5] - as.numeric(Band2)
twentyBB_BBB <- MosDet_TRISOMY_Mit[6,5] - as.numeric(Band1)
twentyTotal <- abs(twentyBB_BBB) + abs(twentyBB_BBA) + abs(twentyAB_ABB) + abs(twentyAB_ABA) + abs(twentyAA_AAB) + abs(twentyAA_AAA)

twentyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,6] - as.numeric(Band6)
twentyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,6] - as.numeric(Band5)
twentyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,6] - as.numeric(Band4)
twentyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,6] - as.numeric(Band3)
twentyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,6] - as.numeric(Band2)
twentyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,6] - as.numeric(Band1)
twentyfiveTotal <- abs(twentyfiveBB_BBB) + abs(twentyfiveBB_BBA) + abs(twentyfiveAB_ABB) + abs(twentyfiveAB_ABA) + abs(twentyfiveAA_AAB) + abs(twentyfiveAA_AAA)

thirtyAA_AAA <- MosDet_TRISOMY_Mit[1,7] - as.numeric(Band6)
thirtyAA_AAB <- MosDet_TRISOMY_Mit[2,7] - as.numeric(Band5)
thirtyAB_ABA <- MosDet_TRISOMY_Mit[3,7] - as.numeric(Band4)
thirtyAB_ABB <- MosDet_TRISOMY_Mit[4,7] - as.numeric(Band3)
thirtyBB_BBA <- MosDet_TRISOMY_Mit[5,7] - as.numeric(Band2)
thirtyBB_BBB <- MosDet_TRISOMY_Mit[6,7] - as.numeric(Band1)
thirtyTotal <- abs(thirtyBB_BBB) + abs(thirtyBB_BBA) + abs(thirtyAB_ABB) + abs(thirtyAB_ABA) + abs(thirtyAA_AAB) + abs(thirtyAA_AAA)

thirtyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,8] - as.numeric(Band6)
thirtyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,8] - as.numeric(Band5)
thirtyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,8] - as.numeric(Band4)
thirtyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,8] - as.numeric(Band3)
thirtyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,8] - as.numeric(Band2)
thirtyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,8] - as.numeric(Band1)
thirtyfiveTotal <- abs(thirtyfiveBB_BBB) + abs(thirtyfiveBB_BBA) + abs(thirtyfiveAB_ABB) + abs(thirtyfiveAB_ABA) + abs(thirtyfiveAA_AAB) + abs(thirtyfiveAA_AAA)

fortyAA_AAA <- MosDet_TRISOMY_Mit[1,9] - as.numeric(Band6)
fortyAA_AAB <- MosDet_TRISOMY_Mit[2,9] - as.numeric(Band5)
fortyAB_ABA <- MosDet_TRISOMY_Mit[3,9] - as.numeric(Band4)
fortyAB_ABB <- MosDet_TRISOMY_Mit[4,9] - as.numeric(Band3)
fortyBB_BBA <- MosDet_TRISOMY_Mit[5,9] - as.numeric(Band2)
fortyBB_BBB <- MosDet_TRISOMY_Mit[6,9] - as.numeric(Band1)
fortyTotal <- abs(fortyBB_BBB) + abs(fortyBB_BBA) + abs(fortyAB_ABB) + abs(fortyAB_ABA) + abs(fortyAA_AAB) + abs(fortyAA_AAA)

fortyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,10] - as.numeric(Band6)
fortyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,10] - as.numeric(Band5)
fortyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,10] - as.numeric(Band4)
fortyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,10] - as.numeric(Band3)
fortyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,10] - as.numeric(Band2)
fortyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,10] - as.numeric(Band1)
fortyfiveTotal <- abs(fortyfiveBB_BBB) + abs(fortyfiveBB_BBA) + abs(fortyfiveAB_ABB) + abs(fortyfiveAB_ABA) + abs(fortyfiveAA_AAB) + abs(fortyfiveAA_AAA)

fiftyAA_AAA <- MosDet_TRISOMY_Mit[1,11] - as.numeric(Band6)
fiftyAA_AAB <- MosDet_TRISOMY_Mit[2,11] - as.numeric(Band5)
fiftyAB_ABA <- MosDet_TRISOMY_Mit[3,11] - as.numeric(Band4)
fiftyAB_ABB <- MosDet_TRISOMY_Mit[4,11] - as.numeric(Band3)
fiftyBB_BBA <- MosDet_TRISOMY_Mit[5,11] - as.numeric(Band2)
fiftyBB_BBB <- MosDet_TRISOMY_Mit[6,11] - as.numeric(Band1)
fiftyTotal <- abs(fiftyBB_BBB) + abs(fiftyBB_BBA) + abs(fiftyAB_ABB) + abs(fiftyAB_ABA) + abs(fiftyAA_AAB) + abs(fiftyAA_AAA)

fiftyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,12] - as.numeric(Band6)
fiftyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,12] - as.numeric(Band5)
fiftyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,12] - as.numeric(Band4)
fiftyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,12] - as.numeric(Band3)
fiftyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,12] - as.numeric(Band2)
fiftyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,12] - as.numeric(Band1)
fiftyfiveTotal <- abs(fiftyfiveBB_BBB) + abs(fiftyfiveBB_BBA) + abs(fiftyfiveAB_ABB) + abs(fiftyfiveAB_ABA) + abs(fiftyfiveAA_AAB) + abs(fiftyfiveAA_AAA)

sixtyAA_AAA <- MosDet_TRISOMY_Mit[1,13] - as.numeric(Band6)
sixtyAA_AAB <- MosDet_TRISOMY_Mit[2,13] - as.numeric(Band5)
sixtyAB_ABA <- MosDet_TRISOMY_Mit[3,13] - as.numeric(Band4)
sixtyAB_ABB <- MosDet_TRISOMY_Mit[4,13] - as.numeric(Band3)
sixtyBB_BBA <- MosDet_TRISOMY_Mit[5,13] - as.numeric(Band2)
sixtyBB_BBB <- MosDet_TRISOMY_Mit[6,13] - as.numeric(Band1)
sixtyTotal <- abs(sixtyBB_BBB) + abs(sixtyBB_BBA) + abs(sixtyAB_ABB) + abs(sixtyAB_ABA) + abs(sixtyAA_AAB) + abs(sixtyAA_AAA)

sixtyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,14] - as.numeric(Band6)
sixtyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,14] - as.numeric(Band5)
sixtyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,14] - as.numeric(Band4)
sixtyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,14] - as.numeric(Band3)
sixtyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,14] - as.numeric(Band2)
sixtyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,14] - as.numeric(Band1)
sixtyfiveTotal <- abs(sixtyfiveBB_BBB) + abs(sixtyfiveBB_BBA) + abs(sixtyfiveAB_ABB) + abs(sixtyfiveAB_ABA) + abs(sixtyfiveAA_AAB) + abs(sixtyfiveAA_AAA)

seventyAA_AAA <- MosDet_TRISOMY_Mit[1,15] - as.numeric(Band6)
seventyAA_AAB <- MosDet_TRISOMY_Mit[2,15] - as.numeric(Band5)
seventyAB_ABA <- MosDet_TRISOMY_Mit[3,15] - as.numeric(Band4)
seventyAB_ABB <- MosDet_TRISOMY_Mit[4,15] - as.numeric(Band3)
seventyBB_BBA <- MosDet_TRISOMY_Mit[5,15] - as.numeric(Band2)
seventyBB_BBB <- MosDet_TRISOMY_Mit[6,15] - as.numeric(Band1)
seventyTotal <- abs(seventyBB_BBB) + abs(seventyBB_BBA) + abs(seventyAB_ABB) + abs(seventyAB_ABA) + abs(seventyAA_AAB) + abs(seventyAA_AAA)

seventyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,16] - as.numeric(Band6)
seventyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,16] - as.numeric(Band5)
seventyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,16] - as.numeric(Band4)
seventyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,16] - as.numeric(Band3)
seventyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,16] - as.numeric(Band2)
seventyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,16] - as.numeric(Band1)
seventyfiveTotal <- abs(seventyfiveBB_BBB) + abs(seventyfiveBB_BBA) + abs(seventyfiveAB_ABB) + abs(seventyfiveAB_ABA) + abs(seventyfiveAA_AAB) + abs(seventyfiveAA_AAA)

eightyAA_AAA <- MosDet_TRISOMY_Mit[1,17] - as.numeric(Band6)
eightyAA_AAB <- MosDet_TRISOMY_Mit[2,17] - as.numeric(Band5)
eightyAB_ABA <- MosDet_TRISOMY_Mit[3,17] - as.numeric(Band4)
eightyAB_ABB <- MosDet_TRISOMY_Mit[4,17] - as.numeric(Band3)
eightyBB_BBA <- MosDet_TRISOMY_Mit[5,17] - as.numeric(Band2)
eightyBB_BBB <- MosDet_TRISOMY_Mit[6,17] - as.numeric(Band1)
eightyTotal <- abs(eightyBB_BBB) + abs(eightyBB_BBA) + abs(eightyAB_ABB) + abs(eightyAB_ABA) + abs(eightyAA_AAB) + abs(eightyAA_AAA)

eightyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,18] - as.numeric(Band6)
eightyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,18] - as.numeric(Band5)
eightyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,18] - as.numeric(Band4)
eightyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,18] - as.numeric(Band3)
eightyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,18] - as.numeric(Band2)
eightyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,18] - as.numeric(Band1)
eightyfiveTotal <- abs(eightyfiveBB_BBB) + abs(eightyfiveBB_BBA) + abs(eightyfiveAB_ABB) + abs(eightyfiveAB_ABA) + abs(eightyfiveAA_AAB) + abs(eightyfiveAA_AAA)

ninetyAA_AAA <- MosDet_TRISOMY_Mit[1,19] - as.numeric(Band6)
ninetyAA_AAB <- MosDet_TRISOMY_Mit[2,19] - as.numeric(Band5)
ninetyAB_ABA <- MosDet_TRISOMY_Mit[3,19] - as.numeric(Band4)
ninetyAB_ABB <- MosDet_TRISOMY_Mit[4,19] - as.numeric(Band3)
ninetyBB_BBA <- MosDet_TRISOMY_Mit[5,19] - as.numeric(Band2)
ninetyBB_BBB <- MosDet_TRISOMY_Mit[6,19] - as.numeric(Band1)
ninetyTotal <- abs(ninetyBB_BBB) + abs(ninetyBB_BBA) + abs(ninetyAB_ABB) + abs(ninetyAB_ABA) + abs(ninetyAA_AAB) + abs(ninetyAA_AAA)

ninetyfiveAA_AAA <- MosDet_TRISOMY_Mit[1,20] - as.numeric(Band6)
ninetyfiveAA_AAB <- MosDet_TRISOMY_Mit[2,20] - as.numeric(Band5)
ninetyfiveAB_ABA <- MosDet_TRISOMY_Mit[3,20] - as.numeric(Band4)
ninetyfiveAB_ABB <- MosDet_TRISOMY_Mit[4,20] - as.numeric(Band3)
ninetyfiveBB_BBA <- MosDet_TRISOMY_Mit[5,20] - as.numeric(Band2)
ninetyfiveBB_BBB <- MosDet_TRISOMY_Mit[6,20] - as.numeric(Band1)
ninetyfiveTotal <- abs(ninetyfiveBB_BBB) + abs(ninetyfiveBB_BBA) + abs(ninetyfiveAB_ABB) + abs(ninetyfiveAB_ABA) + abs(ninetyfiveAA_AAB) + abs(ninetyfiveAA_AAA)

onehundredAA_AAA <- MosDet_TRISOMY_Mit[1,21] - as.numeric(Band6)
onehundredAA_AAB <- MosDet_TRISOMY_Mit[2,21] - as.numeric(Band5)
onehundredAB_ABA <- MosDet_TRISOMY_Mit[3,21] - as.numeric(Band4)
onehundredAB_ABB <- MosDet_TRISOMY_Mit[4,21] - as.numeric(Band3)
onehundredBB_BBA <- MosDet_TRISOMY_Mit[5,21] - as.numeric(Band2)
onehundredBB_BBB <- MosDet_TRISOMY_Mit[6,21] - as.numeric(Band1)
onehundredTotal <- abs(onehundredBB_BBB) + abs(onehundredBB_BBA) + abs(onehundredAB_ABB) + abs(onehundredAB_ABA) + abs(onehundredAA_AAB) + abs(onehundredAA_AAA)

######################

CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))

Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
Percentage
