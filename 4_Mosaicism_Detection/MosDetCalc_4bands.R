###Trisomy MII_Mit (4 BAF bands)###

rm(list=ls(all=T))

##Load mosaicism percentage reference files from Conlin et al. 2010.
MosDet_TRISOMY_MII_Mit <- read.table("/Users/G10039937/Surfdrive/ClinicalGenetics/Miscarriage/MosaicismDetection/MosDet_TRISOMY_MII_Mit.csv",sep="\t",header=T,stringsAsFactors=F,row.names = 1)

options(digits = 16)

###BAF - 4 bands
Band1   <- "0.984285513337574"
Band2   <- "0.645075890176089"
Band3   <- "0.328900497261774"
Band4   <- "0.0111009247790074"

###Trisomy MII_Mit (4 BAF bands)###

zeroAA_AAA <- MosDet_TRISOMY_MII_Mit[1,1] - as.numeric(Band4)
zeroAB_ABA <- MosDet_TRISOMY_MII_Mit[2,1] - as.numeric(Band3)
zeroAB_ABB <- MosDet_TRISOMY_MII_Mit[3,1] - as.numeric(Band2)
zeroBB_BBB <- MosDet_TRISOMY_MII_Mit[4,1] - as.numeric(Band1)
zeroTotal <- abs(zeroBB_BBB) + abs(zeroAB_ABA) + abs(zeroAB_ABB) + abs(zeroAA_AAA)

fiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,2] - as.numeric(Band4)
fiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,2] - as.numeric(Band3)
fiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,2] - as.numeric(Band2)
fiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,2] - as.numeric(Band1)
fiveTotal <- abs(fiveBB_BBB) + abs(fiveAB_ABA) + abs(fiveAB_ABB) + abs(fiveAA_AAA)

tenAA_AAA <- MosDet_TRISOMY_MII_Mit[1,3] - as.numeric(Band4)
tenAB_ABA <- MosDet_TRISOMY_MII_Mit[2,3] - as.numeric(Band3)
tenAB_ABB <- MosDet_TRISOMY_MII_Mit[3,3] - as.numeric(Band2)
tenBB_BBB <- MosDet_TRISOMY_MII_Mit[4,3] - as.numeric(Band1)
tenTotal <- abs(tenBB_BBB) + abs(tenAB_ABA) + abs(tenAB_ABB) + abs(tenAA_AAA)

fifteenAA_AAA <- MosDet_TRISOMY_MII_Mit[1,4] - as.numeric(Band4)
fifteenAB_ABA <- MosDet_TRISOMY_MII_Mit[2,4] - as.numeric(Band3)
fifteenAB_ABB <- MosDet_TRISOMY_MII_Mit[3,4] - as.numeric(Band2)
fifteenBB_BBB <- MosDet_TRISOMY_MII_Mit[4,4] - as.numeric(Band1)
fifteenTotal <- abs(fifteenBB_BBB) + abs(fifteenAB_ABA) + abs(fifteenAB_ABB) + abs(fifteenAA_AAA)

twentyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,5] - as.numeric(Band4)
twentyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,5] - as.numeric(Band3)
twentyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,5] - as.numeric(Band2)
twentyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,5] - as.numeric(Band1)
twentyTotal <- abs(twentyBB_BBB) + abs(twentyAB_ABA) + abs(twentyAB_ABB) + abs(twentyAA_AAA)

twentyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,6] - as.numeric(Band4)
twentyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,6] - as.numeric(Band3)
twentyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,6] - as.numeric(Band2)
twentyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,6] - as.numeric(Band1)
twentyfiveTotal <- abs(twentyfiveBB_BBB) + abs(twentyfiveAB_ABA) + abs(twentyfiveAB_ABB) + abs(twentyfiveAA_AAA)


thirtyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,7] - as.numeric(Band4)
thirtyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,7] - as.numeric(Band3)
thirtyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,7] - as.numeric(Band2)
thirtyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,7] - as.numeric(Band1)
thirtyTotal <- abs(thirtyBB_BBB) + abs(thirtyAB_ABA) + abs(thirtyAB_ABB) + abs(thirtyAA_AAA)

thirtyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,8] - as.numeric(Band4)
thirtyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,8] - as.numeric(Band3)
thirtyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,8] - as.numeric(Band2)
thirtyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,8] - as.numeric(Band1)
thirtyfiveTotal <- abs(thirtyfiveBB_BBB) + abs(thirtyfiveAB_ABA) + abs(thirtyfiveAB_ABB) + abs(thirtyfiveAA_AAA)

fortyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,9] - as.numeric(Band4)
fortyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,9] - as.numeric(Band3)
fortyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,9] - as.numeric(Band2)
fortyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,9] - as.numeric(Band1)
fortyTotal <- abs(fortyBB_BBB) + abs(fortyAB_ABA) + abs(fortyAB_ABB) + abs(fortyAA_AAA)

fortyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,10] - as.numeric(Band4)
fortyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,10] - as.numeric(Band3)
fortyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,10] - as.numeric(Band2)
fortyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,10] - as.numeric(Band1)
fortyfiveTotal <- abs(fortyfiveBB_BBB) + abs(fortyfiveAB_ABA) + abs(fortyfiveAB_ABB) + abs(fortyfiveAA_AAA)

fiftyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,11] - as.numeric(Band4)
fiftyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,11] - as.numeric(Band3)
fiftyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,11] - as.numeric(Band2)
fiftyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,11] - as.numeric(Band1)
fiftyTotal <- abs(fiftyBB_BBB) + abs(fiftyAB_ABA) + abs(fiftyAB_ABB) + abs(fiftyAA_AAA)

fiftyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,12] - as.numeric(Band4)
fiftyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,12] - as.numeric(Band3)
fiftyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,12] - as.numeric(Band2)
fiftyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,12] - as.numeric(Band1)
fiftyfiveTotal <- abs(fiftyfiveBB_BBB) + abs(fiftyfiveAB_ABA) + abs(fiftyfiveAB_ABB) + abs(fiftyfiveAA_AAA)

sixtyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,13] - as.numeric(Band4)
sixtyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,13] - as.numeric(Band3)
sixtyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,13] - as.numeric(Band2)
sixtyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,13] - as.numeric(Band1)
sixtyTotal <- abs(sixtyBB_BBB) + abs(sixtyAB_ABA) + abs(sixtyAB_ABB) + abs(sixtyAA_AAA)

sixtyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,14] - as.numeric(Band4)
sixtyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,14] - as.numeric(Band3)
sixtyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,14] - as.numeric(Band2)
sixtyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,14] - as.numeric(Band1)
sixtyfiveTotal <- abs(sixtyfiveBB_BBB) + abs(sixtyfiveAB_ABA) + abs(sixtyfiveAB_ABB) + abs(sixtyfiveAA_AAA)

seventyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,15] - as.numeric(Band4)
seventyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,15] - as.numeric(Band3)
seventyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,15] - as.numeric(Band2)
seventyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,15] - as.numeric(Band1)
seventyTotal <- abs(seventyBB_BBB) + abs(seventyAB_ABA) + abs(seventyAB_ABB) + abs(seventyAA_AAA)

seventyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,16] - as.numeric(Band4)
seventyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,16] - as.numeric(Band3)
seventyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,16] - as.numeric(Band2)
seventyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,16] - as.numeric(Band1)
seventyfiveTotal <- abs(seventyfiveBB_BBB) + abs(seventyfiveAB_ABA) + abs(seventyfiveAB_ABB) + abs(seventyfiveAA_AAA)

eightyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,17] - as.numeric(Band4)
eightyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,17] - as.numeric(Band3)
eightyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,17] - as.numeric(Band2)
eightyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,17] - as.numeric(Band1)
eightyTotal <- abs(eightyBB_BBB) + abs(eightyAB_ABA) + abs(eightyAB_ABB) + abs(eightyAA_AAA)

eightyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,18] - as.numeric(Band4)
eightyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,18] - as.numeric(Band3)
eightyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,18] - as.numeric(Band2)
eightyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,18] - as.numeric(Band1)
eightyfiveTotal <- abs(eightyfiveBB_BBB) + abs(eightyfiveAB_ABA) + abs(eightyfiveAB_ABB) + abs(eightyfiveAA_AAA)

ninetyAA_AAA <- MosDet_TRISOMY_MII_Mit[1,19] - as.numeric(Band4)
ninetyAB_ABA <- MosDet_TRISOMY_MII_Mit[2,19] - as.numeric(Band3)
ninetyAB_ABB <- MosDet_TRISOMY_MII_Mit[3,19] - as.numeric(Band2)
ninetyBB_BBB <- MosDet_TRISOMY_MII_Mit[4,19] - as.numeric(Band1)
ninetyTotal <- abs(ninetyBB_BBB) + abs(ninetyAB_ABA) + abs(ninetyAB_ABB) + abs(ninetyAA_AAA)

ninetyfiveAA_AAA <- MosDet_TRISOMY_MII_Mit[1,20] - as.numeric(Band4)
ninetyfiveAB_ABA <- MosDet_TRISOMY_MII_Mit[2,20] - as.numeric(Band3)
ninetyfiveAB_ABB <- MosDet_TRISOMY_MII_Mit[3,20] - as.numeric(Band2)
ninetyfiveBB_BBB <- MosDet_TRISOMY_MII_Mit[4,20] - as.numeric(Band1)
ninetyfiveTotal <- abs(ninetyfiveBB_BBB) + abs(ninetyfiveAB_ABA) + abs(ninetyfiveAB_ABB) + abs(ninetyfiveAA_AAA)

onehundredAA_AAA <- MosDet_TRISOMY_MII_Mit[1,21] - as.numeric(Band4)
onehundredAB_ABA <- MosDet_TRISOMY_MII_Mit[2,21] - as.numeric(Band3)
onehundredAB_ABB <- MosDet_TRISOMY_MII_Mit[3,21] - as.numeric(Band2)
onehundredBB_BBB <- MosDet_TRISOMY_MII_Mit[4,21] - as.numeric(Band1)
onehundredTotal <- abs(onehundredBB_BBB) + abs(onehundredAB_ABA) + abs(onehundredAB_ABB) + abs(onehundredAA_AAA)

######################3

CompletePerc <- data.frame(Percentage = c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%"),value = c(as.numeric(zeroTotal),as.numeric(fiveTotal),as.numeric(tenTotal),as.numeric(fifteenTotal),as.numeric(twentyTotal),as.numeric(twentyfiveTotal),as.numeric(thirtyTotal),as.numeric(thirtyfiveTotal),as.numeric(fortyTotal),as.numeric(fortyfiveTotal),as.numeric(fiftyTotal),as.numeric(fiftyfiveTotal),as.numeric(sixtyTotal),as.numeric(sixtyfiveTotal),as.numeric(seventyTotal),as.numeric(seventyfiveTotal),as.numeric(eightyTotal),as.numeric(eightyfiveTotal),as.numeric(ninetyTotal),as.numeric(ninetyfiveTotal),as.numeric(onehundredTotal)))

Percentage <- CompletePerc[which(CompletePerc$value == min(CompletePerc$value)),]
Percentage
