setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
#function for rr in BLCA
rrfunction<- function(x)
{
	(x/505)-(1/(2*505))
}

BLCA<- data.frame(read.table(file="BLCA pvalues and rvalues of log2ratios.txt"))
BLCAorder <- BLCA[order(-BLCA[,c("pvalues")]),, drop=FALSE]
#rid NA rows
BLCAorderp<- BLCAorder[,1]
BLCAorderNA<- BLCAorderp[1:505]
#rank
BLCAorderrank<- data.frame(rank(BLCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(BLCAorderrank)= rownames(BLCAorder[c(1:505),])
#find rr
rrBLCA <- rrfunction(BLCAorderrank)

#find Hknot
Hknot<- function(x)
{
	-2*(sum(log(x)))
}

BLCAH<- Hknot(rrBLCA)

#create an Hknot data frame
Hknotresults=NULL
Hknotresults <-cbind(Hknotresults, BLCAH)

#BRCA
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2ratios.txt"))
BRCAorder <- BRCA[order(-BRCA[,c("pvalues")]),, drop=FALSE]

#rid NA rows
BRCAorderp<- BRCAorder[,1]
BRCAorderNA<- BRCAorderp[1:554]
#rank
BRCAorderrank<- data.frame(rank(BRCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(BRCAorderrank)= rownames(BRCAorder[c(1:554),])

#find rr
#554 is the number of mirna that is not NA 
#function for rr in BRCA
rrfunction<- function(x)
{
	(x/554)-(1/(2*554))
}
rrBRCA <- rrfunction(BRCAorderrank)
BRCAH<- Hknot(rrBRCA)
Hknotresults<- cbind(Hknotresults, BRCAH)


#CESC
CESC<- data.frame(read.table(file="CESC pvalues and rvalues of log2ratios.txt"))
CESCorder <- CESC[order(-CESC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
CESCorderp<- CESCorder[,1]
CESCorderNA<- CESCorderp[1:327]
#rank
CESCorderrank<- data.frame(rank(CESCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(CESCorderrank)= rownames(CESCorder[c(1:327),])

#find rr
#function for rr in CESC
rrfunction<- function(x)
{
	(x/327)-(1/(2*327))
}
rrCESC <- rrfunction(CESCorderrank)
CESCH<- Hknot(rrCESC)
Hknotresults<- cbind(Hknotresults, CESCH)

#CHOL
CHOL<- data.frame(read.table(file="CHOL pvalues and rvalues of log2ratios.txt"))
CHOLorder <- CHOL[order(-CHOL[,c("pvalues")]),, drop=FALSE]

#rid NA rows
CHOLorderp<- CHOLorder[,1]
CHOLorderNA<- CHOLorderp[1:430]
#rank
CHOLorderrank<- data.frame(rank(CHOLorderNA, na.last = TRUE, ties.method=c("min")))
rownames(CHOLorderrank)= rownames(CHOLorder[c(1:430),])

#find rr
#function for rr in CHOL
rrfunction<- function(x)
{
	(x/430)-(1/(2*430))
}
rrCHOL <- rrfunction(CHOLorderrank)
CHOLH<- Hknot(rrCHOL)
Hknotresults<- cbind(Hknotresults, CHOLH)

#ESCA
ESCA<- data.frame(read.table(file="ESCA pvalues and rvalues of log2ratios.txt"))
ESCAorder <- ESCA[order(-ESCA[,c("pvalues")]),, drop=FALSE]

#rid NA rows
ESCAorderp<- ESCAorder[,1]
ESCAorderNA<- ESCAorderp[1:421]
#rank
ESCAorderrank<- data.frame(rank(ESCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(ESCAorderrank)= rownames(ESCAorder[c(1:421),])

#find rr
#function for rr in ESCA
rrfunction<- function(x)
{
	(x/421)-(1/(2*421))
}
rrESCEA <- rrfunction(ESCAorderrank)
ESCAH<- Hknot(rrESCEA)
Hknotresults<- cbind(Hknotresults, ESCAH)

#HNSC
HNSC<- data.frame(read.table(file="HNSC pvalues and rvalues of log2ratios.txt"))
HNSCorder <- HNSC[order(-HNSC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
HNSCorderp<- HNSCorder[,1]
HNSCorderNA<- HNSCorderp[1:585]
#rank
HNSCorderrank<- data.frame(rank(HNSCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(HNSCorderrank)= rownames(HNSCorder[c(1:585),])

#find rr
#function for rr in HNSC
rrfunction<- function(x)
{
	(x/585)-(1/(2*585))
}
rrHNSC <- rrfunction(HNSCorderrank)
HNSCH<- Hknot(rrHNSC)
Hknotresults<- cbind(Hknotresults, HNSCH)

#KICH
KICH<- data.frame(read.table(file="KICH pvalues and rvalues of log2ratios.txt"))
KICHorder <- KICH[order(-KICH[,c("pvalues")]),, drop=FALSE]

#rid NA rows
KICHorderp<- KICHorder[,1]
KICHorderNA<- KICHorderp[1:491]
#rank
KICHorderrank<- data.frame(rank(KICHorderNA, na.last = TRUE, ties.method=c("min")))
rownames(KICHorderrank)= rownames(KICHorder[c(1:491),])

#find rr
#function for rr in KICH
rrfunction<- function(x)
{
	(x/491)-(1/(2*491))
}
rrKICH <- rrfunction(KICHorderrank)
KICHH<- Hknot(rrKICH)
Hknotresults<- cbind(Hknotresults, KICHH)

#KIRC
KIRC<- data.frame(read.table(file="KIRC pvalues and rvalues of log2ratios.txt"))
KIRCorder <- KIRC[order(-KIRC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
KIRCorderp<- KIRCorder[,1]
KIRCorderNA<- KIRCorderp[1:498]
#rank
KIRCorderrank<- data.frame(rank(KIRCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(KIRCorderrank)= rownames(KIRCorder[c(1:498),])

#find rr
#function for rr in KIRC
rrfunction<- function(x)
{
	(x/498)-(1/(2*498))
}
rrKIRC <- rrfunction(KIRCorderrank)
KIRCH<- Hknot(rrKIRC)
Hknotresults<- cbind(Hknotresults, KIRCH)

#KIRP
KIRP<- data.frame(read.table(file="KIRP pvalues and rvalues of log2ratios.txt"))
KIRPorder <- KIRP[order(-KIRP[,c("pvalues")]),, drop=FALSE]

#rid NA rows
KIRPorderp<- KIRPorder[,1]
KIRPorderNA<- KIRPorderp[1:526]
#rank
KIRPorderrank<- data.frame(rank(KIRPorderNA, na.last = TRUE, ties.method=c("min")))
rownames(KIRPorderrank)= rownames(KIRPorder[c(1:526),])

#find rr
#function for rr in KIRP
rrfunction<- function(x)
{
	(x/526)-(1/(2*526))
}
rrKIRP <- rrfunction(KIRPorderrank)
KIRPH<- Hknot(rrKIRP)
Hknotresults<- cbind(Hknotresults, KIRPH)

#LIHC
LIHC<- data.frame(read.table(file="LIHC pvalues and rvalues of log2ratios.txt"))
LIHCorder <- LIHC[order(-LIHC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
LIHCorderp<- LIHCorder[,1]
LIHCorderNA<- LIHCorderp[1:561]
#rank
LIHCorderrank<- data.frame(rank(LIHCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(LIHCorderrank)= rownames(LIHCorder[c(1:561),])

#find rr
#function for rr in LIHC
rrfunction<- function(x)
{
	(x/561)-(1/(2*561))
}
rrLIHC <- rrfunction(LIHCorderrank)
LIHCH<- Hknot(rrLIHC)
Hknotresults<- cbind(Hknotresults, LIHCH)

#LUAD
LUAD<- data.frame(read.table(file="LUAD pvalues and rvalues of log2ratios.txt"))
LUADorder <- LUAD[order(-LUAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
LUADorderp<- LUADorder[,1]
LUADorderNA<- LUADorderp[1:465]
#rank
LUADorderrank<- data.frame(rank(LUADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(LUADorderrank)= rownames(LUADorder[c(1:465),])

#find rr
#function for rr in LUAD
rrfunction<- function(x)
{
	(x/465)-(1/(2*465))
}
rrLUAD <- rrfunction(LUADorderrank)
LUADH<- Hknot(rrLUAD)
Hknotresults<- cbind(Hknotresults, LUADH)

#LUSC
LUSC<- data.frame(read.table(file="LUSC pvalues and rvalues of log2ratios.txt"))
LUSCorder <- LUSC[order(-LUSC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
LUSCorderp<- LUSCorder[,1]
LUSCorderNA<- LUSCorderp[1:569]
#rank
LUSCorderrank<- data.frame(rank(LUSCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(LUSCorderrank)= rownames(LUSCorder[c(1:569),])

#find rr
#function for rr in LUSC
rrfunction<- function(x)
{
	(x/569)-(1/(2*569))
}
rrLUSC <- rrfunction(LUSCorderrank)
LUSCH<- Hknot(rrLUSC)
Hknotresults<- cbind(Hknotresults, LUSCH)

#PAAD
PAAD<- data.frame(read.table(file="PAAD pvalues and rvalues of log2ratios.txt"))
PAADorder <- PAAD[order(-PAAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
PAADorderp<- PAADorder[,1]
PAADorderNA<- PAADorderp[1:395]
#rank
PAADorderrank<- data.frame(rank(PAADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(PAADorderrank)= rownames(PAADorder[c(1:395),])

#find rr
#function for rr in PAAD
rrfunction<- function(x)
{
	(x/395)-(1/(2*395))
}
rrPAAD <- rrfunction(PAADorderrank)
PAADH<- Hknot(rrPAAD)
Hknotresults<- cbind(Hknotresults, PAADH)

#PCPG
PCPG<- data.frame(read.table(file="PCPG pvalues and rvalues of log2ratios.txt"))
PCPGorder <- PCPG[order(-PCPG[,c("pvalues")]),, drop=FALSE]

#rid NA rows
PCPGorderp<- PCPGorder[,1]
PCPGorderNA<- PCPGorderp[1:338]
#rank
PCPGorderrank<- data.frame(rank(PCPGorderNA, na.last = TRUE, ties.method=c("min")))
rownames(PCPGorderrank)= rownames(PCPGorder[c(1:338),])

#find rr
#function for rr in PCPG
rrfunction<- function(x)
{
	(x/338)-(1/(2*338))
}
rrPCPG <- rrfunction(PCPGorderrank)
PCPGH<- Hknot(rrPCPG)
Hknotresults<- cbind(Hknotresults, PCPGH)

#PRAD
PRAD<- data.frame(read.table(file="PRAD pvalues and rvalues of log2ratios.txt"))
PRADorder <- PRAD[order(-PRAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
PRADorderp<- PRADorder[,1]
PRADorderNA<- PRADorderp[1:501]
#rank
PRADorderrank<- data.frame(rank(PRADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(PRADorderrank)= rownames(PRADorder[c(1:501),])

#find rr
#function for rr in PRAD
rrfunction<- function(x)
{
	(x/501)-(1/(2*501))
}
rrPRAD <- rrfunction(PRADorderrank)
PRADH<- Hknot(rrPRAD)
Hknotresults<- cbind(Hknotresults, PRADH)

#STAD
STAD<- data.frame(read.table(file="STAD pvalues and rvalues of log2ratios.txt"))
STADorder <- STAD[order(-STAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
STADorderp<- STADorder[,1]
STADorderNA<- STADorderp[1:453]
#rank
STADorderrank<- data.frame(rank(STADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(STADorderrank)= rownames(STADorder[c(1:453),])

#find rr
#function for rr in STAD
rrfunction<- function(x)
{
	(x/453)-(1/(2*453))
}
rrSTAD <- rrfunction(STADorderrank)
STADH<- Hknot(rrSTAD)
Hknotresults<- cbind(Hknotresults, STADH)

#THCA
THCA<- data.frame(read.table(file="THCA pvalues and rvalues of log2ratios.txt"))
THCAorder <- THCA[order(-THCA[,c("pvalues")]),, drop=FALSE]

#rid NA rows
THCAorderp<- THCAorder[,1]
THCAorderNA<- THCAorderp[1:583]
#rank
THCAorderrank<- data.frame(rank(THCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(THCAorderrank)= rownames(THCAorder[c(1:583),])

#find rr
#function for rr in THCA
rrfunction<- function(x)
{
	(x/583)-(1/(2*583))
}
rrTHCA <- rrfunction(THCAorderrank)
THCAH<- Hknot(rrTHCA)
Hknotresults<- cbind(Hknotresults, THCAH)

#THYM
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2ratios.txt"))
THYMorder <- THYM[order(-THYM[,c("pvalues")]),, drop=FALSE]

#rid NA rows
THYMorderp<- THYMorder[,1]
THYMH=NA
Hknotresults<- cbind(Hknotresults, THYMH)

#UCEC
UCEC<- data.frame(read.table(file="UCEC pvalues and rvalues of log2ratios.txt"))
UCECorder <- UCEC[order(-UCEC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
UCECorderp<- UCECorder[,1]
UCECorderNA<- UCECorderp[1:527]
#rank
UCECorderrank<- data.frame(rank(UCECorderNA, na.last = TRUE, ties.method=c("min")))
rownames(UCECorderrank)= rownames(UCECorder[c(1:527),])

#find rr
#function for rr in UCEC
rrfunction<- function(x)
{
	(x/527)-(1/(2*527))
}
rrUCEC <- rrfunction(UCECorderrank)
UCECH<- Hknot(rrUCEC)
Hknotresults<- cbind(Hknotresults, UCECH)
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/Tests from article")
write.table(Hknotresults, file="Hknot results by rank.txt", sep=" ")








#rrfunction has the base of every microRNA
rrfunction<- function(x)
{
	(x/1046)-(1/(2*1046))
}
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")

#BLCA
BLCA<- data.frame(read.table(file="BLCA pvalues and rvalues of log2ratios.txt"))
BLCAorder <- BLCA[order(-BLCA[,c("pvalues")]),, drop=FALSE]
#rid NA rows
BLCAorderp<- BLCAorder[,1]
BLCAorderNA<- BLCAorderp[1:505]
#rank
BLCAorderrank<- data.frame(rank(BLCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(BLCAorderrank)= rownames(BLCAorder[c(1:505),])
#find rr
rrBLCA <- rrfunction(BLCAorderrank)

#find Hknot
Hknot<- function(x)
{
	-2*(sum(log(x)))
}

BLCAH<- Hknot(rrBLCA)

#create an Hknot data frame
HknotresultsA=NULL
HknotresultsA <-cbind(HknotresultsA, BLCAH)

#BRCA
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2ratios.txt"))
BRCAorder <- BRCA[order(-BRCA[,c("pvalues")]),, drop=FALSE]

#rid NA rows
BRCAorderp<- BRCAorder[,1]
BRCAorderNA<- BRCAorderp[1:554]
#rank
BRCAorderrank<- data.frame(rank(BRCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(BRCAorderrank)= rownames(BRCAorder[c(1:554),])

rrBRCA <- rrfunction(BRCAorderrank)
BRCAH<- Hknot(rrBRCA)
HknotresultsA<- cbind(HknotresultsA, BRCAH)

#CESC
CESC<- data.frame(read.table(file="CESC pvalues and rvalues of log2ratios.txt"))
CESCorder <- CESC[order(-CESC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
CESCorderp<- CESCorder[,1]
CESCorderNA<- CESCorderp[1:327]
#rank
CESCorderrank<- data.frame(rank(CESCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(CESCorderrank)= rownames(CESCorder[c(1:327),])

rrCESC <- rrfunction(CESCorderrank)
CESCH<- Hknot(rrCESC)
HknotresultsA<- cbind(HknotresultsA, CESCH)

#CHOL
CHOL<- data.frame(read.table(file="CHOL pvalues and rvalues of log2ratios.txt"))
CHOLorder <- CHOL[order(-CHOL[,c("pvalues")]),, drop=FALSE]

#rid NA rows
CHOLorderp<- CHOLorder[,1]
CHOLorderNA<- CHOLorderp[1:430]
#rank
CHOLorderrank<- data.frame(rank(CHOLorderNA, na.last = TRUE, ties.method=c("min")))
rownames(CHOLorderrank)= rownames(CHOLorder[c(1:430),])

rrCHOL <- rrfunction(CHOLorderrank)
CHOLH<- Hknot(rrCHOL)
HknotresultsA<- cbind(HknotresultsA, CHOLH)

#ESCA
ESCA<- data.frame(read.table(file="ESCA pvalues and rvalues of log2ratios.txt"))
ESCAorder <- ESCA[order(-ESCA[,c("pvalues")]),, drop=FALSE]

#rid NA rows
ESCAorderp<- ESCAorder[,1]
ESCAorderNA<- ESCAorderp[1:421]
#rank
ESCAorderrank<- data.frame(rank(ESCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(ESCAorderrank)= rownames(ESCAorder[c(1:421),])

#find rr
#function for rr in ESCA
rrfunction<- function(x)

rrESCEA <- rrfunction(ESCAorderrank)
ESCAH<- Hknot(rrESCEA)
HknotresultsA<- cbind(HknotresultsA, ESCAH)

#HNSC
HNSC<- data.frame(read.table(file="HNSC pvalues and rvalues of log2ratios.txt"))
HNSCorder <- HNSC[order(-HNSC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
HNSCorderp<- HNSCorder[,1]
HNSCorderNA<- HNSCorderp[1:585]
#rank
HNSCorderrank<- data.frame(rank(HNSCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(HNSCorderrank)= rownames(HNSCorder[c(1:585),])

rrHNSC <- rrfunction(HNSCorderrank)
HNSCH<- Hknot(rrHNSC)
HknotresultsA<- cbind(HknotresultsA, HNSCH)

#KICH
KICH<- data.frame(read.table(file="KICH pvalues and rvalues of log2ratios.txt"))
KICHorder <- KICH[order(-KICH[,c("pvalues")]),, drop=FALSE]

#rid NA rows
KICHorderp<- KICHorder[,1]
KICHorderNA<- KICHorderp[1:491]
#rank
KICHorderrank<- data.frame(rank(KICHorderNA, na.last = TRUE, ties.method=c("min")))
rownames(KICHorderrank)= rownames(KICHorder[c(1:491),])

rrKICH <- rrfunction(KICHorderrank)
KICHH<- Hknot(rrKICH)
HknotresultsA<- cbind(HknotresultsA, KICHH)

#KIRC
KIRC<- data.frame(read.table(file="KIRC pvalues and rvalues of log2ratios.txt"))
KIRCorder <- KIRC[order(-KIRC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
KIRCorderp<- KIRCorder[,1]
KIRCorderNA<- KIRCorderp[1:498]
#rank
KIRCorderrank<- data.frame(rank(KIRCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(KIRCorderrank)= rownames(KIRCorder[c(1:498),])

rrKIRC <- rrfunction(KIRCorderrank)
KIRCH<- Hknot(rrKIRC)
HknotresultsA<- cbind(HknotresultsA, KIRCH)

#KIRP
KIRP<- data.frame(read.table(file="KIRP pvalues and rvalues of log2ratios.txt"))
KIRPorder <- KIRP[order(-KIRP[,c("pvalues")]),, drop=FALSE]

#rid NA rows
KIRPorderp<- KIRPorder[,1]
KIRPorderNA<- KIRPorderp[1:526]
#rank
KIRPorderrank<- data.frame(rank(KIRPorderNA, na.last = TRUE, ties.method=c("min")))
rownames(KIRPorderrank)= rownames(KIRPorder[c(1:526),])

rrKIRP <- rrfunction(KIRPorderrank)
KIRPH<- Hknot(rrKIRP)
HknotresultsA<- cbind(HknotresultsA, KIRPH)

#LIHC
LIHC<- data.frame(read.table(file="LIHC pvalues and rvalues of log2ratios.txt"))
LIHCorder <- LIHC[order(-LIHC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
LIHCorderp<- LIHCorder[,1]
LIHCorderNA<- LIHCorderp[1:561]
#rank
LIHCorderrank<- data.frame(rank(LIHCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(LIHCorderrank)= rownames(LIHCorder[c(1:561),])

rrLIHC <- rrfunction(LIHCorderrank)
LIHCH<- Hknot(rrLIHC)
HknotresultsA<- cbind(HknotresultsA, LIHCH)

#LUAD
LUAD<- data.frame(read.table(file="LUAD pvalues and rvalues of log2ratios.txt"))
LUADorder <- LUAD[order(-LUAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
LUADorderp<- LUADorder[,1]
LUADorderNA<- LUADorderp[1:465]
#rank
LUADorderrank<- data.frame(rank(LUADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(LUADorderrank)= rownames(LUADorder[c(1:465),])

rrLUAD <- rrfunction(LUADorderrank)
LUADH<- Hknot(rrLUAD)
HknotresultsA<- cbind(HknotresultsA, LUADH)

#LUSC
LUSC<- data.frame(read.table(file="LUSC pvalues and rvalues of log2ratios.txt"))
LUSCorder <- LUSC[order(-LUSC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
LUSCorderp<- LUSCorder[,1]
LUSCorderNA<- LUSCorderp[1:569]
#rank
LUSCorderrank<- data.frame(rank(LUSCorderNA, na.last = TRUE, ties.method=c("min")))
rownames(LUSCorderrank)= rownames(LUSCorder[c(1:569),])

rrLUSC <- rrfunction(LUSCorderrank)
LUSCH<- Hknot(rrLUSC)
HknotresultsA<- cbind(HknotresultsA, LUSCH)

#PAAD
PAAD<- data.frame(read.table(file="PAAD pvalues and rvalues of log2ratios.txt"))
PAADorder <- PAAD[order(-PAAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
PAADorderp<- PAADorder[,1]
PAADorderNA<- PAADorderp[1:395]
#rank
PAADorderrank<- data.frame(rank(PAADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(PAADorderrank)= rownames(PAADorder[c(1:395),])

rrPAAD <- rrfunction(PAADorderrank)
PAADH<- Hknot(rrPAAD)
HknotresultsA<- cbind(HknotresultsA, PAADH)

#PCPG
PCPG<- data.frame(read.table(file="PCPG pvalues and rvalues of log2ratios.txt"))
PCPGorder <- PCPG[order(-PCPG[,c("pvalues")]),, drop=FALSE]

#rid NA rows
PCPGorderp<- PCPGorder[,1]
PCPGorderNA<- PCPGorderp[1:338]
#rank
PCPGorderrank<- data.frame(rank(PCPGorderNA, na.last = TRUE, ties.method=c("min")))
rownames(PCPGorderrank)= rownames(PCPGorder[c(1:338),])

rrPCPG <- rrfunction(PCPGorderrank)
PCPGH<- Hknot(rrPCPG)
HknotresultsA<- cbind(HknotresultsA, PCPGH)

#PRAD
PRAD<- data.frame(read.table(file="PRAD pvalues and rvalues of log2ratios.txt"))
PRADorder <- PRAD[order(-PRAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
PRADorderp<- PRADorder[,1]
PRADorderNA<- PRADorderp[1:501]
#rank
PRADorderrank<- data.frame(rank(PRADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(PRADorderrank)= rownames(PRADorder[c(1:501),])

rrPRAD <- rrfunction(PRADorderrank)
PRADH<- Hknot(rrPRAD)
HknotresultsA<- cbind(HknotresultsA, PRADH)

#STAD
STAD<- data.frame(read.table(file="STAD pvalues and rvalues of log2ratios.txt"))
STADorder <- STAD[order(-STAD[,c("pvalues")]),, drop=FALSE]

#rid NA rows
STADorderp<- STADorder[,1]
STADorderNA<- STADorderp[1:453]
#rank
STADorderrank<- data.frame(rank(STADorderNA, na.last = TRUE, ties.method=c("min")))
rownames(STADorderrank)= rownames(STADorder[c(1:453),])

rrSTAD <- rrfunction(STADorderrank)
STADH<- Hknot(rrSTAD)
HknotresultsA<- cbind(HknotresultsA, STADH)

#THCA
THCA<- data.frame(read.table(file="THCA pvalues and rvalues of log2ratios.txt"))
THCAorder <- THCA[order(-THCA[,c("pvalues")]),, drop=FALSE]

#rid NA rows
THCAorderp<- THCAorder[,1]
THCAorderNA<- THCAorderp[1:583]
#rank
THCAorderrank<- data.frame(rank(THCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(THCAorderrank)= rownames(THCAorder[c(1:583),])

rrTHCA <- rrfunction(THCAorderrank)
THCAH<- Hknot(rrTHCA)
HknotresultsA<- cbind(HknotresultsA, THCAH)

#THYM
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2ratios.txt"))
THYMorder <- THYM[order(-THYM[,c("pvalues")]),, drop=FALSE]

#rid NA rows
THYMorderp<- THYMorder[,1]
THYMH=NA
HknotresultsA<- cbind(HknotresultsA, THYMH)

#UCEC
UCEC<- data.frame(read.table(file="UCEC pvalues and rvalues of log2ratios.txt"))
UCECorder <- UCEC[order(-UCEC[,c("pvalues")]),, drop=FALSE]

#rid NA rows
UCECorderp<- UCECorder[,1]
UCECorderNA<- UCECorderp[1:527]
#rank
UCECorderrank<- data.frame(rank(UCECorderNA, na.last = TRUE, ties.method=c("min")))
rownames(UCECorderrank)= rownames(UCECorder[c(1:527),])

rrUCEC <- rrfunction(UCECorderrank)
UCECH<- Hknot(rrUCEC)
HknotresultsA<- cbind(HknotresultsA, UCECH)
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/Tests from article")
write.table(HknotresultsA, file="Hknot results with a constant Lrr.txt", sep=" ")














#BRCA
lowpBRCA <- singletailslow(tpresultsBRCA)
highpBRCA<- singletailshigh(lowpBRCA)

#multiply by 2
lowpBRCA2<- lowpBRCA*2
highpBRCA2<- highpBRCA*2

lowandhighBRCA<- cbind(lowpBRCA2, highpBRCA2)
colnames(lowandhighBRCA)=c("lowP", "highP")

#calculate REC score
REC=NULL
for (i in 1:nrow(lowandhighBRCA))
{
	y=NULL
	y<- if(i[,c("lowP")] < i[,c("highP")]) 
	{
		log10(2*i[,c("lowP")])
	} 		else{
				-log10(2*i[,c("highP")])
			}
	RECScores<- cbind(REC, y)
}
