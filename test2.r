#create an Hknot data frame
Hknotresults2=NULL


#2nd sided test
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
#find rr BLCA
rr2function<- function(x)
{
	505-x+1
}
rrBLCA <- rr2function(BLCAorderrank)
BLCAH<- Hknot(rrBLCA)
Hknotresults2 <-cbind(Hknotresults2, BLCAH)

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
#function for rr in BRCA
rr2function<- function(x)
{
	554-x+1
}
rrBRCA <- rr2function(BRCAorderrank)
BRCAH<- Hknot(rrBRCA)
Hknotresults2<- cbind(Hknotresults2, BRCAH)

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
rr2function<- function(x)
{
	327-x+1
}
rrCESC <- rr2function(CESCorderrank)
CESCH<- Hknot(rrCESC)
Hknotresults2<- cbind(Hknotresults2, CESCH)

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
rr2function<- function(x)
{
	430-x+1
}
rrCHOL <- rr2function(CHOLorderrank)
CHOLH<- Hknot(rrCHOL)
Hknotresults2<- cbind(Hknotresults2, CHOLH)

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
rr2function<- function(x)
{
	421-x+1
}
rrESCEA <- rr2function(ESCAorderrank)
ESCAH<- Hknot(rrESCEA)
Hknotresults2<- cbind(Hknotresults2, ESCAH)

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
rr2function<- function(x)
{
	585-x+1
}
rrHNSC <- rr2function(HNSCorderrank)
HNSCH<- Hknot(rrHNSC)
Hknotresults2<- cbind(Hknotresults2, HNSCH)

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
rr2function<- function(x)
{
	491-x+1
}
rrKICH <- rr2function(KICHorderrank)
KICHH<- Hknot(rrKICH)
Hknotresults2<- cbind(Hknotresults2, KICHH)

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
rr2function<- function(x)
{
	498-x+1
}
rrKIRC <- rr2function(KIRCorderrank)
KIRCH<- Hknot(rrKIRC)
Hknotresults2<- cbind(Hknotresults2, KIRCH)

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
rr2function<- function(x)
{
	526-x+1
}
rrKIRP <- rr2function(KIRPorderrank)
KIRPH<- Hknot(rrKIRP)
Hknotresults2<- cbind(Hknotresults2, KIRPH)

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
rr2function<- function(x)
{
	561-x+1
}
rrLIHC <- rr2function(LIHCorderrank)
LIHCH<- Hknot(rrLIHC)
Hknotresults2<- cbind(Hknotresults2, LIHCH)

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
rr2function<- function(x)
{
	465-x+1
}
rrLUAD <- rr2function(LUADorderrank)
LUADH<- Hknot(rrLUAD)
Hknotresults2<- cbind(Hknotresults2, LUADH)

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
rr2function<- function(x)
{
	569-x+1
}
rrLUSC <- rr2function(LUSCorderrank)
LUSCH<- Hknot(rrLUSC)
Hknotresults2<- cbind(Hknotresults2, LUSCH)

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
rr2function<- function(x)
{
	395-x+1
}
rrPAAD <- rr2function(PAADorderrank)
PAADH<- Hknot(rrPAAD)
Hknotresults2<- cbind(Hknotresults2, PAADH)

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
rr2function<- function(x)
{
	338-x+1
}
rrPCPG <- rr2function(PCPGorderrank)
PCPGH<- Hknot(rrPCPG)
Hknotresults2<- cbind(Hknotresults2, PCPGH)

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
rr2function<- function(x)
{
	501-x+1
}
rrPRAD <- rr2function(PRADorderrank)
PRADH<- Hknot(rrPRAD)
Hknotresults2<- cbind(Hknotresults2, PRADH)

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
rr2function<- function(x)
{
	453-x+1
}
rrSTAD <- rr2function(STADorderrank)
STADH<- Hknot(rrSTAD)
Hknotresults2<- cbind(Hknotresults2, STADH)

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
rr2function<- function(x)
{
	583-x+1
}
rrTHCA <- rr2function(THCAorderrank)
THCAH<- Hknot(rrTHCA)
Hknotresults2<- cbind(Hknotresults2, THCAH)

#THYM
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2ratios.txt"))
THYMorder <- THYM[order(-THYM[,c("pvalues")]),, drop=FALSE]

#rid NA rows
THYMorderp<- THYMorder[,1]
THYMH=NA
Hknotresults2<- cbind(Hknotresults2, THYMH)

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
rr2function<- function(x)
{
	527-x+1
}
rrUCEC <- rr2function(UCECorderrank)
UCECH<- Hknot(rrUCEC)
Hknotresults2<- cbind(Hknotresults2, UCECH)

#save results
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/Tests from article")
write.table(Hknotresults2, file="Hknot results for the 2nd negative side.txt", sep=" ")

positivetail<- read.table(file="Hknot results by rank.txt", header=TRUE)
negativetail<- read.table(file="Hknot results for the 2nd negative side.txt", header=TRUE)
negativetail<- negativetail[,-c(18)]
positivetail<- positivetail[,-c(18)]
negativetail<- -1*negativetail


RECfunction<- function(x)
{
	-log10(2*(x))
	
}
REC<-RECfunction(positivetail)
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/Tests from article")
write.table(REC, file="REC Scores from the first positive positivetail.txt", sep=" ")



