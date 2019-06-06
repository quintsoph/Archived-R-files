
#read in files
#read in the created files
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results/log2 foldchange ratio by each")
BLCA<- data.frame(read.table(file="BLCA pvalues and rvalues of log2 foldchange ratio.txt"))
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2 foldchange ratio.txt"))
CESC<- data.frame(read.table(file="CESC pvalues and rvalues of log2 foldchange ratio.txt"))
CHOL<- data.frame(read.table(file="CHOL pvalues and rvalues of log2 foldchange ratio.txt"))
ESCA<- data.frame(read.table(file="ESCA pvalues and rvalues of log2 foldchange ratio.txt"))
HNSC<- data.frame(read.table(file="HNSC pvalues and rvalues of log2 foldchange ratio.txt"))
KICH<- data.frame(read.table(file="KICH pvalues and rvalues of log2 foldchange ratio.txt"))
KIRC<- data.frame(read.table(file="KIRC pvalues and rvalues of log2 foldchange ratio.txt"))
KIRP<- data.frame(read.table(file="KIRP pvalues and rvalues of log2 foldchange ratio.txt"))
LIHC<- data.frame(read.table(file="LIHC pvalues and rvalues of log2 foldchange ratio.txt"))
LUAD<- data.frame(read.table(file="LUAD pvalues and rvalues of log2 foldchange ratio.txt"))
LUSC<- data.frame(read.table(file="LUSC pvalues and rvalues of log2 foldchange ratio.txt"))
PAAD<- data.frame(read.table(file="PAAD pvalues and rvalues of log2 foldchange ratio.txt"))
PCPG<- data.frame(read.table(file="PCPG pvalues and rvalues of log2 foldchange ratio.txt"))
PRAD<- data.frame(read.table(file="PRAD pvalues and rvalues of log2 foldchange ratio.txt"))
STAD<- data.frame(read.table(file="STAD pvalues and rvalues of log2 foldchange ratio.txt"))
THCA<- data.frame(read.table(file="THCA pvalues and rvalues of log2 foldchange ratio.txt"))
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2 foldchange ratio.txt"))
UCEC<- data.frame(read.table(file="UCEC pvalues and rvalues of log2 foldchange ratio.txt"))
READ<- data.frame(read.table(file="READ pvalues and rvalues of log2 foldchange ratio.txt"))
FPPP<- data.frame(read.table(file="FPPP pvalues and rvalues of log2 foldchange ratio.txt"))
COAD<- data.frame(read.table(file="COAD pvalues and rvalues of log2 foldchange ratio.txt"))

L=NULL

#BLCA
tBLCA<- t(BLCA)
BLCAp<- tBLCA[-c(2),]
BLCANA<- BLCAp[!is.na(BLCAp)]
BLCANA1<- data.frame(BLCANA)
BLCANA<- data.frame((rank(BLCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(BLCANA)=rownames(BLCANA1)
BLCAs<- data.frame(BLCANA[sample(1:nrow(BLCANA)),])
rownames(BLCAs)=rownames(BLCANA)
BLCAr<- nrow(BLCANA)
L<- cbind(BLCAr, L)

#BRCA
BRCAorder <- BRCA[order(-BRCA[,c("pvalues")]),, drop=FALSE]
tBRCA<- t(BRCAorder)
BRCAp<- tBRCA[-c(2),]
BRCANA<- BRCAp[!is.na(BRCAp)]
BRCANA1<- data.frame(BRCANA)
BRCANA<- data.frame((rank(BRCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(BRCANA)=rownames(BRCANA1)
BRCAs<- data.frame(BRCANA[sample(1:nrow(BRCANA)),])
rownames(BRCAs)=rownames(BRCANA)
BRCAr<- nrow(BRCANA)
L<- cbind(BRCAr, L)

#CESC
CESCorder <- CESC[order(-CESC[,c("pvalues")]),, drop=FALSE]
tCESC<- t(CESCorder)
CESCp<- tCESC[-c(2),]
CESCNA<- CESCp[!is.na(CESCp)]
CESCNA1<- data.frame(CESCNA)
CESCNA<- data.frame((rank(CESCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(CESCNA)=rownames(CESCNA1)
CESCs<- data.frame(CESCNA[sample(1:nrow(CESCNA)),])
rownames(CESCs)=rownames(CESCNA)
CESCr<- nrow(CESCNA)
L<- cbind(CESCr, L)

#CHOL
CHOLorder <- CHOL[order(-CHOL[,c("pvalues")]),, drop=FALSE]
tCHOL<- t(CHOLorder)
CHOLp<- tCHOL[-c(2),]
CHOLNA<- CHOLp[!is.na(CHOLp)]
CHOLNA1<- data.frame(CHOLNA)
CHOLNA<- data.frame((rank(CHOLNA1, na.last = TRUE, ties.method=c("min"))))
rownames(CHOLNA)=rownames(CHOLNA1)
CHOLs<- data.frame(CHOLNA[sample(1:nrow(CHOLNA)),])
rownames(CHOLs)=rownames(CHOLNA)
CHOLr<- nrow(CHOLNA)
L<- cbind(CHOLr, L)

#ESCA
ESCAorder <- ESCA[order(-ESCA[,c("pvalues")]),, drop=FALSE]
tESCA<- t(ESCAorder)
ESCAp<- tESCA[-c(2),]
ESCANA<- ESCAp[!is.na(ESCAp)]
ESCANA1<- data.frame(ESCANA)
ESCANA<- data.frame((rank(ESCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(ESCANA)=rownames(ESCANA1)
ESCAs<- data.frame(ESCANA[sample(1:nrow(ESCANA)),])
rownames(ESCAs)=rownames(ESCANA)
ESCAr<- nrow(ESCANA)
L<- cbind(ESCAr, L)

#HNSC
HNSCorder <- HNSC[order(-HNSC[,c("pvalues")]),, drop=FALSE]
tHNSC<- t(HNSCorder)
HNSCp<- tHNSC[-c(2),]
HNSCNA<- HNSCp[!is.na(HNSCp)]
HNSCNA1<- data.frame(HNSCNA)
HNSCNA<- data.frame((rank(HNSCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(HNSCNA)=rownames(HNSCNA1)
HNSCs<- data.frame(HNSCNA[sample(1:nrow(HNSCNA)),])
rownames(HNSCs)=rownames(HNSCNA)
HNSCr<- nrow(HNSCNA)
L<- cbind(HNSCr, L)

#KICH
KICHorder <- KICH[order(-KICH[,c("pvalues")]),, drop=FALSE]
tKICH<- t(KICHorder)
KICHp<- tKICH[-c(2),]
KICHNA<- KICHp[!is.na(KICHp)]
KICHNA1<- data.frame(KICHNA)
KICHNA<- data.frame((rank(KICHNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KICHNA)=rownames(KICHNA1)
KICHs<- data.frame(KICHNA[sample(1:nrow(KICHNA)),])
rownames(KICHs)=rownames(KICHNA)
KICHr<- nrow(KICHNA)
L<- cbind(KICHr, L)

#KIRC
KIRCorder <- KIRC[order(-KIRC[,c("pvalues")]),, drop=FALSE]
tKIRC<- t(KIRCorder)
KIRCp<- tKIRC[-c(2),]
KIRCNA<- KIRCp[!is.na(KIRCp)]
KIRCNA1<- data.frame(KIRCNA)
KIRCNA<- data.frame((rank(KIRCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KIRCNA)=rownames(KIRCNA1)
KIRCs<- data.frame(KIRCNA[sample(1:nrow(KIRCNA)),])
rownames(KIRCs)=rownames(KIRCNA)
KIRCr<- nrow(KIRCNA)
L<- cbind(KIRCr, L)

#KIRP
KIRPorder <- KIRP[order(-KIRP[,c("pvalues")]),, drop=FALSE]
tKIRP<- t(KIRPorder)
KIRPp<- tKIRP[-c(2),]
KIRPNA<- KIRPp[!is.na(KIRPp)]
KIRPNA1<- data.frame(KIRPNA)
KIRPNA<- data.frame((rank(KIRPNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KIRPNA)=rownames(KIRPNA1)
KIRPs<- data.frame(KIRPNA[sample(1:nrow(KIRPNA)),])
rownames(KIRPs)=rownames(KIRPNA)
KIRPr<- nrow(KIRPNA)
L<- cbind(KIRPr, L)

#LIHC
LIHCorder <- LIHC[order(-LIHC[,c("pvalues")]),, drop=FALSE]
tLIHC<- t(LIHCorder)
LIHCp<- tLIHC[-c(2),]
LIHCNA<- LIHCp[!is.na(LIHCp)]
LIHCNA1<- data.frame(LIHCNA)
LIHCNA<- data.frame((rank(LIHCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LIHCNA)=rownames(LIHCNA1)
LIHCs<- data.frame(LIHCNA[sample(1:nrow(LIHCNA)),])
rownames(LIHCs)=rownames(LIHCNA)
LIHCr<- nrow(LIHCNA)
L<- cbind(LIHCr, L)

#LUAD
LUADorder <- LUAD[order(-LUAD[,c("pvalues")]),, drop=FALSE]
tLUAD<- t(LUADorder)
LUADp<- tLUAD[-c(2),]
LUADNA<- LUADp[!is.na(LUADp)]
LUADNA1<- data.frame(LUADNA)
LUADNA<- data.frame((rank(LUADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LUADNA)=rownames(LUADNA1)
LUADs<- data.frame(LUADNA[sample(1:nrow(LUADNA)),])
rownames(LUADs)=rownames(LUADNA)
LUADr<- nrow(LUADNA)
L<- cbind(LUADr, L)

#LUSC
LUSCorder <- LUSC[order(-LUSC[,c("pvalues")]),, drop=FALSE]
tLUSC<- t(LUSCorder)
LUSCp<- tLUSC[-c(2),]
LUSCNA<- LUSCp[!is.na(LUSCp)]
LUSCNA1<- data.frame(LUSCNA)
LUSCNA<- data.frame((rank(LUSCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LUSCNA)=rownames(LUSCNA1)
LUSCs<- data.frame(LUSCNA[sample(1:nrow(LUSCNA)),])
rownames(LUSCs)=rownames(LUSCNA)
LUSCr<- nrow(LUSCNA)
L<- cbind(LUSCr, L)

#PAAD
PAADorder <- PAAD[order(-PAAD[,c("pvalues")]),, drop=FALSE]
tPAAD<- t(PAADorder)
PAADp<- tPAAD[-c(2),]
PAADNA<- PAADp[!is.na(PAADp)]
PAADNA1<- data.frame(PAADNA)
PAADNA<- data.frame((rank(PAADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(PAADNA)=rownames(PAADNA1)
PAADs<- data.frame(PAADNA[sample(1:nrow(PAADNA)),])
rownames(PAADs)=rownames(PAADNA)
PAADr<- nrow(PAADNA)
L<- cbind(PAADr, L)

#PCPG
PCPGorder <- PCPG[order(-PCPG[,c("pvalues")]),, drop=FALSE]
tPCPG<- t(PCPGorder)
PCPGp<- tPCPG[-c(2),]
PCPGNA<- PCPGp[!is.na(PCPGp)]
PCPGNA1<- data.frame(PCPGNA)
PCPGNA<- data.frame((rank(PCPGNA1, na.last = TRUE, ties.method=c("min"))))
rownames(PCPGNA)=rownames(PCPGNA1)
PCPGs<- data.frame(PCPGNA[sample(1:nrow(PCPGNA)),])
rownames(PCPGs)=rownames(PCPGNA)
PCPGr<- nrow(PCPGNA)
L<- cbind(PCPGr, L)

#PRAD
PRADorder <- PRAD[order(-PRAD[,c("pvalues")]),, drop=FALSE]
tPRAD<- t(PRADorder)
PRADp<- tPRAD[-c(2),]
PRADNA<- PRADp[!is.na(PRADp)]
PRADNA1<- data.frame(PRADNA)
PRADNA<- data.frame((rank(PRADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(PRADNA)=rownames(PRADNA1)
PRADs<- data.frame(PRADNA[sample(1:nrow(PRADNA)),])
rownames(PRADs)=rownames(PRADNA)
PRADr<- nrow(PRADNA)
L<- cbind(PRADr, L)

#STAD
STADorder <- STAD[order(-STAD[,c("pvalues")]),, drop=FALSE]
tSTAD<- t(STADorder)
STADp<- tSTAD[-c(2),]
STADNA<- STADp[!is.na(STADp)]
STADNA1<- data.frame(STADNA)
STADNA<- data.frame((rank(STADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(STADNA)=rownames(STADNA1)
STADs<- data.frame(STADNA[sample(1:nrow(STADNA)),])
rownames(STADs)=rownames(STADNA)
STADr<- nrow(STADNA)
L<- cbind(STADr, L)

#THCA
THCAorder <- THCA[order(-THCA[,c("pvalues")]),, drop=FALSE]
tTHCA<- t(THCAorder)
THCAp<- tTHCA[-c(2),]
THCANA<- THCAp[!is.na(THCAp)]
THCANA1<- data.frame(THCANA)
THCANA<- data.frame((rank(THCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(THCANA)=rownames(THCANA1)
THCAs<- data.frame(THCANA[sample(1:nrow(THCANA)),])
rownames(THCAs)=rownames(THCANA)
THCAr<- nrow(THCANA)
L<- cbind(THCAr, L)

#UCEC
UCECorder <- UCEC[order(-UCEC[,c("pvalues")]),, drop=FALSE]
tUCEC<- t(UCECorder)
UCECp<- tUCEC[-c(2),]
UCECNA<- UCECp[!is.na(UCECp)]
UCECNA1<- data.frame(UCECNA)
UCECNA<- data.frame((rank(UCECNA1, na.last = TRUE, ties.method=c("min"))))
rownames(UCECNA)=rownames(UCECNA1)
UCECs<- data.frame(UCECNA[sample(1:nrow(UCECNA)),])
rownames(UCECs)=rownames(UCECNA)
UCECr<- nrow(UCECNA)
L<- cbind(UCECr, L)

#FPPP
FPPPorder <- FPPP[order(-FPPP[,c("pvalues")]),, drop=FALSE]
tFPPP<- t(FPPPorder)
FPPPp<- tFPPP[-c(2),]
FPPPNA<- FPPPp[!is.na(FPPPp)]
FPPPNA1<- data.frame(FPPPNA)
FPPPNA<- data.frame((rank(FPPPNA1, na.last = TRUE, ties.method=c("min"))))
rownames(FPPPNA)=rownames(FPPPNA1)
FPPPs<- data.frame(FPPPNA[sample(1:nrow(FPPPNA)),])
rownames(FPPPs)=rownames(FPPPNA)
FPPPr<- nrow(FPPPNA)
L<- cbind(FPPPr, L)

#COAD
COADorder <- COAD[order(-COAD[,c("pvalues")]),, drop=FALSE]
tCOAD<- t(COADorder)
COADp<- tCOAD[-c(2),]
COADNA<- COADp[!is.na(COADp)]
COADNA1<- data.frame(COADNA)
COADNA<- data.frame((rank(COADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(COADNA)=rownames(COADNA1)
COADs<- data.frame(COADNA[sample(1:nrow(COADNA)),])
rownames(COADs)=rownames(COADNA)
COADr<- nrow(COADNA)
L<- cbind(COADr, L)


microcut<- function(x)
{
	mirna=NULL
	BLCA<- BLCAs[c(x),]
	BRCA<- cbind(BRCAs[c(x),], BLCA)
	CESC<- cbind(CESCs[c(x),], BRCA)
	CHOL<- cbind(CHOLs[c(x),], CESC)
	ESCA<- cbind(ESCAs[c(x),], CHOL)
	HNSC<- cbind(HNSCs[c(x),], ESCA)
	KICH<- cbind(KICHs[c(x),], HNSC)
	KIRC<- cbind(KIRCs[c(x),], KICH)
	KIRP<- cbind(KIRPs[c(x),], KIRC)
	LIHC<- cbind(LIHCs[c(x),], KIRP)
	LUAD<- cbind(LUADs[c(x),], LIHC)
	LUSC<- cbind(LUSCs[c(x),], LUAD)
	PAAD<- cbind(PAADs[c(x),], LUSC)
	PCPG<- cbind(PCPGs[c(x),], PAAD)
	PRAD<- cbind(PRADs[c(x),], PCPG)
	STAD<- cbind(STADs[c(x),], PRAD)
	THCA<- cbind(THCAs[c(x),], STAD)
	UCEC<- cbind(UCECs[c(x),], THCA)
	FPPP<- cbind(FPPPs[c(x),], UCEC)
	COAD<- cbind(COADs[c(x),], FPPP)

	return(COAD)
}
	
#rank function
	rrfunction<- function(x)
{
	(x/abs(L))-(1/(2*abs(L)))
}

#inverted rank
invrrfunction<- function(x)
{
	abs(L)-x+1
}

#find Hknot
Hknot<- function(x)
{
	y<-t(x)
	y<-na.omit(y)
	-2*(sum(log(y)))
}

#function for ranking and Hknots

Hknots=NULL
Hknotsinv=NULL
testpos<- function(y)
{
	w <- microcut(y)
	colnames(w)<- c("COAD", "FPPP", "UCEC", "THCA", "STAD", "PRAD", "PCPG", "PAAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "ESCA", "CHOL", "CESC", "BRCA", "BLCA")
	#rank miRNA 
	rrw <- rrfunction(w)
	Hknotresult <- Hknot(rrw)
	return(Hknotresult)
}

testneg<- function(y)
{
	w <- microcut(y)
	colnames(w)<- c("COAD", "FPPP", "UCEC", "THCA", "STAD", "PRAD", "PCPG", "PAAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "ESCA", "CHOL", "CESC", "BRCA", "BLCA")
	invrrw<- rrfunction(invrrfunction(w))
	Hknotinvresult<- Hknot(invrrw)
}

#STOP
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
normalmirnanames<- read.csv(file="normal mirna names.csv")
names<- apply(normalmirnanames, 2, gsub, pat="-", replace=".")
tnames<- t(names)
colnames(tnames)=names[,1]

positiveHknot=NULL
for (i in colnames(tnames))
{
	k=NULL
	k<-testpos(i)
	positiveHknot<- cbind(positiveHknot, k)
	
}
colnames(positiveHknot)=colnames(tnames)

negativeHknot=NULL
for (i in colnames(tnames))
{
	k=NULL
	k<-testneg(i)
	negativeHknot<- cbind(negativeHknot, k)
	
}
colnames(negativeHknot)=colnames(tnames)

tailedH<-rbind(positiveHknot, negativeHknot)
rownames(tailedH)=c("positiveHknot", "negativeHknot")

#chi square find df
df=NULL
for (i in colnames(tnames))
{
	slice<- microcut(i)
	colnames(slice)<- c("COAD", "FPPP", "UCEC", "THCA", "STAD", "PRAD", "PCPG", "PAAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "ESCA", "CHOL", "CESC", "BRCA", "BLCA")
	tslice<- t(slice)
	NAslice<- na.omit(tslice)
	y<- nrow(NAslice)
	df<- cbind(df, y)
	}
colnames(df)=colnames(tnames)

#chi square tails
poschi<-1-pchisq(positiveHknot,(2*df))
negchi<-1-pchisq(negativeHknot,(2*df))
chi<- rbind(poschi, negchi)

#recscores

REC=NULL
for (i in colnames(chi))
{
	y=NULL
	y<- if(chi[2, i] < chi[1, i]) 
	{
		log10(2*(chi[2, i]))
	} 		else{
				-log10(2*(chi[1, i]))
			}
	REC<- cbind(REC, y)
}

colnames(REC)=colnames(tnames)

#Final results
results<- rbind(REC, chi, df, tailedH)
rownames(results)=c("REC", "positivechi", "negativechi", "no Na cancers", "positiveHknot", "negativeHknot")
tresults<- t(results)
#STOP
is.na(tresults) <- do.call(cbind,lapply(tresults, is.infinite))
nonaresults<- na.omit(tresults)
Resultrank <- nonaresults[order(-nonaresults[,c("REC")]),, drop=FALSE]
hist(Resultrank[,1])
lines(density(Resultrank[,1]), col="red")
