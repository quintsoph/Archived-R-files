#REC score genes with new data NEGATIVES of GDG
#read in files
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/BLCA")
BLCA<- read.table("POS coef of BLCA.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/BRCA")
BRCA<- read.table("POS coef of BRCA.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/HNSC")
HNSC<- read.table("POS coef of HNSC.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/KICH")
KICH<- read.table("POS coef of KICH.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/KIRC")
KIRC<- read.table("POS coef of KIRC.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/KIRP")
KIRP<- read.table("POS coef of KIRP.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/LIHC")
LIHC<- read.table("POS coef of LIHC.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/LUAD")
LUAD<- read.table("POS coef of LUAD.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/LUSC")
LUSC<- read.table("POS coef of LUSC.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/PRAD")
PRAD<- read.table("POS coef of PRAD.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/THCA")
THCA<- read.table("POS coef of THCA.txt", header=TRUE, row.names=1)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2/UCEC")
UCEC<- read.table("POS coef of UCEC.txt", header=TRUE, row.names=1)
#create ranks
L=NULL
#BLCA
BLCAorder <- BLCA[order(-BLCA[,c("generatio_pval")]),, drop=FALSE]
tBLCA<- t(BLCAorder)
BLCAp<- tBLCA[-c(1),]
BLCANA<- BLCAp[!is.na(BLCAp)]
BLCANA1<- data.frame(BLCANA)
BLCANA<- data.frame((rank(BLCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(BLCANA)=rownames(BLCANA1)
BLCAr<- nrow(BLCANA)
L<- cbind(BLCAr, L)

#BRCA
BRCAorder <- BRCA[order(-BRCA[,c("generatio_pval")]),, drop=FALSE]
tBRCA<- t(BRCAorder)
BRCAp<- tBRCA[-c(1),]
BRCANA<- BRCAp[!is.na(BRCAp)]
BRCANA1<- data.frame(BRCANA)
BRCANA<- data.frame((rank(BRCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(BRCANA)=rownames(BRCANA1)
BRCAr<- nrow(BRCANA)
L<- cbind(BRCAr, L)

#HNSC
HNSCorder <- HNSC[order(-HNSC[,c("generatio_pval")]),, drop=FALSE]
tHNSC<- t(HNSCorder)
HNSCp<- tHNSC[-c(1),]
HNSCNA<- HNSCp[!is.na(HNSCp)]
HNSCNA1<- data.frame(HNSCNA)
HNSCNA<- data.frame((rank(HNSCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(HNSCNA)=rownames(HNSCNA1)
HNSCr<- nrow(HNSCNA)
L<- cbind(HNSCr, L)

#KICH
KICHorder <- KICH[order(-KICH[,c("generatio_pval")]),, drop=FALSE]
tKICH<- t(KICHorder)
KICHp<- tKICH[-c(1),]
KICHNA<- KICHp[!is.na(KICHp)]
KICHNA1<- data.frame(KICHNA)
KICHNA<- data.frame((rank(KICHNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KICHNA)=rownames(KICHNA1)
KICHr<- nrow(KICHNA)
L<- cbind(KICHr, L)

#KIRC
KIRCorder <- KIRC[order(-KIRC[,c("generatio_pval")]),, drop=FALSE]
tKIRC<- t(KIRCorder)
KIRCp<- tKIRC[-c(1),]
KIRCNA<- KIRCp[!is.na(KIRCp)]
KIRCNA1<- data.frame(KIRCNA)
KIRCNA<- data.frame((rank(KIRCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KIRCNA)=rownames(KIRCNA1)
KIRCr<- nrow(KIRCNA)
L<- cbind(KIRCr, L)

#KIRP
KIRPorder <- KIRP[order(-KIRP[,c("generatio_pval")]),, drop=FALSE]
tKIRP<- t(KIRPorder)
KIRPp<- tKIRP[-c(1),]
KIRPNA<- KIRPp[!is.na(KIRPp)]
KIRPNA1<- data.frame(KIRPNA)
KIRPNA<- data.frame((rank(KIRPNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KIRPNA)=rownames(KIRPNA1)
KIRPr<- nrow(KIRPNA)
L<- cbind(KIRPr, L)

#LIHC
LIHCorder <- LIHC[order(-LIHC[,c("generatio_pval")]),, drop=FALSE]
tLIHC<- t(LIHCorder)
LIHCp<- tLIHC[-c(1),]
LIHCNA<- LIHCp[!is.na(LIHCp)]
LIHCNA1<- data.frame(LIHCNA)
LIHCNA<- data.frame((rank(LIHCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LIHCNA)=rownames(LIHCNA1)
LIHCr<- nrow(LIHCNA)
L<- cbind(LIHCr, L)

#LUAD
LUADorder <- LUAD[order(-LUAD[,c("generatio_pval")]),, drop=FALSE]
tLUAD<- t(LUADorder)
LUADp<- tLUAD[-c(1),]
LUADNA<- LUADp[!is.na(LUADp)]
LUADNA1<- data.frame(LUADNA)
LUADNA<- data.frame((rank(LUADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LUADNA)=rownames(LUADNA1)
LUADr<- nrow(LUADNA)
L<- cbind(LUADr, L)

#LUSC
LUSCorder <- LUSC[order(-LUSC[,c("generatio_pval")]),, drop=FALSE]
tLUSC<- t(LUSCorder)
LUSCp<- tLUSC[-c(1),]
LUSCNA<- LUSCp[!is.na(LUSCp)]
LUSCNA1<- data.frame(LUSCNA)
LUSCNA<- data.frame((rank(LUSCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LUSCNA)=rownames(LUSCNA1)
LUSCr<- nrow(LUSCNA)
L<- cbind(LUSCr, L)

#PRAD
PRADorder <- PRAD[order(-PRAD[,c("generatio_pval")]),, drop=FALSE]
tPRAD<- t(PRADorder)
PRADp<- tPRAD[-c(1),]
PRADNA<- PRADp[!is.na(PRADp)]
PRADNA1<- data.frame(PRADNA)
PRADNA<- data.frame((rank(PRADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(PRADNA)=rownames(PRADNA1)
PRADr<- nrow(PRADNA)
L<- cbind(PRADr, L)

#THCA
THCAorder <- THCA[order(-THCA[,c("generatio_pval")]),, drop=FALSE]
tTHCA<- t(THCAorder)
THCAp<- tTHCA[-c(1),]
THCANA<- THCAp[!is.na(THCAp)]
THCANA1<- data.frame(THCANA)
THCANA<- data.frame((rank(THCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(THCANA)=rownames(THCANA1)
THCAr<- nrow(THCANA)
L<- cbind(THCAr, L)

#UCEC
UCECorder <- UCEC[order(-UCEC[,c("generatio_pval")]),, drop=FALSE]
tUCEC<- t(UCECorder)
UCECp<- tUCEC[-c(1),]
UCECNA<- UCECp[!is.na(UCECp)]
UCECNA1<- data.frame(UCECNA)
UCECNA<- data.frame((rank(UCECNA1, na.last = TRUE, ties.method=c("min"))))
rownames(UCECNA)=rownames(UCECNA1)
UCECr<- nrow(UCECNA)
L<- cbind(UCECr, L)

#Gene cutting
microcut<- function(x)
{
	microcut=NULL
	BLCA<- BLCANA[c(x),]
	BRCA<- cbind(BRCANA[c(x),], BLCA)
	HNSC<- cbind(HNSCNA[c(x),], BRCA)
	KICH<- cbind(KICHNA[c(x),], HNSC)
	KIRC<- cbind(KIRCNA[c(x),], KICH)
	KIRP<- cbind(KIRPNA[c(x),], KIRC)
	LIHC<- cbind(LIHCNA[c(x),], KIRP)
	LUAD<- cbind(LUADNA[c(x),], LIHC)
	LUSC<- cbind(LUSCNA[c(x),], LUAD)
	PRAD<- cbind(PRADNA[c(x),], LUSC)
	THCA<- cbind(THCANA[c(x),], PRAD)
	UCEC<- cbind(UCECNA[c(x),], THCA)
	return(UCEC)
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
Hknottest<- function(x)
{
	y<-t(x)
	y<-na.omit(y)
	-2*(sum(log(y)))
}

#function for ranking and Hknots

Hknots=NULL
Hknotsinv=NULL
#mirna 1
x="hsa-mir-651"
test<- function(y)
{
	w <- microcut(y)
	colnames(w)<- c("UCEC", "THCA", "PRAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "BRCA", "BLCA")
	#rank miRNA 
	rrw <- rrfunction(w)
	Hknotresult <- Hknottest(rrw)
	return(Hknotresult)
}

invrtest<- function(y)
{
	w <- microcut(y)
	colnames(w)<- c("UCEC", "THCA", "PRAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "BRCA", "BLCA")
	invrrw<- rrfunction(invrrfunction(w))
	Hknotinvresult<- Hknottest(invrrw)
}

#STOP
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
normalmirnanames<- read.csv(file="normal mirna names.csv")
names<- normalmirnanames
tnames<- t(names)
colnames(tnames)=names[,1]

Hknot=NULL
for (i in colnames(tnames))
{
	k=NULL
	k<-test(i)
	Hknot<- cbind(Hknot, k)
	
}
colnames(Hknot)=colnames(tnames)

invrHknot=NULL
for (i in colnames(tnames))
{
	k=NULL
	k<-invrtest(i)
	invrHknot<- cbind(invrHknot, k)
	
}
colnames(invrHknot)=colnames(tnames)

tailedH<-rbind(Hknot, invrHknot)
rownames(tailedH)=c("Hknot", "invrHknot")

#chi square find df
df=NULL
for (i in colnames(tnames))
{
	slice<- microcut(i)
	colnames(slice)<- c("UCEC", "THCA", "PRAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "BRCA", "BLCA")
	tslice<- t(slice)
	NAslice<- na.omit(tslice)
	y<- nrow(NAslice)
	df<- cbind(df, y)
	}
colnames(df)=colnames(tnames)

#chi square tails
pneg<-1-pchisq(Hknot,(2*df))
pplus<-1-pchisq(invrHknot,(2*df))
chi<- rbind(pneg, pplus)
rownames(chi)= c("pneg", "pplus")
#recscores
#THIS WAS WHERE I HAD TO SWITCH THE COLUMNS
REC=NULL
for (i in colnames(chi))
{
	y=NULL
	y<- if(chi[1, i] < chi[2, i]) 
	{
		log((2*(chi[1, i])), base=10)
	} 		else{
				-log((2*(chi[2, i])), base=10)
			}
	REC<- cbind(REC, y)
}

colnames(REC)=colnames(tnames)
REC2<- cbind(colnames(REC), t(REC))
finalREC<-data.frame(ifelse(chi[1,] == chi[2,], 0, REC[1,]))
finalREC2<- cbind(rownames(finalREC), finalREC)
#Final results
results<- rbind(t(finalREC), chi, df, tailedH)
rownames(results)=c("REC", "leftchi", "rightchi", "no Na cancers", "leftHknot", "rightHknot")
tresults<- t(results)
orderedREC<- tresults[order(tresults[,c("REC")]),, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.miRNA.log2")
write.table(orderedREC, "OrderedREC miRNA pos new data2.txt", quote=FALSE, sep="\t")
