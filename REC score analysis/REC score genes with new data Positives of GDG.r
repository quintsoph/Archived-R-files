#REC score genes with new data POSITIVES of GDG
#read in files
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/BLCA")
BLCA<- read.table("POS coef of BLCA.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/BRCA")
BRCA<- read.table("POS coef of BRCA.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/COAD")
COAD<- read.table("POS coef of COAD.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/HNSC")
HNSC<- read.table("POS coef of HNSC.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/KICH")
KICH<- read.table("POS coef of KICH.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/KIRC")
KIRC<- read.table("POS coef of KIRC.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/KIRP")
KIRP<- read.table("POS coef of KIRP.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/LIHC")
LIHC<- read.table("POS coef of LIHC.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/LUAD")
LUAD<- read.table("POS coef of LUAD.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/LUSC")
LUSC<- read.table("POS coef of LUSC.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/PRAD")
PRAD<- read.table("POS coef of PRAD.txt", header=TRUE, row.names=1)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/THCA")
THCA<- read.table("POS coef of THCA.txt", header=TRUE, row.names=1)

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

#COAD
COADorder <- COAD[order(-COAD[,c("generatio_pval")]),, drop=FALSE]
tCOAD<- t(COADorder)
COADp<- tCOAD[-c(1),]
COADNA<- COADp[!is.na(COADp)]
COADNA1<- data.frame(COADNA)
COADNA<- data.frame((rank(COADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(COADNA)=rownames(COADNA1)
COADr<- nrow(COADNA)
L<- cbind(COADr, L)

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

#Gene cutting
genecut<- function(x)
{
	gene=NULL
	BLCA<- BLCANA[c(x),]
	BRCA<- cbind(BRCANA[c(x),], BLCA)
	COAD<- cbind(COADNA[c(x),], BRCA)
	HNSC<- cbind(HNSCNA[c(x),], COAD)
	KICH<- cbind(KICHNA[c(x),], HNSC)
	KIRC<- cbind(KIRCNA[c(x),], KICH)
	KIRP<- cbind(KIRPNA[c(x),], KIRC)
	LIHC<- cbind(LIHCNA[c(x),], KIRP)
	LUAD<- cbind(LUADNA[c(x),], LIHC)
	LUSC<- cbind(LUSCNA[c(x),], LUAD)
	PRAD<- cbind(PRADNA[c(x),], LUSC)
	THCA<- cbind(THCANA[c(x),], PRAD)
	return(THCA)
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

#HKNOTS
Hknots=NULL
Hknotsinv=NULL
testpos<- function(y)
{
	w <- genecut(y)
	colnames(w)<- c("THCA", "PRAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "COAD", "BRCA", "BLCA")
	#rank miRNA 
	rrw <- rrfunction(w)
	Hknotresult <- Hknot(rrw)
	return(Hknotresult)
}

testneg<- function(y)
{
	w <- genecut(y)
	colnames(w)<- c("THCA", "PRAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "COAD", "BRCA", "BLCA")
	invrrw<- rrfunction(invrrfunction(w))
	Hknotinvresult<- Hknot(invrrw)
	return(Hknotinvresult)
}

#create the name files
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis")
genes<- read.table("Allgenes.txt")
rownames(genes)=genes[,1]
positiveHknot=NULL
for (i in rownames(genes))
{
	k=NULL
	k<-testpos(i)
	positiveHknot<- cbind(positiveHknot, k)
	
}
colnames(positiveHknot)=rownames(genes)

negativeHknot=NULL
for (i in rownames(genes))
{
	k=NULL
	k<-testneg(i)
	negativeHknot<- cbind(negativeHknot, k)
	
}
colnames(negativeHknot)=rownames(genes)

tailedH<-rbind(positiveHknot, negativeHknot)
rownames(tailedH)=c("positiveHknot", "negativeHknot")

#chi square find df
df=NULL
for (i in rownames(genes))
{
	slice<- genecut(i)
	colnames(slice)<- c("THCA", "PRAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "COAD", "BRCA", "BLCA")
	tslice<- t(slice)
	NAslice<- na.omit(tslice)
	y<- nrow(NAslice)
	df<- cbind(df, y)
	}
colnames(df)=rownames(genes)

#chi square tails
poschi<-1-pchisq(positiveHknot,(2*df))
negchi<-1-pchisq(negativeHknot,(2*df))
chi<- rbind(poschi, negchi)

#recscores
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

colnames(REC)=rownames(genes)

#Final results
results<- rbind(REC, chi, df, tailedH)
rownames(results)=c("REC", "positivechi", "negativechi", "no Na cancers", "positiveHknot", "negativeHknot")
tresults<- t(results)
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG")
orderedREC<- tresults[order(tresults[,c("REC")]),, drop=FALSE]
write.table(orderedREC, "Gene REC Scores POS.txt", sep="\t", quote=FALSE)

#Find the ranks in each cancer

#Gene cutting
genecut<- function(x)
{
	gene=NULL
	BLCA<- BLCANA[c(x),]
	BRCA<- cbind(BRCANA[c(x),], BLCA)
	COAD<- cbind(COADNA[c(x),], BRCA)
	HNSC<- cbind(HNSCNA[c(x),], COAD)
	KICH<- cbind(KICHNA[c(x),], HNSC)
	KIRC<- cbind(KIRCNA[c(x),], KICH)
	KIRP<- cbind(KIRPNA[c(x),], KIRC)
	LIHC<- cbind(LIHCNA[c(x),], KIRP)
	LUAD<- cbind(LUADNA[c(x),], LIHC)
	LUSC<- cbind(LUSCNA[c(x),], LUAD)
	PRAD<- cbind(PRADNA[c(x),], LUSC)
	THCA<- cbind(THCANA[c(x),], PRAD)
	return(THCA)
}
 totalranks=NULL
 for (i in rownames(genes))
 {
	ranks=NULL
	ranks<- genecut(i)
	colnames(ranks)<- c("THCA", "PRAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "COAD", "BRCA", "BLCA")
	rownames(ranks)<- c(i)
	totalranks<- rbind(totalranks, ranks)
 }
 
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG")
write.table(totalranks, "ranks in each cancer of each gene FIX POS.txt", sep="\t", quote=FALSE)
