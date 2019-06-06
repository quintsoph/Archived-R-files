#Read in files
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
input<- data.frame(read.table("geneQvals.txt"))
input<- input[-c(1),]
input2<- input[,-c(4)]
#seperate by cancer type
BLCA<- subset(input2, input2[,2]=="BLCA")
rownames(BLCA)=BLCA[,1]
BLCA<- data.frame(BLCA[,-c(1,2)])
BRCA<- subset(input2, input2[,2]=="BRCA")
rownames(BRCA)=BRCA[,1]
BRCA<- data.frame(BRCA[,-c(1,2)])
CHOL<- subset(input2, input2[,2]=="CHOL")
rownames(CHOL)=CHOL[,1]
CHOL<- data.frame(CHOL[,-c(1,2)])
COAD<- subset(input2, input2[,2]=="COAD")
rownames(COAD)=COAD[,1]
COAD<- data.frame(COAD[,-c(1,2)])
ESCA<- subset(input2, input2[,2]=="ESCA")
rownames(ESCA)=ESCA[,1]
ESCA<- data.frame(ESCA[,-c(1,2)])
HNSC<- subset(input2, input2[,2]=="HNSC")
rownames(HNSC)=HNSC[,1]
HNSC<- data.frame(HNSC[,-c(1,2)])
KIRC<- subset(input2, input2[,2]=="KIRC")
rownames(KIRC)=KIRC[,1]
KIRC<- data.frame(KIRC[,-c(1,2)])
KIRP<- subset(input2, input2[,2]=="KIRP")
rownames(KIRP)=KIRP[,1]
KIRP<- data.frame(KIRP[,-c(1,2)])
KICH<- subset(input2, input2[,2]=="KICH")
rownames(KICH)=KICH[,1]
KICH<- data.frame(KICH[,-c(1,2)])
LIHC<- subset(input2, input2[,2]=="LIHC")
rownames(LIHC)=LIHC[,1]
LIHC<- data.frame(LIHC[,-c(1,2)])
LUAD<- subset(input2, input2[,2]=="LUAD")
rownames(LUAD)=LUAD[,1]
LUAD<- data.frame(LUAD[,-c(1,2)])
PRAD<- subset(input2, input2[,2]=="PRAD")
rownames(PRAD)=PRAD[,1]
PRAD<- data.frame(PRAD[,-c(1,2)])
THCA<- subset(input2, input2[,2]=="THCA")
rownames(THCA)=THCA[,1]
THCA<- data.frame(THCA[,-c(1,2)])
#create ranks
L=NULL
#BLCA
BLCAorder <- BLCA
tBLCA<- t(BLCAorder)
BLCAp<- tBLCA[-c(2),]
BLCANA<- BLCAp[!is.na(BLCAp)]
BLCANA1<- data.frame(BLCANA)
BLCANA<- data.frame((rank(BLCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(BLCANA)=rownames(BLCANA1)
BLCAr<- nrow(BLCANA)
L<- cbind(BLCAr, L)

#BRCA
BRCAorder <- BRCA
tBRCA<- t(BRCAorder)
BRCAp<- tBRCA[-c(2),]
BRCANA<- BRCAp[!is.na(BRCAp)]
BRCANA1<- data.frame(BRCANA)
BRCANA<- data.frame((rank(BRCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(BRCANA)=rownames(BRCANA1)
BRCAr<- nrow(BRCANA)
L<- cbind(BRCAr, L)

#CHOL
CHOLorder <- CHOL
tCHOL<- t(CHOLorder)
CHOLp<- tCHOL[-c(2),]
CHOLNA<- CHOLp[!is.na(CHOLp)]
CHOLNA1<- data.frame(CHOLNA)
CHOLNA<- data.frame((rank(CHOLNA1, na.last = TRUE, ties.method=c("min"))))
rownames(CHOLNA)=rownames(CHOLNA1)
CHOLr<- nrow(CHOLNA)
L<- cbind(CHOLr, L)

#COAD
COADorder <- COAD
tCOAD<- t(COADorder)
COADp<- tCOAD[-c(2),]
COADNA<- COADp[!is.na(COADp)]
COADNA1<- data.frame(COADNA)
COADNA<- data.frame((rank(COADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(COADNA)=rownames(COADNA1)
COADr<- nrow(COADNA)
L<- cbind(COADr, L)

#ESCA
ESCAorder <- ESCA
tESCA<- t(ESCAorder)
ESCAp<- tESCA[-c(2),]
ESCANA<- ESCAp[!is.na(ESCAp)]
ESCANA1<- data.frame(ESCANA)
ESCANA<- data.frame((rank(ESCANA1, na.last = TRUE, ties.method=c("min"))))
rownames(ESCANA)=rownames(ESCANA1)
ESCAr<- nrow(ESCANA)
L<- cbind(ESCAr, L)

#HNSC
HNSCorder <- HNSC
tHNSC<- t(HNSCorder)
HNSCp<- tHNSC[-c(2),]
HNSCNA<- HNSCp[!is.na(HNSCp)]
HNSCNA1<- data.frame(HNSCNA)
HNSCNA<- data.frame((rank(HNSCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(HNSCNA)=rownames(HNSCNA1)
HNSCr<- nrow(HNSCNA)
L<- cbind(HNSCr, L)

#KIRC
KIRCorder <- KIRC
tKIRC<- t(KIRCorder)
KIRCp<- tKIRC[-c(2),]
KIRCNA<- KIRCp[!is.na(KIRCp)]
KIRCNA1<- data.frame(KIRCNA)
KIRCNA<- data.frame((rank(KIRCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KIRCNA)=rownames(KIRCNA1)
KIRCr<- nrow(KIRCNA)
L<- cbind(KIRCr, L)

#KIRP
KIRPorder <- KIRP
tKIRP<- t(KIRPorder)
KIRPp<- tKIRP[-c(2),]
KIRPNA<- KIRPp[!is.na(KIRPp)]
KIRPNA1<- data.frame(KIRPNA)
KIRPNA<- data.frame((rank(KIRPNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KIRPNA)=rownames(KIRPNA1)
KIRPr<- nrow(KIRPNA)
L<- cbind(KIRPr, L)

#KICH
KICHorder <- KICH
tKICH<- t(KICHorder)
KICHp<- tKICH[-c(2),]
KICHNA<- KICHp[!is.na(KICHp)]
KICHNA1<- data.frame(KICHNA)
KICHNA<- data.frame((rank(KICHNA1, na.last = TRUE, ties.method=c("min"))))
rownames(KICHNA)=rownames(KICHNA1)
KICHr<- nrow(KICHNA)
L<- cbind(KICHr, L)

#LIHC
LIHCorder <- LIHC
tLIHC<- t(LIHCorder)
LIHCp<- tLIHC[-c(2),]
LIHCNA<- LIHCp[!is.na(LIHCp)]
LIHCNA1<- data.frame(LIHCNA)
LIHCNA<- data.frame((rank(LIHCNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LIHCNA)=rownames(LIHCNA1)
LIHCr<- nrow(LIHCNA)
L<- cbind(LIHCr, L)

#LUAD
LUADorder <- LUAD
tLUAD<- t(LUADorder)
LUADp<- tLUAD[-c(2),]
LUADNA<- LUADp[!is.na(LUADp)]
LUADNA1<- data.frame(LUADNA)
LUADNA<- data.frame((rank(LUADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(LUADNA)=rownames(LUADNA1)
LUADr<- nrow(LUADNA)
L<- cbind(LUADr, L)

#PRAD
PRADorder <- PRAD
tPRAD<- t(PRADorder)
PRADp<- tPRAD[-c(2),]
PRADNA<- PRADp[!is.na(PRADp)]
PRADNA1<- data.frame(PRADNA)
PRADNA<- data.frame((rank(PRADNA1, na.last = TRUE, ties.method=c("min"))))
rownames(PRADNA)=rownames(PRADNA1)
PRADr<- nrow(PRADNA)
L<- cbind(PRADr, L)

#THCA
THCAorder <- THCA
tTHCA<- t(THCAorder)
THCAp<- tTHCA[-c(2),]
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
	CHOL<- cbind(CHOLNA[c(x),], BRCA)
	COAD<- cbind(COADNA[c(x),], CHOL)
	ESCA<- cbind(ESCANA[c(x),], COAD)
	HNSC<- cbind(HNSCNA[c(x),], ESCA)
	KIRC<- cbind(KIRCNA[c(x),], HNSC)
	KIRP<- cbind(KIRPNA[c(x),], KIRC)
	KICH<- cbind(KICHNA[c(x),], KIRP)
	LIHC<- cbind(LIHCNA[c(x),], KICH)
	LUAD<- cbind(LUADNA[c(x),], LIHC)
	PRAD<- cbind(PRADNA[c(x),], LUAD)
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
	colnames(w)<- c("THCA", "PRAD", "LUAD", "LIHC", "KICH", "KIRP", "KIRC", "HNSC", "ESCA", "COAD", "CHOL", "BRCA", "BLCA")
	#rank miRNA 
	rrw <- rrfunction(w)
	Hknotresult <- Hknot(rrw)
	return(Hknotresult)
}

testneg<- function(y)
{
	w <- genecut(y)
	colnames(w)<- c("THCA", "PRAD", "LUAD", "LIHC", "KICH", "KIRP", "KIRC", "HNSC", "ESCA", "COAD", "CHOL", "BRCA", "BLCA")
	invrrw<- rrfunction(invrrfunction(w))
	Hknotinvresult<- Hknot(invrrw)
	return(Hknotinvresult)
}

#create the name files
tests<- data.frame(table(input[,1]))
genes<- data.frame(tests[,1])
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
	colnames(slice)<- c("THCA", "PRAD", "LUAD", "LIHC", "KICH", "KIRP", "KIRC", "HNSC", "ESCA", "COAD", "CHOL", "BRCA", "BLCA")
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

colnames(REC)=rownames(genes)

#Final results
results<- rbind(REC, chi, df, tailedH)
rownames(results)=c("REC", "positivechi", "negativechi", "no Na cancers", "positiveHknot", "negativeHknot")
tresults<- t(results)
setwd("E:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
write.table(tresults, "Gene REC Scores.txt", sep="\t", quote=FALSE)

#Find the ranks in each cancer

#Gene cutting
genecut<- function(x)
{
	gene=NULL
	BLCA<- BLCANA[c(x),]
	BRCA<- cbind(BRCANA[c(x),], BLCA)
	CHOL<- cbind(CHOLNA[c(x),], BRCA)
	COAD<- cbind(COADNA[c(x),], CHOL)
	ESCA<- cbind(ESCANA[c(x),], COAD)
	HNSC<- cbind(HNSCNA[c(x),], ESCA)
	KIRC<- cbind(KIRCNA[c(x),], HNSC)
	KIRP<- cbind(KIRPNA[c(x),], KIRC)
	KICH<- cbind(KICHNA[c(x),], KIRP)
	LIHC<- cbind(LIHCNA[c(x),], KICH)
	LUAD<- cbind(LUADNA[c(x),], LIHC)
	PRAD<- cbind(PRADNA[c(x),], LUAD)
	THCA<- cbind(THCANA[c(x),], PRAD)
	return(THCA)
}
 
 totalranks=NULL
 for (i in rownames(genes))
 {
	ranks=NULL
	ranks<- genecut(i)
	colnames(ranks)<- c("THCA", "PRAD", "LUAD", "LIHC", "KICH", "KIRP", "KIRC", "HNSC", "ESCA", "COAD", "CHOL", "BRCA", "BLCA")
	rownames(ranks)<- c(i)
	totalranks<- rbind(totalranks, ranks)
 }
 
 setwd("E:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
write.table(totalranks, "ranks in each cancer of each gene.txt", sep="\t", quote=FALSE)

 