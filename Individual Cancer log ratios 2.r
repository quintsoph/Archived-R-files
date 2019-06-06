setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
AllCancers <-data.frame(read.table(file="Cancer Types Sorted"))
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/formated L1HS and mirna")
L1HS<- data.frame(read.csv("L1HS.csv"))

#L1HS cancer
L1HScancer<- L1HS[,c(1,4)]
colnames(L1HScancer)=c("id", "BaseMeanB Cancer")
rownames(L1HScancer)= L1HScancer[,1]
L1HScancerrevised<- gsub("\\..*","", rownames(L1HScancer))
L1HScancerF<- L1HScancer[,2]
L1HScancerF<- data.frame(L1HScancerF)
rownames(L1HScancerF)=L1HScancerrevised
L1HScancert<- t(L1HScancerF)

#L1HS normal
L1HSnormal<- L1HS[,c(1,3)]
colnames(L1HSnormal)=c("id", "BaseMeanA Normal")
L1HSnormalF<- L1HSnormal[,2]
L1HSnormalF<- data.frame(L1HSnormalF)
rownames(L1HSnormalF)=L1HScancerrevised
L1HSnormalt<- t(L1HSnormalF)

#Mirna cancer
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/MiRNA files")
mirnafiles <- dir(pattern = "*.miRna", full.names = TRUE)
lst <- lapply(mirnafiles, read.table, header=TRUE, sep='')
mirna<- sapply(lst, '[[',3)
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/formated L1HS and mirna")
names<- read.table("miRNA names.txt")
tnames=t(names)
rownames(mirna)=tnames
tmirna=t(mirna)
rownames(tmirna)=L1HScancerrevised
cancermirna<- t(tmirna)

#mirnanormal
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/Log2 ratio")
normalmirna<- read.table(file="normalMiRNA data ALL.ReadCountPerMillion", header=TRUE)
framenormalmirna<- data.matrix(normalmirna)
framecancermirna<- data.matrix(mirna)
framenormalmirna<- framenormalmirna[,-c(1)]


Allraw <- subset(AllCancers, V1=="BLCA")
AllcancerL1HS<- L1HScancert[,c(Allraw[,2])]
AllcancerL1HS<- t(AllcancerL1HS)
AllnormalL1HS<- L1HSnormalt[,c(Allraw[,2])]
AllnormalL1HS<- t(AllnormalL1HS)

#BRCA MiRNA
Allcancermirna<- cancermirna[,c(Allraw[,2])]
Allcancermirna<- t(Allcancermirna)
Allnormalmirna<- framenormalmirna[,c(Allraw[,2])]
Allnormalmirna<- t(Allnormalmirna)

#All log2ratios
#find the dividen and log
DividensofL1HSAll<- AllcancerL1HS/AllnormalL1HS
logofL1HSAll <- log2(DividensofL1HSAll)
logofL1HSAll <- t(logofL1HSAll)
#find the dividensofmirna and log
DividensofmirnaAll<- Allcancermirna/Allnormalmirna
logofmirnaAll<- log2(DividensofmirnaAll)
#Change inf to NAs.
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.infinite))
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.nan))

 #L1HS and mirna
ratiopvaluesAll<- function(x)
{
	idea <- lm(logofL1HSAll[,c(1)] ~ logofmirnaAll[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaAll <- !apply (is.na(logofmirnaAll), 2, all)
framefalsenaAll<- data.frame(falsenaAll)
NONAAll <- subset(framefalsenaAll, falsenaAll==TRUE)
sublogmirnaAll<- logofmirnaAll[,c(rownames(NONAAll))]

#iterate
presultsAll=NULL
for (i in colnames(sublogmirnaAll))
{
	q = NULL
	
		q <-ratiopvaluesAll(i)

	presultsAll <- cbind(presultsAll, q) 
}
 #organize pvalues
 colnames(presultsAll)=colnames(sublogmirnaAll)
tpresultsAll<- t(presultsAll)
tpresultsAllBLCA<- tpresults
#organize rvalues
rvaluesAll=NULL
for (i in colnames(sublogmirnaAll))
{
	r = NULL
	
		r <-cor(logofL1HSAll[,c(1)], logofmirnaAll[,i])

	rvaluesAll <- cbind(rvaluesAll, r) 
}

trvaluesAll<- t(rvaluesAll)
colnames(trvaluesAll)=c("data")
rownames(trvaluesAll)=colnames(sublogmirnaAll)
resultsAll<- cbind(tpresultsAll, trvaluesAll)
colnames(resultsAll)=c("pvalues", "rvalues")
sortedbypvaluesAll<- resultsAll[order(-resultsAll[,c("pvalues")]),, drop=FALSE]

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
 write.table(sortedbypvaluesAll, file="ESCA pvalues and rvalues of log2ratios.txt", sep=" ")

 #combine them all
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
BLCA<- read.table(file="BLCA pvalues and rvalues of log2ratios.txt")
BRCA<- read.table(file="BRCA pvalues and rvalues of log2ratios.txt")
CESC<- read.table(file="CESC pvalues and rvalues of log2ratios.txt")
CHOL<- read.table(file="CHOL pvalues and rvalues of log2ratios.txt")
ESCA<- read.table(file="ESCA pvalues and rvalues of log2ratios.txt")
HNSC<- read.table(file="HNSC pvalues and rvalues of log2ratios.txt")
KICH<- read.table(file="KICH pvalues and rvalues of log2ratios.txt")
KIRC<- read.table(file="KIRC pvalues and rvalues of log2ratios.txt")
KIRP<- read.table(file="KIRP pvalues and rvalues of log2ratios.txt")
LIHC<- read.table(file="LIHC pvalues and rvalues of log2ratios.txt")
LUAD<- read.table(file="LUAD pvalues and rvalues of log2ratios.txt")
LUSC<- read.table(file="LUSC pvalues and rvalues of log2ratios.txt")
PAAD<- read.table(file="PAAD pvalues and rvalues of log2ratios.txt")
PCPG<- read.table(file="PCPG pvalues and rvalues of log2ratios.txt")
PRAD<- read.table(file="PRAD pvalues and rvalues of log2ratios.txt")
STAD<- read.table(file="STAD pvalues and rvalues of log2ratios.txt")
THCA<- read.table(file="THCA pvalues and rvalues of log2ratios.txt")
THYM<- read.table(file="THYM pvalues and rvalues of log2ratios.txt")
UCEC<- read.table(file="UCEC pvalues and rvalues of log2ratios.txt")


rank()
rank(x, na.last = TRUE, ties.method = c("average", "first", "random", "max", "min"))


#p.adjust(KICH[,c(1)])
KICHqnone<- p.adjust(KICH[,c(1)], method="none")
KICHqBH<- p.adjust(KICH[,c(1)], method="BH")
KICHqbon<- p.adjust(KICH[,c(1)], method="bonferroni")
KICHq<- cbind(KICH, KICHqnone, KICHqBH, KICHqbon)
hsa-mir-153-1

BLCAqnone<- p.adjust(BLCA[,c(1)], method="none")
BLCAqBH<- p.adjust(BLCA[,c(1)], method="BH")
BLCAqbon<- p.adjust(BLCA[,c(1)], method="bonferroni")
BLCAq<- cbind(BLCA, BLCAqnone, BLCAqBH, BLCAqbon)
hsa-mir-651

BRCAqnone<- p.adjust(BRCA[,c(1)], method="none")
BRCAqBH<- p.adjust(BRCA[,c(1)], method="BH")
BRCAqbon<- p.adjust(BRCA[,c(1)], method="bonferroni")
BRCAq<- cbind(BRCA, BRCAqnone, BRCAqBH, BRCAqbon)
hsa-mir-375

HNSCqnone<- p.adjust(HNSC[,c(1)], method="none")
HNSCqBH<- p.adjust(HNSC[,c(1)], method="BH")
HNSCqbon<- p.adjust(HNSC[,c(1)], method="bonferroni")
HNSCq<- cbind(HNSC, HNSCqnone, HNSCqBH, HNSCqbon)
hsa-mir-767

KIRCqnone<- p.adjust(KIRC[,c(1)], method="none")
KIRCqBH<- p.adjust(KIRC[,c(1)], method="BH")
KIRCqbon<- p.adjust(KIRC[,c(1)], method="bonferroni")
KIRCq<- cbind(KIRC, KIRCqnone, KIRCqBH, KIRCqbon)
hsa-mir-628

KIRPqnone<- p.adjust(KIRP[,c(1)], method="none")
KIRPqBH<- p.adjust(KIRP[,c(1)], method="BH")
KIRPqbon<- p.adjust(KIRP[,c(1)], method="bonferroni")
KIRPq<- cbind(KIRP, KIRPqnone, KIRPqBH, KIRPqbon)
hsa-mir-9-1

LIHCqnone<- p.adjust(LIHC[,c(1)], method="none")
LIHCqBH<- p.adjust(LIHC[,c(1)], method="BH")
LIHCqbon<- p.adjust(LIHC[,c(1)], method="bonferroni")
LIHCq<- cbind(LIHC, LIHCqnone, LIHCqBH, LIHCqbon)
hsa-mir-599

LUSCqnone<- p.adjust(LUSC[,c(1)], method="none")
LUSCqBH<- p.adjust(LUSC[,c(1)], method="BH")
LUSCqbon<- p.adjust(LUSC[,c(1)], method="bonferroni")
LUSCq<- cbind(LUSC, LUSCqnone, LUSCqBH, LUSCqbon)
hsa-let-7e

PRADqnone<- p.adjust(PRAD[,c(1)], method="none")
PRADqBH<- p.adjust(PRAD[,c(1)], method="BH")
PRADqbon<- p.adjust(PRAD[,c(1)], method="bonferroni")
PRADq<- cbind(PRAD, PRADqnone, PRADqBH, PRADqbon)
hsa-mir-204

THCAqnone<- p.adjust(THCA[,c(1)], method="none")
THCAqBH<- p.adjust(THCA[,c(1)], method="BH")
THCAqbon<- p.adjust(THCA[,c(1)], method="bonferroni")
THCAq<- cbind(THCA, THCAqnone, THCAqBH, THCAqbon)
hsa-mir-30b

UCECqnone<- p.adjust(UCEC[,c(1)], method="none")
UCECqBH<- p.adjust(UCEC[,c(1)], method="BH")
UCECqbon<- p.adjust(UCEC[,c(1)], method="bonferroni")
UCECq<- cbind(UCEC, UCECqnone, UCECqBH, UCECqbon)
hsa-mir-380

