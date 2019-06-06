setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
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

#Mirna normal
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/Log2 ratio")
normalmirna<- read.table(file="normalMiRNA data ALL.ReadCountPerMillion", header=TRUE)
framenormalmirna<- data.matrix(normalmirna)
framenormalmirna<- framenormalmirna[,-c(1)]
rownames(framenormalmirna)=rownames(cancermirna)

#BLCA L1HS
Bladderraw <- subset(AllCancers, V1=="BLCA")
subset(L1HScancert, colnames==Bladderraw[,2])
BladdercancerL1HS<- L1HScancert[,c(Bladderraw[,2])]
BladdercancerL1HS<- t(BladdercancerL1HS)
BladdernormalL1HS<- L1HSnormalt[,c(Bladderraw[,2])]
BladdernormalL1HS<- t(BladdernormalL1HS)

#BLCA MiRNA
Bladdercancermirna<- cancermirna[,c(Bladderraw[,2])]
Bladdercancermirna<- t(Bladdercancermirna)
Bladdernormalmirna<- framenormalmirna[,c(Bladderraw[,2])]
Bladdernormalmirna<- t(Bladdernormalmirna)

#BLCA log2ratios
#find the dividen and log
DividensofL1HSBladder<- BladdercancerL1HS/BladdernormalL1HS
logofL1HSBladder <- log2(DividensofL1HSBladder)
logofL1HSBladder<- t(logofL1HSBladder)
#find the dividensofmirna and log
DividensofmirnaBladder<- Bladdercancermirna/Bladdernormalmirna
logofmirnaBladder<- log2(DividensofmirnaBladder)
#Change inf to NAs.
 is.na(logofmirnaBladder) <- do.call(cbind,lapply(logofmirnaBladder, is.infinite))
 is.na(logofmirnaBladder) <- do.call(cbind,lapply(logofmirnaBladder, is.nan))

 #L1HS and mirna bladder pvalues
ratiopvaluesBLCA<- function(x)
{
	idea <- lm(logofL1HSBladder[,c(1)] ~ logofmirnaBladder[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaBLCA <- !apply (is.na(logofmirnaBladder), 2, all)
framefalsenaBLCA<- data.frame(falsenaBLCA)
NONABLCA <- subset(framefalsenaBLCA, falsenaBLCA==TRUE)
sublogmirnaBLCA<- logofmirnaBladder[,c(rownames(NONABLCA))]

#iterate
presultsBLCA=NULL
for (i in colnames(sublogmirnaBLCA))
{
	q = NULL
	
		q <-ratiopvaluesBLCA(i)

	presultsBLCA <- cbind(presultsBLCA, q) 
}
 #organize pvalues
 colnames(presultsBLCA)=colnames(sublogmirnaBLCA)
tpresultsBLCA<- t(presultsBLCA)

#organize rvalues
rvaluesBLCA=NULL
for (i in colnames(sublogmirnaBLCA))
{
	r = NULL
	
		r <-cor(logofL1HSBladder[,c(1)], logofmirnaBladder[,i])

	rvaluesBLCA <- cbind(rvaluesBLCA, r) 
}

trvaluesBLCA<- t(rvaluesBLCA)
colnames(trvaluesBLCA)=c("data")
rownames(trvaluesBLCA)=colnames(sublogmirnaBLCA)
resultsBLCA<- cbind(tpresultsBLCA, trvaluesBLCA)
colnames(resultsBLCA)=c("pvalues", "rvalues")
sortedbypvaluesBLCA<- resultsBLCA[order(-resultsBLCA[,c("pvalues")]),, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
write.table(sortedbypvaluesBLCA, file="BLCA pvalues and rvalues of log2ratios.txt", sep=" ")




 
 
 BRCAraw <- subset(AllCancers, V1=="BRCA")
BRCAcancerL1HS<- L1HScancert[,c(BRCAraw[,2])]
BRCAcancerL1HS<- t(BRCAcancerL1HS)
BRCAnormalL1HS<- L1HSnormalt[,c(BRCAraw[,2])]
BRCAnormalL1HS<- t(BRCAnormalL1HS)

#BRCA MiRNA
BRCAcancermirna<- cancermirna[,c(BRCAraw[,2])]
BRCAcancermirna<- t(BRCAcancermirna)
BRCAnormalmirna<- framenormalmirna[,c(BRCAraw[,2])]
BRCAnormalmirna<- t(BRCAnormalmirna)

#BRCA log2ratios
#find the dividen and log
DividensofL1HSBRCA<- BRCAcancerL1HS/BRCAnormalL1HS
logofL1HSBRCA <- log2(DividensofL1HSBRCA)
logofL1HSBRCA<- t(logofL1HSBRCA)
#find the dividensofmirna and log
DividensofmirnaBRCA<- BRCAcancermirna/BRCAnormalmirna
logofmirnaBRCA<- log2(DividensofmirnaBRCA)
#Change inf to NAs.
 is.na(logofmirnaBRCA) <- do.call(cbind,lapply(logofmirnaBRCA, is.infinite))
 is.na(logofmirnaBRCA) <- do.call(cbind,lapply(logofmirnaBRCA, is.nan))

 #L1HS and mirna bladder pvalues
ratiopvaluesBRCA<- function(x)
{
	idea <- lm(logofL1HSBRCA[,c(1)] ~ logofmirnaBRCA[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaBRCA <- !apply (is.na(logofmirnaBRCA), 2, all)
framefalsenaBRCA<- data.frame(falsenaBRCA)
NONABRCA <- subset(framefalsenaBRCA, falsenaBRCA==TRUE)
sublogmirnaBRCA<- logofmirnaBRCA[,c(rownames(NONABRCA))]

#iterate
presultsBRCA=NULL
for (i in colnames(sublogmirnaBRCA))
{
	q = NULL
	
		q <-ratiopvaluesBRCA(i)

	presultsBRCA <- cbind(presultsBRCA, q) 
}
 #organize pvalues
 colnames(presultsBRCA)=colnames(sublogmirnaBRCA)
tpresultsBRCA<- t(presultsBRCA)

#organize rvalues
rvaluesBRCA=NULL
for (i in colnames(sublogmirnaBRCA))
{
	r = NULL
	
		r <-cor(logofL1HSBRCA[,c(1)], logofmirnaBRCA[,i])

	rvaluesBRCA <- cbind(rvaluesBRCA, r) 
}

trvaluesBRCA<- t(rvaluesBRCA)
colnames(trvaluesBRCA)=c("data")
rownames(trvaluesBRCA)=colnames(sublogmirnaBRCA)
resultsBRCA<- cbind(tpresultsBRCA, trvaluesBRCA)
colnames(resultsBRCA)=c("pvalues", "rvalues")
sortedbypvaluesBRCA<- resultsBRCA[order(-resultsBRCA[,c("pvalues")]),, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
 write.table(sortedbypvaluesBRCA, file="BRCA pvalues and rvalues of log2ratios.txt", sep=" ")






#All
#All
#All
#All L1HS
Allraw <- subset(AllCancers, V1=="ESCA")
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

 #L1HS and mirna bladder pvalues
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
BLCA<- data.frame(read.table(file="BLCA pvalues and rvalues of log2ratios.txt"))
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2ratios.txt"))
CESC<- data.frame(read.table(file="CESC pvalues and rvalues of log2ratios.txt"))
CHOL<- data.frame(read.table(file="CHOL pvalues and rvalues of log2ratios.txt"))
ESCA<- data.frame(read.table(file="ESCA pvalues and rvalues of log2ratios.txt"))
HNSC<- data.frame(read.table(file="HNSC pvalues and rvalues of log2ratios.txt"))
KICH<- data.frame(read.table(file="KICH pvalues and rvalues of log2ratios.txt"))
KIRC<- data.frame(read.table(file="KIRC pvalues and rvalues of log2ratios.txt"))
KIRP<- data.frame(read.table(file="KIRP pvalues and rvalues of log2ratios.txt"))
LIHC<- data.frame(read.table(file="LIHC pvalues and rvalues of log2ratios.txt"))
LUAD<- data.frame(read.table(file="LUAD pvalues and rvalues of log2ratios.txt"))
LUSC<- data.frame(read.table(file="LUSC pvalues and rvalues of log2ratios.txt"))
PAAD<- data.frame(read.table(file="PAAD pvalues and rvalues of log2ratios.txt"))
PCPG<- data.frame(read.table(file="PCPG pvalues and rvalues of log2ratios.txt"))
PRAD<- data.frame(read.table(file="PRAD pvalues and rvalues of log2ratios.txt"))
STAD<- data.frame(read.table(file="STAD pvalues and rvalues of log2ratios.txt"))
THCA<- data.frame(read.table(file="THCA pvalues and rvalues of log2ratios.txt"))
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2ratios.txt"))
UCEC<- data.frame(read.table(file="UCEC pvalues and rvalues of log2ratios.txt"))
 
 
 #REC Scores
 my.scores<- data.frame()
 
 BLCArank<- rank(BLCA, na.last=TRUE, ties.method=c("min"))
 
 
 
 
 
 
 
 
 
