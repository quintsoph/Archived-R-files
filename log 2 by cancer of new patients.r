setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
NewNames<- data.frame(read.table(file="SophiaPatCanTypes"), header=FALSE)
NewNames<- t(NewNames)
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
NewNames<- data.frame(read.table(file="SophiapatCanTypes"))
NewNames<- t(NewNames)
#test
LUSCraw<- subset(NewNames, NewNames[,2]=="LUSC")

#format cancer
newL1HScancer<- data.frame(read.table(file="L1HS.Cancer.BaseMeansB", header=TRUE))
newL1HScancer<- data.matrix(newL1HScancer)
tnewL1HScancer<- t(newL1HScancer)
newmirnacancer<- read.csv(file="Mirna.Cancer.csv")
newmirnacancer<- data.matrix(newmirnacancer)
cancermirnanames<- read.csv(file="cancer mirna names.csv")
rownames(newmirnacancer)=cancermirnanames[,1]
tnewmirnacancer<- t(newmirnacancer)
tnewmirnacancer<- tnewmirnacancer[-1,]
tnewL1HScancer<- tnewL1HScancer[-1,]
cancer<- cbind(tnewL1HScancer, tnewmirnacancer)
cancermirna<- cancer[,-1]

#format normal
newL1HSnormal<- data.frame(read.table(file="L1HS.Norma.BaseMeansA", header=TRUE))
newL1HSnormal<- data.matrix(newL1HSnormal)
tnewL1HSnormal<- t(newL1HSnormal)
newmirnanormal<- read.csv(file="Mirna.Normal.csv")
newmirnanormal<- data.matrix(newmirnanormal)
normalmirnanames<- read.csv(file="normal mirna names.csv")
rownames(newmirnanormal)=normalmirnanames[,1]
tnewmirnanormal<- t(newmirnanormal)
tnewmirnanormal<- tnewmirnanormal[-1,]
tnewL1HSnormal<- tnewL1HSnormal[-1,]
normal<- cbind(tnewL1HSnormal, tnewmirnanormal)
normalmirna<- normal[,-1]
dividenL1HS<- newL1HScancer/newL1HSnormal
foldedL1HS<- log2(dividenL1HS)
dividenmirna<- newmirnacancer/newmirnanormal
foldedmirna<- t(log2(dividenmirna))
foldedmirna<- foldedmirna[-1,]
#
#new stuff
L1HScancer<- data.frame((cancer[,1]))
L1HSnormal<- data.frame(normal[,1])
mirnacancer<- data.frame(cancermirna)
mirnanormal<- data.frame(normalmirna)

#TEST WITH ONE
BLCAraw<- subset(NewNames, NewNames[,2]=="BLCA")
dfBLCAraw<- data.frame(BLCAraw)

#cancer L1HS
matched<- match(dfBLCAraw[,1], rownames(L1HScancer))
BLCAL1HScancer <-L1HScancer[matched,]
BLCAL1HScancer<- data.frame(BLCAL1HScancer)
rownames(BLCAL1HScancer)=dfBLCAraw[,1]
logofL1HSBLCA<- log2(BLCAL1HScancer)

#cancermir
matchedmirna<- match(dfBLCAraw[,1], rownames(mirnacancer))
BLCAmirnacancer<- mirnacancer[matched,]
BLCAmirnacancer<- data.frame(BLCAmirnacancer)
rownames(BLCAmirnacancer)=dfBLCAraw[,1]
logofmirnaBLCA<- log2(BLCAmirnacancer)
is.na(logofmirnaBLCA) <- do.call(cbind,lapply(logofmirnaBLCA, is.infinite))

#pvalues
#Function to calculate p-values.
testpvalues<- function(x)
{
	print(x)
	value <-lm(logofL1HSBLCA[,c("BLCAL1HScancer")] ~ logofmirnaBLCA[,x], sep="")
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(logofmirnaBLCA))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(logofmirnaBLCA[i]))){
		p <-testpvalues(i)
	}
	else{
		print(paste("Skipping", i))
	}

	pvalues <- cbind(pvalues, p) 
}
pvalues<- data.frame(pvalues)
tpvalues<- t(pvalues)
skip=NULL
for (i in colnames(logofmirnaBLCA))
{
	skipped=NULL
	if(!all(is.na(logofmirnaBLCA[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}

skip<- t(skip)
rownames(tpvalues)= skip[,1]

#rvalues
rvalues=NULL
for (i in rownames(tpvalues))
{
	r = NULL
	
		r <-cor(logofmirnaBLCA[,i], logofL1HSBLCA[,c(1)])

	rvalues <- cbind(rvalues, r) 
}
trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)
colnames(randpvalues)=c("pvalues", "rvalues")
sortedbypvaluesBLCA<- randpvalues[order(-randpvalues[,c("pvalues")]),, drop=FALSE]
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results/log2 cancer by each")
write.table(sortedbypvaluesBLCA, file="BLCA pvalues and rvalues of log2 cancer reverse cor.txt", sep=" ")



#ALL LOG2 CANCER 
Allraw<- subset(NewNames, NewNames[,2]=="UCEC")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]
logofL1HSall<- log2(allL1HScancer)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]
logofmirnaall<- log2(allmirnacancer)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))

#pvalues
#Function to calculate p-values.
testpvalues<- function(x)
{
	print(x)
	value <-lm(logofL1HSall[,c("allL1HScancer")] ~ logofmirnaall[,x], sep="")
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(logofmirnaall))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(logofmirnaall[i]))){
		p <-testpvalues(i)
	}
	else{
		print(paste("Skipping", i))
	}

	pvalues <- cbind(pvalues, p) 
}
pvalues<- data.frame(pvalues)
tpvalues<- t(pvalues)
skip=NULL
for (i in colnames(logofmirnaall))
{
	skipped=NULL
	if(!all(is.na(logofmirnaall[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}

skip<- t(skip)
rownames(tpvalues)= skip[,1]

#rvalues
rvalues=NULL
for (i in rownames(tpvalues))
{
	r = NULL
	
		r <-cor(logofL1HSall[,c(1)], logofmirnaall[,i])

	rvalues <- cbind(rvalues, r) 
}
trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)
colnames(randpvalues)=c("pvalues", "rvalues")
sortedbypvaluesall<- randpvalues[order(-randpvalues[,c("pvalues")]),, drop=FALSE]

#create files
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results/log2 cancer by each")
write.table(sortedbypvaluesall, file="BRCA pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="CESC pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="CHOL pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="ESCA pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="HNSC pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="KICH pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="KIRC pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="KIRP pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="LIHC pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="LUAD pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="LUSC pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="PAAD pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="PCPG pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="PRAD pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="STAD pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="THCA pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="THYM pvalues and rvalues of log2 cancer.txt", sep=" ")
write.table(sortedbypvaluesall, file="UCEC pvalues and rvalues of log2 cancer.txt", sep=" ")

#ALL CANCER FOLDCHANGE RATIO
Allraw<- subset(NewNames, NewNames[,2]=="BLCA")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))
 
 
 falsenaBLCA <- !apply (is.na(logofmirnaall), 2, all)
framefalsenaBLCA<- data.frame(falsenaBLCA)
NONABLCA <- subset(framefalsenaBLCA, falsenaBLCA==TRUE)
sublogmirnaBLCA<- logofmirnaall[,c(rownames(NONABLCA))]


#pvalues
#Function to calculate p-values.
testpvalues<- function(x)
{
	value <-lm(logofL1HSall[,c("allL1HScancer")] ~ logofmirnaall[,x])
	anova(value)$Pr[1]
}
value <-lm(logofL1HSall[,c(1)] ~ logofmirnaall[,c("hsa.let.7a.1")])


#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(logofmirnaall))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(logofmirnaall[i]))){
		p <-testpvalues(i)
	}
	else{
		print(paste("Skipping", i))
	}

	pvalues <- cbind(pvalues, p) 
}
pvalues<- data.frame(pvalues)
tpvalues<- t(pvalues)
skip=NULL
for (i in colnames(logofmirnaall))
{
	skipped=NULL
	if(!all(is.na(logofmirnaall[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}

skip<- t(skip)
rownames(tpvalues)= skip[,1]

#rvalues
rvalues=NULL
for (i in rownames(tpvalues))
{
	r = NULL
	
		r <-cor(logofmirnaall[,i], logofL1HSall[,c(1)])

	rvalues <- cbind(rvalues, r) 
}
trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)
colnames(randpvalues)=c("pvalues", "rvalues")
sortedbypvaluesall<- randpvalues[order(-randpvalues[,c("pvalues")]),, drop=FALSE]

#create files
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results/log2 foldchange ratio by each")
write.table(sortedbypvaluesall, file="BLCA pvalues and rvalues of log2 foldchange ratio reverse.txt", sep=" ")
write.table(sortedbypvaluesall, file="BRCA pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="CESC pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="CHOL pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="ESCA pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="HNSC pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="KICH pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="KIRC pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="KIRP pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="LIHC pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="LUAD pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="LUSC pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="PAAD pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="PCPG pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="PRAD pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="STAD pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="THCA pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="THYM pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="UCEC pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="READ pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="FPPP pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")
write.table(sortedbypvaluesall, file="COAD pvalues and rvalues of log2 foldchange ratio.txt", sep=" ")


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



qKICH<- p.adjust(KICH[,c("pvalues")], method="BH")
finalqKICH<- cbind(KICH, qKICH)