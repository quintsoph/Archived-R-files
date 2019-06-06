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

#STOP
#All Cancer pvalues log R setup
falsenacancer <- !apply (is.na(foldedmirna), 2, all)
framefalsenacancer<- data.frame(falsenacancer)
NONAcancer <- subset(framefalsenacancer, falsenacancer==TRUE)
sublogmirna<- foldedmirna[,c(rownames(NONAcancer))]

#START
is.na(foldedmirna) <- do.call(cbind,lapply(foldedmirna, is.infinite))
is.na(foldedmirna) <- do.call(cbind,lapply(foldedmirna, is.nan))
#Function to calculate p-values.
testpvalues<- function(x)
{
	value <-lm(foldedL1HS[,1] ~ foldedmirna[,x], sep="")
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(foldedmirna))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(foldedmirna[i]))){
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
for (i in colnames(foldedmirna))
{
	skipped=NULL
	if(!all(is.na(foldedmirna[i]))){
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
	
		r <-cor(t(foldedL1HS[,c(1)]), foldedmirna[,i])

	rvalues <- cbind(rvalues, r) 
}
trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)

trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)
colnames(randpvalues)=c("pvalues", "rvalues")
sortedbypvaluesall<- randpvalues[order(-randpvalues[,c("pvalues")]),, drop=FALSE]

#qvalues
qvalues<- p.adjust(tpvalues, method="BH")
final<- cbind(randpvalues, qvalues)

setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results")
write.table(final, file="log2 cancer only.txt", sep=" ")











#log2 foldchange
ratio<- cancer/normal
logratio<- log2(ratio)
logratiomirna<- logratio[,-1]
is.na(logratio) <- sapply(logratio, is.infinite)
is.na(logratiomirna) <- sapply(logratiomirna, is.infinite)
logratio<- data.frame(logratio)
logratiomirna<- data.frame(logratiomirna)
testpvalues<- function(x)
{
	print(x)
	value <-lm(logratio[,c("tnewL1HScancer")] ~ logratiomirna[,x], sep="")
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(logratiomirna))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(logratiomirna[i]))){
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
for (i in colnames(logratiomirna))
{
	skipped=NULL
	if(!all(is.na(logratiomirna[i]))){
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
	
		r <-cor(logratio[,c(1)], logratiomirna[,i])

	rvalues <- cbind(rvalues, r) 
}
trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)

#qvalues
qvalues<- p.adjust(tpvalues, method="BH")
final<- cbind(randpvalues, qvalues)

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results")
write.table(final, file="log2 foldchange ratio only.txt", sep=" ")


setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
MedMut<- data.frame(read.table("Mutationsnoout.txt", header=TRUE, row.names=1))
library(stringi)
edfold<- stri_sub(rownames(foldedmirna),9)
rownames(foldedmirna)=edfold



value <-lm(foldedL1HS[,1] ~ foldedmirna[,c("hsa-mir-99a")] +, sep="")
	


