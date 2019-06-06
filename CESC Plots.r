setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
NewNames<- data.frame(read.table(file="SophiapatCanTypes"))
NewNames<- t(NewNames)
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
#log2ratio
L1HSratio<- newL1HScancer/newL1HSnormal
L1HS<- data.frame(log2(L1HSratio))
L1HS<- data.frame(L1HS[,-c(1)])
mirnaratio<- newmirnacancer/newmirnanormal
mirna<- data.frame(log2(mirnaratio))
mirna<- data.frame(mirna[,-c(1)])
#Select CESC cancer type
CESC<- subset(NewNames, NewNames[,2]=="CESC")
matched<- match(CESC[,1], colnames(L1HS))
CESCL1HS <-L1HS[,matched]
tL1HS<- t(CESCL1HS)
colnames(tL1HS)= c("L1HS")
matched<- match(CESC[,1], colnames(mirna))
CESCmirna <-mirna[,matched]
num1<- CESCmirna[c("hsa-mir-128-1"),]
tnum1<- t(num1)
num2<- CESCmirna[c("hsa-mir-128-2"),]
tnum2<- t(num2)
CESCmirna128<- cbind(tL1HS, tnum1, tnum2)
#Time to plot
plot(CESCmirna128[,3], CESCmirna128[,1], xlab="log2 fold change mirna", ylab="log2 fold change L1HS", main="CESC mir-128-2", col="dark green")
value<- lm(CESCmirna128[,1] ~ CESCmirna128[,3])
abline(value)
#pvalues
p<- anova(value)$Pr[1]

#cancer only correlation
mirnacancer<- data.frame(newmirnacancer[,-c(1)])
mir1281<- mirnacancer[c("hsa-mir-128-1"),]
mir1282<- mirnacancer[c("hsa-mir-128-2"),]
mirnaf<- rbind(mir1281, mir1282)
L1HScancer<- data.frame(newL1HScancer[,-c(1)])
#Select CESC cancer type
CESC<- subset(NewNames, NewNames[,2]=="CESC")
matched<- match(CESC[,1], rownames(L1HScancer))
tL1HScancer<- t(L1HScancer)
L1HS <-data.frame(tL1HScancer[,matched])
colnames(L1HS)= c("L1HS")
matched<- match(CESC[,1], colnames(mirnaf))
mirna <-mirnaf[,matched]
tmirna<- data.frame(t(mirna))
#pvalues
#Function to calculate p-values.
testpvalues<- function(x)
{
	print(x)
	value <-lm(L1HS[,c("L1HS")] ~ tmirna[,x])
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(tmirna))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(tmirna[i]))){
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
for (i in colnames(tmirna))
{
	skipped=NULL
	if(!all(is.na(tmirna[i]))){
	skipped<- print(paste(i))
	}
	skip<- data.frame(cbind(skip, skipped))
	}

skip<- t(skip)
rownames(tpvalues)= skip[,1]

#rvalues
rvalues=NULL
for (i in rownames(tpvalues))
{
	r = NULL
	
		r <-cor(L1HS[,c(1)], tmirna[,i])

	rvalues <- cbind(rvalues, r) 
}
trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)
colnames(randpvalues)=c("pvalues", "rvalues")
sortedbypvalues<- randpvalues[order(-randpvalues[,c("pvalues")]),, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/mir 128 p and r/Cancer log2")
write.table(sortedbypvalues, file="COAD pvalues and rvalues mir 128 cancer.txt", sep="\t", quote=FALSE)



#cancer log2 correlation
mirnacancer<- data.frame(log2(newmirnacancer[,-c(1)]))
mir1281<- mirnacancer[c("hsa-mir-128-1"),]
mir1282<- mirnacancer[c("hsa-mir-128-2"),]
mirnaf<- rbind(mir1281, mir1282)
L1HScancer<- data.frame(log2(newL1HScancer[,-c(1)]))
#Select CESC cancer type
COAD<- subset(NewNames, NewNames[,2]=="COAD")
matched<- match(COAD[,1], rownames(L1HScancer))
tL1HScancer<- t(L1HScancer)
L1HS <-data.frame(tL1HScancer[,matched])
colnames(L1HS)= c("L1HS")
matched<- match(COAD[,1], colnames(mirnaf))
mirna <-mirnaf[,matched]
tmirna<- data.frame(t(mirna))

#pvalues
#Function to calculate p-values.
testpvalues<- function(x)
{
	print(x)
	value <-lm(L1HS[,c("L1HS")] ~ tmirna[,x])
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(tmirna))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(tmirna[i]))){
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
for (i in colnames(tmirna))
{
	skipped=NULL
	if(!all(is.na(tmirna[i]))){
	skipped<- print(paste(i))
	}
	skip<- data.frame(cbind(skip, skipped))
	}

skip<- t(skip)
rownames(tpvalues)= skip[,1]

#rvalues
rvalues=NULL
for (i in rownames(tpvalues))
{
	r = NULL
	
		r <-cor(L1HS[,c(1)], tmirna[,i])

	rvalues <- cbind(rvalues, r) 
}
trvalues<- t(rvalues)
randpvalues<- cbind(tpvalues, trvalues)
colnames(randpvalues)=c("pvalues", "rvalues")
sortedbypvalues<- randpvalues[order(-randpvalues[,c("pvalues")]),, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/mir 128 p and r/Cancer log2")
write.table(sortedbypvalues, file="COAD pvalues and rvalues mir 128 cancer.txt", sep="\t", quote=FALSE)

