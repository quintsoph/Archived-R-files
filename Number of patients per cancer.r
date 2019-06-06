#Find patients per cancer patients
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
Names<- read.table("SophiaPatCanTypes", header=FALSE)
tNames<- t(Names)
FPPP<- subset(tNames, tNames[,2] == "FPPP", drop=FALSE)


setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
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
#ratio<- cancer/normal

#BLCA
BLCA<- subset(tNames, tNames[,2] == "BLCA", drop=FALSE)
matched<- match(BLCA[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(BLCA[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.195")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.143")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.497")])
third<- cbind(L1HS, mirna3)
plot(third[,2], third[,1], xlab="hsa.mir.497 in BLCA", ylab="L1HS", main="MiRNA 497 in BLCA", col="blue")
abline(lm(third[,1] ~ third[,2]))
cor(third[,2], third[,1])

#BRCA
BRCA<- subset(tNames, tNames[,2] == "BRCA", drop=FALSE)
matched<- match(BRCA[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(BRCA[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.22")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.486")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.378")])
third<- cbind(L1HS, mirna3)
plot(first[,2], first[,1], xlab="hsa.mir.22 in BRCA", ylab="L1HS", main="MiRNA 22 in BRCA", col="blue")
abline(lm(first[,1] ~ first[,2]))
cor(first[,2], first[,1])

#HNSC
HNSC<- subset(tNames, tNames[,2] == "HNSC", drop=FALSE)
matched<- match(HNSC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(HNSC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.29a")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.424")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.450a.2")])
third<- cbind(L1HS, mirna3)
plot(third[,2], third[,1], xlab="hsa.mir.450a.2 in HNSC", ylab="L1HS", main="MiRNA 450a.2 in HNSC", col="blue")
abline(lm(third[,1] ~ third[,2]))
cor(third[,2], third[,1])

#KICH
KICH<- subset(tNames, tNames[,2] == "KICH", drop=FALSE)
matched<- match(KICH[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(KICH[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.374a")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.503")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.505")])
third<- cbind(L1HS, mirna3)
plot(first[,2], first[,1], xlab="hsa.mir.374a in KICH", ylab="L1HS", main="MiRNA 374a in KICH", col="blue")
abline(lm(first[,1] ~ first[,2]))
cor(first[,2], first[,1])

#KIRC
KIRC<- subset(tNames, tNames[,2] == "KIRC", drop=FALSE)
matched<- match(KIRC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(KIRC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.379")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.127")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.21")])
third<- cbind(L1HS, mirna3)
plot(third[,2], third[,1], xlab="hsa.mir.21 in KIRC", ylab="L1HS", main="MiRNA 21 in KIRC", col="blue")
abline(lm(third[,1] ~ third[,2]))
cor(third[,2], third[,1])

#KIRP
KIRP<- subset(tNames, tNames[,2] == "KIRP", drop=FALSE)
matched<- match(KIRP[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(KIRP[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.let.7d")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.181b.2")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.148a")])
third<- cbind(L1HS, mirna3)
plot(first[,2], first[,1], xlab="hsa.let.7d in KIRP", ylab="L1HS", main="MiRNA let.7d in KIRP", col="blue")
abline(lm(first[,1] ~ first[,2]))
cor(first[,2], first[,1])

#LIHC
LIHC<- subset(tNames, tNames[,2] == "LIHC", drop=FALSE)
matched<- match(LIHC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(LIHC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.29a")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.29b.2")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.29b.1")])
third<- cbind(L1HS, mirna3)
plot(third[,2], third[,1], xlab="hsa.mir.29b.1 in LIHC", ylab="L1HS", main="MiRNA 29b.1 in LIHC", col="blue")
abline(lm(third[,1] ~ third[,2]))
cor(third[,2], third[,1])

#LUAD
LUAD<- subset(tNames, tNames[,2] == "LUAD", drop=FALSE)
matched<- match(LUAD[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(LUAD[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.542")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.10b")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.101.1")])
third<- cbind(L1HS, mirna3)
plot(third[,2], third[,1], xlab="hsa.mir.101.1 in LUAD", ylab="L1HS", main="MiRNA 101.1 in LUAD", col="blue")
abline(lm(third[,1] ~ third[,2]))
cor(third[,2], third[,1])

#LUSC
LUSC<- subset(tNames, tNames[,2] == "LUSC", drop=FALSE)
matched<- match(LUSC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(LUSC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.let.7e")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.181b.1")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.153.2")])
third<- cbind(L1HS, mirna3)
plot(third[,2], third[,1], xlab="hsa.mir.153.2 in LUSC", ylab="L1HS", main="MiRNA 153.2 in LUSC", col="blue")
abline(lm(third[,1] ~ third[,2]))
cor(third[,2], third[,1])

#PRAD
PRAD<- subset(tNames, tNames[,2] == "PRAD", drop=FALSE)
matched<- match(PRAD[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(PRAD[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.204")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.378c")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.487b")])
third<- cbind(L1HS, mirna3)
plot(first[,2], first[,1], xlab="hsa.mir.204 in PRAD", ylab="L1HS", main="MiRNA 204 in PRAD", col="blue")
abline(lm(first[,1] ~ first[,2]))
cor(first[,2], first[,1])

#THCA
THCA<- subset(tNames, tNames[,2] == "THCA", drop=FALSE)
matched<- match(THCA[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(THCA[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.514.3")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.508")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.514.1")])
third<- cbind(L1HS, mirna3)
plot(third[,2], third[,1], xlab="hsa.mir.514.1 in THCA", ylab="L1HS", main="MiRNA 514.1 in THCA", col="blue")
abline(lm(third[,1] ~ third[,2]))
cor(third[,2], third[,1])

#UCEC
UCEC<- subset(tNames, tNames[,2] == "UCEC", drop=FALSE)
matched<- match(UCEC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(UCEC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.133a.1")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.1.2")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.145")])
third<- cbind(L1HS, mirna3)
plot(first[,2], first[,1], xlab="hsa.mir.133a.1 in UCEC", ylab="L1HS", main="MiRNA 133a.1 in UCEC", col="blue")
abline(lm(first[,1] ~ first[,2]))
cor(first[,2], first[,1])
