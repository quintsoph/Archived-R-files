#####NEW miRNA L1HS correlation with Median Mutations
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
NewNames<- data.frame(read.table(file="SophiapatCanTypes"))
NewNames<- t(NewNames)

#read in libraries
library(methods)
library(stringi)
library(matrixcalc)

#read in L1HS matches
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
L1HSfold<- read.table("L1HS_Gene_Matches.txt", header=FALSE, row.names=1)
tL1HS<- t(L1HSfold)
colnames(tL1HS)=rownames(L1HSfold)
#read in the Mutation file
Mutations<- data.frame(read.table("Mutation_Matches.txt", header=FALSE, row.names=1))
tMutations<- data.frame(t(Mutations))
colnames(tMutations)=rownames(Mutations)

#TEST WITH ONE
BLCAraw<- subset(NewNames, NewNames[,2]=="BLCA")
dfBLCAraw<- data.frame(BLCAraw)

#create log2 foldchange of miRNA
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
newmirnacancer<- read.csv(file="Mirna.Cancer.csv")
newmirnacancer<- data.matrix(newmirnacancer)
cancermirnanames<- read.csv(file="cancer mirna names.csv")
rownames(newmirnacancer)=cancermirnanames[,1]
tnewmirnacancer<- t(newmirnacancer)
tnewmirnacancer<- tnewmirnacancer[-1,]
newmirnanormal<- read.csv(file="Mirna.Normal.csv")
newmirnanormal<- data.matrix(newmirnanormal)
normalmirnanames<- read.csv(file="normal mirna names.csv")
rownames(newmirnanormal)=normalmirnanames[,1]
tnewmirnanormal<- t(newmirnanormal)
tnewmirnanormal<- as.matrix(tnewmirnanormal[-1,])
#remove patient A0DB
removename<- "TCGA.A7.A0DB"
tnewmirnacancernop<- tnewmirnacancer[!rownames(tnewmirnacancer) %in% removename, ]
tnewmirnanormalnop<- tnewmirnanormal[!rownames(tnewmirnanormal) %in% removename, ]
#create log2 foldchange
dividenmirna<- tnewmirnacancernop/tnewmirnanormalnop
foldedmirna<- t(log2(dividenmirna))
colnames(foldedmirna)=stri_sub(colnames(foldedmirna),9)

#upload mutations
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
Mutations<- data.frame(read.table("Mutation_Matches.txt", header=FALSE, row.names=1))
tMutations<- data.frame(t(Mutations))
colnames(tMutations)=rownames(Mutations)
#combine L1HS and mirna
tL1HS<- data.frame(tL1HS, check.names=FALSE)
L1HS<- data.frame(t(tL1HS), check.names=FALSE)
foldedmirna<- data.frame(foldedmirna, check.names=FALSE)
tfoldedmirna<- data.frame(t(foldedmirna), check.names=FALSE)
inter<-intersect(rownames(tfoldedmirna), rownames(L1HS))
testm<-match(inter, rownames(tfoldedmirna))
testm<- testm[!is.na(testm)]
comb<- cbind(rownames(tfoldedmirna), tfoldedmirna)
resultnamesmirna<- data.frame(comb[testm,])
match(rownames(resultnamesmirna), rownames(L1HS))
resultnamesmirnaedit<- resultnamesmirna[,-1]
mirnacut<- data.frame(t(resultnamesmirnaedit), check.names=FALSE)
library(dplyr)
comp<- bind_rows(tL1HS, mirnacut)
mirnas<- rownames(mirnacut)
rownames(comp)= c("L1HS", mirnas)
comp2<- bind_rows(tMutations, comp)
rownames(comp2)=c("Mutations", "L1HS", mirnas)
tcomp<- data.frame(t(comp2))
detach("package:dplyr")

##BLCA
is.na(tcomp)<- sapply(tcomp, is.infinite)
is.na(tcomp)<- sapply(tcomp, is.nan)
mirnaset<- tcomp[,-c(1,2)]
edits<- cbind(NewNames, stri_sub(NewNames[,1], 9))
CancerPatients<- edits[,-c(1)]
test<-intersect(CancerPatients[,2], rownames(mirnaset))
match(test, rownames(mirnaset))
testm<- match(test, CancerPatients[,2])
testm<- testm[!is.na(testm)]
comb<- CancerPatients
rownames(comb)=CancerPatients[,2]
CPatients<- data.frame(comb[testm,])
BLCA<- subset(CPatients, CPatients[,1]=="BLCA")
BLCAcut<- tcomp[which(row.names(tcomp) %in% BLCA[,2]),]
ALLL1HSandMUT<- BLCAcut
ALLGENES<- data.frame(BLCAcut[,-c(1,2)])
