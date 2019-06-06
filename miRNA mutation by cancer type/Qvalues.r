###QVALUES for each cancer type

#BLCA
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
BLCA<- read.table("BLCAresults.txt", header=TRUE, row.names=1)
BLCApvals<- BLCA[,1]
BLCAQ<-  data.frame(p.adjust(BLCApvals, method="BH"))
rownames(BLCAQ)= rownames(BLCA)
colnames(BLCAQ)= c("BLCA Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(BLCAQ, file="BLCA Q-Values.txt", sep="\t", quote=FALSE)

#BRCA
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
BRCA<- read.table("BRCAresults.txt", header=TRUE, row.names=1)
BRCApvals<- BRCA[,1]
BRCAQ<-  data.frame(p.adjust(BRCApvals, method="BH"))
rownames(BRCAQ)= rownames(BRCA)
colnames(BRCAQ)= c("BRCA Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(BRCAQ, file="BRCA Q-Values.txt", sep="\t", quote=FALSE)

#CHOL
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
CHOL<- read.table("CHOLresults.txt", header=TRUE, row.names=1)
CHOLpvals<- CHOL[,1]
CHOLQ<-  data.frame(p.adjust(CHOLpvals, method="BH"))
rownames(CHOLQ)= rownames(CHOL)
colnames(CHOLQ)= c("CHOL Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(CHOLQ, file="CHOL Q-Values.txt", sep="\t", quote=FALSE)

#COAD
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
COAD<- read.table("COADresults.txt", header=TRUE, row.names=1)
COADpvals<- COAD[,1]
COADQ<-  data.frame(p.adjust(COADpvals, method="BH"))
rownames(COADQ)= rownames(COAD)
colnames(COADQ)= c("COAD Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(COADQ, file="COAD Q-Values.txt", sep="\t", quote=FALSE)

#ESCA
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
ESCA<- read.table("ESCAresults.txt", header=TRUE, row.names=1)
ESCApvals<- ESCA[,1]
ESCAQ<-  data.frame(p.adjust(ESCApvals, method="BH"))
rownames(ESCAQ)= rownames(ESCA)
colnames(ESCAQ)= c("ESCA Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(ESCAQ, file="ESCA Q-Values.txt", sep="\t", quote=FALSE)

#KICH
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
KICH<- read.table("KICHresults.txt", header=TRUE, row.names=1)
KICHpvals<- KICH[,1]
KICHQ<-  data.frame(p.adjust(KICHpvals, method="BH"))
rownames(KICHQ)= rownames(KICH)
colnames(KICHQ)= c("KICH Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(KICHQ, file="KICH Q-Values.txt", sep="\t", quote=FALSE)

#KIRC
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
KIRC<- read.table("KIRCresults.txt", header=TRUE, row.names=1)
KIRCpvals<- KIRC[,1]
KIRCQ<-  data.frame(p.adjust(KIRCpvals, method="BH"))
rownames(KIRCQ)= rownames(KIRC)
colnames(KIRCQ)= c("KIRC Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(KIRCQ, file="KIRC Q-Values.txt", sep="\t", quote=FALSE)

#KIRP
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
KIRP<- read.table("KIRPresults.txt", header=TRUE, row.names=1)
KIRPpvals<- KIRP[,1]
KIRPQ<-  data.frame(p.adjust(KIRPpvals, method="BH"))
rownames(KIRPQ)= rownames(KIRP)
colnames(KIRPQ)= c("KIRP Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(KIRPQ, file="KIRP Q-Values.txt", sep="\t", quote=FALSE)

#LIHC
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
LIHC<- read.table("LIHCresults.txt", header=TRUE, row.names=1)
LIHCpvals<- LIHC[,1]
LIHCQ<-  data.frame(p.adjust(LIHCpvals, method="BH"))
rownames(LIHCQ)= rownames(LIHC)
colnames(LIHCQ)= c("LIHC Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(LIHCQ, file="LIHC Q-Values.txt", sep="\t", quote=FALSE)

#LUAD
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
LUAD<- read.table("LUADresults.txt", header=TRUE, row.names=1)
LUADpvals<- LUAD[,1]
LUADQ<-  data.frame(p.adjust(LUADpvals, method="BH"))
rownames(LUADQ)= rownames(LUAD)
colnames(LUADQ)= c("LUAD Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(LUADQ, file="LUAD Q-Values.txt", sep="\t", quote=FALSE)

#PAAD
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
PAAD<- read.table("PAADresults.txt", header=TRUE, row.names=1)
PAADpvals<- PAAD[,1]
PAADQ<-  data.frame(p.adjust(PAADpvals, method="BH"))
rownames(PAADQ)= rownames(PAAD)
colnames(PAADQ)= c("PAAD Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(PAADQ, file="PAAD Q-Values.txt", sep="\t", quote=FALSE)

#PRAD
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
PRAD<- read.table("PRADresults.txt", header=TRUE, row.names=1)
PRADpvals<- PRAD[,1]
PRADQ<-  data.frame(p.adjust(PRADpvals, method="BH"))
rownames(PRADQ)= rownames(PRAD)
colnames(PRADQ)= c("PRAD Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(PRADQ, file="PRAD Q-Values.txt", sep="\t", quote=FALSE)

#THCA
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
THCA<- read.table("THCAresults.txt", header=TRUE, row.names=1)
THCApvals<- THCA[,1]
THCAQ<-  data.frame(p.adjust(THCApvals, method="BH"))
rownames(THCAQ)= rownames(THCA)
colnames(THCAQ)= c("THCA Q-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov/QVals")
write.table(THCAQ, file="THCA Q-Values.txt", sep="\t", quote=FALSE)
