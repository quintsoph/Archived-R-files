#### Create a patient file with cancer and the TCGA data
setwd("E:/Bioinformatics Lab/Cancer Data/Methylation Comparison")
GDCpatients<- read.table("GDC.patient_ids.paired.txt")
L1HSpatients<- read.table("L1HS normal for Graphs.txt", header=TRUE)
#creating a fixed file
library(stringi)
GDCsmall<-stri_sub(GDCpatients[,1],9,-20)
GDCtogether<- cbind(GDCpatients, GDCsmall)
GDCcut<- GDCtogether[-c(643, 663, 760, 770, 797, 809, 853, 879, 951, 984, 1035, 1093, 1144, 1196, 1203, 1236, 1294),]
GDC1<- GDCcut[c(1:642),]
GDC<- stri_sub(GDC1[,2],0,-2)
GDCfix<- cbind(GDC1, GDC)
GDCfix<- GDCfix[,-2]
colnames(GDCfix)=c("full list", "short list")
GDCremoved<- GDCcut[-c(1:642),]
colnames(GDCremoved)=c("full list", "short list")
CombinedGDC<- rbind(GDCfix, GDCremoved)
matches<- CombinedGDC[which(L1HSpatients[,1] %in% CombinedGDC[,2]),]
matches2<- matches[order(matches[,2]),, drop=FALSE]
library(dplyr)
rownames(matches2)=matches2[,2]
rownames(L1HSpatients)=L1HSpatients[,1]
Finalized<- merge(matches2, L1HSpatients, by="row.names", all=TRUE)
Finalized1<- na.omit(Finalized)
Finalized2<- Finalized1[,-c(3,4)]
Finalized3<- Finalized2[,c(2,3)]
cuts<-stri_sub(Finalized3[,1], 0, -21)
Resultout<- cbind(cuts, Finalized3)
Resultout<- Resultout[,c(1,3)]
#seperating my cancer type
setwd("E:/Bioinformatics Lab/Cancer Data/Methylation Comparison/By cancer L1HS matches")
BLCA<- subset(Resultout, Resultout[,2]=="BLCA")
write.table(BLCA[,1], file="BLCA patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
BRCA<- subset(Resultout, Resultout[,2]=="BRCA")
write.table(BRCA[,1], file="BRCA patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
COAD<- subset(Resultout, Resultout[,2]=="COAD")
write.table(COAD[,1], file="COAD patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
ESCA<- subset(Resultout, Resultout[,2]=="ESCA")
write.table(ESCA[,1], file="ESCA patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
HNSC<- subset(Resultout, Resultout[,2]=="HNSC")
write.table(HNSC[,1], file="HNSC patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
KIRC<- subset(Resultout, Resultout[,2]=="KIRC")
write.table(KIRC[,1], file="KIRC patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
KIRP<- subset(Resultout, Resultout[,2]=="KIRP")
write.table(KIRP[,1], file="KIRP patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
LIHC<- subset(Resultout, Resultout[,2]=="LIHC")
write.table(LIHC[,1], file="LIHC patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
LUAD<- subset(Resultout, Resultout[,2]=="LUAD")
write.table(LUAD[,1], file="LUAD patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
LUSC<- subset(Resultout, Resultout[,2]=="LUSC")
write.table(LUSC[,1], file="LUSC patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
PRAD<- subset(Resultout, Resultout[,2]=="PRAD")
write.table(PRAD[,1], file="PRAD patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
READ<- subset(Resultout, Resultout[,2]=="READ")
write.table(READ[,1], file="READ patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
STAD<- subset(Resultout, Resultout[,2]=="STAD")
write.table(STAD[,1], file="STAD patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
THCA<- subset(Resultout, Resultout[,2]=="THCA")
write.table(THCA[,1], file="THCA patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
UCEC<- subset(Resultout, Resultout[,2]=="UCEC")
write.table(UCEC[,1], file="UCEC patients.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)




