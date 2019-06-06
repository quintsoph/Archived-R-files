###Seperate the positive from the negative
###GDG
#BLCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#COAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/COAD")
COADinput<- read.table("results_lm.txt", header=TRUE)
COADALL<- COADinput[,c(1,4,5)]
LessThanzero<- ifelse(COADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(COADALL, LessThanzero)
negCOAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCOAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCOAD, "NEG coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCOAD, "POS coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")

#Regular
#BLCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#COAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/COAD")
COADinput<- read.table("results_lm.txt", header=TRUE)
COADALL<- COADinput[,c(1,4,5)]
LessThanzero<- ifelse(COADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(COADALL, LessThanzero)
negCOAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCOAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCOAD, "NEG coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCOAD, "POS coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#UCEC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular/UCEC")
UCECinput<- read.table("results_lm.txt", header=TRUE)
UCECALL<- UCECinput[,c(1,4,5)]
LessThanzero<- ifelse(UCECALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(UCECALL, LessThanzero)
negUCEC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posUCEC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negUCEC, "NEG coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posUCEC, "POS coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")

#miRNA
#BLCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/MiRNA/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#COAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/COAD")
COADinput<- read.table("results_lm.txt", header=TRUE)
COADALL<- COADinput[,c(1,4,5)]
LessThanzero<- ifelse(COADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(COADALL, LessThanzero)
negCOAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCOAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCOAD, "NEG coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCOAD, "POS coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#UCEC
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/miRNA/UCEC")
UCECinput<- read.table("results_lm.txt", header=TRUE)
UCECALL<- UCECinput[,c(1,4,5)]
LessThanzero<- ifelse(UCECALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(UCECALL, LessThanzero)
negUCEC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posUCEC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negUCEC, "NEG coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posUCEC, "POS coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")
