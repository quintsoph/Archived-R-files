###Seperate the positive from the negative
#Not GDC
#BLCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#CHOL
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/CHOL")
CHOLinput<- read.table("results_lm.txt", header=TRUE)
CHOLALL<- CHOLinput[,c(1,4,5)]
LessThanzero<- ifelse(CHOLALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(CHOLALL, LessThanzero)
negCHOL<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCHOL<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCHOL, "NEG coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCHOL, "POS coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
#COADREAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/COADREAD")
COADREADinput<- read.table("results_lm.txt", header=TRUE)
COADREADALL<- COADREADinput[,c(1,4,5)]
LessThanzero<- ifelse(COADREADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(COADREADALL, LessThanzero)
negCOADREAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCOADREAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCOADREAD, "NEG coef of COADREAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCOADREAD, "POS coef of COADREAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#ESCASTAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/ESCASTAD")
ESCASTADinput<- read.table("results_lm.txt", header=TRUE)
ESCASTADALL<- ESCASTADinput[,c(1,4,5)]
LessThanzero<- ifelse(ESCASTADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(ESCASTADALL, LessThanzero)
negESCASTAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posESCASTAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negESCASTAD, "NEG coef of ESCASTAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posESCASTAD, "POS coef of ESCASTAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")

#GDC
#BLCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#CHOL
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/CHOL")
CHOLinput<- read.table("results_lm.txt", header=TRUE)
CHOLALL<- CHOLinput[,c(1,4,5)]
LessThanzero<- ifelse(CHOLALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(CHOLALL, LessThanzero)
negCHOL<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCHOL<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCHOL, "NEG coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCHOL, "POS coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
#COADREAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/COADREAD")
COADREADinput<- read.table("results_lm.txt", header=TRUE)
COADREADALL<- COADREADinput[,c(1,4,5)]
LessThanzero<- ifelse(COADREADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(COADREADALL, LessThanzero)
negCOADREAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCOADREAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCOADREAD, "NEG coef of COADREAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCOADREAD, "POS coef of COADREAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#ESCASTAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/ESCASTAD")
ESCASTADinput<- read.table("results_lm.txt", header=TRUE)
ESCASTADALL<- ESCASTADinput[,c(1,4,5)]
LessThanzero<- ifelse(ESCASTADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(ESCASTADALL, LessThanzero)
negESCASTAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posESCASTAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negESCASTAD, "NEG coef of ESCASTAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posESCASTAD, "POS coef of ESCASTAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")

#miRNA
#BLCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#CHOL
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/CHOL")
CHOLinput<- read.table("results_lm.txt", header=TRUE)
CHOLALL<- CHOLinput[,c(1,4,5)]
LessThanzero<- ifelse(CHOLALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(CHOLALL, LessThanzero)
negCHOL<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCHOL<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCHOL, "NEG coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCHOL, "POS coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
#ESCASTAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/ESCASTAD")
ESCASTADinput<- read.table("results_lm.txt", header=TRUE)
ESCASTADALL<- ESCASTADinput[,c(1,4,5)]
LessThanzero<- ifelse(ESCASTADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(ESCASTADALL, LessThanzero)
negESCASTAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posESCASTAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negESCASTAD, "NEG coef of ESCASTAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posESCASTAD, "POS coef of ESCASTAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")


####GDC check
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4")
check<- read.table("geneGDC check.txt", header=TRUE)
check$matches<- ifelse(check[,1] == check[,2], "TRUE", "FALSE")



