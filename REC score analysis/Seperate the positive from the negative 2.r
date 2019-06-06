###Seperate the positive from the negative
###GDG
#BLCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#CHOL
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/CHOL")
CHOLinput<- read.table("results_lm.txt", header=TRUE)
CHOLALL<- CHOLinput[,c(1,4,5)]
LessThanzero<- ifelse(CHOLALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(CHOLALL, LessThanzero)
negCHOL<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCHOL<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of CHOL.txt", quote=FALSE, row.names=FALSE, sep="\t")
#COAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/COAD")
COADinput<- read.table("results_lm.txt", header=TRUE)
COADALL<- COADinput[,c(1,4,5)]
LessThanzero<- ifelse(COADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(COADALL, LessThanzero)
negCOAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCOAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCOAD, "NEG coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCOAD, "POS coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#ESCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/ESCA")
ESCAinput<- read.table("results_lm.txt", header=TRUE)
ESCAALL<- ESCAinput[,c(1,4,5)]
LessThanzero<- ifelse(ESCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(ESCAALL, LessThanzero)
negESCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posESCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negESCA, "NEG coef of ESCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posESCA, "POS coef of ESCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#READ
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/READ")
READinput<- read.table("results_lm.txt", header=TRUE)
READALL<- READinput[,c(1,4,5)]
LessThanzero<- ifelse(READALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(READALL, LessThanzero)
negREAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posREAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negREAD, "NEG coef of READ.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posREAD, "POS coef of READ.txt", quote=FALSE, row.names=FALSE, sep="\t")
#STAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/STAD")
STADinput<- read.table("results_lm.txt", header=TRUE)
STADALL<- STADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negSTAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posSTAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negSTAD, "NEG coef of STAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posSTAD, "POS coef of STAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#UCEC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/UCEC")
UCECinput<- read.table("results_lm.txt", header=TRUE)
UCECALL<- UCECinput[,c(1,4,5)]
LessThanzero<- ifelse(UCECALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(UCECALL, LessThanzero)
negUCEC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posUCEC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negUCEC, "NEG coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posUCEC, "POS coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")

###miRNA
#BLCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/BLCA")
BLCAinput<- read.table("results_lm.txt", header=TRUE)
BLCAALL<- BLCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BLCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BLCAALL, LessThanzero)
negBLCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBLCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBLCA, "NEG coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBLCA, "POS coef of BLCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/BRCA")
BRCAinput<- read.table("results_lm.txt", header=TRUE)
BRCAALL<- BRCAinput[,c(1,4,5)]
LessThanzero<- ifelse(BRCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(BRCAALL, LessThanzero)
negBRCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posBRCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negBRCA, "NEG coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posBRCA, "POS coef of BRCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#COAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/COAD")
COADinput<- read.table("results_lm.txt", header=TRUE)
COADALL<- COADinput[,c(1,4,5)]
LessThanzero<- ifelse(COADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(COADALL, LessThanzero)
negCOAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posCOAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negCOAD, "NEG coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posCOAD, "POS coef of COAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#HNSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/HNSC")
HNSCinput<- read.table("results_lm.txt", header=TRUE)
HNSCALL<- HNSCinput[,c(1,4,5)]
LessThanzero<- ifelse(HNSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(HNSCALL, LessThanzero)
negHNSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posHNSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negHNSC, "NEG coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posHNSC, "POS coef of HNSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KICH
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/KICH")
KICHinput<- read.table("results_lm.txt", header=TRUE)
KICHALL<- KICHinput[,c(1,4,5)]
LessThanzero<- ifelse(KICHALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KICHALL, LessThanzero)
negKICH<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKICH<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKICH, "NEG coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKICH, "POS coef of KICH.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/KIRC")
KIRCinput<- read.table("results_lm.txt", header=TRUE)
KIRCALL<- KIRCinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRCALL, LessThanzero)
negKIRC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRC, "NEG coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRC, "POS coef of KIRC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/KIRP")
KIRPinput<- read.table("results_lm.txt", header=TRUE)
KIRPALL<- KIRPinput[,c(1,4,5)]
LessThanzero<- ifelse(KIRPALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(KIRPALL, LessThanzero)
negKIRP<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posKIRP<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negKIRP, "NEG coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posKIRP, "POS coef of KIRP.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/LIHC")
LIHCinput<- read.table("results_lm.txt", header=TRUE)
LIHCALL<- LIHCinput[,c(1,4,5)]
LessThanzero<- ifelse(LIHCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LIHCALL, LessThanzero)
negLIHC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLIHC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLIHC, "NEG coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLIHC, "POS coef of LIHC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/LUAD")
LUADinput<- read.table("results_lm.txt", header=TRUE)
LUADALL<- LUADinput[,c(1,4,5)]
LessThanzero<- ifelse(LUADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUADALL, LessThanzero)
negLUAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUAD, "NEG coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUAD, "POS coef of LUAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#LUSC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/LUSC")
LUSCinput<- read.table("results_lm.txt", header=TRUE)
LUSCALL<- LUSCinput[,c(1,4,5)]
LessThanzero<- ifelse(LUSCALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(LUSCALL, LessThanzero)
negLUSC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posLUSC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negLUSC, "NEG coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posLUSC, "POS coef of LUSC.txt", quote=FALSE, row.names=FALSE, sep="\t")
#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/PRAD")
PRADinput<- read.table("results_lm.txt", header=TRUE)
PRADALL<- PRADinput[,c(1,4,5)]
LessThanzero<- ifelse(PRADALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(PRADALL, LessThanzero)
negPRAD<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posPRAD<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negPRAD, "NEG coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posPRAD, "POS coef of PRAD.txt", quote=FALSE, row.names=FALSE, sep="\t")
#THCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/THCA")
THCAinput<- read.table("results_lm.txt", header=TRUE)
THCAALL<- THCAinput[,c(1,4,5)]
LessThanzero<- ifelse(THCAALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(THCAALL, LessThanzero)
negTHCA<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posTHCA<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negTHCA, "NEG coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posTHCA, "POS coef of THCA.txt", quote=FALSE, row.names=FALSE, sep="\t")
#UCEC
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/REC Score Stuff/UCEC")
UCECinput<- read.table("results_lm.txt", header=TRUE)
UCECALL<- UCECinput[,c(1,4,5)]
LessThanzero<- ifelse(UCECALL[,2] < 0 , TRUE, FALSE)
combined<- cbind(UCECALL, LessThanzero)
negUCEC<- subset(combined, LessThanzero=="TRUE", select=c(genenames, generatio_coef, generatio_pval))
posUCEC<- subset(combined, LessThanzero=="FALSE", select=c(genenames, generatio_coef, generatio_pval))
write.table(negUCEC, "NEG coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(posUCEC, "POS coef of UCEC.txt", quote=FALSE, row.names=FALSE, sep="\t")
