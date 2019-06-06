#Upload new L1HS data
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/L1HS expression/New Data")
L1HSupload<- read.table(file="New L1HS.txt", header=TRUE)
rownames(L1HSupload)= L1HSupload[,1]
L1HSnames<- sub(pattern="-01A_gene_TE_analysis.txt:L1HS:L1:LINE*", replacement="", x=rownames(L1HSupload), perl=TRUE)
L1HSnames2<- gsub(pattern="^TCGA-*-", replacement="", L1HSnames)
L1HSnames3<- gsub(pattern="^..-", replacement="", L1HSnames2)
rownames(L1HSupload)= L1HSnames3
L1HSALL<- L1HSupload[,-1]
ratioL1HS<- data.frame(L1HSALL[,4])
rownames(ratioL1HS)= rownames(L1HSALL)
tratioL1HS<- t(ratioL1HS)

#create averages by cancer/normal
#BRCA/BLCA
m1<- merge(BLCAmerge, BRCAmerge)
notna<- colSums(!is.na(m1))
sum<- colSums(m1, na.rm=TRUE)
avmutBLCABRCA<- data.frame(sum/notna)
#CESC
notna<- colSums(!is.na(CESCmerge))
sum<- colSums(CESCmerge, na.rm=TRUE)
avmutCESC<- data.frame(sum/notna)
#CHOL
notna<- colSums(!is.na(CHOLmerge))
sum<- colSums(CHOLmerge, na.rm=TRUE)
avmutCHOL<- data.frame(sum/notna)
#COAD
notna<- colSums(!is.na(COADmerge))
sum<- colSums(Coadmerge, na.rm=TRUE)
avmutCOAD<- data.frame(sum/notna)
#ESCA
notna<- colSums(!is.na(ESCAmerge))
sum<- colSums(ESCAmerge, na.rm=TRUE)
avmutESCA<- data.frame(sum/notna)
#HNSCm4
notna<- colSums(!is.na(HNSCm4))
sum<- colSums(HNSCm4, na.rm=TRUE)
avmutHNSC<- data.frame(sum/notna)
#KICH
notna<- colSums(!is.na(KICHmerge))
sum<- colSums(KICHmerge, na.rm=TRUE)
avmutKICH<- data.frame(sum/notna)
#KIRC
notna<- colSums(!is.na(KIRCmerge))
sum<- colSums(KIRCmerge, na.rm=TRUE)
avmutKIRC<- data.frame(sum/notna)
#KIRP
notna<- colSums(!is.na(KIRPm6))
sum<- colSums(KIRPm6, na.rm=TRUE)
avmutKIRP<- data.frame(sum/notna)
#LIHC
notna<- colSums(!is.na(LIHCm5))
sum<- colSums(LIHCm5, na.rm=TRUE)
avmutLIHC<- data.frame(sum/notna)
#LUAD
notna<- colSums(!is.na(LUADmerge))
sum<- colSums(LIHCm5, na.rm=TRUE)
avmutLUAD<- data.frame(sum/notna)
#PAAD
notna<- colSums(!is.na(PAADm5))
sum<- colSums(PAADm5, na.rm=TRUE)
avmutPAAD<- data.frame(sum/notna)
#PCPG
notna<- colSums(!is.na(PCPGm3))
sum<- colSums(PCPGm3, na.rm=TRUE)
avmutPCPG<- data.frame(sum/notna)
#PRAD
notna<- colSums(!is.na(PRADm4))
sum<- colSums(PRADm4, na.rm=TRUE)
avmutPRAD<- data.frame(sum/notna)
#SARCm3
notna<- colSums(!is.na(SARCm3))
sum<- colSums(SARCm3, na.rm=TRUE)
avmutSARC<- data.frame(sum/notna)
#THCA
notna<- colSums(!is.na(THCAm4))
sum<- colSums(THCAm4, na.rm=TRUE)
avmutTHCA<- data.frame(sum/notna)
#THYM
notna<- colSums(!is.na(THYMm4))
sum<- colSums(THYMm4, na.rm=TRUE)
avmutTHYM<- data.frame(sum/notna)
#UCEC
notna<- colSums(!is.na(UCECm3))
sum<- colSums(UCECm3, na.rm=TRUE)
avmutUCEC<- data.frame(sum/notna)
#time to rbind!
mergemut1<- rbind(avmutBLCABRCA, avmutCESC)
mergemut2<- rbind(mergemut1, avmutCHOL)
mergemut3<- rbind(mergemut2, avmutESCA)
mergemut5<- rbind(mergemut3, avmutHNSC)
mergemut6<- rbind(mergemut5, avmutKIRP)
mergemut7<- rbind(mergemut6, avmutLIHC)
mergemut8<- rbind(mergemut7, avmutPAAD)
mergemut9<- rbind(mergemut8, avmutPCPG)
mergemut10<- rbind(mergemut9, avmutPRAD)
mergemut11<- rbind(mergemut10, avmutSARC)
mergemut12<- rbind(mergemut11, avmutTHCA)
mergemut13<- rbind(mergemut12, avmutTHYM)
mergemut14<- rbind(mergemut13, avmutUCEC)
mergemut15<- rbind(mergemut14, avmutCOAD)
mergemut16<- rbind(mergemut15, avmutKICH)
mergemut17<- rbind(mergemut16, avmutKIRC)
mergemut18<- rbind(mergemut17, avmutLUAD)
mergemut18<- data.frame(mergemut18)
#match L1HS
tmergemut18<- t(mergemut18)
matchedL<- match(rownames(mergemut18), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HSmatched<- data.frame(tratioL1HS[,matchedL])
tL1HSmatched<- t(L1HSmatched)
combineddata<- merge(tmergemut18, tL1HSmatched, all=TRUE)
combineddata<- t(combineddata)
colnames(combineddata)= c("L1HS", "mutations")
logcombineddata<- log2(combineddata)
cor(logcombineddata[,2], logcombineddata[,1])
value<- lm(logcombineddata[,1] ~ logcombineddata[,2])
plot(logcombineddata[,2], logcombineddata[,1], xlab="log average mutation frequency", ylab="log fold change L1HS", main="New L1HS Log", col="magenta")
abline(value)
#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")
r2<- summary(value)$r.squared



#median
#BLCA and BRCA
m1<- merge(BLCAmerge, BRCAmerge)
BLCABRCAmed<- data.frame(apply(m1, 2, median, na.rm=TRUE))
colnames(BLCABRCAmed)=c("mut")
BLCAmed<- data.frame(apply(BLCAmerge, 2, median, na.rm=TRUE))
colnames(BLCAmed)= c("mut")
BRCAmed<- data.frame(apply(BRCAmerge, 2, median, na.rm=TRUE))
#CESC
CESCmed<- data.frame(apply(CESCmerge, 2, median, na.rm=TRUE))
colnames(CESCmed)=c("mut")
#CHOL
CHOLmed<- data.frame(apply(CHOLmerge, 2, median, na.rm=TRUE))
colnames(CHOLmed)=c("mut")
#COAD
COADmed<- data.frame(apply(COADmerge, 2, median, na.rm=TRUE))
colnames(COADmed)=c("mut")
#ESCA
ESCAmed<- data.frame(apply(ESCAmerge, 2, median, na.rm=TRUE))
colnames(ESCAmed)=c("mut")
#HNSCm4
HNSCmed<- data.frame(apply(HNSCm4, 2, median, na.rm=TRUE))
colnames(HNSCmed)=c("mut")
#KICH
KICHmed<- data.frame(apply(KICHmerge, 2, median, na.rm=TRUE))
colnames(KICHmed)=c("mut")
#KIRC
KIRCmed<- data.frame(apply(KIRCmerge, 2, median, na.rm=TRUE))
colnames(KIRCmed)=c("mut")
#KIRP
KIRPmed<- data.frame(apply(KIRPm6, 2, median, na.rm=TRUE))
colnames(KIRPmed)=c("mut")
#LIHC
LIHCmed<- data.frame(apply(LIHCm5, 2, median, na.rm=TRUE))
colnames(LIHCmed)=c("mut")
#LUAD
LUADmed<- data.frame(apply(LUADmerge, 2, median, na.rm=TRUE))
colnames(LUADmed)=c("mut")
#PAAD
PAADmed<- data.frame(apply(PAADm5, 2, median, na.rm=TRUE))
colnames(PAADmed)=c("mut")
#PCPG
PCPGmed<- data.frame(apply(PCPGm3, 2, median, na.rm=TRUE))
colnames(PCPGmed)=c("mut")
#PRAD
PRADmed<- data.frame(apply(PRADm4, 2, median, na.rm=TRUE))
colnames(PRADmed)=c("mut")
#SARCm3
SARCmed<- data.frame(apply(SARCm3, 2, median, na.rm=TRUE))
colnames(SARCmed)=c("mut")
#THCA
THCAmed<- data.frame(apply(THCAm4, 2, median, na.rm=TRUE))
colnames(THCAmed)=c("mut")
#THYM
THYMmed<- data.frame(apply(THYMm4, 2, median, na.rm=TRUE))
colnames(THYMmed)=c("mut")
#UCEC
UCECmed<- data.frame(apply(UCECm3, 2, median, na.rm=TRUE))
colnames(UCECmed)=c("mut")

#time to rbind!
mergemut1<- rbind(BLCABRCAmed, CESCmed)
mergemut2<- rbind(mergemut1, CHOLmed)
mergemut3<- rbind(mergemut2, ESCAmed)
mergemut5<- rbind(mergemut3, HNSCmed)
mergemut6<- rbind(mergemut5, KIRPmed)
mergemut7<- rbind(mergemut6, LIHCmed)
mergemut8<- rbind(mergemut7, PAADmed)
mergemut9<- rbind(mergemut8, PCPGmed)
mergemut10<- rbind(mergemut9, PRADmed)
mergemut11<- rbind(mergemut10, SARCmed)
mergemut12<- rbind(mergemut11, THCAmed)
mergemut13<- rbind(mergemut12, THYMmed)
mergemut14<- rbind(mergemut13, UCECmed)
mergemut15<- rbind(mergemut14, COADmed)
mergemut16<- rbind(mergemut15, KICHmed)
mergemut17<- rbind(mergemut16, KIRCmed)
mergemut18<- rbind(mergemut17, LUADmed)
mergemut18<- data.frame(mergemut18)
#match L1HS
tmergemut18<- t(mergemut18)
matchedL<- match(rownames(mergemut18), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HSmatched<- data.frame(tratioL1HS[,matchedL])
tL1HSmatched<- t(L1HSmatched)
combineddata<- merge(tmergemut18, tL1HSmatched, all=TRUE)
combineddata<- t(combineddata)
colnames(combineddata)= c("L1HS", "mutations")
logcombineddata<- log2(combineddata)
cor(logcombineddata[,2], logcombineddata[,1])
value<- lm(logcombineddata[,1] ~ logcombineddata[,2])
plot(logcombineddata[,2], logcombineddata[,1], xlab="log median mutation frequency", ylab="log fold change L1HS", main="New L1HS Log", col="red")
abline(value)
#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")
r2<- summary(value)$r.squared
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
write.table(logcombineddata, file="logcombinedfile.txt", sep="\t")
L1HS<- data.frame(logcombineddata[,1])
Mutations<- data.frame(logcombineddata[,2])
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
write.table(L1HS, file="L1HSlog.txt", sep="\t", row.names=TRUE)
write.table(Mutations, file="Mutations.txt", sep="\t", row.names=TRUE)

#Remove Outlier!
#read in the file
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
oldmutation<- read.table("Mutations.txt", header=TRUE, row.names=1)
toldmutation<- t(oldmutation)
newmutation<- subset(toldmutation, select = -c(A0DB) )
tnewmutation<- t(newmutation)
#Upload New L1HS Data
#match with L1HS data
matchedL<- match(rownames(tnewmutation), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HSmatched<- data.frame(tratioL1HS[,matchedL])
tL1HSmatched<- t(L1HSmatched)
logL1HS<- log2(L1HSmatched)
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
write.table(logL1HS, file="L1HSnoout.txt", sep="\t", row.names=TRUE, quote=FALSE)
write.table(tnewmutation, file="Mutationsnoout.txt", sep="\t", row.names=TRUE, quote=FALSE)

#create cancer type covariate
BLCAnames<- colnames(BLCAmerge)
BLCA<- c("BLCA")
BLCApatients<- cbind(BLCAnames, BLCA)
BRCAnames<- colnames(BRCAmerge)
BRCA<- c("BRCA")
BRCApatients<- cbind(BRCAnames, BRCA)
CESCnames<- colnames(CESCmerge)
CESC<- c("CESC")
CESCpatients<- cbind(CESCnames, CESC)
CHOLnames<- colnames(CHOLmerge)
CHOL<- c("CHOL")
CHOLpatients<- cbind(CHOLnames, CHOL)
ESCAnames<- colnames(ESCAmerge)
ESCA<- c("ESCA")
ESCApatients<- cbind(ESCAnames, ESCA)
HNSCnames<- colnames(HNSCm4)
HNSC<- c("HNSC")
HNSCpatients<- cbind(HNSCnames, HNSC)
KIRPnames<- colnames(KIRPm6)
KIRP<- c("KIRP")
KIRPpatients<- cbind(KIRPnames, KIRP)
LIHCnames<- colnames(LIHCm5)
LIHC<- c("LIHC")
LIHCpatients<- cbind(LIHCnames, LIHC)
PAADnames<- colnames(PAADm5)
PAAD<- c("PAAD")
PAADpatients<- cbind(PAADnames, PAAD)
PCPGnames<- colnames(PCPGm3)
PCPG<- c("PCPG")
PCPGpatients<- cbind(PCPGnames, PCPG)
PRADnames<- colnames(PRADm4)
PRAD<- c("PRAD")
PRADpatients<- cbind(PRADnames, PRAD)
SARCnames<- colnames(SARCm3)
SARC<- c("SARC")
SARCpatients<- cbind(SARCnames, SARC)
THCAnames<- colnames(THCAm4)
THCA<- c("THCA")
THCApatients<- cbind(THCAnames, THCA)
THYMnames<- colnames(THYMm4)
THYM<- c("THYM")
THYMpatients<- cbind(THYMnames, THYM)
UCECnames<- colnames(UCECm3)
UCEC<- c("UCEC")
UCECpatients<- cbind(UCECnames, UCEC)
COADnames<- colnames(COADmerge)
COAD<- c("COAD")
COADpatients<- cbind(COADnames, COAD)
KICHnames<- colnames(KICHmerge)
KICH<- c("KICH")
KICHpatients<- cbind(KICHnames, KICH)
KIRCnames<- colnames(KIRCmerge)
KIRC<- c("KIRC")
KIRCpatients<- cbind(KIRCnames, KIRC)
LUADnames<- colnames(LUADmerge)
LUAD<- c("LUAD")
LUADpatients<- cbind(LUADnames, LUAD)
cancertype<- rbind(BLCApatients, BRCApatients, CESCpatients, CHOLpatients, ESCApatients, HNSCpatients, KIRPpatients, LIHCpatients, PAADpatients, PCPGpatients, PRADpatients, SARCpatients, THCApatients, THYMpatients, UCECpatients, COADpatients, KICHpatients, KIRCpatients, LUADpatients)
rownames(cancertype)= cancertype[,1]
tcancertype<- t(cancertype)
logL1Hsmatched<- log2(L1HSmatched)
L1HS<- cbind(rownames(logL1HSmatched), logL1HSmatched)
cancertest<- cbind(cancertype, L1HS)
cancerfile<- cancertest[-c(1,3)]
tcancerfile<- t(cancerfile)
tcancerfile2<- subset(tcancerfile, select = -c(A0DB) )
cancerfile2<- t(tcancerfile2)
colnames(cancerfile)= c("CancerType", "TEInsertions")
colnames(cancerfile2)= c("CancerType", "TEInsertions")
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
write.table(cancerfile, file="cancertype.txt", sep="\t", row.names=TRUE, quote=FALSE)
write.table(cancerfile2, file="cancertypenoout.txt", sep="\t", row.names=TRUE, quote=FALSE)
