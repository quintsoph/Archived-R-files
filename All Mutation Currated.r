#ignore the warning
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/L1HS expression")
L1HSnormal<- read.table(file="L1HS.Normal.Final", header=FALSE)
L1HScancer<- read.table(file="L1HS.Cancer.Final", header=FALSE)
#create ratio of cancer/normal
ratioL1HS<- data.frame(L1HScancer[,2]/L1HSnormal[,2])
rownames(ratioL1HS)=L1HScancer[,1]
colnames(ratioL1HS)<- c("ratio")
tratioL1HS<- t(ratioL1HS)
#BLCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/BLCA")
BLCA1<- read.table(file="BLCA130_somatic_updated.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
BLCAnames1<- sub(pattern=".01", replacement="", x=colnames(BLCA1), perl=TRUE)
colnames(BLCA1)=BLCAnames1

#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/BRCA")
BRCA1<- data.frame(read.table(file="genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.1.1.0.curated.somatic.TOTALCOUNT.tsv", header=TRUE))
BRCAnames1<- sub(pattern=".01", replacement="", x=colnames(BRCA1), perl=TRUE)
colnames(BRCA1)=BRCAnames1

#CESC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/CESC")
CESC1<- read.table(file="genome.wustl.edu_CESC.IlluminaGA_DNASeq_curated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
CESCnames1<- sub(pattern=".01", replacement="", x=colnames(CESC1), perl=TRUE)
colnames(CESC1)=CESCnames1

#FPPP
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/FPPP")
FPPP<- read.table(file="hgsc.bcm.edu_FPPP.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
FPPPnames1<- sub(pattern=".01", replacement="", x=colnames(FPPP), perl=TRUE)
colnames(FPPP)=FPPPnames1
matched<- match(FPPPnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
FPPPfinal1<- data.frame(tratioL1HS[,matched])
FPPPcut1<- FPPP[,rownames(FPPPfinal1)]

#HNSC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/HNSC")
HNSC1<- read.table(file="pair_set_279_freeze_Mar262013.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSCnames1<- sub(pattern=".01", replacement="", x=colnames(HNSC1), perl=TRUE)
colnames(HNSC1)=HNSCnames1

#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/KIRP")
KIRP1<- read.table(file="An_TCGA_KIRP_MultiCenterCalling_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTA", header=TRUE)
KIRPnames1<- sub(pattern=".01", replacement="", x=colnames(KIRP1), perl=TRUE)
colnames(KIRP1)=KIRPnames1

#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LIHC")
LIHC1<- read.table(file="An_TCGA_LIHC_External_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHCnames1<- sub(pattern=".01", replacement="", x=colnames(LIHC1), perl=TRUE)
colnames(LIHC1)=LIHCnames1

#PAAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PAAD")
PAAD1<- read.table(file="freeze3.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
PAADnames1<- sub(pattern=".01", replacement="", x=colnames(PAAD1), perl=TRUE)
colnames(PAAD1)=PAADnames1

#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PRAD")
PRAD1<- read.table(file="PRAD_Capture_All_Pairs_QCPASS_v6_Nikki_Nov_25.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOU", header=TRUE)
PRADnames1<- sub(pattern=".01", replacement="", x=colnames(PRAD1), perl=TRUE)
colnames(PRAD1)=PRADnames1

#THYM
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/THYM")
THYM1<- read.table(file="genome.wustl.edu_THYM.IlluminaGA_DNASeq_curated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
THYMnames1<- sub(pattern=".01", replacement="", x=colnames(THYM1), perl=TRUE)
colnames(THYM1)=THYMnames1

#combind
allcur<- cbind(BLCA1, BRCA1, CESC1, HNSC1, KIRP1, LIHC1, PAAD1, PRAD1, THYM1)
allcurt<- t(allcur)
matchedL<- match(rownames(allcurt), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HS<- data.frame(tratioL1HS[,matchedL])
tL1HS<- t(L1HS)
matchedcur<- match(rownames(L1HS), rownames(allcurt))
matchedcur<- matchedcur[!is.na(matchedcur)]
cutcur<- data.frame(allcur[,matchedcur])
cutcurt<- t(cutcur)
logcutcurt<- log2(cutcurt)
logL1HS<- log2(L1HS)
L1HSt<- t(L1HS)
finaldata<- merge(L1HSt, cutcur, all=TRUE)
tfinaldata<- t(finaldata)
colnames(tfinaldata)= c("L1HS", "mutations")
logfinal<- log2(tfinaldata)

#covaraiate
#BLCA
tBLCA1<- t(BLCA1)
matchedBLCA<- match(rownames(tBLCA1), rownames(cutcurt))
matchedBLCA<- matchedBLCA[!is.na(matchedBLCA)]
BLCAcut1<- data.frame(cutcur[,matchedBLCA])

#BRCA
tBRCA1<- t(BRCA1)
matchedBRCA<- match(rownames(tBRCA1), rownames(cutcurt))
matchedBRCA<- matchedBRCA[!is.na(matchedBRCA)]
BRCAcut1<- data.frame(cutcur[,matchedBRCA])

#CESC
tCESC1<- t(CESC1)
matchedCESC<- match(rownames(tCESC1), rownames(cutcurt))
matchedCESC<- matchedCESC[!is.na(matchedCESC)]
CESCcut1<- data.frame(cutcur[,matchedCESC])
#HNSC
tHNSC1<- t(HNSC1)
matchedHNSC<- match(rownames(tHNSC1), rownames(cutcurt))
matchedHNSC<- matchedHNSC[!is.na(matchedHNSC)]
HNSCcut1<- data.frame(cutcur[,matchedHNSC])
#KIRP
tKIRP1<- t(KIRP1)
matchedKIRP<- match(rownames(tKIRP1), rownames(cutcurt))
matchedKIRP<- matchedKIRP[!is.na(matchedKIRP)]
KIRPcut1<- data.frame(cutcur[,matchedKIRP])
#LIHC
tLIHC1<- t(LIHC1)
matchedLIHC<- match(rownames(tLIHC1), rownames(cutcurt))
matchedLIHC<- matchedLIHC[!is.na(matchedLIHC)]
LIHCcut1<- data.frame(cutcur[,matchedLIHC])
#PAAD
tPAAD1<- t(PAAD1)
matchedPAAD<- match(rownames(tPAAD1), rownames(cutcurt))
matchedPAAD<- matchedPAAD[!is.na(matchedPAAD)]
PAADcut1<- data.frame(cutcur[,matchedPAAD])
#PRAD
tPRAD1<- t(PRAD1)
matchedPRAD<- match(rownames(tPRAD1), rownames(cutcurt))
matchedPRAD<- matchedPRAD[!is.na(matchedPRAD)]
PRADcut1<- data.frame(cutcur[,matchedPRAD])
#THYM
tTHYM1<- t(THYM1)
matchedTHYM<- match(rownames(tTHYM1), rownames(cutcurt))
matchedTHYM<- matchedTHYM[!is.na(matchedTHYM)]
THYMcut1<- data.frame(cutcur[,matchedTHYM])

m1<- merge(tL1HS, BLCAcut1, all=TRUE)
m2<- merge(m1, BRCAcut1, all=TRUE)
m3<- merge(m2, CESCcut1, all=TRUE)
m4<- merge(m3, HNSCcut1, all=TRUE)
m5<- merge(m4, KIRPcut1, all=TRUE)
m6<- merge(m5, LIHCcut1, all=TRUE)
m7<- merge(m6, PAADcut1, all=TRUE)
m8<- merge(m7, PRADcut1, all=TRUE)
m9<- merge(m8, THYMcut1, all=TRUE)
rownames(m9)= c("L1HS", "THYM", "PAAD", "LIHC", "KIRP", "CESC", "BRCA", "BLCA")
tm9<- t(m9)
logtm9<- log2(tm9)
test<- lm(tm9[,1] ~ tm9[,2] + tm9[,3] + tm9[,4] + tm9[,5] + tm9[,6] + tm9[,7] + tm9[,8], na.action= na.exclude)
summary(test)
set<- logtm9[,-1]

#loop
test=NULL
for (i in colnames(set))
{
	l<- lm(logtm9[,1] ~ set[,i])
	test<- cbind(test, l)
}

#graph for curated only
#graph data points
plot(logfinal[,2], logfinal[,1], xlab="curated mutations", ylab="L1HS levels", main="curated mutations", col="darkgreen")
value<- lm(logfinal[,1] ~ logfinal[,2])
abline(value)
cor(logfinal[,2], logfinal[,1])
p<- anova(value)$Pr[1]
r2<- summary(value)$r.squared

# addition of other cancers
#CHOL
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/CHOL")
CHOL1<- read.table(file="bcgsc.ca_CHOL.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
CHOL2<- read.table(file="CHOL_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
CHOL3<- read.table(file="hgsc.bcm.edu_CHOL.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
CHOL4<- read.table(file="hgsc.bcm.edu_CHOL.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
CHOL5<- read.table(file="ucsc.edu_CHOL.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
CHOLnames1<- sub(pattern=".01", replacement="", x=colnames(CHOL1), perl=TRUE)
CHOLnames2<- sub(pattern=".01", replacement="", x=colnames(CHOL2), perl=TRUE)
CHOLnames3<- sub(pattern=".01", replacement="", x=colnames(CHOL3), perl=TRUE)
CHOLnames4<- sub(pattern=".01", replacement="", x=colnames(CHOL4), perl=TRUE)
CHOLnames5<- sub(pattern=".01", replacement="", x=colnames(CHOL5), perl=TRUE)
colnames(CHOL1)=CHOLnames1
colnames(CHOL2)=CHOLnames2
colnames(CHOL3)=CHOLnames3
colnames(CHOL4)=CHOLnames4
colnames(CHOL5)=CHOLnames5
matched<- match(CHOLnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
CHOLfinal1<- data.frame(tratioL1HS[,matched])
CHOLcut1<- CHOL1[,rownames(CHOLfinal1)]
matched2<- match(CHOLnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
CHOLfinal2<- data.frame(tratioL1HS[,matched2])
CHOLcut2<- CHOL2[,rownames(CHOLfinal2)]
matched3<- match(CHOLnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
CHOLfinal3<- data.frame(tratioL1HS[,matched3])
CHOLcut3<- CHOL3[,rownames(CHOLfinal3)]
matched4<- match(CHOLnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
CHOLfinal4<- data.frame(tratioL1HS[,matched4])
CHOLcut4<- CHOL4[,rownames(CHOLfinal4)]
matched5<- match(CHOLnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
CHOLfinal5<- data.frame(tratioL1HS[,matched5])
CHOLcut5<- CHOL5[,rownames(CHOLfinal5)]
CHOLm1<- merge(CHOLcut5, CHOLcut4, all=TRUE)
CHOLm2<- merge(CHOLm1, CHOLcut3, all=TRUE)
CHOLm3<- merge(CHOLm2, CHOLcut2, all=TRUE)
CHOLm4<- merge(CHOLm3, CHOLcut1, all=TRUE)
#CHOL
notna<- colSums(!is.na(CHOLm4))
sum<- colSums(CHOLm4, na.rm=TRUE)
avmutCHOL<- data.frame(sum/notna)
#ESCA
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/ESCA")
ESCA1<- read.table(file="An_TCGA_ESCA_External_capture_All_Pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.t", header=TRUE)
ESCA2<- read.table(file="bcgsc.ca_ESCA.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCA3<- read.table(file="genome.wustl.edu_ESCA.IlluminaHiSeq_DNASeq_automated.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCA4<- read.table(file="hgsc.bcm.edu_ESCA.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCA5<- read.table(file="ucsc.edu_ESCA.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCAnames1<- sub(pattern=".01", replacement="", x=colnames(ESCA1), perl=TRUE)
ESCAnames2<- sub(pattern=".01", replacement="", x=colnames(ESCA2), perl=TRUE)
ESCAnames3<- sub(pattern=".01", replacement="", x=colnames(ESCA3), perl=TRUE)
ESCAnames4<- sub(pattern=".01", replacement="", x=colnames(ESCA4), perl=TRUE)
ESCAnames5<- sub(pattern=".01", replacement="", x=colnames(ESCA5), perl=TRUE)
colnames(ESCA1)=ESCAnames1
colnames(ESCA2)=ESCAnames2
colnames(ESCA3)=ESCAnames3
colnames(ESCA4)=ESCAnames4
colnames(ESCA5)=ESCAnames5
matched<- match(ESCAnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
ESCAfinal1<- data.frame(tratioL1HS[,matched])
ESCAcut1<- ESCA1[,rownames(ESCAfinal1)]
matched2<- match(ESCAnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
ESCAfinal2<- data.frame(tratioL1HS[,matched2])
ESCAcut2<- ESCA2[,rownames(ESCAfinal2)]
matched3<- match(ESCAnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
ESCAfinal3<- data.frame(tratioL1HS[,matched3])
ESCAcut3<- ESCA3[,rownames(ESCAfinal3)]
matched4<- match(ESCAnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
ESCAfinal4<- data.frame(tratioL1HS[,matched4])
ESCAcut4<- ESCA4[,rownames(ESCAfinal4)]
matched5<- match(ESCAnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
ESCAfinal5<- data.frame(tratioL1HS[,matched5])
ESCAcut5<- ESCA5[,rownames(ESCAfinal5)]
ESCAm1<- merge(ESCAcut5, ESCAcut4, all=TRUE)
ESCAm2<- merge(ESCAm1, ESCAcut3, all=TRUE)
ESCAm3<- merge(ESCAm2, ESCAcut2, all=TRUE)
ESCAm4<- merge(ESCAm3, ESCAcut1, all=TRUE)
notna<- colSums(!is.na(ESCAm4))
sum<- colSums(ESCAm4, na.rm=TRUE)
avmutESCA<- data.frame(sum/notna)
#PCPG
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PCPG")
PCPG1<- read.table(file="bcgsc.ca_PCPG.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPG2<- read.table(file="hgsc.bcm.edu_PCPG.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPG3<- read.table(file="PCPG_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPG4<- read.table(file="ucsc.edu_PCPG.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPGnames1<- sub(pattern=".01", replacement="", x=colnames(PCPG1), perl=TRUE)
PCPGnames2<- sub(pattern=".01", replacement="", x=colnames(PCPG2), perl=TRUE)
PCPGnames3<- sub(pattern=".01", replacement="", x=colnames(PCPG3), perl=TRUE)
PCPGnames4<- sub(pattern=".01", replacement="", x=colnames(PCPG4), perl=TRUE)
colnames(PCPG1)=PCPGnames1
colnames(PCPG2)=PCPGnames2
colnames(PCPG3)=PCPGnames3
colnames(PCPG4)=PCPGnames4
matched<- match(PCPGnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PCPGfinal1<- data.frame(tratioL1HS[,matched])
PCPGcut1<- PCPG1[,rownames(PCPGfinal1)]
matched2<- match(PCPGnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
PCPGfinal2<- data.frame(tratioL1HS[,matched2])
PCPGcut2<- PCPG2[,rownames(PCPGfinal2)]
matched3<- match(PCPGnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
PCPGfinal3<- data.frame(tratioL1HS[,matched3])
PCPGcut3<- PCPG3[,rownames(PCPGfinal3)]
matched4<- match(PCPGnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
PCPGfinal4<- data.frame(tratioL1HS[,matched4])
PCPGcut4<- PCPG4[,rownames(PCPGfinal4)]
PCPGm1<- merge(PCPGcut4, PCPGcut3, all=TRUE)
PCPGm2<- merge(PCPGm1, PCPGcut2, all=TRUE)
PCPGm3<- merge(PCPGm2, PCPGcut1, all=TRUE)
notna<- colSums(!is.na(PCPGm3))
sum<- colSums(PCPGm3, na.rm=TRUE)
avmutPCPG<- data.frame(sum/notna)
#SARC
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/SARC")
SARC1<- read.table(file="bcgsc.ca_SARC.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
SARC2<- read.table(file="genome.wustl.edu_SARC.IlluminaHiSeq_DNASeq_automated.1.5.0.somatic.TOTALCOUNT.tsv", header=TRUE)
SARC3<- read.table(file="SARC_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
SARC4<- read.table(file="ucsc.edu_SARC.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
SARCnames1<- sub(pattern=".01", replacement="", x=colnames(SARC1), perl=TRUE)
SARCnames2<- sub(pattern=".01", replacement="", x=colnames(SARC2), perl=TRUE)
SARCnames3<- sub(pattern=".01", replacement="", x=colnames(SARC3), perl=TRUE)
SARCnames4<- sub(pattern=".01", replacement="", x=colnames(SARC4), perl=TRUE)
colnames(SARC1)=SARCnames1
colnames(SARC2)=SARCnames2
colnames(SARC3)=SARCnames3
colnames(SARC4)=SARCnames4
matched<- match(SARCnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
SARCfinal1<- data.frame(tratioL1HS[,matched])
SARCcut1<- SARC1[,rownames(SARCfinal1)]
matched2<- match(SARCnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
SARCfinal2<- data.frame(tratioL1HS[,matched2])
SARCcut2<- SARC2[,rownames(SARCfinal2)]
matched3<- match(SARCnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
SARCfinal3<- data.frame(tratioL1HS[,matched3])
SARCcut3<- SARC3[,rownames(SARCfinal3)]
matched4<- match(SARCnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
SARCfinal4<- data.frame(tratioL1HS[,matched4])
SARCcut4<- SARC4[,rownames(SARCfinal4)]
SARCm1<- merge(SARCcut4, SARCcut3, all=TRUE)
SARCm2<- merge(SARCm1, SARCcut2, all=TRUE)
SARCm3<- merge(SARCm2, SARCcut1, all=TRUE)
notna<- colSums(!is.na(SARCm3))
sum<- colSums(SARCm3, na.rm=TRUE)
avmutSARC<- data.frame(sum/notna)
#THCA
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/THCA")
THCA1<- read.table(file="AN_TCGA_THCA_PAIR_Capture_ALLQC_14Aug2013_429.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA2<- read.table(file="bcgsc.ca_THCA.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA3<- read.table(file="hgsc.bcm.edu_THCA.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA4<- read.table(file="PR_TCGA_THCA_PAIR_Capture_All_Pairs_QCPASS_v2.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA5<- read.table(file="THCA_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
THCAnames1<- sub(pattern=".01", replacement="", x=colnames(THCA1), perl=TRUE)
THCAnames2<- sub(pattern=".01", replacement="", x=colnames(THCA2), perl=TRUE)
THCAnames3<- sub(pattern=".01", replacement="", x=colnames(THCA3), perl=TRUE)
THCAnames4<- sub(pattern=".01", replacement="", x=colnames(THCA4), perl=TRUE)
THCAnames5<- sub(pattern=".01", replacement="", x=colnames(THCA5), perl=TRUE)
colnames(THCA1)=THCAnames1
colnames(THCA2)=THCAnames2
colnames(THCA3)=THCAnames3
colnames(THCA4)=THCAnames4
colnames(THCA5)=THCAnames5
matched<- match(THCAnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
THCAfinal1<- data.frame(tratioL1HS[,matched])
THCAcut1<- THCA1[,rownames(THCAfinal1)]
matched2<- match(THCAnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
THCAfinal2<- data.frame(tratioL1HS[,matched2])
THCAcut2<- THCA2[,rownames(THCAfinal2)]
matched3<- match(THCAnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
THCAfinal3<- data.frame(tratioL1HS[,matched3])
THCAcut3<- THCA3[,rownames(THCAfinal3)]
matched4<- match(THCAnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
THCAfinal4<- data.frame(tratioL1HS[,matched4])
THCAcut4<- THCA4[,rownames(THCAfinal4)]
matched5<- match(THCAnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
THCAfinal5<- data.frame(tratioL1HS[,matched5])
THCAcut5<- THCA5[,rownames(THCAfinal5)]
THCAm1<- merge(THCAcut5, THCAcut4, all=TRUE)
THCAm2<- merge(THCAm1, THCAcut3, all=TRUE)
THCAm3<- merge(THCAm2, THCAcut2, all=TRUE)
THCAm4<- merge(THCAm3, THCAcut1, all=TRUE)
notna<- colSums(!is.na(THCAm4))
sum<- colSums(THCAm4, na.rm=TRUE)
avmutTHCA<- data.frame(sum/notna)

#Covariate cancer type





















#Check
#ignore the warning
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/L1HS expression")
L1HSnormal<- read.table(file="L1HS.Normal.Final", header=FALSE)
L1HScancer<- read.table(file="L1HS.Cancer.Final", header=FALSE)
#create ratio of cancer/normal
ratioL1HS<- data.frame(L1HScancer[,2]/L1HSnormal[,2])
rownames(ratioL1HS)=L1HScancer[,1]
colnames(ratioL1HS)<- c("ratio")
tratioL1HS<- t(ratioL1HS)
#BLCA
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/BLCA")
BLCA1<- read.table(file="BLCA130_somatic_updated.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
BLCAnames1<- sub(pattern=".01", replacement="", x=colnames(BLCA1), perl=TRUE)
colnames(BLCA1)=BLCAnames1
matched<- match(BLCAnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
BLCAfinal1<- data.frame(tratioL1HS[,matched])
BLCAcut1<- BLCA1[,rownames(BLCAfinal1)]
#BRCA
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/BRCA")
BRCA1<- data.frame(read.table(file="genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.1.1.0.curated.somatic.TOTALCOUNT.tsv", header=TRUE))
BRCAnames1<- sub(pattern=".01", replacement="", x=colnames(BRCA1), perl=TRUE)
colnames(BRCA1)=BRCAnames1
matched<- match(BRCAnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
BRCAfinal1<- data.frame(tratioL1HS[,matched])
BRCAcut1<- BRCA1[,rownames(BRCAfinal1)]
#CESC
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/CESC")
CESC1<- read.table(file="genome.wustl.edu_CESC.IlluminaGA_DNASeq_curated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
CESCnames1<- sub(pattern=".01", replacement="", x=colnames(CESC1), perl=TRUE)
colnames(CESC1)=CESCnames1
matched<- match(CESCnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
CESCfinal1<- data.frame(tratioL1HS[,matched])
CESCcut1<- CESC1[,rownames(CESCfinal1)]
#FPPP
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/FPPP")
FPPP<- read.table(file="hgsc.bcm.edu_FPPP.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
FPPPnames1<- sub(pattern=".01", replacement="", x=colnames(FPPP), perl=TRUE)
colnames(FPPP)=FPPPnames1
matched<- match(FPPPnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
FPPPfinal1<- data.frame(tratioL1HS[,matched])
FPPPcut1<- FPPP[,rownames(FPPPfinal1)]
#HNSC
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/HNSC")
HNSC1<- read.table(file="pair_set_279_freeze_Mar262013.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSCnames1<- sub(pattern=".01", replacement="", x=colnames(HNSC1), perl=TRUE)
colnames(HNSC1)=HNSCnames1
matched<- match(HNSCnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
HNSCfinal1<- data.frame(tratioL1HS[,matched])
HNSCcut1<- HNSC1[,rownames(HNSCfinal1)]
#KIRP
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/KIRP")
KIRP1<- read.table(file="An_TCGA_KIRP_MultiCenterCalling_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTA", header=TRUE)
KIRPnames1<- sub(pattern=".01", replacement="", x=colnames(KIRP1), perl=TRUE)
colnames(KIRP1)=KIRPnames1
matched<- match(KIRPnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
KIRPfinal1<- data.frame(tratioL1HS[,matched])
KIRPcut1<- KIRP1[,rownames(KIRPfinal1)]
#LIHC
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LIHC")
LIHC1<- read.table(file="An_TCGA_LIHC_External_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHCnames1<- sub(pattern=".01", replacement="", x=colnames(LIHC1), perl=TRUE)
colnames(LIHC1)=LIHCnames1
matched<- match(LIHCnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
LIHCfinal1<- data.frame(tratioL1HS[,matched])
LIHCcut1<- LIHC1[,rownames(LIHCfinal1)]
tLIHCcut1<- t(LIHCcut1)
#PAAD
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PAAD")
PAAD1<- read.table(file="freeze3.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
PAADnames1<- sub(pattern=".01", replacement="", x=colnames(PAAD1), perl=TRUE)
colnames(PAAD1)=PAADnames1
matched<- match(PAADnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PAADfinal1<- data.frame(tratioL1HS[,matched])
PAADcut1<- PAAD1[,rownames(PAADfinal1)]
#PRAD
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PRAD")
PRAD1<- read.table(file="PRAD_Capture_All_Pairs_QCPASS_v6_Nikki_Nov_25.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOU", header=TRUE)
PRADnames1<- sub(pattern=".01", replacement="", x=colnames(PRAD1), perl=TRUE)
colnames(PRAD1)=PRADnames1
matched<- match(PRADnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PRADfinal1<- data.frame(tratioL1HS[,matched])
PRADcut1<- PRAD1[,rownames(PRADfinal1)]
#THYM
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/THYM")
THYM1<- read.table(file="genome.wustl.edu_THYM.IlluminaGA_DNASeq_curated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
THYMnames1<- sub(pattern=".01", replacement="", x=colnames(THYM1), perl=TRUE)
colnames(THYM1)=THYMnames1
matched<- match(THYMnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
THYMfinal1<- data.frame(tratioL1HS[,matched])
THYMcut1<- THYM1[,rownames(THYMfinal1)]
#merge
mutcheck<- cbind(BLCAcut1, BRCAcut1, CESCcut1, HNSCcut1, KIRPcut1, LIHCcut1, PAADcut1, PRADcut1, THYMcut1)
tmutcheck<- t(mutcheck)
match(rownames(tmutcheck), rownames(cutcurt))


#merge
mutcheck<- cbind(BLCAcut1, BRCAcut1, CESCcut1, HNSCcut1, KIRPcut1, LIHCcut1, PAADcut1, PRADcut1, THYMcut1)
tmutcheck<- t(mutcheck)
match(rownames(tmutcheck), rownames(cutcurt))
#combind
merge1<- rbind(avmutTHCA, avmutSARC)
merge2<- rbind(merge1, avmutPCPG)
merge3<- rbind(merge2, avmutESCA)
merge4<- rbind(merge3, avmutCHOL)
tFPPP<- t(FPPPcut1)
colnames(tFPPP)= c("sum.notna")
merge5<- rbind(merge4, tFPPP)
colnames(tmutcheck)= c("sum.notna")
ALLmut<- rbind(tmutcheck, merge5)
ALLmutt<- t(ALLmut)
matchedL<- match(rownames(ALLmut), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HS2<- data.frame(tratioL1HS[,matchedL])
matchedcur<- match(rownames(L1HS2), rownames(ALLmut))
matchedcur<- matchedcur[!is.na(matchedcur)]
cutALL<- data.frame(ALLmutt[,matchedcur])
L1HSt2<- t(L1HS2)
tcutALL<- t(cutALL)
finalalldata<- merge(L1HSt2, tcutALL, all=TRUE)
tfinalalldata<- t(finalalldata)
colnames(tfinalalldata)= c("L1HS", "mutations")
logfinalall<- log2(tfinalalldata)

#graph
#graph data points
plot(logfinalall[,2], logfinalall[,1], xlab="curated mutations and other cancers", ylab="L1HS levels", main="mutations", col="darkorange")
value<- lm(logfinalall[,1] ~ logfinalall[,2])
abline(value)
cor(logfinalall[,2], logfinalall[,1])
p<- anova(value)$Pr[1]
r2<- summary(value)$r.squared

setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Testing")
write.table(logfinal, file="Curated Version ONLY.txt", sep=" ")
write.table(logfinalall, file="Curated and average mutations.txt", sep=" 0")
write.table(tfinaldata, file="Not Log, Currated Version.txt", sep=" ")

