setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations")
ACC1straw <-read.table(file.choose(), header=TRUE)
ACC2ndraw<-read.table(file.choose(),header=TRUE)
ACC3rdraw<- read.table(file.choose(), header=TRUE)
ACCnames1<- sub(pattern=".01", replacement="", x=colnames(ACC1straw), perl=TRUE)
ACCnames2<- sub(pattern=".01*", replacement="", x=colnames(ACC2ndraw), perl=TRUE)
ACCnames3<- sub(pattern=".01*", replacement="", x=colnames(ACC3rdraw), perl=TRUE)
colnames(ACC1straw)=ACCnames1
colnames(ACC2ndraw)=ACCnames2
colnames(ACC3rdraw)=ACCnames1
#ignore the warning DO NOT USE THIS
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/L1HS expression")
L1HSnormal<- read.table(file="L1HS.Normal.Final", header=FALSE)
L1HScancer<- read.table(file="L1HS.Cancer.Final", header=FALSE)
ratioL1HS<- data.frame(L1HScancer[,2]/L1HSnormal[,2])
rownames(ratioL1HS)=L1HScancer[,1]
colnames(ratioL1HS)<- c("ratio")
tratioL1HS<- t(ratioL1HS)

#ACC matches to L1HS
matched<- match(ACCnames1, rownames(ratioL1HS))
matched2<- match(ACCnames2, rownames(ratioL1HS))
matched3<- match(ACCnames3, rownames(ratioL1HS))
#test... not as many patients in list than in L1HS expressions
names1<-read.table(file.choose())
names<-sub(pattern="-01*", replacement="", x=names1[,1], perl=TRUE)
matchnam<- match(rownames(ratioL1HS), names)
#BLCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations")
BLCA1straw <-data.frame(read.table(file.choose(), header=TRUE))
BLCA2ndraw<- data.frame(read.table(file.choose(),header=TRUE))
BLCA3rdraw<- data.frame(read.table(file.choose(), header=TRUE))
BLCA4rdraw<- data.frame(read.table(file.choose(), header=TRUE))
BLCA5rdraw<- data.frame(read.table(file.choose(), header=TRUE))
BLCAnames1<- sub(pattern=".01", replacement="", x=colnames(BLCA1straw), perl=TRUE)
BLCAnames2<- sub(pattern=".01", replacement="", x=colnames(BLCA2ndraw), perl=TRUE)
BLCAnames3<- sub(pattern=".01", replacement="", x=colnames(BLCA3rdraw), perl=TRUE)
BLCAnames4<- sub(pattern=".01", replacement="", x=colnames(BLCA4rdraw), perl=TRUE)
BLCAnames5<- sub(pattern=".01", replacement="", x=colnames(BLCA5rdraw), perl=TRUE)
colnames(BLCA1straw)=BLCAnames1
colnames(BLCA2ndraw)=BLCAnames2
colnames(BLCA3rdraw)=BLCAnames3
colnames(BLCA4rdraw)=BLCAnames4
colnames(BLCA5rdraw)=BLCAnames5
matched<- match(BLCAnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
BLCAfinal1<- data.frame(tratioL1HS[,matched])
BLCAcut1<- BLCA1straw[,rownames(BLCAfinal1)]
tBLCAcut1<- t(BLCAcut1)
matched2<- match(BLCAnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
BLCAfinal2<- data.frame(tratioL1HS[,matched2])
BLCAcut2<- BLCA2ndraw[,rownames(BLCAfinal2)]
tBLCAcut2<- t(BLCAcut2)
matched3<- match(BLCAnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
BLCAfinal3<- data.frame(tratioL1HS[,matched3])
BLCAcut3<- BLCA3rdraw[,rownames(BLCAfinal3)]
tBLCAcut3<- t(BLCAcut3)
matched4<- match(BLCAnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
BLCAfinal4<- data.frame(tratioL1HS[,matched4])
BLCAcut4<- BLCA4rdraw[,rownames(BLCAfinal4)]
tBLCAcut4<- t(BLCAcut4)
matched5<- match(BLCAnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
BLCAfinal5<- data.frame(tratioL1HS[,matched5])
BLCAcut5<- BLCA5rdraw[,rownames(BLCAfinal5)]
tBLCAcut5<- t(BLCAcut5)
m1<- merge(BLCAcut5, BLCAcut4, all=TRUE)
m2<- merge(m1, BLCAcut3, all=TRUE)
m3<- merge(m2, BLCAcut2, all=TRUE)
BLCAmerge<- merge(m3, BLCAcut1, all=TRUE)

#BRCA
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/BRCA")
BRCA1<- data.frame(read.table(file="genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.1.1.0.curated.somatic.TOTALCOUNT.tsv", header=TRUE))
BRCA2<- data.frame(read.table(file="genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.2.0.0.TOTALCOUNT.tsv", header=TRUE))
BRCA3<- data.frame(read.table(file="genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.TOTALCOUNT.tsv", header=TRUE))
BRCA4<- data.frame(read.table(file="genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.3.0.somatic.TOTALCOUNT.tsv", header=TRUE))
BRCA5<- data.frame(read.table(file="genome.wustl.edu_BRCA.IlluminaHiSeq_DNASeq_automated.1.4.0.somatic.TOTALCOUNT.tsv", header=TRUE))
BRCAnames1<- sub(pattern=".01", replacement="", x=colnames(BRCA1), perl=TRUE)
BRCAnames2<- sub(pattern=".01", replacement="", x=colnames(BRCA2), perl=TRUE)
BRCAnames3<- sub(pattern=".01", replacement="", x=colnames(BRCA3), perl=TRUE)
BRCAnames4<- sub(pattern=".01", replacement="", x=colnames(BRCA4), perl=TRUE)
BRCAnames5<- sub(pattern=".01", replacement="", x=colnames(BRCA5), perl=TRUE)
colnames(BRCA1)=BRCAnames1
colnames(BRCA2)=BRCAnames2
colnames(BRCA3)=BRCAnames3
colnames(BRCA4)=BRCAnames4
colnames(BRCA5)=BRCAnames5
matched<- match(BRCAnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
BRCAfinal1<- data.frame(tratioL1HS[,matched])
BRCAcut1<- BRCA1[,rownames(BRCAfinal1)]
tBRCAcut1<- t(BRCAcut1)
matched2<- match(BRCAnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
BRCAfinal2<- data.frame(tratioL1HS[,matched2])
BRCAcut2<- BRCA2[,rownames(BRCAfinal2)]
tBRCAcut2<- t(BRCAcut2)
matched3<- match(BRCAnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
BRCAfinal3<- data.frame(tratioL1HS[,matched3])
BRCAcut3<- BRCA3[,rownames(BRCAfinal3)]
tBRCAcut3<- t(BRCAcut3)
matched4<- match(BRCAnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
BRCAfinal4<- data.frame(tratioL1HS[,matched4])
BRCAcut4<- BRCA4[,rownames(BRCAfinal4)]
tBRCAcut4<- t(BRCAcut4)
matched5<- match(BRCAnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
BRCAfinal5<- data.frame(tratioL1HS[,matched5])
BRCAcut5<- BRCA5[,rownames(BRCAfinal5)]
tBRCAcut5<- t(BRCAcut5)
m1<- merge(BRCAcut5, BRCAcut4)
m2<- merge(m1, BRCAcut3, all=TRUE)
m3<- merge(m2, BRCAcut2, all=TRUE)
BRCAmerge<- merge(m3, BRCAcut1, all=TRUE)

#CESC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/CESC")
CESC1<- read.table(file="bcgsc.ca_CESC.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
CESC2<- read.table(file="genome.wustl.edu_CESC.IlluminaGA_DNASeq_curated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
CESC3<- read.table(file="genome.wustl.edu_CESC.IlluminaHiSeq_DNASeq_automated.1.3.0.somatic.TOTALCOUNT.tsv", header=TRUE)
CESC4<- read.table(file="PR_TCGA_CESC_PAIR_Capture_All_Pairs.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
CESC5<- read.table(file="PR_TCGA_CESC_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
CESC6<- read.table(file="PR_TCGA_CESC_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.automated.somatic.TOTALC", header=TRUE)
CESC7<- read.table(file="ucsc.edu_CESC.IlluminaGA_DNASeq_automated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
CESCnames1<- sub(pattern=".01", replacement="", x=colnames(CESC1), perl=TRUE)
CESCnames2<- sub(pattern=".01", replacement="", x=colnames(CESC2), perl=TRUE)
CESCnames3<- sub(pattern=".01", replacement="", x=colnames(CESC3), perl=TRUE)
CESCnames4<- sub(pattern=".01", replacement="", x=colnames(CESC4), perl=TRUE)
CESCnames5<- sub(pattern=".01", replacement="", x=colnames(CESC5), perl=TRUE)
CESCnames6<- sub(pattern=".01", replacement="", x=colnames(CESC6), perl=TRUE)
CESCnames7<- sub(pattern=".01", replacement="", x=colnames(CESC7), perl=TRUE)
colnames(CESC1)=CESCnames1
colnames(CESC2)=CESCnames2
colnames(CESC3)=CESCnames3
colnames(CESC4)=CESCnames4
colnames(CESC5)=CESCnames5
colnames(CESC6)=CESCnames6
colnames(CESC7)=CESCnames7
matched<- match(CESCnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
CESCfinal1<- data.frame(tratioL1HS[,matched])
CESCcut1<- CESC1[,rownames(CESCfinal1)]
matched2<- match(CESCnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
CESCfinal2<- data.frame(tratioL1HS[,matched2])
CESCcut2<- CESC2[,rownames(CESCfinal2)]
matched3<- match(CESCnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
CESCfinal3<- data.frame(tratioL1HS[,matched3])
CESCcut3<- CESC3[,rownames(CESCfinal3)]
matched4<- match(CESCnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
CESCfinal4<- data.frame(tratioL1HS[,matched4])
CESCcut4<- CESC4[,rownames(CESCfinal4)]
matched5<- match(CESCnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
CESCfinal5<- data.frame(tratioL1HS[,matched5])
CESCcut5<- CESC5[,rownames(CESCfinal5)]
matched6<- match(CESCnames6, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
CESCfinal6<- data.frame(tratioL1HS[,matched6])
CESCcut6<- CESC6[,rownames(CESCfinal6)]
matched7<- match(CESCnames7, rownames(ratioL1HS))
matched7<- matched7[!is.na(matched7)]
CESCfinal7<- data.frame(tratioL1HS[,matched7])
CESCcut7<- CESC7[,rownames(CESCfinal7)]
m1<- merge(CESCcut1, CESCcut2, all=TRUE)
m2<- merge(m1, CESCcut3, all=TRUE)
m3<- merge(m2, CESCcut4, all=TRUE)
m4<- merge(m3, CESCcut5, all=TRUE)
m5<- merge(m4, CESCcut6, all=TRUE)
CESCmerge<- merge(m5, CESCcut7, all=TRUE)

#CHOL
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/CHOL")
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
m1<- merge(CHOLcut1, CHOLcut2, all=TRUE)
m2<- merge(m1, CHOLcut3, all=TRUE)
m3<- merge(m2, CHOLcut4, all=TRUE)
CHOLmerge<- merge(m3, CHOLcut5, all=TRUE)

#COAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/COAD")
COAD1<- read.table(file="hgsc.bcm.edu_COAD.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
COAD2<- read.table(file="hgsc.bcm.edu_COAD.SOLiD_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
COADnames1<- gsub(pattern="\\.01", replacement="", x=colnames(COAD1), perl=TRUE)
COADnames2<- gsub(pattern="\\.01", replacement="", x=colnames(COAD2), perl=TRUE)
COADnames1.2<- gsub(pattern="^X*", replacement="", x=COADnames1)
COADnames2.2<- gsub(pattern="^X*", replacement="", x=COADnames2)
colnames(COAD1)=COADnames1.2
colnames(COAD2)=COADnames2.2
matched<- match(COADnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
COADfinal1<- data.frame(tratioL1HS[,matched])
COADcut1<- COAD1[,rownames(COADfinal1)]
matched2<- match(COADnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
COADfinal2<- data.frame(tratioL1HS[,matched2])
COADcut2<- COAD2[,rownames(COADfinal2)]
COADmerge<- merge(COADcut1, COADcut2, all=TRUE)

#DLBC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/DLBC")
DLBC1<- read.table(file="hgsc.bcm.edu_DLBC.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
DLBCnames1<- gsub(pattern="\\.01", replacement="", x=colnames(DLBC1), perl=TRUE)
DLBCnames1.2<- gsub(pattern="^X*", replacement="", x=DLBCnames1, perl=TRUE)
colnames(DLBC1)=DLBCnames1.2
matched<- match(DLBCnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
DLBCfinal1<- tratioL1HS[,matched]
#ESCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/ESCA")
ESCA1<- read.table(file="An_TCGA_ESCA_External_capture_All_Pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.t", header=TRUE)
ESCA2<- read.table(file="bcgsc.ca_ESCA.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCA3<- read.table(file="genome.wustl.edu_ESCA.IlluminaHiSeq_DNASeq_automated.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCA4<- read.table(file="hgsc.bcm.edu_ESCA.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCA5<- read.table(file="ucsc.edu_ESCA.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
ESCAnames1<- gsub(pattern="\\.01", replacement="", x=colnames(ESCA1), perl=TRUE)
ESCAnames2<- gsub(pattern="\\.01", replacement="", x=colnames(ESCA2), perl=TRUE)
ESCAnames3<- gsub(pattern="\\.01", replacement="", x=colnames(ESCA3), perl=TRUE)
ESCAnames4<- gsub(pattern="\\.01", replacement="", x=colnames(ESCA4), perl=TRUE)
ESCAnames5<- gsub(pattern="\\.01", replacement="", x=colnames(ESCA5), perl=TRUE)
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
m1<- merge(ESCAcut1, ESCAcut2, all=TRUE)
m2<- merge(m1, ESCAcut3, all=TRUE)
m3<- merge(m2, ESCAcut4, all=TRUE)
ESCAmerge<- merge(m3, ESCAcut5, all=TRUE)

#FPPP
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/FPPP")
FPPP<- read.table(file="hgsc.bcm.edu_FPPP.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
FPPPnames1<- gsub(pattern="\\.01", replacement="", x=colnames(FPPP), perl=TRUE)
FPPPnames1.2<- gsub(pattern="^X*", replacement="", x=FPPPnames1)
colnames(FPPP)=FPPPnames1.2
matched<- match(FPPPnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
FPPPfinal1<- data.frame(tratioL1HS[,matched])
FPPPcut<- FPPP[,rownames(FPPPfinal1)]
#GBM
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/GBM")
GBM1<- read.table(file="gbm_liftover.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
GBM2<- read.table(file="mdanderson.org_GBM.IlluminaGA_DNASeq.Level_2.1.6.somatic.TOTALCOUNT.tsv", header=TRUE)
GBM3<- read.table(file="mdanderson.org_GBM.IlluminaGA_DNASeq.Level_2.1.somatic.TOTALCOUNT.tsv", header=TRUE)
GBM4<- read.table(file="step4_gbm_liftover.aggregated.capture.tcga.uuid.maf2.4.migrated.somatic.TOTALCOUNT.tsv", header=TRUE)
GBM5<- read.table(file="ucsc.edu_GBM.IlluminaGA_DNASeq_automated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
GBMnames1<- gsub(pattern="\\.01", replacement="", x=colnames(GBM1), perl=TRUE)
GBMnames1.2<- gsub(pattern="^X*", replacement="", x=GBMnames1, perl=TRUE)
GBMnames2<- gsub(pattern="\\.01", replacement="", x=colnames(GBM2), perl=TRUE)
GBMnames2.2<- gsub(pattern="^X*", replacement="", x=GBMnames2, perl=TRUE)
GBMnames3<- gsub(pattern="\\.01", replacement="", x=colnames(GBM3), perl=TRUE)
GBMnames3.2<- gsub(pattern="^X*", replacement="", x=GBMnames3, perl=TRUE)
GBMnames4<- gsub(pattern="\\.01", replacement="", x=colnames(GBM4), perl=TRUE)
GBMnames4.2<- gsub(pattern="^X*", replacement="", x=GBMnames4, perl=TRUE)
GBMnames5<- gsub(pattern="\\.01*", replacement="", x=colnames(GBM5), perl=TRUE)
GBMnames5.2<- gsub(pattern="^X*", replacement="", x=GBMnames5, perl=TRUE)
colnames(GBM1)=GBMnames1.2
colnames(GBM2)=GBMnames2.2
colnames(GBM3)=GBMnames3.2
colnames(GBM4)=GBMnames4.2
colnames(GBM5)=GBMnames5.2
matched<- match(GBMnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
GBMfinal1<- tratioL1HS[,matched]
matched2<- match(GBMnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
GBMfinal2<- tratioL1HS[,matched2]
matched3<- match(GBMnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
GBMfinal3<- tratioL1HS[,matched3]
matched4<- match(GBMnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
GBMfinal4<- tratioL1HS[,matched4]
matched5<- match(GBMnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
GBMfinal5<- tratioL1HS[,matched5]
#HNSC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/HNSC")
HNSC1<- read.table(file="bcgsc.ca_HNSC.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC2<- read.table(file="HNSC_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC3<- read.table(file="pair_set_279_freeze_Mar262013.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC4<- read.table(file="PR_TCGA_HNSC_PAIR_Capture_All_Pairs_QCPASS_v2.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC5<- read.table(file="PR_TCGA_HNSC_PAIR_Capture_TP-NT_TP-NB.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSCnames1<- gsub(pattern="\\.01", replacement="", x=colnames(HNSC1), perl=TRUE)
HNSCnames1.2<- gsub("^X*", "", x=HNSCnames1)
HNSCnames2<- sub(pattern="\\.01", replacement="", x=colnames(HNSC2), perl=TRUE)
HNSCnames2.2<- gsub("^X*", "", x=HNSCnames2)
HNSCnames3<- sub(pattern="\\.01", replacement="", x=colnames(HNSC3), perl=TRUE)
HNSCnames3.2<- gsub("^X*", "", x=HNSCnames3)
HNSCnames4<- sub(pattern="\\.01", replacement="", x=colnames(HNSC4), perl=TRUE)
HNSCnames4.2<- gsub("^X*", "", x=HNSCnames4)
HNSCnames5<- gsub(pattern="\\.01*", replacement="", x=colnames(HNSC5), perl=TRUE)
HNSCnames5.2<- gsub("^X*", "", x=HNSCnames5)
colnames(HNSC1)=HNSCnames1.2
colnames(HNSC2)=HNSCnames2.2
colnames(HNSC3)=HNSCnames3.2
colnames(HNSC4)=HNSCnames4.2
colnames(HNSC5)=HNSCnames5.2
matched<- match(HNSCnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
HNSCfinal1<- data.frame(tratioL1HS[,matched])
HNSCcut1<- HNSC1[,rownames(HNSCfinal1)]
matched2<- match(HNSCnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
HNSCfinal2<- data.frame(tratioL1HS[,matched2])
HNSCcut2<- HNSC2[,rownames(HNSCfinal2)]
matched3<- match(HNSCnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
HNSCfinal3<- data.frame(tratioL1HS[,matched3])
HNSCcut3<- HNSC3[,rownames(HNSCfinal3)]
matched4<- match(HNSCnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
HNSCfinal4<- data.frame(tratioL1HS[,matched4])
HNSCcut4<- HNSC4[,rownames(HNSCfinal4)]
matched5<- match(HNSCnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
HNSCfinal5<- data.frame(tratioL1HS[,matched5])
HNSCcut5<- HNSC5[,rownames(HNSCfinal5)]
HNSCm1<- merge(HNSCcut5, HNSCcut4, all=TRUE)
HNSCm2<- merge(HNSCm1, HNSCcut3, all=TRUE)
HNSCm3<- merge(HNSCm2, HNSCcut2, all=TRUE)
HNSCm4<- merge(HNSCm3, HNSCcut1, all=TRUE)

#KICH
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/KICH")
KICH1<- read.table(file="bcgsc.ca_KICH.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
KICH3<- read.table(file="hgsc.bcm.edu_KICH.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
KICH2<- read.table(file="BCM-KICH-TCGA.broad_calls.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
KICH4<- read.table(file="hgsc.bcm.edu_KICH.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
KICH5<- read.table(file="hgsc.bcm.edu_KICH.IlluminaGA_DNASeq.1.somatic_2.TOTALCOUNT.tsv", header=TRUE)
KICH6<- read.table(file="hgsc.bcm.edu_KICH.IlluminaGA_DNASeq.mitochondria.1.somatic.TOTALCOUNT.tsv", header=TRUE)
KICHnames1<- gsub(pattern="\\.01", replacement="", x=colnames(KICH1), perl=TRUE)
KICHnames1.2<- gsub("^X*", "", x=KICHnames1)
KICHnames2<- gsub(pattern="\\.01", replacement="", x=colnames(KICH2), perl=TRUE)
KICHnames2.2<- gsub("^X*", "", x=KICHnames2)
KICHnames3<- gsub(pattern="\\.01", replacement="", x=colnames(KICH3), perl=TRUE)
KICHnames3.2<- gsub("^X*", "", x=KICHnames3)
KICHnames4<- gsub(pattern="\\.01", replacement="", x=colnames(KICH4), perl=TRUE)
KICHnames4.2<- gsub("^X*", "", x=KICHnames4)
KICHnames5<- gsub(pattern="\\.01", replacement="", x=colnames(KICH5), perl=TRUE)
KICHnames5.2<- gsub("^X*", "", x=KICHnames5)
KICHnames6<- gsub(pattern="\\.01", replacement="", x=colnames(KICH6), perl=TRUE)
KICHnames6.2<- gsub("^X*", "", x=KICHnames6)
colnames(KICH1)=KICHnames1.2
colnames(KICH2)=KICHnames2.2
colnames(KICH3)=KICHnames3.2
colnames(KICH4)=KICHnames4.2
colnames(KICH5)=KICHnames5.2
colnames(KICH6)=KICHnames6.2
matched<- match(KICHnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
KICHfinal1<- data.frame(tratioL1HS[,matched])
KICHcut1<- KICH1[,rownames(KICHfinal1)]
matched2<- match(KICHnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
KICHfinal2<- data.frame(tratioL1HS[,matched2])
KICHcut2<- KICH2[,rownames(KICHfinal2)]
matched3<- match(KICHnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
KICHfinal3<- data.frame(tratioL1HS[,matched3])
KICHcut3<- KICH3[,rownames(KICHfinal3)]
matched4<- match(KICHnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
KICHfinal4<- data.frame(tratioL1HS[,matched4])
KICHcut4<- KICH4[,rownames(KICHfinal4)]
matched5<- match(KICHnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
KICHfinal5<- data.frame(tratioL1HS[,matched5])
KICHcut5<- KICH5[,rownames(KICHfinal5)]
matched6<- match(KICHnames6.2, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
KICHfinal6<- data.frame(tratioL1HS[,matched6])
KICHcut6<- KICH6[,rownames(KICHfinal6)]
KICHm1<- merge(KICHcut6, KICHcut5, all=TRUE)
KICHm2<- merge(KICHm1, KICHcut4, all=TRUE)
KICHm3<- merge(KICHm2, KICHcut3, all=TRUE)
KICHm4<- merge(KICHm3, KICHcut2, all=TRUE)
KICHmerge<- merge(KICHm4, KICHcut1, all=TRUE)

#KIRC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/KIRC")
KIRC1<- read.table(file="BI_and_BCM_1.4.aggregated.tcga.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRC2<- read.table(file="hgsc.bcm.edu_KIRC.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRC3<- read.table(file="hgsc.bcm.edu_KIRC.Mixed_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRC4<- read.table(file="PR_TCGA_KIRC_PAIR_Capture_All_Pairs_QCPASS_v2.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRC5<- read.table(file="PR_TCGA_KIRC_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.automated.somatic.TOTALC", header=TRUE)
KIRCnames1<- gsub(pattern="\\.01", replacement="", x=colnames(KIRC1), perl=TRUE)
KIRCnames1.2<- gsub("^X*", "", x=KIRCnames1)
KIRCnames2<- gsub(pattern="\\.01", replacement="", x=colnames(KIRC2), perl=TRUE)
KIRCnames2.2<- gsub("^X*", "", x=KIRCnames2)
KIRCnames3<- gsub(pattern="\\.01", replacement="", x=colnames(KIRC3), perl=TRUE)
KIRCnames3.2<- gsub("^X*", "", x=KIRCnames3)
KIRCnames4<- gsub(pattern="\\.01", replacement="", x=colnames(KIRC4), perl=TRUE)
KIRCnames4.2<- gsub("^X*", "", x=KIRCnames4)
KIRCnames5<- gsub(pattern="\\.01", replacement="", x=colnames(KIRC5), perl=TRUE)
KIRCnames5.2<- gsub("^X*", "", x=KIRCnames5)
colnames(KIRC1)=KIRCnames1.2
colnames(KIRC2)=KIRCnames2.2
colnames(KIRC3)=KIRCnames3.2
colnames(KIRC4)=KIRCnames4.2
colnames(KIRC5)=KIRCnames5.2
matched<- match(KIRCnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
KIRCfinal1<- data.frame(tratioL1HS[,matched])
KIRCcut1<- KIRC1[,rownames(KIRCfinal1)]
matched2<- match(KIRCnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
KIRCfinal2<- data.frame(tratioL1HS[,matched2])
KIRCcut2<- KIRC2[,rownames(KIRCfinal2)]
matched3<- match(KIRCnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
KIRCfinal3<- data.frame(tratioL1HS[,matched3])
KIRCcut3<- KIRC3[,rownames(KIRCfinal3)]
matched4<- match(KIRCnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
KIRCfinal4<- tratioL1HS[,matched4]
KIRCcut4<- KIRC4[,rownames(KIRCfinal4)]
matched5<- match(KIRCnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
KIRCfinal5<- data.frame(tratioL1HS[,matched5])
KIRCcut5<- KIRC5[,rownames(KIRCfinal5)]
KIRCm1<- merge(KIRCcut5, KIRCcut4, all=TRUE)
KIRCm2<- merge(KIRCm1, KIRCcut3, all=TRUE)
KIRCm3<- merge(KIRCm2, KIRCcut2, all=TRUE)
KIRCmerge<- merge(KIRCm3, KIRCcut1, all=TRUE)

#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/KIRP")
KIRP1<- read.table(file="An_TCGA_KIRP_MultiCenterCalling_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTA", header=TRUE)
KIRP2<- read.table(file="hgsc.bcm.edu_KIRP.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRP3<- read.table(file="hgsc.bcm.edu_KIRP.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
KIRP4<- read.table(file="hgsc.bcm.edu_KIRP.IlluminaGA_DNASeq.1.somatic_2.TOTALCOUNT.tsv", header=TRUE)
KIRP5<- read.table(file="PR_TCGA_KIRP_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRP6<- read.table(file="PR_TCGA_KIRP_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.automated.somatic.TOTALC", header=TRUE)
KIRP7<- read.table(file="ucsc.edu_KIRP.IlluminaGA_DNASeq_automated.Level_2.1.2.0.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRPnames1<- gsub(pattern="\\.01", replacement="", x=colnames(KIRP1), perl=TRUE)
KIRPnames1.2<- gsub("^X*", "", x=KIRPnames1)
KIRPnames2<- gsub(pattern="\\.01", replacement="", x=colnames(KIRP2), perl=TRUE)
KIRPnames2.2<- gsub("^X*", "", x=KIRPnames2)
KIRPnames3<- gsub(pattern="\\.01", replacement="", x=colnames(KIRP3), perl=TRUE)
KIRPnames3.2<- gsub("^X*", "", x=KIRPnames3)
KIRPnames4<- gsub(pattern="\\.01", replacement="", x=colnames(KIRP4), perl=TRUE)
KIRPnames4.2<- gsub("^X*", "", x=KIRPnames4)
KIRPnames5<- gsub(pattern="\\.01", replacement="", x=colnames(KIRP5), perl=TRUE)
KIRPnames5.2<- gsub("^X*", "", x=KIRPnames5)
KIRPnames6<- gsub(pattern="\\.01", replacement="", x=colnames(KIRP6), perl=TRUE)
KIRPnames6.2<- gsub("^X*", "", x=KIRPnames6)
KIRPnames7<- gsub(pattern="\\.01", replacement="", x=colnames(KIRP7), perl=TRUE)
KIRPnames7.2<- gsub("^X*", "", x=KIRPnames7)
colnames(KIRP1)=KIRPnames1.2
colnames(KIRP2)=KIRPnames2.2
colnames(KIRP3)=KIRPnames3.2
colnames(KIRP4)=KIRPnames4.2
colnames(KIRP5)=KIRPnames5.2
colnames(KIRP6)=KIRPnames6.2
colnames(KIRP7)=KIRPnames7.2
matched<- match(KIRPnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
KIRPfinal1<- data.frame(tratioL1HS[,matched])
KIRPcut1<- KIRP1[,rownames(KIRPfinal1)]
matched2<- match(KIRPnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
KIRPfinal2<- data.frame(tratioL1HS[,matched2])
KIRPcut2<- KIRP2[,rownames(KIRPfinal2)]
matched3<- match(KIRPnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
KIRPfinal3<- data.frame(tratioL1HS[,matched3])
KIRPcut3<- KIRP3[,rownames(KIRPfinal3)]
matched4<- match(KIRPnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
KIRPfinal4<- data.frame(tratioL1HS[,matched4])
KIRPcut4<- KIRP4[,rownames(KIRPfinal4)]
matched5<- match(KIRPnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
KIRPfinal5<- data.frame(tratioL1HS[,matched5])
KIRPcut5<- KIRP5[,rownames(KIRPfinal5)]
matched6<- match(KIRPnames6.2, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
KIRPfinal6<- data.frame(tratioL1HS[,matched6])
KIRPcut6<- KIRP6[,rownames(KIRPfinal6)]
matched7<- match(KIRPnames7.2, rownames(ratioL1HS))
matched7<- matched7[!is.na(matched7)]
KIRPfinal7<- data.frame(tratioL1HS[,matched7])
KIRPcut7<- KIRP7[,rownames(KIRPfinal7)]
KIRPm1<- merge(KIRPcut7, KIRPcut6, all=TRUE)
KIRPm2<- merge(KIRPm1, KIRPcut5, all=TRUE)
KIRPm3<- merge(KIRPm2, KIRPcut4, all=TRUE)
KIRPm4<- merge(KIRPm3, KIRPcut3, all=TRUE)
KIRPm5<- merge(KIRPm4, KIRPcut2, all=TRUE)
KIRPm6<- merge(KIRPm5, KIRPcut1, all=TRUE)

#LAML
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LAML")
LAML1<- read.table(file="genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.13.0.somatic.TOTALCOUNT.tsv", header=TRUE)
LAML2<- read.table(file="genome.wustl.edu_LAML.IlluminaHiSeq_DNASeq_automated.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
LAMLnames1<- gsub(pattern="\\.01", replacement="", x=colnames(LAML1), perl=TRUE)
LAMLnames1.2<- gsub("^X*", "", x=LAMLnames1)
LAMLnames1.3<- gsub("\\.03", "", x=LAMLnames1.2)
LAMLnames2<- gsub(pattern="\\.01", replacement="", x=colnames(LAML2), perl=TRUE)
LAMLnames2.2<- gsub("^X*", "", x=LAMLnames2)
LAMLnames2.3<- gsub("\\.03", "", x=LAMLnames2.2)
colnames(LAML1)=LAMLnames1.3
colnames(LAML2)=LAMLnames2.3
matched<- match(LAMLnames1.3, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
LAMLfinal1<- tratioL1HS[,matched]
matched2<- match(LAMLnames2.3, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
LAMLfinal2<- tratioL1HS[,matched2]
#LGG
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LGG")
LGG1<- read.table(file="hgsc.bcm.edu_LGG.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
LGG2<- read.table(file="LGG_FINAL_ANALYSIS.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
LGG3<- read.table(file="LGG_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
LGG4<- read.table(file="mdanderson.org_LGG.IlluminaGA_DNASeq.Level_2.1.6.somatic.TOTALCOUNT.tsv", header=TRUE)
LGG5<- read.table(file="mdanderson.org_LGG.IlluminaGA_DNASeq.Level_2.1.somatic.TOTALCOUNT.tsv", header=TRUE)
LGG6<- read.table(file="PR_TCGA_LGG_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
LGG7<- read.table(file="ucsc.edu_LGG.IlluminaGA_DNASeq_automated.Level_2.1.4.0.somatic.TOTALCOUNT.tsv", header=TRUE)
LGGnames1<- gsub(pattern="\\.01", replacement="", x=colnames(LGG1), perl=TRUE)
LGGnames1.2<- gsub("^X*", "", x=LGGnames1)
LGGnames2<- gsub(pattern="\\.01", replacement="", x=colnames(LGG2), perl=TRUE)
LGGnames2.2<- gsub("^X*", "", x=LGGnames2)
LGGnames3<- gsub(pattern="\\.01", replacement="", x=colnames(LGG3), perl=TRUE)
LGGnames3.2<- gsub("^X*", "", x=LGGnames3)
LGGnames4<- gsub(pattern="\\.01", replacement="", x=colnames(LGG4), perl=TRUE)
LGGnames4.2<- gsub("^X*", "", x=LGGnames4)
LGGnames5<- gsub(pattern="\\.01", replacement="", x=colnames(LGG5), perl=TRUE)
LGGnames5.2<- gsub("^X*", "", x=LGGnames5)
LGGnames6<- gsub(pattern="\\.01", replacement="", x=colnames(LGG6), perl=TRUE)
LGGnames6.2<- gsub("^X*", "", x=LGGnames6)
LGGnames7<- gsub(pattern="\\.01", replacement="", x=colnames(LGG7), perl=TRUE)
LGGnames7.2<- gsub("^X*", "", x=LGGnames7)
colnames(LGG1)=LGGnames1.2
colnames(LGG2)=LGGnames2.2
colnames(LGG3)=LGGnames3.2
colnames(LGG4)=LGGnames4.2
colnames(LGG5)=LGGnames5.2
colnames(LGG6)=LGGnames6.2
colnames(LGG7)=LGGnames7.2
matched<- match(LGGnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
LGGfinal1<- tratioL1HS[,matched]
matched2<- match(LGGnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
LGGfinal2<- tratioL1HS[,matched2]
matched3<- match(LGGnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
LGGfinal3<- tratioL1HS[,matched3]
matched4<- match(LGGnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
LGGfinal4<- tratioL1HS[,matched4]
matched5<- match(LGGnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
LGGfinal5<- tratioL1HS[,matched5]
matched6<- match(LGGnames6.2, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
LGGfinal6<- tratioL1HS[,matched6]
matched7<- match(LGGnames7.2, rownames(ratioL1HS))
matched7<- matched7[!is.na(matched7)]
LGGfinal7<- tratioL1HS[,matched7]
#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LIHC")
LIHC1<- read.table(file="An_TCGA_LIHC_External_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHC2<- read.table(file="bcgsc.ca_LIHC.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHC3<- read.table(file="hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHC4<- read.table(file="hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
LIHC5<- read.table(file="hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic_2.TOTALCOUNT.tsv", header=TRUE)
LIHC6<- read.table(file="ucsc.edu_LIHC.IlluminaGA_DNASeq_automated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHCnames1<- gsub(pattern="\\.01", replacement="", x=colnames(LIHC1), perl=TRUE)
LIHCnames1.2<- gsub("^X*", "", x=LIHCnames1)
LIHCnames2<- gsub(pattern="\\.01", replacement="", x=colnames(LIHC2), perl=TRUE)
LIHCnames2.2<- gsub("^X*", "", x=LIHCnames2)
LIHCnames3<- gsub(pattern="\\.01", replacement="", x=colnames(LIHC3), perl=TRUE)
LIHCnames3.2<- gsub("^X*", "", x=LIHCnames3)
LIHCnames4<- gsub(pattern="\\.01", replacement="", x=colnames(LIHC4), perl=TRUE)
LIHCnames4.2<- gsub("^X*", "", x=LIHCnames4)
LIHCnames5<- gsub(pattern="\\.01", replacement="", x=colnames(LIHC5), perl=TRUE)
LIHCnames5.2<- gsub("^X*", "", x=LIHCnames5)
LIHCnames6<- gsub(pattern="\\.01", replacement="", x=colnames(LIHC6), perl=TRUE)
LIHCnames6.2<- gsub("^X*", "", x=LIHCnames6)
colnames(LIHC1)=LIHCnames1.2
colnames(LIHC2)=LIHCnames2.2
colnames(LIHC3)=LIHCnames3.2
colnames(LIHC4)=LIHCnames4.2
colnames(LIHC5)=LIHCnames5.2
colnames(LIHC6)=LIHCnames6.2
matched<- match(LIHCnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
LIHCfinal1<- data.frame(tratioL1HS[,matched])
LIHCcut1<- LIHC1[,rownames(LIHCfinal1)]
tLIHCcut1<- t(LIHCcut1)
matched2<- match(LIHCnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
LIHCfinal2<- data.frame(tratioL1HS[,matched2])
LIHCcut2<- LIHC2[,rownames(LIHCfinal2)]
tLIHCcut2<- t(LIHCcut2)
matched3<- match(LIHCnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
LIHCfinal3<- data.frame(tratioL1HS[,matched3])
LIHCcut3<- LIHC3[,rownames(LIHCfinal3)]
tLIHCcut3<- t(LIHCcut3)
matched4<- match(LIHCnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
LIHCfinal4<- data.frame(tratioL1HS[,matched4])
LIHCcut4<- LIHC4[,rownames(LIHCfinal4)]
tLIHCcut4<- t(LIHCcut4)
matched5<- match(LIHCnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
LIHCfinal5<- data.frame(tratioL1HS[,matched5])
LIHCcut5<- LIHC5[,rownames(LIHCfinal5)]
tLIHCcut5<- t(LIHCcut5)
matched6<- match(LIHCnames6.2, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
LIHCfinal6<- data.frame(tratioL1HS[,matched6])
LIHCcut6<- LIHC6[,rownames(LIHCfinal6)]
tLIHCcut6<- t(LIHCcut6)
LIHCm1<- merge(LIHCcut6, LIHCcut5, all=TRUE)
LIHCm2<- merge(LIHCm1, LIHCcut4, all=TRUE)
LIHCm3<- merge(LIHCm2, LIHCcut3, all=TRUE)
LIHCm4<- merge(LIHCm3, LIHCcut2, all=TRUE)
LIHCm5<- merge(LIHCm4, LIHCcut1, all=TRUE)

#LUAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LUAD")
LUAD1<- read.table(file="AN_TCGA_LUAD_PAIR_capture_freeze_FINAL_230.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT", header=TRUE)
LUAD2<- read.table(file="PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
LUAD3<- read.table(file="PR_TCGA_LUAD_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.automated.somatic.TOTALC", header=TRUE)
LUADnames1<- gsub(pattern="\\.01", replacement="", x=colnames(LUAD1), perl=TRUE)
LUADnames1.2<- gsub("^X*", "", x=LUADnames1)
LUADnames2<- gsub(pattern="\\.01", replacement="", x=colnames(LUAD2), perl=TRUE)
LUADnames2.2<- gsub("^X*", "", x=LUADnames2)
LUADnames3<- gsub(pattern="\\.01", replacement="", x=colnames(LUAD3), perl=TRUE)
LUADnames3.2<- gsub("^X*", "", x=LUADnames3)
colnames(LUAD1)=LUADnames1.2
colnames(LUAD2)=LUADnames2.2
colnames(LUAD3)=LUADnames3.2
matched<- match(LUADnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
LUADfinal1<- data.frame(tratioL1HS[,matched])
LUADcut1<- LUAD1[,rownames(LUADfinal1)]
matched2<- match(LUADnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
LUADfinal2<- data.frame(tratioL1HS[,matched2])
LUADcut2<- LUAD2[,rownames(LUADfinal2)]
matched3<- match(LUADnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
LUADfinal3<- data.frame(tratioL1HS[,matched3])
LUADcut3<- LUAD3[,rownames(LUADfinal3)]
LUADm1<- merge(LUADcut3, LUADcut2, all=TRUE)
LUADmerge<- merge(LUADm1, LUADcut1, all=TRUE)

#LUSC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LUSC")
LUSC1<- read.table(file="LUSC_Paper_v8.aggregated.tcga.somatic.TOTALCOUNT.tsv", header=TRUE)
LUSC2<- read.table(file="step4_LUSC_Paper_v8.aggregated.tcga.maf2.4.migrated.somatic.TOTALCOUNT.tsv", header=TRUE)
LUSCnames1<- gsub(pattern="\\.01", replacement="", x=colnames(LUSC1), perl=TRUE)
LUSCnames1.2<- gsub("^X*", "", x=LUSCnames1)
LUSCnames2<- gsub(pattern="\\.01", replacement="", x=colnames(LUSC2), perl=TRUE)
LUSCnames2.2<- gsub("^X*", "", x=LUSCnames2)
colnames(LUSC1)=LUSCnames1.2
colnames(LUSC2)=LUSCnames2.2
matched<- match(LUSCnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
LUSCfinal1<- data.frame(tratioL1HS[,matched])
LUSCcut1<- LUSC1[,rownames(LUSCfinal1)]
matched2<- match(LUSCnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
LUSCfinal2<- data.frame(tratioL1HS[,matched2])
LUSCcut2<- LUSC2[,rownames(LUSCfinal2)]
LUSCmerge<- merge(LUSCcut2, LUSCcut1, all=TRUE)

#MESO
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/MESO")
MESO1<- read.table(file="bcgsc.ca_MESO.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
MESO2<- read.table(file="MESO_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
MESO3<- read.table(file="sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
MESO4<- read.table(file="sanger.ac.uk_MESO.IlluminaHiSeq_DNASeq_automated.Level_2.1.2.0.somatic.TOTALCOUNT.tsv", header=TRUE)
MESOnames1<- gsub(pattern="\\.01", replacement="", x=colnames(MESO1), perl=TRUE)
MESOnames1.2<- gsub("^X*", "", x=MESOnames1)
MESOnames2<- gsub(pattern="\\.01", replacement="", x=colnames(MESO2), perl=TRUE)
MESOnames2.2<- gsub("^X*", "", x=MESOnames2)
MESOnames3<- gsub(pattern="\\.01", replacement="", x=colnames(MESO3), perl=TRUE)
MESOnames3.2<- gsub("^X*", "", x=MESOnames3)
MESOnames4<- gsub(pattern="\\.01", replacement="", x=colnames(MESO4), perl=TRUE)
MESOnames4.2<- gsub("^X*", "", x=MESOnames4)
colnames(MESO1)=MESOnames1.2
colnames(MESO2)=MESOnames2.2
colnames(MESO3)=MESOnames3.2
colnames(MESO4)=MESOnames4.2
matched<- match(MESOnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
MESOfinal1<- tratioL1HS[,matched]
matched2<- match(MESOnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
MESOfinal2<- tratioL1HS[,matched2]
matched3<- match(MESOnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
MESOfinal3<- tratioL1HS[,matched3]
matched4<- match(MESOnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
MESOfinal4<- tratioL1HS[,matched4]
#OV
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/OV")
OV1<- read.table(file="broad.mit.edu_OV.IlluminaGA_DNASeq.Level_2.7.somatic.TOTALCOUNT.tsv", header=TRUE)
OV2<- read.table(file="genome.wustl.edu_OV.IlluminaGA_DNASeq.1.3.somatic.TOTALCOUNT.tsv", header=TRUE)
OV3<- read.table(file="genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.2.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
OV4<- read.table(file="genome.wustl.edu_OV.IlluminaHiSeq_DNASeq_automated.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
OV5<- read.table(file="hgsc.bcm.edu_OV.SOLiD_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
OV6<- read.table(file="hgsc.bcm.edu_OV.SOLiD_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
OV7<- read.table(file="ov_liftover.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
OV8<- read.table(file="step4_ov_liftover.aggregated.capture.tcga.uuid.maf2.4.migrated.somatic.TOTALCOUNT.tsv", header=TRUE)
OVnames1<- gsub(pattern="\\.01", replacement="", x=colnames(OV1), perl=TRUE)
OVnames1.2<- gsub("^X*", "", x=OVnames1)
OVnames2<- gsub(pattern="\\.01", replacement="", x=colnames(OV2), perl=TRUE)
OVnames2.2<- gsub("^X*", "", x=OVnames2)
OVnames3<- gsub(pattern="\\.01", replacement="", x=colnames(OV3), perl=TRUE)
OVnames3.2<- gsub("^X*", "", x=OVnames3)
OVnames4<- gsub(pattern="\\.01", replacement="", x=colnames(OV4), perl=TRUE)
OVnames4.2<- gsub("^X*", "", x=OVnames4)
OVnames5<- gsub(pattern="\\.01", replacement="", x=colnames(OV5), perl=TRUE)
OVnames5.2<- gsub("^X*", "", x=OVnames5)
OVnames6<- gsub(pattern="\\.01", replacement="", x=colnames(OV6), perl=TRUE)
OVnames6.2<- gsub("^X*", "", x=OVnames6)
OVnames7<- gsub(pattern="\\.01", replacement="", x=colnames(OV7), perl=TRUE)
OVnames7.2<- gsub("^X*", "", x=OVnames7)
OVnames8<- gsub(pattern="\\.01", replacement="", x=colnames(OV8), perl=TRUE)
OVnames8.2<- gsub("^X*", "", x=OVnames8)
colnames(OV1)=OVnames1.2
colnames(OV2)=OVnames2.2
colnames(OV3)=OVnames3.2
colnames(OV4)=OVnames4.2
colnames(OV5)=OVnames5.2
colnames(OV6)=OVnames6.2
colnames(OV7)=OVnames7.2
colnames(OV8)=OVnames8.2
matched<- match(OVnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
OVfinal1<- tratioL1HS[,matched]
matched2<- match(OVnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
OVfinal2<- tratioL1HS[,matched2]
matched3<- match(OVnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
OVfinal3<- tratioL1HS[,matched3]
matched4<- match(OVnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
OVfinal4<- tratioL1HS[,matched4]
matched5<- match(OVnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
OVfinal5<- tratioL1HS[,matched]
matched6<- match(OVnames6.2, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
OVfinal6<- tratioL1HS[,matched6]
matched7<- match(OVnames7.2, rownames(ratioL1HS))
matched7<- matched7[!is.na(matched7)]
OVfinal7<- tratioL1HS[,matched7]
matched8<- match(OVnames8.2, rownames(ratioL1HS))
matched8<- matched8[!is.na(matched8)]
OVfinal8<- tratioL1HS[,matched8]
#PAAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PAAD")
PAAD1<- read.table(file="bcgsc.ca_PAAD.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD2<- read.table(file="freeze3.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD3<- read.table(file="hgsc.bcm.edu_PAAD.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD4<- read.table(file="PAAD_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD5<- read.table(file="PR_TCGA_PAAD_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD6<- read.table(file="ucsc.edu_PAAD.IlluminaGA_DNASeq_automated.Level_2.1.2.0.somatic.TOTALCOUNT.tsv", header=TRUE)
PAADnames1<- gsub(pattern="\\.01", replacement="", x=colnames(PAAD1), perl=TRUE)
PAADnames1.2<- gsub("^X*", "", x=PAADnames1)
PAADnames2<- gsub(pattern="\\.01", replacement="", x=colnames(PAAD2), perl=TRUE)
PAADnames2.2<- gsub("^X*", "", x=PAADnames2)
PAADnames3<- gsub(pattern="\\.01", replacement="", x=colnames(PAAD3), perl=TRUE)
PAADnames3.2<- gsub("^X*", "", x=PAADnames3)
PAADnames4<- gsub(pattern="\\.01", replacement="", x=colnames(PAAD4), perl=TRUE)
PAADnames4.2<- gsub("^X*", "", x=PAADnames4)
PAADnames5<- gsub(pattern="\\.01", replacement="", x=colnames(PAAD5), perl=TRUE)
PAADnames5.2<- gsub("^X*", "", x=PAADnames5)
PAADnames6<- gsub(pattern="\\.01", replacement="", x=colnames(PAAD6), perl=TRUE)
PAADnames6.2<- gsub("^X*", "", x=PAADnames6)
colnames(PAAD1)=PAADnames1.2
colnames(PAAD2)=PAADnames2.2
colnames(PAAD3)=PAADnames3.2
colnames(PAAD4)=PAADnames4.2
colnames(PAAD5)=PAADnames5.2
colnames(PAAD6)=PAADnames6.2
matched<- match(PAADnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PAADfinal1<- data.frame(tratioL1HS[,matched])
PAADcut1<- PAAD1[,rownames(PAADfinal1)]
matched2<- match(PAADnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
PAADfinal2<- data.frame(tratioL1HS[,matched2])
PAADcut2<- PAAD2[,rownames(PAADfinal2)]
matched3<- match(PAADnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
PAADfinal3<- data.frame(tratioL1HS[,matched3])
PAADcut3<- PAAD3[,rownames(PAADfinal3)]
matched4<- match(PAADnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
PAADfinal4<- data.frame(tratioL1HS[,matched4])
PAADcut4<- PAAD4[,rownames(PAADfinal4)]
matched5<- match(PAADnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
PAADfinal5<- data.frame(tratioL1HS[,matched5])
PAADcut5<- PAAD5[,rownames(PAADfinal5)]
matched6<- match(PAADnames6.2, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
PAADfinal6<- data.frame(tratioL1HS[,matched6])
PAADcut6<- PAAD6[,rownames(PAADfinal6)]
PAADm1<- merge(PAADcut6, PAADcut5, all=TRUE)
PAADm2<- merge(PAADm1, PAADcut4, all=TRUE)
PAADm3<- merge(PAADm2, PAADcut3, all=TRUE)
PAADm4<- merge(PAADm3, PAADcut2, all=TRUE)
PAADm5<- merge(PAADm4, PAADcut1, all=TRUE)

#PCPG
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PCPG")
PCPG1<- read.table(file="bcgsc.ca_PCPG.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPG2<- read.table(file="hgsc.bcm.edu_PCPG.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPG3<- read.table(file="PCPG_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPG4<- read.table(file="ucsc.edu_PCPG.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
PCPGnames1<- gsub(pattern="\\.01", replacement="", x=colnames(PCPG1), perl=TRUE)
PCPGnames1.2<- gsub("^X*", "", x=PCPGnames1)
PCPGnames2<- gsub(pattern="\\.01", replacement="", x=colnames(PCPG2), perl=TRUE)
PCPGnames2.2<- gsub("^X*", "", x=PCPGnames2)
PCPGnames3<- gsub(pattern="\\.01", replacement="", x=colnames(PCPG3), perl=TRUE)
PCPGnames3.2<- gsub("^X*", "", x=PCPGnames3)
PCPGnames4<- gsub(pattern="\\.01", replacement="", x=colnames(PCPG4), perl=TRUE)
PCPGnames4.2<- gsub("^X*", "", x=PCPGnames4)
colnames(PCPG1)=PCPGnames1.2
colnames(PCPG2)=PCPGnames2.2
colnames(PCPG3)=PCPGnames3.2
colnames(PCPG4)=PCPGnames4.2
matched<- match(PCPGnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PCPGfinal1<- data.frame(tratioL1HS[,matched])
PCPGcut1<- PCPG1[,rownames(PCPGfinal1)]
matched2<- match(PCPGnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
PCPGfinal2<- data.frame(tratioL1HS[,matched2])
PCPGcut2<- PCPG2[,rownames(PCPGfinal2)]
matched3<- match(PCPGnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
PCPGfinal3<- data.frame(tratioL1HS[,matched3])
PCPGcut3<- PCPG3[,rownames(PCPGfinal3)]
matched4<- match(PCPGnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
PCPGfinal4<- data.frame(tratioL1HS[,matched4])
PCPGcut4<- PCPG4[,rownames(PCPGfinal4)]
PCPGm1<- merge(PCPGcut4, PCPGcut3, all=TRUE)
PCPGm2<- merge(PCPGm1, PCPGcut2, all=TRUE)
PCPGm3<- merge(PCPGm2, PCPGcut1, all=TRUE)

#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PRAD")
PRAD1<- read.table(file="hgsc.bcm.edu_PRAD.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PRAD2<- read.table(file="PR_TCGA_PRAD_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
PRAD3<- read.table(file="PRAD_Capture_All_Pairs_QCPASS_v6_Nikki_Nov_25.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOU", header=TRUE)
PRAD4<- read.table(file="PRAD_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
PRAD5<- read.table(file="PR-TCGA-Analysis_set.aggregated.tcga.somatic.TOTALCOUNT.tsv", header=TRUE)
PRADnames1<- gsub(pattern="\\.01", replacement="", x=colnames(PRAD1), perl=TRUE)
PRADnames1.2<- gsub("^X*", "", x=PRADnames1)
PRADnames2<- gsub(pattern="\\.01", replacement="", x=colnames(PRAD2), perl=TRUE)
PRADnames2.2<- gsub("^X*", "", x=PRADnames2)
PRADnames3<- gsub(pattern="\\.01", replacement="", x=colnames(PRAD3), perl=TRUE)
PRADnames3.2<- gsub("^X*", "", x=PRADnames3)
PRADnames4<- gsub(pattern="\\.01", replacement="", x=colnames(PRAD4), perl=TRUE)
PRADnames4.2<- gsub("^X*", "", x=PRADnames4)
PRADnames5<- gsub(pattern="\\.01", replacement="", x=colnames(PRAD5), perl=TRUE)
PRADnames5.2<- gsub("^X*", "", x=PRADnames5)
colnames(PRAD1)=PRADnames1.2
colnames(PRAD2)=PRADnames2.2
colnames(PRAD3)=PRADnames3.2
colnames(PRAD4)=PRADnames4.2
colnames(PRAD5)=PRADnames5.2
matched<- match(PRADnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PRADfinal1<- data.frame(tratioL1HS[,matched])
PRADcut1<- PRAD1[,rownames(PRADfinal1)]
matched2<- match(PRADnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
PRADfinal2<- data.frame(tratioL1HS[,matched2])
PRADcut2<- PRAD2[,rownames(PRADfinal2)]
matched3<- match(PRADnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
PRADfinal3<- data.frame(tratioL1HS[,matched3])
PRADcut3<- PRAD3[,rownames(PRADfinal3)]
matched4<- match(PRADnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
PRADfinal4<- data.frame(tratioL1HS[,matched4])
PRADcut4<- PRAD4[,rownames(PRADfinal4)]
matched5<- match(PRADnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
PRADfinal5<- data.frame(tratioL1HS[,matched5])
PRADcut5<- PRAD5[,rownames(PRADfinal5)]
PRADm1<- merge(PRADcut5, PRADcut4, all=TRUE)
PRADm2<- merge(PRADm1, PRADcut3, all=TRUE)
PRADm3<- merge(PRADm2, PRADcut2, all=TRUE)
PRADm4<- merge(PRADm3, PRADcut1, all=TRUE)

#READ
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/READ")
READ1<- read.table(file="hgsc.bcm.edu_READ.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
READ2<- read.table(file="hgsc.bcm.edu_READ.SOLiD_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
READnames1<- gsub(pattern="\\.01", replacement="", x=colnames(READ1), perl=TRUE)
READnames1.2<- gsub("^X*", "", x=READnames1)
READnames2<- gsub(pattern="\\.01", replacement="", x=colnames(READ2), perl=TRUE)
READnames2.2<- gsub("^X*", "", x=READnames2)
colnames(READ1)=READnames1.2
colnames(READ2)=READnames2.2
matched<- match(READnames1.2, rownames(ratioL1HS))
matched<- data.frame(matched[!is.na(matched)])
READfinal1<- data.frame(tratioL1HS[,matched])
READcut1<- READ1[,rownames(READfinal1)]
matched2<- match(READnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
READfinal2<- data.frame(tratioL1HS[,matched2])
READcut2<- READ2[,rownames(READfinal2)]
READmerge<- merge(READcut1, READcut2, all=TRUE)

#SARC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/SARC")
SARC1<- read.table(file="bcgsc.ca_SARC.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
SARC2<- read.table(file="genome.wustl.edu_SARC.IlluminaHiSeq_DNASeq_automated.1.5.0.somatic.TOTALCOUNT.tsv", header=TRUE)
SARC3<- read.table(file="SARC_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
SARC4<- read.table(file="ucsc.edu_SARC.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
SARCnames1<- gsub(pattern="\\.01", replacement="", x=colnames(SARC1), perl=TRUE)
SARCnames1.2<- gsub("^X*", "", x=SARCnames1)
SARCnames2<- gsub(pattern="\\.01", replacement="", x=colnames(SARC2), perl=TRUE)
SARCnames2.2<- gsub("^X*", "", x=SARCnames2)
SARCnames3<- gsub(pattern="\\.01", replacement="", x=colnames(SARC3), perl=TRUE)
SARCnames3.2<- gsub("^X*", "", x=SARCnames3)
SARCnames4<- gsub(pattern="\\.01", replacement="", x=colnames(SARC4), perl=TRUE)
SARCnames4.2<- gsub("^X*", "", x=SARCnames4)
colnames(SARC1)=SARCnames1.2
colnames(SARC2)=SARCnames2.2
colnames(SARC3)=SARCnames3.2
colnames(SARC4)=SARCnames4.2
matched<- match(SARCnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
SARCfinal1<- data.frame(tratioL1HS[,matched])
SARCcut1<- SARC1[,rownames(SARCfinal1)]
matched2<- match(SARCnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
SARCfinal2<- data.frame(tratioL1HS[,matched2])
SARCcut2<- SARC2[,rownames(SARCfinal2)]
matched3<- match(SARCnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
SARCfinal3<- data.frame(tratioL1HS[,matched3])
SARCcut3<- SARC3[,rownames(SARCfinal3)]
matched4<- match(SARCnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
SARCfinal4<- data.frame(tratioL1HS[,matched4])
SARCcut4<- SARC4[,rownames(SARCfinal4)]
SARCm1<- merge(SARCcut4, SARCcut3, all=TRUE)
SARCm2<- merge(SARCm1, SARCcut2, all=TRUE)
SARCm3<- merge(SARCm2, SARCcut1, all=TRUE)

#SKCM
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/SKCM")
SKCM1<- read.table(file="hgsc.bcm.edu_SKCM.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
SKCM2<- read.table(file="PR_TCGA_SKCM_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
SKCM3<- read.table(file="skcm_clean_pairs.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
SKCM4<- read.table(file="SKCM_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
SKCMnames1<- gsub(pattern="\\.01", replacement="", x=colnames(SKCM1), perl=TRUE)
SKCMnames1.2<- gsub("^X*", "", x=SKCMnames1)
SKCMnames1.3<- gsub("\\.06", "", x=SKCMnames1.2)
SKCMnames2<- gsub(pattern="\\.01", replacement="", x=colnames(SKCM2), perl=TRUE)
SKCMnames2.2<- gsub("^X*", "", x=SKCMnames2)
SKCMnames2.3<- gsub("\\.06", "", x=SKCMnames2.2)
SKCMnames3<- gsub(pattern="\\.01", replacement="", x=colnames(SKCM3), perl=TRUE)
SKCMnames3.2<- gsub("^X*", "", x=SKCMnames3)
SKCMnames3.3<- gsub("\\.06", "", x=SKCMnames3.2)
SKCMnames4<- gsub(pattern="\\.01", replacement="", x=colnames(SKCM4), perl=TRUE)
SKCMnames4.2<- gsub("^X*", "", x=SKCMnames4)
SKCMnames4.3<- gsub("\\.06", "", x=SKCMnames4.2)
colnames(SKCM1)=SKCMnames1.3
colnames(SKCM2)=SKCMnames2.3
colnames(SKCM3)=SKCMnames3.3
colnames(SKCM4)=SKCMnames4.3
matched<- match(SKCMnames1.3, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
SKCMfinal1<- tratioL1HS[,matched]
matched2<- match(SKCMnames2.3, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
SKCMfinal2<- tratioL1HS[,matched2]
matched3<- match(SKCMnames3.3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
SKCMfinal3<- tratioL1HS[,matched3]
matched4<- match(SKCMnames4.3, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
SKCMfinal4<- tratioL1HS[,matched4]
#STAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/STAD")
STAD1<- read.table(file="An_TCGA_STAD_External_capture_All_Pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.t", header=TRUE)
STAD2<- read.table(file="bcgsc.ca_STAD.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
STAD3<- read.table(file="genome.wustl.edu_STAD.IlluminaHiSeq_DNASeq_automated.1.3.0.somatic.TOTALCOUNT.tsv", header=TRUE)
STAD4<- read.table(file="hgsc.bcm.edu_STAD.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
STAD5<- read.table(file="hgsc.bcm.edu_STAD.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
STAD6<- read.table(file="PR_TCGA_STAD_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
STAD7<- read.table(file="QCv5_blacklist_Pass.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
STAD8<- read.table(file="ucsc.edu_STAD.IlluminaGA_DNASeq_automated.Level_2.2.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
STADnames1<- gsub(pattern="\\.01", replacement="", x=colnames(STAD1), perl=TRUE)
STADnames1.2<- gsub("^X*", "", x=STADnames1)
STADnames2<- gsub(pattern="\\.01", replacement="", x=colnames(STAD2), perl=TRUE)
STADnames2.2<- gsub("^X*", "", x=STADnames2)
STADnames3<- gsub(pattern="\\.01", replacement="", x=colnames(STAD3), perl=TRUE)
STADnames3.2<- gsub("^X*", "", x=STADnames3)
STADnames4<- gsub(pattern="\\.01", replacement="", x=colnames(STAD4), perl=TRUE)
STADnames4.2<- gsub("^X*", "", x=STADnames4)
STADnames5<- gsub(pattern="\\.01", replacement="", x=colnames(STAD5), perl=TRUE)
STADnames5.2<- gsub("^X*", "", x=STADnames5)
STADnames6<- gsub(pattern="\\.01", replacement="", x=colnames(STAD6), perl=TRUE)
STADnames6.2<- gsub("^X*", "", x=STADnames6)
STADnames7<- gsub(pattern="\\.01", replacement="", x=colnames(STAD7), perl=TRUE)
STADnames7.2<- gsub("^X*", "", x=STADnames7)
STADnames8<- gsub(pattern="\\.01", replacement="", x=colnames(STAD8), perl=TRUE)
STADnames8.2<- gsub("^X*", "", x=STADnames8)
colnames(STAD1)=STADnames1.2
colnames(STAD2)=STADnames2.2
colnames(STAD3)=STADnames3.2
colnames(STAD4)=STADnames4.2
colnames(STAD5)=STADnames5.2
colnames(STAD6)=STADnames6.2
colnames(STAD7)=STADnames7.2
colnames(STAD8)=STADnames8.2
matched<- match(STADnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
STADfinal1<- data.frame(tratioL1HS[,matched])
STADcut1<- STAD1[,rownames(STADfinal1)]
matched2<- match(STADnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
STADfinal2<- data.frame(tratioL1HS[,matched2])
STADcut2<- STAD2[,rownames(STADfinal2)]
matched3<- match(STADnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
STADfinal3<- data.frame(tratioL1HS[,matched3])
STADcut3<- STAD3[,rownames(STADfinal3)]
matched4<- match(STADnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
STADfinal4<- data.frame(tratioL1HS[,matched4])
STADcut4<- STAD4[,rownames(STADfinal4)]
matched5<- match(STADnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
STADfinal5<- data.frame(tratioL1HS[,matched])
STADcut5<- STAD5[,rownames(STADfinal5)]
matched6<- match(STADnames6.2, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
STADfinal6<- data.frame(tratioL1HS[,matched6])
STADcut6<- STAD6[,rownames(STADfinal6)]
matched7<- match(STADnames7.2, rownames(ratioL1HS))
matched7<- matched7[!is.na(matched7)]
STADfinal7<- data.frame(tratioL1HS[,matched7])
STADcut7<- STAD7[,rownames(STADfinal7)]
matched8<- match(STADnames8.2, rownames(ratioL1HS))
matched8<- matched8[!is.na(matched8)]
STADfinal8<- data.frame(tratioL1HS[,matched8])
STADcut8<- STAD8[,rownames(STADfinal8)]
STADm1<- merge(STADcut8, STADcut7, all=TRUE)
STADm2<- merge(STADm1, STADcut6, all=TRUE)
STADm3<- merge(STADm2, STADcut5, all=TRUE)
STADm4<- merge(STADm3, STADcut4, all=TRUE)
STADm5<- merge(STADm4, STADcut3, all=TRUE)
STADm6<- merge(STADm5, STADcut2, all=TRUE)
STADmerge<- merge(STADm6, STADcut1, all=TRUE)

#TGCT
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/TGCT")
TGCT1<- read.table(file="bcgsc.ca_TGCT.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
TGCT2<- read.table(file="hgsc.bcm.edu_TGCT.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
TGCT3<- read.table(file="hgsc.bcm.edu_TGCT.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
TGCT4<- read.table(file="TGCT_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
TGCT5<- read.table(file="ucsc.edu_TGCT.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
TGCTnames1<- gsub(pattern="\\.01", replacement="", x=colnames(TGCT1), perl=TRUE)
TGCTnames1.2<- gsub(pattern="\\.05", replacement="", x=TGCTnames1, perl=TRUE)
TGCTnames2<- gsub(pattern="\\.01", replacement="", x=colnames(TGCT2), perl=TRUE)
TGCTnames2.2<- gsub(pattern="\\.05", replacement="", x=TGCTnames2, perl=TRUE)
TGCTnames3<- gsub(pattern="\\.01", replacement="", x=colnames(TGCT3), perl=TRUE)
TGCTnames3.2<- gsub(pattern="\\.05", replacement="", x=TGCTnames3, perl=TRUE)
TGCTnames4<- gsub(pattern="\\.01", replacement="", x=colnames(TGCT4), perl=TRUE)
TGCTnames4.2<- gsub(pattern="\\.05", replacement="", x=TGCTnames4, perl=TRUE)
TGCTnames5<- gsub(pattern="\\.01", replacement="", x=colnames(TGCT5), perl=TRUE)
TGCTnames5.2<- gsub(pattern="\\.05", replacement="", x=TGCTnames5, perl=TRUE)
colnames(TGCT1)=TGCTnames1.2
colnames(TGCT2)=TGCTnames2.2
colnames(TGCT3)=TGCTnames3.2
colnames(TGCT4)=TGCTnames4.2
colnames(TGCT5)=TGCTnames5.2
matched<- match(TGCTnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
TGCTfinal1<- tratioL1HS[,matched]
matched2<- match(TGCTnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
TGCTfinal2<- tratioL1HS[,matched2]
matched3<- match(TGCTnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
TGCTfinal3<- tratioL1HS[,matched3]
matched4<- match(TGCTnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
TGCTfinal4<- tratioL1HS[,matched4]
matched5<- match(TGCTnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
TGCTfinal5<- tratioL1HS[,matched5]
#THCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/THCA")
THCA1<- read.table(file="AN_TCGA_THCA_PAIR_Capture_ALLQC_14Aug2013_429.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA2<- read.table(file="bcgsc.ca_THCA.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA3<- read.table(file="hgsc.bcm.edu_THCA.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA4<- read.table(file="PR_TCGA_THCA_PAIR_Capture_All_Pairs_QCPASS_v2.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
THCA5<- read.table(file="THCA_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
THCAnames1<- gsub(pattern="\\.01", replacement="", x=colnames(THCA1), perl=TRUE)
THCAnames1.2<- gsub("\\.06", "", THCAnames1)
THCAnames2<- gsub(pattern="\\.01", replacement="", x=colnames(THCA2), perl=TRUE)
THCAnames2.2<- gsub("\\.06", "", THCAnames2)
THCAnames3<- gsub(pattern="\\.01", replacement="", x=colnames(THCA3), perl=TRUE)
THCAnames3.2<- gsub("\\.06", "", THCAnames3)
THCAnames4<- gsub(pattern="\\.01", replacement="", x=colnames(THCA4), perl=TRUE)
THCAnames4.2<- gsub("\\.06", "", THCAnames4)
THCAnames5<- gsub(pattern="\\.01", replacement="", x=colnames(THCA5), perl=TRUE)
THCAnames5.2<- gsub("\\.06", "", THCAnames5)
colnames(THCA1)=THCAnames1.2
colnames(THCA2)=THCAnames2.2
colnames(THCA3)=THCAnames3.2
colnames(THCA4)=THCAnames4.2
colnames(THCA5)=THCAnames5.2
matched<- match(THCAnames1.2, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
THCAfinal1<- data.frame(tratioL1HS[,matched])
THCAcut1<- THCA1[,rownames(THCAfinal1)]
matched2<- match(THCAnames2.2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
THCAfinal2<- data.frame(tratioL1HS[,matched2])
THCAcut2<- THCA2[,rownames(THCAfinal2)]
matched3<- match(THCAnames3.2, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
THCAfinal3<- data.frame(tratioL1HS[,matched3])
THCAcut3<- THCA3[,rownames(THCAfinal3)]
matched4<- match(THCAnames4.2, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
THCAfinal4<- data.frame(tratioL1HS[,matched4])
THCAcut4<- THCA4[,rownames(THCAfinal4)]
matched5<- match(THCAnames5.2, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
THCAfinal5<- data.frame(tratioL1HS[,matched5])
THCAcut5<- THCA5[,rownames(THCAfinal5)]
THCAm1<- merge(THCAcut5, THCAcut4, all=TRUE)
THCAm2<- merge(THCAm1, THCAcut3, all=TRUE)
THCAm3<- merge(THCAm2, THCAcut2, all=TRUE)
THCAm4<- merge(THCAm3, THCAcut1, all=TRUE)

#THYM
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/THYM")
THYM1<- read.table(file="bcgsc.ca_THYM.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM2<- read.table(file="genome.wustl.edu_THYM.IlluminaGA_DNASeq_curated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM3<- read.table(file="hgsc.bcm.edu_THYM.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM4<- read.table(file="THYM_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM5<- read.table(file="ucsc.edu_THYM.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
THYMnames1<- gsub(pattern="\\.01", replacement="", x=colnames(THYM1), perl=TRUE)
THYMnames2<- gsub(pattern="\\.01", replacement="", x=colnames(THYM2), perl=TRUE)
THYMnames3<- gsub(pattern="\\.01", replacement="", x=colnames(THYM3), perl=TRUE)
THYMnames4<- gsub(pattern="\\.01", replacement="", x=colnames(THYM4), perl=TRUE)
THYMnames5<- gsub(pattern="\\.01", replacement="", x=colnames(THYM5), perl=TRUE)
colnames(THYM1)=THYMnames1
colnames(THYM2)=THYMnames2
colnames(THYM3)=THYMnames3
colnames(THYM4)=THYMnames4
colnames(THYM5)=THYMnames5
matched<- match(THYMnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
THYMfinal1<- data.frame(tratioL1HS[,matched])
THYMcut1<- THYM1[,rownames(THYMfinal1)]
matched2<- match(THYMnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
THYMfinal2<- data.frame(tratioL1HS[,matched2])
THYMcut2<- THYM2[,rownames(THYMfinal2)]
matched3<- match(THYMnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
THYMfinal3<- data.frame(tratioL1HS[,matched3])
THYMcut3<- THYM3[,rownames(THYMfinal3)]
matched4<- match(THYMnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
THYMfinal4<- data.frame(tratioL1HS[,matched4])
THYMcut4<- THYM4[,rownames(THYMfinal4)]
matched5<- match(THYMnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
THYMfinal5<- data.frame(tratioL1HS[,matched5])
THYMcut5<- THYM5[,rownames(THYMfinal5)]
THYMm1<- merge(THYMcut5, THYMcut4, all=TRUE)
THYMm2<- merge(THYMm1, THYMcut3, all=TRUE)
THYMm3<- merge(THYMm2, THYMcut2, all=TRUE)
THYMm4<- merge(THYMm3, THYMcut1, all=TRUE)

#UCEC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/UCEC")
UCEC1<- read.table(file="An_UCEC_194.aggregated.tcga.somatic.TOTALCOUNT.tsv", header=TRUE)
UCEC2<- read.table(file="genome.wustl.edu_UCEC.IlluminaGA_DNASeq.Level_2.1.7.somatic.TOTALCOUNT.tsv", header=TRUE)
UCEC3<- read.table(file="genome.wustl.edu_UCEC.IlluminaHiSeq_DNASeq_automated.1.8.0.somatic.TOTALCOUNT.tsv", header=TRUE)
UCEC4<- read.table(file="step4_An_UCEC_194.aggregated.tcga.maf2.4.migrated.somatic.TOTALCOUNT.tsv", header=TRUE)
UCECnames1<- gsub(pattern="\\.01", replacement="", x=colnames(UCEC1), perl=TRUE)
UCECnames2<- gsub(pattern="\\.01", replacement="", x=colnames(UCEC2), perl=TRUE)
UCECnames3<- gsub(pattern="\\.01", replacement="", x=colnames(UCEC3), perl=TRUE)
UCECnames4<- gsub(pattern="\\.01", replacement="", x=colnames(UCEC4), perl=TRUE)
colnames(UCEC1)=UCECnames1
colnames(UCEC2)=UCECnames2
colnames(UCEC3)=UCECnames3
colnames(UCEC4)=UCECnames4
matched<- match(UCECnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
UCECfinal1<- data.frame(tratioL1HS[,matched])
UCECcut1<- UCEC1[,rownames(UCECfinal1)]
matched2<- match(UCECnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
UCECfinal2<- data.frame(tratioL1HS[,matched2])
UCECcut2<- UCEC2[,rownames(UCECfinal2)]
matched3<- match(UCECnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
UCECfinal3<- data.frame(tratioL1HS[,matched3])
UCECcut3<- UCEC3[,rownames(UCECfinal3)]
matched4<- match(UCECnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
UCECfinal4<- data.frame(tratioL1HS[,matched4])
UCECcut4<- UCEC4[,rownames(UCECfinal4)]
UCECm1<- rbind(UCECcut4, UCECcut1)
UCECm2<- rbind(UCECcut3, UCECcut2)
UCECm3<- merge(UCECm2, UCECm1, all=TRUE)

#UCS
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/UCS")
UCS1<- read.table(file="bcgsc.ca_UCS.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
UCS2<- read.table(file="hgsc.bcm.edu_UCS.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
UCS3<- read.table(file="PR_TCGA_UCS_PAIR_Capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
UCS4<- read.table(file="UCS_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
UCSnames1<- gsub(pattern="\\.01", replacement="", x=colnames(UCS1), perl=TRUE)
UCSnames2<- gsub(pattern="\\.01", replacement="", x=colnames(UCS2), perl=TRUE)
UCSnames3<- gsub(pattern="\\.01", replacement="", x=colnames(UCS3), perl=TRUE)
UCSnames4<- gsub(pattern="\\.01", replacement="", x=colnames(UCS4), perl=TRUE)
colnames(UCS1)=UCSnames1
colnames(UCS2)=UCSnames2
colnames(UCS3)=UCSnames3
colnames(UCS4)=UCSnames4
matched<- match(UCSnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
UCSfinal1<- tratioL1HS[,matched]
matched2<- match(UCSnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
UCSfinal2<- tratioL1HS[,matched2]
matched3<- match(UCSnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
UCSfinal3<- tratioL1HS[,matched3]
matched4<- match(UCSnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
UCSfinal4<- tratioL1HS[,matched4]
#UVM
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/UVM")
UVM1<- read.table(file="bcgsc.ca_UVM.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
UVM2<- read.table(file="hgsc.bcm.edu_UVM.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
UVM3<- read.table(file="PR_TCGA_UVM_PAIR_Capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
UVM4<- read.table(file="ucsc.edu_UVM.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
UVM5<- read.table(file="UVM_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
UVMnames1<- gsub(pattern="\\.01", replacement="", x=colnames(UVM1), perl=TRUE)
UVMnames2<- gsub(pattern="\\.01", replacement="", x=colnames(UVM2), perl=TRUE)
UVMnames3<- gsub(pattern="\\.01", replacement="", x=colnames(UVM3), perl=TRUE)
UVMnames4<- gsub(pattern="\\.01", replacement="", x=colnames(UVM4), perl=TRUE)
UVMnames5<- gsub(pattern="\\.01", replacement="", x=colnames(UVM5), perl=TRUE)
colnames(UVM1)=UVMnames1
colnames(UVM2)=UVMnames2
colnames(UVM3)=UVMnames3
colnames(UVM4)=UVMnames4
colnames(UVM5)=UVMnames5
matched<- match(UVMnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
UVMfinal1<- tratioL1HS[,matched]
matched2<- match(UVMnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
UVMfinal2<- tratioL1HS[,matched2]
matched3<- match(UVMnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
UVMfinal3<- tratioL1HS[,matched3]
matched4<- match(UVMnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
UVMfinal4<- tratioL1HS[,matched4]
matched5<- match(UVMnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
UVMfinal5<- tratioL1HS[,matched5]

#merge
m1<- merge(BRCAcut5, BRCAcut4)
m2<- merge(m1, BRCAcut3, all=TRUE)
m3<- merge(m2, BRCAcut2, all=TRUE)
m4<- merge(m3, BRCAcut1, all=TRUE)
notna<- colSums(!is.na(m4))
sum<- colSums(m4, na.rm=TRUE)
avmutBRCA<- sum/notna

#BRCA
#log fold it
logBRCAmut<- log2(avmutBRCA)
logBRCAmut<- data.frame(logBRCAmut)
matchedL<- match(rownames(logBRCAmut), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HSBRCA<- data.frame(tratioL1HS[,matchedL])
L1HSBRCA<- data.frame(log2(L1HSBRCA))
cor(logBRCAmut, L1HSBRCA)
value<- lm(L1HSBRCA[,1] ~ logBRCAmut[,1])
plot(logBRCAmut[,1], L1HSBRCA[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="BRCA Correlation no outlier", col="magenta")
#correlation determination
value<- lm(L1HSBRCA[,1] ~ logBRCAmut[,1])
r2<- summary(value)$r.squared
#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")
#get rid of outlier
L1HSBRCA2<- t(L1HSBRCA)
L1HSBRCA2<- subset(L1HSBRCA2, select = -c(A0DB))
logBRCAmut2<- t(logBRCAmut)
logBRCAmut2<- subset(logBRCAmut2, select = -c(A0DB))
L1HS<- t(L1HSBRCA2)
mut<- t(logBRCAmut2)
#rvalues
cor(mut, L1HS)
plot(mut[,1], L1HS[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="BRCA Correlation no outlier", col="magenta")
abline(lm(L1HS ~ mut))
#without second outlier
L1HS2<- t(L1HS)
L1HS2<- data.frame(L1HS2)
L1HS2<- subset(L1HS2, select = -c(A0H9))
mut2<- t(mut)
mut2<- data.frame(mut2)
mut2<- subset(mut2, select = -c(A0H9))
L1HS<- t(L1HS2)
mut<- t(mut2)
plot(mut[,1], L1HS[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="BRCA Correlation no outlier", col="magenta")
#correlation determination
value<- lm(L1HS[,1] ~ mut[,1])
r2<- summary(value)$r.squared
#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")

#LIHC
#merge
m1<- merge(LIHCcut6, LIHCcut5, all=TRUE)
m2<- merge(m1, LIHCcut4, all=TRUE)
m3<- merge(m2, LIHCcut3, all=TRUE)
m4<- merge(m3, LIHCcut2, all=TRUE)
m5<- merge(m4, LIHCcut1, all=TRUE)
notna<- colSums(!is.na(m5))
sum<- colSums(m5, na.rm=TRUE)
avmutLIHC<- sum/notna
logLIHCmut<- log2(avmutLIHC)
logLIHCmut<- data.frame(logLIHCmut)
matchedL<- match(rownames(logLIHCmut), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HSLIHC<- data.frame(tratioL1HS[,matchedL])
L1HSLIHC<- data.frame(log2(L1HSLIHC))
cor(logLIHCmut, L1HSLIHC)
value<- lm(L1HSLIHC[,1] ~ logLIHCmut[,1])
plot(logLIHCmut[,1], L1HSLIHC[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="LIHC Correlation all", col="dark green")
abline(lm(logLIHCmut[,1] ~ L1HSLIHC[,1]))
#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")
r2<- summary(value)$r.squared

#THCA
#merge
m1<- merge(THCAcut5, THCAcut4, all=TRUE)
m2<- merge(m1, THCAcut3, all=TRUE)
m3<- merge(m2, THCAcut2, all=TRUE)
m4<- merge(m3, THCAcut1, all=TRUE)
notna<- colSums(!is.na(m4))
sum<- colSums(m4, na.rm=TRUE)
avmutTHCA<- sum/notna
logTHCAmut<- log2(avmutTHCA)
logTHCAmut<- data.frame(logTHCAmut)
matchedL<- match(rownames(logTHCAmut), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HSTHCA<- data.frame(tratioL1HS[,matchedL])
L1HSTHCA<- data.frame(log2(L1HSTHCA))
cor(logTHCAmut, L1HSTHCA)
value<- lm(L1HSTHCA[,1] ~ logTHCAmut[,1])
plot(logTHCAmut[,1], L1HSTHCA[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="THCA Correlation all", col="orange")
abline(lm(logTHCAmut[,1] ~ L1HSTHCA[,1]))
#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")
r2<- summary(value)$r.squared

#BLCA
#merge
m1<- merge(BLCAcut5, BLCAcut4, all=TRUE)
m2<- merge(m1, BLCAcut3, all=TRUE)
m3<- merge(m2, BLCAcut2, all=TRUE)
m4<- merge(m3, BLCAcut1, all=TRUE)
notna<- colSums(!is.na(m4))
sum<- colSums(m4, na.rm=TRUE)
avmutBLCA<- sum/notna
logBLCAmut<- log2(avmutBLCA)
logBLCAmut<- data.frame(logBLCAmut)
matchedL<- match(rownames(logBLCAmut), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HSBLCA<- data.frame(tratioL1HS[,matchedL])
L1HSBLCA<- data.frame(log2(L1HSBLCA))
cor(logBLCAmut, L1HSBLCA)
value<- lm(L1HSBLCA[,1] ~ logBLCAmut[,1])
plot(logBLCAmut[,1], L1HSBLCA[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="BLCA Correlation all", col="red")
abline(lm(logBLCAmut[,1] ~ L1HSBLCA[,1]))
#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")
r2<- summary(value)$r.squared

setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
write.table(combineddata, file="Median Mutations", sep="	")

