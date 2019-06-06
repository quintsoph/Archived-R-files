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
BLCA1straw <-read.table(file="bcgsc.ca_BLCA.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
BLCA2ndraw<- read.table(file="BLCA-28-original.aggregated.tcga.somatic.TOTALCOUNT.tsv", header=TRUE)
BLCA3rdraw<- read.table(file="BLCA130_somatic_updated.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
BLCA4rdraw<- read.table(file="PR_TCGA_BLCA_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
BLCA5rdraw<- read.table(file="PR_TCGA_BLCA_PAIR_Capture_All_Pairs_QCPASS_v6.aggregated.capture.tcga.uuid.automated.somatic.TOTALC", header=TRUE)
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
BLCAm1<- merge(BLCAcut5, BLCAcut4, all=TRUE)
BLCAm2<- merge(BLCAm1, BLCAcut3, all=TRUE)
BLCAm3<- merge(BLCAm2, BLCAcut2, all=TRUE)
BLCAm4<- data.frame(merge(BLCAm3, BLCAcut1, all=TRUE))

#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/BRCA")
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
BRCAm1<- rbind(BRCAcut5, BRCAcut4)
BRCAm2<- merge(BRCAm1, BRCAcut3, all=TRUE)
BRCAm3<- merge(BRCAm2, BRCAcut2, all=TRUE)
BRCAm4<- data.frame(merge(BRCAm3, BRCAcut1, all=TRUE))

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
CESCm1<- merge(CESCcut7, CESCcut6, all=TRUE)
CESCm2<- merge(CESCm1, CESCcut5, all=TRUE)
CESCm3<- merge(CESCm2, CESCcut4, all=TRUE)
CESCm4<- merge(CESCm3, CESCcut3, all=TRUE)
CESCm5<- merge(CESCm4, CESCcut2, all=TRUE)
CESCm6<- merge(CESCm5, CESCcut1, all=TRUE)

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
CHOLm1<- merge(CHOLcut5, CHOLcut4, all=TRUE)
CHOLm2<- merge(CHOLm1, CHOLcut3, all=TRUE)
CHOLm3<- merge(CHOLm2, CHOLcut2, all=TRUE)
CHOLm4<- merge(CHOLm3, CHOLcut1, all=TRUE)

#ESCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/ESCA")
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
HNSC1<- read.table(file="bcgsc.ca_HNSC.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC2<- read.table(file="HNSC_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC3<- read.table(file="pair_set_279_freeze_Mar262013.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC4<- read.table(file="PR_TCGA_HNSC_PAIR_Capture_All_Pairs_QCPASS_v2.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSC5<- read.table(file="PR_TCGA_HNSC_PAIR_Capture_TP-NT_TP-NB.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
HNSCnames1<- sub(pattern=".01", replacement="", x=colnames(HNSC1), perl=TRUE)
HNSCnames2<- sub(pattern=".01", replacement="", x=colnames(HNSC2), perl=TRUE)
HNSCnames3<- sub(pattern=".01", replacement="", x=colnames(HNSC3), perl=TRUE)
HNSCnames4<- sub(pattern=".01", replacement="", x=colnames(HNSC4), perl=TRUE)
HNSCnames5<- sub(pattern=".01", replacement="", x=colnames(HNSC5), perl=TRUE)
colnames(HNSC1)=HNSCnames1
colnames(HNSC2)=HNSCnames2
colnames(HNSC3)=HNSCnames3
colnames(HNSC4)=HNSCnames4
colnames(HNSC5)=HNSCnames5
matched<- match(HNSCnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
HNSCfinal1<- data.frame(tratioL1HS[,matched])
HNSCcut1<- HNSC1[,rownames(HNSCfinal1)]
matched2<- match(HNSCnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
HNSCfinal2<- data.frame(tratioL1HS[,matched2])
HNSCcut2<- HNSC2[,rownames(HNSCfinal2)]
matched3<- match(HNSCnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
HNSCfinal3<- data.frame(tratioL1HS[,matched3])
HNSCcut3<- HNSC3[,rownames(HNSCfinal3)]
matched4<- match(HNSCnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
HNSCfinal4<- data.frame(tratioL1HS[,matched4])
HNSCcut4<- HNSC4[,rownames(HNSCfinal4)]
matched5<- match(HNSCnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
HNSCfinal5<- data.frame(tratioL1HS[,matched5])
HNSCcut5<- HNSC5[,rownames(HNSCfinal5)]
HNSCm1<- merge(HNSCcut5, HNSCcut4, all=TRUE)
HNSCm2<- merge(HNSCm1, HNSCcut3, all=TRUE)
HNSCm3<- merge(HNSCm2, HNSCcut2, all=TRUE)
HNSCm4<- merge(HNSCm3, HNSCcut1, all=TRUE)

#KIRP
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/KIRP")
KIRP1<- read.table(file="An_TCGA_KIRP_MultiCenterCalling_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTA", header=TRUE)
KIRP2<- read.table(file="hgsc.bcm.edu_KIRP.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRP3<- read.table(file="hgsc.bcm.edu_KIRP.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
KIRP4<- read.table(file="hgsc.bcm.edu_KIRP.IlluminaGA_DNASeq.1.somatic_2.TOTALCOUNT.tsv", header=TRUE)
KIRP5<- read.table(file="PR_TCGA_KIRP_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRP6<- read.table(file="PR_TCGA_KIRP_PAIR_Capture_All_Pairs_QCPASS_v4.aggregated.capture.tcga.uuid.automated.somatic.TOTALC", header=TRUE)
KIRP7<- read.table(file="ucsc.edu_KIRP.IlluminaGA_DNASeq_automated.Level_2.1.2.0.somatic.TOTALCOUNT.tsv", header=TRUE)
KIRPnames1<- sub(pattern=".01", replacement="", x=colnames(KIRP1), perl=TRUE)
KIRPnames2<- sub(pattern=".01", replacement="", x=colnames(KIRP2), perl=TRUE)
KIRPnames3<- sub(pattern=".01", replacement="", x=colnames(KIRP3), perl=TRUE)
KIRPnames4<- sub(pattern=".01", replacement="", x=colnames(KIRP4), perl=TRUE)
KIRPnames5<- sub(pattern=".01", replacement="", x=colnames(KIRP5), perl=TRUE)
KIRPnames6<- sub(pattern=".01", replacement="", x=colnames(KIRP6), perl=TRUE)
KIRPnames7<- sub(pattern=".01", replacement="", x=colnames(KIRP7), perl=TRUE)
colnames(KIRP1)=KIRPnames1
colnames(KIRP2)=KIRPnames2
colnames(KIRP3)=KIRPnames3
colnames(KIRP4)=KIRPnames4
colnames(KIRP5)=KIRPnames5
colnames(KIRP6)=KIRPnames6
colnames(KIRP7)=KIRPnames7
matched<- match(KIRPnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
KIRPfinal1<- data.frame(tratioL1HS[,matched])
KIRPcut1<- KIRP1[,rownames(KIRPfinal1)]
matched2<- match(KIRPnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
KIRPfinal2<- data.frame(tratioL1HS[,matched2])
KIRPcut2<- KIRP2[,rownames(KIRPfinal2)]
matched3<- match(KIRPnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
KIRPfinal3<- data.frame(tratioL1HS[,matched3])
KIRPcut3<- KIRP3[,rownames(KIRPfinal3)]
matched4<- match(KIRPnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
KIRPfinal4<- data.frame(tratioL1HS[,matched4])
KIRPcut4<- KIRP4[,rownames(KIRPfinal4)]
matched5<- match(KIRPnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
KIRPfinal5<- data.frame(tratioL1HS[,matched5])
KIRPcut5<- KIRP5[,rownames(KIRPfinal5)]
matched6<- match(KIRPnames6, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
KIRPfinal6<- data.frame(tratioL1HS[,matched6])
KIRPcut6<- KIRP6[,rownames(KIRPfinal6)]
matched7<- match(KIRPnames7, rownames(ratioL1HS))
matched7<- matched7[!is.na(matched7)]
KIRPfinal7<- data.frame(tratioL1HS[,matched7])
KIRPcut7<- KIRP7[,rownames(KIRPfinal7)]
KIRPm1<- merge(KIRPcut7, KIRPcut6, all=TRUE)
KIRPm2<- merge(KIRPm1, KIRPcut5, all=TRUE)
KIRPm3<- merge(KIRPm2, KIRPcut4, all=TRUE)
KIRPm4<- merge(KIRPm3, KIRPcut3, all=TRUE)
KIRPm5<- merge(KIRPm4, KIRPcut2, all=TRUE)
KIRPm6<- merge(KIRPm5, KIRPcut1, all=TRUE)


#LIHC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/LIHC")
LIHC1<- read.table(file="An_TCGA_LIHC_External_capture_All_Pairs.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHC2<- read.table(file="bcgsc.ca_LIHC.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHC3<- read.table(file="hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHC4<- read.table(file="hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic_1.TOTALCOUNT.tsv", header=TRUE)
LIHC5<- read.table(file="hgsc.bcm.edu_LIHC.IlluminaGA_DNASeq.1.somatic_2.TOTALCOUNT.tsv", header=TRUE)
LIHC6<- read.table(file="ucsc.edu_LIHC.IlluminaGA_DNASeq_automated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
LIHCnames1<- sub(pattern=".01", replacement="", x=colnames(LIHC1), perl=TRUE)
LIHCnames2<- sub(pattern=".01", replacement="", x=colnames(LIHC2), perl=TRUE)
LIHCnames3<- sub(pattern=".01", replacement="", x=colnames(LIHC3), perl=TRUE)
LIHCnames4<- sub(pattern=".01", replacement="", x=colnames(LIHC4), perl=TRUE)
LIHCnames5<- sub(pattern=".01", replacement="", x=colnames(LIHC5), perl=TRUE)
LIHCnames6<- sub(pattern=".01", replacement="", x=colnames(LIHC6), perl=TRUE)
colnames(LIHC1)=LIHCnames1
colnames(LIHC2)=LIHCnames2
colnames(LIHC3)=LIHCnames3
colnames(LIHC4)=LIHCnames4
colnames(LIHC5)=LIHCnames5
colnames(LIHC6)=LIHCnames6
matched<- match(LIHCnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
LIHCfinal1<- data.frame(tratioL1HS[,matched])
LIHCcut1<- LIHC1[,rownames(LIHCfinal1)]
tLIHCcut1<- t(LIHCcut1)
matched2<- match(LIHCnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
LIHCfinal2<- data.frame(tratioL1HS[,matched2])
LIHCcut2<- LIHC2[,rownames(LIHCfinal2)]
tLIHCcut2<- t(LIHCcut2)
matched3<- match(LIHCnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
LIHCfinal3<- data.frame(tratioL1HS[,matched3])
LIHCcut3<- LIHC3[,rownames(LIHCfinal3)]
tLIHCcut3<- t(LIHCcut3)
matched4<- match(LIHCnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
LIHCfinal4<- data.frame(tratioL1HS[,matched4])
LIHCcut4<- LIHC4[,rownames(LIHCfinal4)]
tLIHCcut4<- t(LIHCcut4)
matched5<- match(LIHCnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
LIHCfinal5<- data.frame(tratioL1HS[,matched5])
LIHCcut5<- LIHC5[,rownames(LIHCfinal5)]
tLIHCcut5<- t(LIHCcut5)
matched6<- match(LIHCnames6, rownames(ratioL1HS))
matched6<- matched6[!is.na(matched6)]
LIHCfinal6<- data.frame(tratioL1HS[,matched6])
LIHCcut6<- LIHC6[,rownames(LIHCfinal6)]
tLIHCcut6<- t(LIHCcut6)
LIHCm1<- merge(LIHCcut6, LIHCcut5, all=TRUE)
LIHCm2<- merge(LIHCm1, LIHCcut4, all=TRUE)
LIHCm3<- merge(LIHCm2, LIHCcut3, all=TRUE)
LIHCm4<- merge(LIHCm3, LIHCcut2, all=TRUE)
LIHCm5<- merge(LIHCm4, LIHCcut1, all=TRUE)

#PAAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PAAD")
PAAD1<- read.table(file="bcgsc.ca_PAAD.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD2<- read.table(file="freeze3.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD3<- read.table(file="hgsc.bcm.edu_PAAD.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD4<- read.table(file="PAAD_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD5<- read.table(file="PR_TCGA_PAAD_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
PAAD6<- read.table(file="ucsc.edu_PAAD.IlluminaGA_DNASeq_automated.Level_2.1.2.0.somatic.TOTALCOUNT.tsv", header=TRUE)
PAADnames1<- sub(pattern=".01", replacement="", x=colnames(PAAD1), perl=TRUE)
PAADnames2<- sub(pattern=".01", replacement="", x=colnames(PAAD2), perl=TRUE)
PAADnames3<- sub(pattern=".01", replacement="", x=colnames(PAAD3), perl=TRUE)
PAADnames4<- sub(pattern=".01", replacement="", x=colnames(PAAD4), perl=TRUE)
PAADnames5<- sub(pattern=".01", replacement="", x=colnames(PAAD5), perl=TRUE)
PAADnames6<- sub(pattern=".01", replacement="", x=colnames(PAAD6), perl=TRUE)
colnames(PAAD1)=PAADnames1
colnames(PAAD2)=PAADnames2
colnames(PAAD3)=PAADnames3
colnames(PAAD4)=PAADnames4
colnames(PAAD5)=PAADnames5
colnames(PAAD6)=PAADnames6
matched<- match(PAADnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PAADfinal1<- data.frame(tratioL1HS[,matched])
PAADcut1<- PAAD1[,rownames(PAADfinal1)]
matched2<- match(PAADnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
PAADfinal2<- data.frame(tratioL1HS[,matched2])
PAADcut2<- PAAD2[,rownames(PAADfinal2)]
matched3<- match(PAADnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
PAADfinal3<- data.frame(tratioL1HS[,matched3])
PAADcut3<- PAAD3[,rownames(PAADfinal3)]
matched4<- match(PAADnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
PAADfinal4<- data.frame(tratioL1HS[,matched4])
PAADcut4<- PAAD4[,rownames(PAADfinal4)]
matched5<- match(PAADnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
PAADfinal5<- data.frame(tratioL1HS[,matched5])
PAADcut5<- PAAD5[,rownames(PAADfinal5)]
matched6<- match(PAADnames6, rownames(ratioL1HS))
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

#PRAD
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/PRAD")
PRAD1<- read.table(file="hgsc.bcm.edu_PRAD.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
PRAD2<- read.table(file="PR_TCGA_PRAD_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.TOTALCOUNT.tsv", header=TRUE)
PRAD3<- read.table(file="PRAD_Capture_All_Pairs_QCPASS_v6_Nikki_Nov_25.aggregated.capture.tcga.uuid.curated.somatic.TOTALCOU", header=TRUE)
PRAD4<- read.table(file="PRAD_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
PRAD5<- read.table(file="PR-TCGA-Analysis_set.aggregated.tcga.somatic.TOTALCOUNT.tsv", header=TRUE)
PRADnames1<- sub(pattern=".01", replacement="", x=colnames(PRAD1), perl=TRUE)
PRADnames2<- sub(pattern=".01", replacement="", x=colnames(PRAD2), perl=TRUE)
PRADnames3<- sub(pattern=".01", replacement="", x=colnames(PRAD3), perl=TRUE)
PRADnames4<- sub(pattern=".01", replacement="", x=colnames(PRAD4), perl=TRUE)
PRADnames5<- sub(pattern=".01", replacement="", x=colnames(PRAD5), perl=TRUE)
colnames(PRAD1)=PRADnames1
colnames(PRAD2)=PRADnames2
colnames(PRAD3)=PRADnames3
colnames(PRAD4)=PRADnames4
colnames(PRAD5)=PRADnames5
matched<- match(PRADnames1, rownames(ratioL1HS))
matched<- matched[!is.na(matched)]
PRADfinal1<- data.frame(tratioL1HS[,matched])
PRADcut1<- PRAD1[,rownames(PRADfinal1)]
matched2<- match(PRADnames2, rownames(ratioL1HS))
matched2<- matched2[!is.na(matched2)]
PRADfinal2<- data.frame(tratioL1HS[,matched2])
PRADcut2<- PRAD2[,rownames(PRADfinal2)]
matched3<- match(PRADnames3, rownames(ratioL1HS))
matched3<- matched3[!is.na(matched3)]
PRADfinal3<- data.frame(tratioL1HS[,matched3])
PRADcut3<- PRAD3[,rownames(PRADfinal3)]
matched4<- match(PRADnames4, rownames(ratioL1HS))
matched4<- matched4[!is.na(matched4)]
PRADfinal4<- data.frame(tratioL1HS[,matched4])
PRADcut4<- PRAD4[,rownames(PRADfinal4)]
matched5<- match(PRADnames5, rownames(ratioL1HS))
matched5<- matched5[!is.na(matched5)]
PRADfinal5<- data.frame(tratioL1HS[,matched5])
PRADcut5<- PRAD5[,rownames(PRADfinal5)]
PRADm1<- merge(PRADcut5, PRADcut4, all=TRUE)
PRADm2<- merge(PRADm1, PRADcut3, all=TRUE)
PRADm3<- merge(PRADm2, PRADcut2, all=TRUE)
PRADm4<- merge(PRADm3, PRADcut1, all=TRUE)

#SARC
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/SARC")
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

#THCA
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/THCA")
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

#THYM
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Counts of Mutations/THYM")
THYM1<- read.table(file="bcgsc.ca_THYM.IlluminaHiSeq_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM2<- read.table(file="genome.wustl.edu_THYM.IlluminaGA_DNASeq_curated.Level_2.1.1.0.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM3<- read.table(file="hgsc.bcm.edu_THYM.IlluminaGA_DNASeq.1.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM4<- read.table(file="THYM_pairs.aggregated.capture.tcga.uuid.automated.somatic.TOTALCOUNT.tsv", header=TRUE)
THYM5<- read.table(file="ucsc.edu_THYM.IlluminaGA_DNASeq_automated.Level_2.1.0.0.somatic.TOTALCOUNT.tsv", header=TRUE)
THYMnames1<- sub(pattern=".01", replacement="", x=colnames(THYM1), perl=TRUE)
THYMnames2<- sub(pattern=".01", replacement="", x=colnames(THYM2), perl=TRUE)
THYMnames3<- sub(pattern=".01", replacement="", x=colnames(THYM3), perl=TRUE)
THYMnames4<- sub(pattern=".01", replacement="", x=colnames(THYM4), perl=TRUE)
THYMnames5<- sub(pattern=".01", replacement="", x=colnames(THYM5), perl=TRUE)
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
UCECnames1<- sub(pattern=".01", replacement="", x=colnames(UCEC1), perl=TRUE)
UCECnames2<- sub(pattern=".01", replacement="", x=colnames(UCEC2), perl=TRUE)
UCECnames3<- sub(pattern=".01", replacement="", x=colnames(UCEC3), perl=TRUE)
UCECnames4<- sub(pattern=".01", replacement="", x=colnames(UCEC4), perl=TRUE)
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
#final merge
m1<- cbind(BLCAm4, BRCAm4)
#check for matches
match(colnames(CESCm6), colnames(m1))
match(colnames(m1), colnames(CHOLm4))
match(colnames(m1), colnames(ESCAm4))
match(colnames(m1), colnames(FPPPcut1))
match(colnames(m1), colnames(HNSCm4))
match(colnames(m1), colnames(KIRPm6))
match(colnames(m1), colnames(LIHCm5))
match(colnames(m1), colnames(PAADm5))
match(colnames(m1), colnames(PCPGm3))
match(colnames(m1), colnames(PRADm4))
match(colnames(m1), colnames(SARCm3))
match(colnames(m1), colnames(THCAm4))
match(colnames(m1), colnames(THYMm4))
match(colnames(m1), colnames(UCECm3))

match(colnames(CESCm6), colnames(CHOLm4))
match(colnames(CESCm6), colnames(ESCAm4))
match(colnames(CESCm6), colnames(FPPPcut1))
match(colnames(CESCm6), colnames(HNSCm4))
match(colnames(CESCm6), colnames(KIRPm6))
match(colnames(CESCm6), colnames(LIHCm5))
match(colnames(CESCm6), colnames(PAADm5))
match(colnames(CESCm6), colnames(PCPGm3))
match(colnames(CESCm6), colnames(PRADm4))
match(colnames(CESCm6), colnames(SARCm3))
match(colnames(CESCm6), colnames(THCAm4))
match(colnames(CESCm6), colnames(THYMm4))
match(colnames(CESCm6), colnames(UCECm3))

match(colnames(CHOLm4), colnames(ESCAm4))
match(colnames(CHOLm4), colnames(FPPPcut1))
match(colnames(CHOLm4), colnames(HNSCm4))
match(colnames(CHOLm4), colnames(KIRPm6))
match(colnames(CHOLm4), colnames(LIHCm5))
match(colnames(CHOLm4), colnames(PAADm5))
match(colnames(CHOLm4), colnames(PCPGm3))
match(colnames(CHOLm4), colnames(PRADm4))
match(colnames(CHOLm4), colnames(SARCm3))
match(colnames(CHOLm4), colnames(THCAm4))
match(colnames(CHOLm4), colnames(THYMm4))
match(colnames(CHOLm4), colnames(UCECm3))

match(colnames(ESCAm4), colnames(FPPPcut1))
match(colnames(ESCAm4), colnames(HNSCm4))
match(colnames(ESCAm4), colnames(KIRPm6))
match(colnames(ESCAm4), colnames(LIHCm5))
match(colnames(ESCAm4), colnames(PAADm5))
match(colnames(ESCAm4), colnames(PCPGm3))
match(colnames(ESCAm4), colnames(PRADm4))
match(colnames(ESCAm4), colnames(SARCm3))
match(colnames(ESCAm4), colnames(THCAm4))
match(colnames(ESCAm4), colnames(THYMm4))
match(colnames(ESCAm4), colnames(UCECm3))

match(colnames(FPPPcut1), colnames(HNSCm4))
match(colnames(FPPPcut1), colnames(KIRPm6))
match(colnames(FPPPcut1), colnames(LIHCm5))
match(colnames(FPPPcut1), colnames(PAADm5))
match(colnames(FPPPcut1), colnames(PCPGm3))
match(colnames(FPPPcut1), colnames(PRADm4))
match(colnames(FPPPcut1), colnames(SARCm3))
match(colnames(FPPPcut1), colnames(THCAm4))
match(colnames(FPPPcut1), colnames(THYMm4))
match(colnames(FPPPcut1), colnames(UCECm3))

match(colnames(HNSCm4), colnames(KIRPm6))
match(colnames(HNSCm4), colnames(LIHCm5))
match(colnames(HNSCm4), colnames(PAADm5))
match(colnames(HNSCm4), colnames(PCPGm3))
match(colnames(HNSCm4), colnames(PRADm4))
match(colnames(HNSCm4), colnames(SARCm3))
match(colnames(HNSCm4), colnames(THCAm4))
match(colnames(HNSCm4), colnames(THYMm4))
match(colnames(HNSCm4), colnames(UCECm3))

match(colnames(KIRPm6), colnames(LIHCm5))
match(colnames(KIRPm6), colnames(PAADm5))
match(colnames(KIRPm6), colnames(PCPGm3))
match(colnames(KIRPm6), colnames(PRADm4))
match(colnames(KIRPm6), colnames(SARCm3))
match(colnames(KIRPm6), colnames(THCAm4))
match(colnames(KIRPm6), colnames(THYMm4))
match(colnames(KIRPm6), colnames(UCECm3))

match(colnames(LIHCm5), colnames(PAADm5))
match(colnames(LIHCm5), colnames(PCPGm3))
match(colnames(LIHCm5), colnames(PRADm4))
match(colnames(LIHCm5), colnames(SARCm3))
match(colnames(LIHCm5), colnames(THCAm4))
match(colnames(LIHCm5), colnames(THYMm4))
match(colnames(LIHCm5), colnames(UCECm3))

match(colnames(PAADm5), colnames(PCPGm3))
match(colnames(PAADm5), colnames(PRADm4))
match(colnames(PAADm5), colnames(SARCm3))
match(colnames(PAADm5), colnames(THCAm4))
match(colnames(PAADm5), colnames(THYMm4))
match(colnames(PAADm5), colnames(UCECm3))

match(colnames(PCPGm3), colnames(PRADm4))
match(colnames(PCPGm3), colnames(SARCm3))
match(colnames(PCPGm3), colnames(THCAm4))
match(colnames(PCPGm3), colnames(THYMm4))
match(colnames(PCPGm3), colnames(UCECm3))

match(colnames(PRADm4), colnames(SARCm3))
match(colnames(PRADm4), colnames(THCAm4))
match(colnames(PRADm4), colnames(THYMm4))
match(colnames(PRADm4), colnames(UCECm3))

match(colnames(SARCm3), colnames(THCAm4))
match(colnames(SARCm3), colnames(THYMm4))
match(colnames(SARCm3), colnames(UCECm3))

match(colnames(THCAm4), colnames(THYMm4))
match(colnames(THCAm4), colnames(UCECm3))

match(colnames(THYMm4), colnames(UCECm3))
#create averages by cancer/normal
#BRCA/BLCA
notna<- colSums(!is.na(m1))
sum<- colSums(m1, na.rm=TRUE)
avmutBLCABRCA<- data.frame(sum/notna)
#CESC
notna<- colSums(!is.na(CESCm6))
sum<- colSums(CESCm6, na.rm=TRUE)
avmutCESC<- data.frame(sum/notna)
#CHOL
notna<- colSums(!is.na(CHOLm4))
sum<- colSums(CHOLm4, na.rm=TRUE)
avmutCHOL<- data.frame(sum/notna)
#ESCA
notna<- colSums(!is.na(ESCAm4))
sum<- colSums(ESCAm4, na.rm=TRUE)
avmutESCA<- data.frame(sum/notna)
#FPPPcut1
notna<- colSums(!is.na(FPPPcut1))
sum<- colSums(FPPPcut1, na.rm=TRUE)
avmutFPPP<- data.frame(sum/notna)
#HNSCm4
notna<- colSums(!is.na(HNSCm4))
sum<- colSums(HNSCm4, na.rm=TRUE)
avmutHNSC<- data.frame(sum/notna)
#KIRP
notna<- colSums(!is.na(KIRPm6))
sum<- colSums(KIRPm6, na.rm=TRUE)
avmutKIRP<- data.frame(sum/notna)
#LIHC
notna<- colSums(!is.na(LIHCm5))
sum<- colSums(LIHCm5, na.rm=TRUE)
avmutLIHC<- data.frame(sum/notna)
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
mergemut4<- rbind(mergemut3, avmutFPPP)
mergemut5<- rbind(mergemut4, avmutHNSC)
mergemut6<- rbind(mergemut5, avmutKIRP)
mergemut7<- rbind(mergemut6, avmutLIHC)
mergemut8<- rbind(mergemut7, avmutPAAD)
mergemut9<- rbind(mergemut8, avmutPCPG)
mergemut10<- rbind(mergemut9, avmutPRAD)
mergemut11<- rbind(mergemut10, avmutSARC)
mergemut12<- rbind(mergemut11, avmutTHCA)
mergemut13<- rbind(mergemut12, avmutTHYM)
mergemut14<- rbind(mergemut13, avmutUCEC)
#To create the L1HS file of all the patients
logmut<- data.frame(mergemut14)
tlogmut<- t(logmut)
matchedL<- match(rownames(logmut), rownames(ratioL1HS))
matchedL<- matchedL[!is.na(matchedL)]
L1HS<- data.frame(tratioL1HS[,matchedL])
L1HSALL<- data.frame(log2(L1HS))
matchedmut<- match(rownames(ratioL1HS), rownames(logmut))
matchedmut<- matchedmut[!is.na(matchedmut)]
final<- data.frame(tlogmut[,matchedmut])
finallogmut<- data.frame(log2(final))
cor(L1HSALL[,1], finallogmut[,1])
value<- lm(L1HSALL[,1] ~ finallogmut[,1])
plot(finallogmut[,1], L1HSALL[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="ALL Cancer Correlation", col="magenta")
abline(lm(L1HSALL[,1] ~ finallogmut[,1]))

#pvalues
p<- anova(value)$Pr[1]
#qvalues
q<- p.adjust(p, method="BH")
r2<- summary(value)$r.squared
test<- resid(lm(finallogmut[,1] ~ L1HSALL[,1]))
plot(finallogmut[,1], test, main="Residuals")
abline(0,0)
plot(L1HSALL[,1], test, main="Residuals")
abline(0,0)
qqnorm(finallogmut[,1], main="Q-Q Plot of mutations")
qqnorm(L1HSALL[,1], main="Q-Q Plot of L1HS")

#cut data down from the extremities
tL1HSALL<- t(L1HSALL)
cutL1HS<- subset(L1HSALL, L1HSALL[,1] > -2 & L1HSALL[,1] < 3)
matched<- match(rownames(cutL1HS), rownames(L1HSALL))
matched<- matched[!is.na(matched)]
finalcutL1HS<- data.frame(tL1HSALL[,matched])
matched2<- match(rownames(finalcutL1HS), rownames(finallogmut))
matched2<- matched2[!is.na(matched2)]
tfinallogmut<- t(finallogmut)
cutmut<- data.frame(tfinallogmut[,matched2])
cut2<- cbind(cutL1HS, cutmut)
colnames(cut2)=c("L1HS", "mutations")
plot(cut2[,2], cut2[,1])
abline(lm(cut2[,1] ~ cut2[,2]))
cutmut2<- subset(cutmut, cutmut[,1] > 3 & cutmut[,1] < 10)
matched3<- match(rownames(cutmut2), rownames(cutL1HS))
matched3<- matched3[!is.na(matched3)] 
tcutL1HS<- t(cutL1HS)
tcutmut<- t(cutmut)
cutL1HS2<- data.frame(tcutL1HS[,matched3])
cutfinal<- cbind(cutL1HS2, cutmut2)
colnames(cutfinal)= c("L1HS", "mutations")
plot(cutfinal[,2], cutfinal[,1], xlab="mutations", ylab="L1HS", main="Correlation w/out extremities", col="magenta")
abline(lm(cutfinal[,1] ~ cutfinal[,2]))
value<- lm(cutfinal[,1] ~ cutfinal[,2])
p<- anova(value)$Pr[1]
r2<- summary(value)$r.squared
qqnorm(cutfinal[,2])
#Fixing by combing the files
combined<- merge(tL1HSALL, tfinallogmut, all=TRUE)
tcombined<- t(combined)
colnames(tcombined)= c("L1HS", "mutation average")
plot(tcombined[,2], tcombined[,1], xlab="log2 average mutation frequency", ylab="log2 fold change L1HS", main="ALL Cancer Correlation", col="magenta")
abline(lm(tcombined[,1] ~ tcombined[,2]))
value<- lm(tcombined[,1] ~ tcombined[,2])
p<- anova(value)$Pr[1]
r2<- summary(value)$r.squared
test<- resid(lm(tcombined[,2] ~ tcombined[,1]))
plot(tcombined[,2], test, main="Residuals", xlab="mutation")
abline(0,0)
plot(tcombined[,1], test, main="Residuals", xlab="L1HS")
abline(0,0)
qqnorm(tcombined[,2], main="Q-Q Plot of mutations")
qqnorm(tcombined[,1], main="Q-Q Plot of L1HS")


















m2<- cbind(m1, CESCm6)
m3<- merge(m2, CHOLm4, all=TRUE)
m4<- merge(m3, ESCAm4, all=TRUE)
m5<- merge(m4, FPPPcut1, all=TRUE)
m6<- merge(m5, HNSCm4, all=TRUE)
m7<- merge(m6, KIRPm6, all=TRUE)
m8<- merge(m7, LIHCm5, all=TRUE)
m9<- merge(m8, PAADm5, all=TRUE)
m10<- merge(m9, PCPGm3, all=TRUE)
m11<- merge(m10, PRADm4, all=TRUE)
m12<- merge(m11, SARCm3, all=TRUE)
m13<- merge(m12, THCAm4, all=TRUE)
m14<- merge(m13, THYMm4, all=TRUE)
m15<- merge(m14, UCECm3, all=TRUE)

bind.pad <- function(l, side="r", len=max(sapply(l,length)))
{
  if (side %in% c("b", "r")) {
    out <- sapply(l, 'length<-', value=len)
  } else {
    out <- sapply(sapply(sapply(l, rev), 'length<-', value=len, simplify=F), rev)}
  if (side %in% c("r", "l")) out <- t(out)
  out
}

setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Testing/merged files")
write.table(BLCAm4, file="BLCA merged.txt", sep="	")
write.table(BRCAm4, file="BRCA merged.txt", sep="	")
write.table(CESCm6, file="CESC merged.txt", sep="	")
write.table(CHOLm4, file="CHOL merged.txt", sep="	")
write.table(ESCAm4, file="ESCA merged.txt", sep="	")
write.table(FPPPcut1, file="FPPP merged.txt", sep="	")
write.table(HNSCm4, file="HNSC merged.txt", sep="	")
write.table(KIRPm6, file="KIRP merged.txt", sep="	")
write.table(LIHCm5, file="LIHC merged.txt", sep="	")
write.table(PAADm5, file="PAAD merged.txt", sep="	")
write.table(PCPGm3, file="PCPG merged.txt", sep="	")
write.table(PRADm4, file="PRAD merged.txt", sep="	")
write.table(SARCm3, file="SARC merged.txt", sep="	")
write.table(THCAm4, file="THCA merged.txt", sep="	")
write.table(THYMm4, file="THYM merged.txt", sep="	")
write.table(UCECm3, file="UCEC merged.txt", sep="	")
write.table(tcombined, file="L1HS and Mutation Averages Logged", sep=" ")

setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Testing/merged files")
mutationfiles <- dir(pattern = "*.txt", full.names = TRUE)
lst <- lapply(mutationfiles, read.table, header=TRUE, sep='', stringsAsFactors=FALSE)
mirna<- sapply(lst, '[[',3)

















