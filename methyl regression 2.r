####methylation regression second run
####first find a set of coordinates that line up and can work
setwd("/home/quints1/Patient_Names/BRCAfiles/")
level2<- read.table("TCGA-E2-A1L7.normal.txt", header=TRUE, fill=TRUE, sep="\t", skip=1)
setwd("/home/quints1/Patient_Names/BRCA/raw_patients/")
level3<- read.table("jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0C3-11A-23D-A12R-05.txt", header=TRUE, fill=TRUE, skip=1, sep="\t")
test<- level2[which(level2[,1] %in% level3[,1]),]
test2<- level3[which(level2[,1] %in% level3[,1]),]
cuttest2<- test2[,c(1,4,5)]
together<- data.frame(cbind(test, cuttest2))
#Beta= methylated allele intensity (M)/ (Unmethylated allele intensity (U) + Methylated allele intensity (M) + 100)
Betameth<- together[,2]/(together[,3] + together[,2] + 100)
together2<- cbind(together, data.frame(Betameth))
finalcut<- together2[,c(1,2,3,8,6,7)]
patientcut<- finalcut[,c(5, 6, 6, 4)]
library(stringi)
x<- "chr"
stri_sub(patientcut[,1], 0, 0)<- x
setwd("/home/quints1/BRCAfiles/editedfiles")
write.table(patientcut, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#####create a file in the server:
setwd("/home/quints1/Patient_Names/BRCAfiles/")
args = commandArgs(trailingOnly= TRUE)
inputfile= read.table(args[1], header=TRUE, fill=TRUE, sep="\t", skip=1)
setwd("/home/quints1/Patient_Names/BRCA/raw_patients/")
level3<- read.table("jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-BH-A0C3-11A-23D-A12R-05.txt", header=TRUE, fill=TRUE, skip=1, sep="\t")
test<- level2[which(level2[,1] %in% level3[,1]),]
test2<- level3[which(level2[,1] %in% level3[,1]),]
cuttest2<- test2[,c(1,4,5)]
together<- cbind(test, cuttest2)
#Beta= methylated allele intensity (M)/ (Unmethylated allele intensity (U) + Methylated allele intensity (M) + 100)
Betameth<- together[,2]/(together[,3] + together[,2] + 100)
together2<- cbind(together, data.frame(Betameth))
finalcut<- together2[,c(1,2,3,8,6,7)]
patientcut<- finalcut[,c(5, 6, 6, 4)]
library(stringi)
x<- "chr"
stri_sub(patientcut[,1], 0, 0)<- x
write.table(patientcut, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
