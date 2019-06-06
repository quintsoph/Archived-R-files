####To compare the patient Beta Values to TE information
###I am assuming that we will have all the patients that cross over
#BLCA
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/bedresults/")
BLCAnames<- read.table("BLCAbedpatients.txt", header=TRUE)
library(stringi)
newBLCA<- as.matrix(stri_sub(BLCAnames[,1], 1, -28))
allBLCA=NULL
for (i in newBLCA)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allBLCA<- rbind(allBLCA, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
BLCAsub<- as.matrix(combined[(combined[,1] %in% allBLCA[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSBLCA=NULL
for (i in BLCAsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSBLCA<- data.frame(cbind(L1HSBLCA, dupsout))
	print(L1HSBLCA)
}

rownames(L1HSBLCA)=L1HSBLCA[,1]
L1HSBLCA2<- L1HSBLCA[,c(2, 4, 6 ,8, 10, 12, 14, 16, 18, 20, 22, 24)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
write.table(L1HSBLCA2, "BLCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#read in exon 105312870
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon1<- read.table("105312870 exon.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% BLCAsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSBLCA2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
ttL1HS52<- cbind(rownames(tL1HS52), tL1HS52)
exon52<- cbind(tL1HS52, outbetasexon52)
finalexon52<- exon52[,c(1,3)]
colnames(finalexon52)= c("expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
write.table(finalexon52, "BLCAexon52.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
###plot
library("ggplot2")
legend_title= "BLCA"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="BLCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=2) +
   labs(x = "expression levels", y = "beta values", title="L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon dup1134
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon2<- read.table("90699440 exon.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% BLCAsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSBLCA2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
ttL1HS1134<- cbind(rownames(tL1HS1134), tL1HS1134)
exon1134<- cbind(tL1HS1134, outbetasexon1134)
finalexon1134<- exon1134[,c(1,3)]
colnames(finalexon1134)= c("expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
write.table(finalexon1134, "BLCAexon1134.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
finalexon1134<- read.table("BLCAexon1134.txt")
library("ggplot2")
legend_title= "BLCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="BLCA")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BLCA L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon3<- read.table("73787793 exon.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% BLCAsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSBLCA2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
ttL1HS795<- cbind(rownames(tL1HS795), tL1HS795)
exon795<- cbind(tL1HS795, outbetasexon795)
finalexon795<- exon795[,c(1,3)]
colnames(finalexon795)= c("expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
write.table(finalexon795, "BLCAexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
finalexon795<- read.table("BLCAexon795.txt")
library("ggplot2")
legend_title= "BLCA"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="BLCA")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BLCA L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon4<- read.table("16943900 exon.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% BLCAsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSBLCA2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
write.table(exon355, "BLCAexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
finalexon355<- read.table("BLCAexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "BLCA"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="BLCA")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BLCA L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon5<- read.table("13258876 exon.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% BLCAsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSBLCA2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
write.table(exon1090, "BLCAexon1090.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
finalexon1090<- read.table("BLCAexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "BLCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="BLCA")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BLCA L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

####BRCA
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/bedresults/")
BRCAnames<- read.table("BRCAbedpatients.txt")
library(stringi)
newBRCA<- as.matrix(stri_sub(BRCAnames[,1], 1, -28))
allBRCA=NULL
for (i in newBRCA)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allBRCA<- rbind(allBRCA, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
BRCAsub<- as.matrix(combined[(combined[,1] %in% allBRCA[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSBRCA=NULL
for (i in BRCAsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSBRCA<- data.frame(cbind(L1HSBRCA, dupsout))
	print(L1HSBRCA)
}

rownames(L1HSBRCA)=L1HSBRCA[,1]
L1HSBRCA2<- data.frame(L1HSBRCA[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
62, 64, 66, 68, 70, 72, 74, 76, 78, 80,
82, 84, 86, 88, 90, 92, 94, 96, 98, 100,
102, 104, 106, 108, 110, 112, 114, 116, 118, 120,
122, 124, 126, 128, 130, 132, 134, 136, 138, 140,
142, 144, 146, 148, 150, 152, 154, 156, 158, 160,
162, 164, 166, 168, 170, 172)])
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
write.table(L1HSBRCA2, "BRCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#read in exon 105312870
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% BRCAsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSBRCA2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
write.table(finalexon52, "BRCAexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
finalexon52<- read.table("BRCAexon52.txt", header=TRUE)
legend_title= "BRCA"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="BRCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BRCA L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% BRCAsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSBRCA2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
write.table(finalexon1134, "BRCAexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
finalexon1134<- read.table("BRCAexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "BRCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="BRCA")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BRCA L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% BRCAsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSBRCA2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
write.table(finalexon795, "BRCAexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
finalexon795<- read.table("BRCAexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "BRCA"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="BRCA")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BRCA L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% BRCAsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSBRCA2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
write.table(exon355, "BRCAexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
finalexon355<- read.table("BRCAexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "BRCA"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="BRCA")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BRCA L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% BRCAsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSBRCA2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
write.table(exon1090, "BRCAexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
finalexon1090<- read.table("BRCAexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "BRCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="BRCA")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="BRCA L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

####COAD
setwd("/home/quints1/Patient_Names/COAD/edited_patients/bedresults/normal_patients")
COADnames<- read.table("COADbedpatientsn.txt")
library(stringi)
newCOAD<- as.matrix(stri_sub(COADnames[,1], 1, -28))
allCOAD=NULL
for (i in newCOAD)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allCOAD<- rbind(allCOAD, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
COADsub<- as.matrix(combined[(combined[,1] %in% allCOAD[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSCOAD=NULL
for (i in COADsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSCOAD<- data.frame(cbind(L1HSCOAD, dupsout))
	print(L1HSCOAD)
}

rownames(L1HSCOAD)=L1HSCOAD[,1]
L1HSCOAD2<- L1HSCOAD[,c(2, 4, 6 ,8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
write.table(L1HSCOAD2, "COADmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% COADsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSCOAD2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
write.table(finalexon52, "COADexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
finalexon52<- read.table("COADexon52.txt", header=TRUE)
legend_title= "COAD"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="COAD")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="COAD L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% COADsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSCOAD2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
write.table(finalexon1134, "COADexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
finalexon1134<- read.table("COADexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "COAD"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="COAD")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="COAD L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% COADsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSCOAD2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
write.table(finalexon795, "COADexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
finalexon795<- read.table("COADexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "COAD"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="COAD")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="COAD L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% COADsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSCOAD2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
write.table(exon355, "COADexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
finalexon355<- read.table("COADexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "COAD"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="COAD")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="COAD L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% COADsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSCOAD2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
write.table(exon1090, "COADexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
finalexon1090<- read.table("COADexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "COAD"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="COAD")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="COAD L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

######ESCA
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/bedresults/normal_patients")
ESCAnames<- read.table("ESCAbedpatientsn.txt")
library(stringi)
newESCA<- as.matrix(stri_sub(ESCAnames[,1], 1, -28))
allESCA=NULL
for (i in newESCA)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allESCA<- rbind(allESCA, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
ESCAsub<- as.matrix(combined[(combined[,1] %in% allESCA[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSESCA=NULL
for (i in ESCAsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSESCA<- data.frame(cbind(L1HSESCA, dupsout))
	print(L1HSESCA)
}

rownames(L1HSESCA)=L1HSESCA[,1]
L1HSESCA2<- L1HSESCA[,c(2, 4, 6)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
write.table(L1HSESCA2, "ESCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% ESCAsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSESCA2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
write.table(finalexon52, "ESCAexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
finalexon52<- read.table("ESCAexon52.txt", header=TRUE)
legend_title= "ESCA"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="ESCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ESCA L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% ESCAsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSESCA2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
write.table(finalexon1134, "ESCAexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
finalexon1134<- read.table("ESCAexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "ESCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="ESCA")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ESCA L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% ESCAsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSESCA2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
write.table(finalexon795, "ESCAexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
finalexon795<- read.table("ESCAexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "ESCA"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="ESCA")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ESCA L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% ESCAsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSESCA2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
write.table(exon355, "ESCAexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
finalexon355<- read.table("ESCAexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "ESCA"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="ESCA")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ESCA L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% ESCAsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSESCA2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
write.table(exon1090, "ESCAexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ESCA")
finalexon1090<- read.table("ESCAexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "ESCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="ESCA")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ESCA L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

######HNSC
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/bedresults/normal_patients")
HNSCnames<- read.table("HNSCbedpatientsn.txt")
library(stringi)
newHNSC<- as.matrix(stri_sub(HNSCnames[,1], 1, -28))
allHNSC=NULL
for (i in newHNSC)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allHNSC<- rbind(allHNSC, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
HNSCsub<- as.matrix(combined[(combined[,1] %in% allHNSC[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSHNSC=NULL
for (i in HNSCsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSHNSC<- data.frame(cbind(L1HSHNSC, dupsout))
	print(L1HSHNSC)
}

rownames(L1HSHNSC)=L1HSHNSC[,1]
L1HSHNSC2<- L1HSHNSC[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
write.table(L1HSHNSC2, "HNSCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% HNSCsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSHNSC2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
write.table(finalexon52, "HNSCexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
finalexon52<- read.table("HNSCexon52.txt", header=TRUE)
legend_title= "HNSC"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="HNSC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="HNSC L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% HNSCsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSHNSC2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
write.table(finalexon1134, "HNSCexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
finalexon1134<- read.table("HNSCexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "HNSC"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="HNSC")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="HNSC L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% HNSCsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSHNSC2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
write.table(finalexon795, "HNSCexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
finalexon795<- read.table("HNSCexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "HNSC"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="HNSC")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="HNSC L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% HNSCsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSHNSC2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
write.table(exon355, "HNSCexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
finalexon355<- read.table("HNSCexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "HNSC"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="HNSC")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="HNSC L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% HNSCsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSHNSC2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
write.table(exon1090, "HNSCexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
finalexon1090<- read.table("HNSCexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "HNSC"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="HNSC")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="HNSC L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

########KIRC
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/bedresults/normal_patients")
KIRCnames<- read.table("KIRCbedpatientsn.txt")
library(stringi)
newKIRC<- as.matrix(stri_sub(KIRCnames[,1], 1, -28))
allKIRC=NULL
for (i in newKIRC)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allKIRC<- rbind(allKIRC, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
KIRCsub<- as.matrix(combined[(combined[,1] %in% allKIRC[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSKIRC=NULL
for (i in KIRCsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSKIRC<- data.frame(cbind(L1HSKIRC, dupsout))
	print(L1HSKIRC)
}

rownames(L1HSKIRC)=L1HSKIRC[,1]
L1HSKIRC2<- L1HSKIRC[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
write.table(L1HSKIRC2, "KIRCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% KIRCsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSKIRC2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
write.table(finalexon52, "KIRCexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
finalexon52<- read.table("KIRCexon52.txt", header=TRUE)
legend_title= "KIRC"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="KIRC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRC L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% KIRCsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSKIRC2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
write.table(finalexon1134, "KIRCexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
finalexon1134<- read.table("KIRCexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRC"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="KIRC")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRC L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% KIRCsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSKIRC2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
write.table(finalexon795, "KIRCexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
finalexon795<- read.table("KIRCexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRC"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="KIRC")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRC L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% KIRCsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSKIRC2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
write.table(exon355, "KIRCexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
finalexon355<- read.table("KIRCexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRC"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="KIRC")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRC L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% KIRCsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSKIRC2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
write.table(exon1090, "KIRCexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
finalexon1090<- read.table("KIRCexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRC"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="KIRC")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRC L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

############KIRP
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/bedresults/normal_patients")
KIRPnames<- read.table("KIRPbedpatientsn.txt")
library(stringi)
newKIRP<- as.matrix(stri_sub(KIRPnames[,1], 1, -28))
allKIRP=NULL
for (i in newKIRP)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allKIRP<- rbind(allKIRP, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
KIRPsub<- as.matrix(combined[(combined[,1] %in% allKIRP[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSKIRP=NULL
for (i in KIRPsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSKIRP<- data.frame(cbind(L1HSKIRP, dupsout))
	print(L1HSKIRP)
}

rownames(L1HSKIRP)=L1HSKIRP[,1]
L1HSKIRP2<- L1HSKIRP[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
write.table(L1HSKIRP2, "KIRPmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% KIRPsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSKIRP2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
write.table(finalexon52, "KIRPexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
finalexon52<- read.table("KIRPexon52.txt", header=TRUE)
legend_title= "KIRP"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="KIRP")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRP L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% KIRPsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSKIRP2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
write.table(finalexon1134, "KIRPexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
finalexon1134<- read.table("KIRPexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRP"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="KIRP")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRP L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% KIRPsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSKIRP2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
write.table(finalexon795, "KIRPexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
finalexon795<- read.table("KIRPexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRP"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="KIRP")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRP L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% KIRPsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSKIRP2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
write.table(exon355, "KIRPexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
finalexon355<- read.table("KIRPexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRP"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="KIRP")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRP L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% KIRPsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSKIRP2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
write.table(exon1090, "KIRPexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
finalexon1090<- read.table("KIRPexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "KIRP"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="KIRP")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="KIRP L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

##########LIHC
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/bedresults/normal_patients")
LIHCnames<- read.table("LIHCbedpatientsn.txt")
library(stringi)
newLIHC<- as.matrix(stri_sub(LIHCnames[,1], 1, -28))
allLIHC=NULL
for (i in newLIHC)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allLIHC<- rbind(allLIHC, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
LIHCsub<- as.matrix(combined[(combined[,1] %in% allLIHC[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSLIHC=NULL
for (i in LIHCsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSLIHC<- data.frame(cbind(L1HSLIHC, dupsout))
	print(L1HSLIHC)
}

rownames(L1HSLIHC)=L1HSLIHC[,1]
L1HSLIHC2<- L1HSLIHC[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
62, 64, 66, 68, 70, 72, 74, 76, 78, 80,
82)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
write.table(L1HSLIHC2, "LIHCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% LIHCsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSLIHC2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
write.table(finalexon52, "LIHCexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
finalexon52<- read.table("LIHCexon52.txt", header=TRUE)
legend_title= "LIHC"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="LIHC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LIHC L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% LIHCsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSLIHC2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
write.table(finalexon1134, "LIHCexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
finalexon1134<- read.table("LIHCexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "LIHC"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="LIHC")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LIHC L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% LIHCsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSLIHC2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
write.table(finalexon795, "LIHCexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
finalexon795<- read.table("LIHCexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "LIHC"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="LIHC")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LIHC L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% LIHCsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSLIHC2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
write.table(exon355, "LIHCexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
finalexon355<- read.table("LIHCexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "LIHC"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="LIHC")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LIHC L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% LIHCsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSLIHC2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
write.table(exon1090, "LIHCexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LIHC")
finalexon1090<- read.table("LIHCexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "LIHC"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="LIHC")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LIHC L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#############LUAD
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/bedresults/normal_patients")
LUADnames<- read.table("LUADbedpatientsn.txt")
library(stringi)
newLUAD<- as.matrix(stri_sub(LUADnames[,1], 1, -28))
allLUAD=NULL
for (i in newLUAD)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allLUAD<- rbind(allLUAD, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
LUADsub<- as.matrix(combined[(combined[,1] %in% allLUAD[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSLUAD=NULL
for (i in LUADsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSLUAD<- data.frame(cbind(L1HSLUAD, dupsout))
	print(L1HSLUAD)
}

rownames(L1HSLUAD)=L1HSLUAD[,1]
L1HSLUAD2<- L1HSLUAD[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
write.table(L1HSLUAD2, "LUADmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% LUADsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSLUAD2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
write.table(finalexon52, "LUADexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
finalexon52<- read.table("LUADexon52.txt", header=TRUE)
legend_title= "LUAD"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="LUAD")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUAD L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% LUADsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSLUAD2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
write.table(finalexon1134, "LUADexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
finalexon1134<- read.table("LUADexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "LUAD"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="LUAD")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUAD L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% LUADsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSLUAD2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
write.table(finalexon795, "LUADexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
finalexon795<- read.table("LUADexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "LUAD"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="LUAD")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUAD L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% LUADsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSLUAD2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
write.table(exon355, "LUADexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
finalexon355<- read.table("LUADexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "LUAD"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="LUAD")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUAD L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% LUADsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSLUAD2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
write.table(exon1090, "LUADexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
finalexon1090<- read.table("LUADexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "LUAD"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="LUAD")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUAD L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#############LUSC
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/bedresults/normal_patients")
LUSCnames<- read.table("LUSCbedpatientsn.txt")
library(stringi)
newLUSC<- as.matrix(stri_sub(LUSCnames[,1], 1, -28))
allLUSC=NULL
for (i in newLUSC)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allLUSC<- rbind(allLUSC, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
LUSCsub<- as.matrix(combined[(combined[,1] %in% allLUSC[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSLUSC=NULL
for (i in LUSCsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSLUSC<- data.frame(cbind(L1HSLUSC, dupsout))
	print(L1HSLUSC)
}

rownames(L1HSLUSC)=L1HSLUSC[,1]
L1HSLUSC2<- L1HSLUSC[,c(2, 4, 6, 8, 10, 12, 14, 16)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
write.table(L1HSLUSC2, "LUSCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% LUSCsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSLUSC2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
write.table(finalexon52, "LUSCexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
finalexon52<- read.table("LUSCexon52.txt", header=TRUE)
legend_title= "LUSC"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="LUSC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUSC L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% LUSCsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSLUSC2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
write.table(finalexon1134, "LUSCexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
finalexon1134<- read.table("LUSCexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "LUSC"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="LUSC")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUSC L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% LUSCsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSLUSC2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
write.table(finalexon795, "LUSCexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
finalexon795<- read.table("LUSCexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "LUSC"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="LUSC")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUSC L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% LUSCsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSLUSC2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
write.table(exon355, "LUSCexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
finalexon355<- read.table("LUSCexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "LUSC"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="LUSC")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUSC L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% LUSCsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSLUSC2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
write.table(exon1090, "LUSCexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUSC")
finalexon1090<- read.table("LUSCexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "LUSC"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="LUSC")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="LUSC L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#############READ
setwd("/home/quints1/Patient_Names/READ/edited_patients/bedresults/normal_patients")
READnames<- read.table("READpatientsn.txt")
library(stringi)
newREAD<- as.matrix(stri_sub(READnames[,1], 1, -28))
allREAD=NULL
for (i in newREAD)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allREAD<- rbind(allREAD, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
READsub<- as.matrix(combined[(combined[,1] %in% allREAD[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSREAD=NULL
for (i in READsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSREAD<- data.frame(cbind(L1HSREAD, dupsout))
	print(L1HSREAD)
}

rownames(L1HSREAD)=L1HSREAD[,1]
L1HSREAD2<- L1HSREAD[,c(2, 4)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
write.table(L1HSREAD2, "READmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% READsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSREAD2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
write.table(finalexon52, "READexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
finalexon52<- read.table("READexon52.txt", header=TRUE)
legend_title= "READ"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="READ")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="READ L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% READsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSREAD2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
write.table(finalexon1134, "READexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
finalexon1134<- read.table("READexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "READ"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="READ")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="READ L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% READsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSREAD2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
write.table(finalexon795, "READexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
finalexon795<- read.table("READexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "READ"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="READ")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="READ L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% READsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSREAD2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
write.table(exon355, "READexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
finalexon355<- read.table("READexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "READ"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="READ")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="READ L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% READsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSREAD2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
write.table(exon1090, "READexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/READ")
finalexon1090<- read.table("READexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "READ"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="READ")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="READ L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

############THCA
setwd("/home/quints1/Patient_Names/THCA/edited_patients/bedresults/normal_patients")
THCAnames<- read.table("THCAbedpatientsn.txt")
library(stringi)
newTHCA<- as.matrix(stri_sub(THCAnames[,1], 1, -28))
allTHCA=NULL
for (i in newTHCA)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	allTHCA<- rbind(allTHCA, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
THCAsub<- as.matrix(combined[(combined[,1] %in% allTHCA[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSTHCA=NULL
for (i in THCAsub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSTHCA<- data.frame(cbind(L1HSTHCA, dupsout))
	print(L1HSTHCA)
}

rownames(L1HSTHCA)=L1HSTHCA[,1]
L1HSTHCA2<- L1HSTHCA[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
62, 64, 66, 68, 70, 72, 74, 76, 78, 80,
82, 84, 86, 88, 90, 92, 94, 96)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
write.table(L1HSTHCA2, "THCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon1<- read.table("105312870 exon2.txt")
exon52betas<- exon1[,c(1,5)]
patients<- stri_sub(exon1[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% THCAsub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSTHCA2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
write.table(finalexon52, "THCAexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
finalexon52<- read.table("THCAexon52.txt", header=TRUE)
legend_title= "THCA"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="THCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="THCA L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% THCAsub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSTHCA2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
write.table(finalexon1134, "THCAexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
finalexon1134<- read.table("THCAexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "THCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="THCA")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="THCA L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
exon795betas<- exon3[,c(1,5)]
patients<- stri_sub(exon3[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% THCAsub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSTHCA2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
write.table(finalexon795, "THCAexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
finalexon795<- read.table("THCAexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "THCA"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="THCA")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="THCA L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
exon355betas<- exon4[,c(1,5)]
patients<- stri_sub(exon4[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% THCAsub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSTHCA2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
write.table(exon355, "THCAexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
finalexon355<- read.table("THCAexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "THCA"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="THCA")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="THCA L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
exon1090betas<- exon5[,c(1,5)]
patients<- stri_sub(exon5[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% THCAsub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSTHCA2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
write.table(exon1090, "THCAexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
finalexon1090<- read.table("THCAexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "THCA"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="THCA")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="THCA L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
########ALL cancer types together
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/bedresults/")
BLCAnames<- read.table("BLCAbedpatients.txt")
setwd("/home/quints1/Patient_Names/THCA/edited_patients/bedresults/normal_patients")
THCAnames<- read.table("THCAbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/bedresults/")
BRCAnames<- read.table("BRCAbedpatients.txt")
setwd("/home/quints1/Patient_Names/COAD/edited_patients/bedresults/normal_patients")
COADnames<- read.table("COADbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/bedresults/normal_patients")
ESCAnames<- read.table("ESCAbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/bedresults/normal_patients")
HNSCnames<- read.table("HNSCbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/bedresults/normal_patients")
KIRCnames<- read.table("KIRCbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/bedresults/normal_patients")
KIRPnames<- read.table("KIRPbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/bedresults/normal_patients")
LIHCnames<- read.table("LIHCbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/bedresults/normal_patients")
LUADnames<- read.table("LUADbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/bedresults/normal_patients")
LUSCnames<- read.table("LUSCbedpatientsn.txt")
setwd("/home/quints1/Patient_Names/READ/edited_patients/bedresults/normal_patients")
READnames<- read.table("READpatientsn.txt")
Allpats<- rbind(BLCAnames, BRCAnames, COADnames, ESCAnames, HNSCnames, 
KIRCnames, KIRPnames, LIHCnames, LUADnames, LUSCnames, READnames, THCAnames)
library(stringi)
new<- as.matrix(stri_sub(Allpats[,1], 1, -28))
alls=NULL
for (i in new)
{
	y=NULL
	y<- data.frame(grep("*11A", i, value=TRUE))
	alls<- rbind(alls, y)
} 
#L1HS readin
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
allfiles<- as.matrix(dir(path="/DataDrives/dd7/pre-mrna/RawCnts"))

Instantfiles=NULL
for (i in allfiles[,1])
{
	y=NULL
	y<- data.frame(grep("*discount.instance*", i, value=TRUE))
	Instantfiles<- rbind(Instantfiles, y)
}

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 16))
combined<- cbind(newfiles, Instantfiles)
ssub<- as.matrix(combined[(combined[,1] %in% alls[,1]),])
setwd("/home/quints1/Patient_Names")
dupall<- read.table("dupsall.txt")
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSall=NULL
for (i in ssub[,2])
{
	name<- stri_sub(i, 0, 16)
	k<-read.table(i)
	outL1HS=NULL
	for (zz in k)
	{
		Y=NULL
		Y<- data.frame(grep("*L1HS*", zz, value=TRUE))
		outL1HS<- rbind(outL1HS, Y)
	}
	finalL1HS<- k[which(k[,1] %in% outL1HS[,1]),]
	dupsout<- finalL1HS[which(finalL1HS[,1] %in% dupall[,2]),]
	colnames(dupsout)=c("L1HS", name)
	dupsout<- as.matrix(dupsout)
	L1HSall<- data.frame(cbind(L1HSall, dupsout))
	print(L1HSall)
}

rownames(L1HSall)=L1HSall[,1]
L1HSall2<- L1HSall[,c(2, 4,	6,	8,	10,	12,	14,	16,	18, 
20,	22,	24,	26,	28,	30,	32,	34,	36, 
38,	40,	42,	44,	46,	48,	50,	52,	54, 
56,	58,	60,	62,	64,	66,	68,	70,	72, 
74,	76,	78,	80,	82,	84,	86,	88,	90, 
92,	94,	96,	98,	100,	102,	104,	106,	108, 
110,	112,	114,	116,	118,	120,	122,	124,	126, 
128,	130,	132,	134,	136,	138,	140,	142,	144, 
146,	148,	150,	152,	154,	156,	158,	160,	162 ,
164,	166,	168,	170,	172,	174,	176,	178,	180, 
182,	184,	186,	188,	190,	192,	194,	196,	198, 
200,	202,	204,	206, 	208,	210,	212,	214,	216, 
218,	220,	222,	224,	226,	228,	230,	232,	234, 
236,	238,	240,	242,	244,	246,	248,	250,	252, 
254,	256,	258,	260,	262,	264,	266,	268,	270, 
272,	274,	276,	278,	280,	282,	284,	286,	288, 
290,	292,	294,	296,	298,	300,	302,	304,	306, 
308,	310,	312,	314,	316,	318,	320,	322,	324, 
326,	328,	330,	332,	334,	336,	338,	340,	342, 
344,	346,	348,	350,	352,	354,	356,	358,	360, 
362,	364,	366,	368,	370,	372,	374,	376,	378, 
380,	382,	384,	386,	388,	390,	392,	394,	396, 
398,	400,	402,	404,	406,	408,	410,	412,	414, 
416,	418,	420,	422,	424,	426,	428,	430,	432, 
434,	436,	438,	440,	442,	444,	446,	448,	450, 
452,	454,	456,	458,	460,	462,	464,	466,	468, 
470,	472,	474,	476,	478,	480,	482,	484,	486, 
488,	490,	492,	494,	496,	498,	500,	502,	504, 
506,	508,	510,	512,	514,	516,	518,	520,	522, 
524,	526,	528,	530,	532,	534,	536,	538,	540, 
542,	544,	546,	548,	550,	552,	554,	556,	558, 
560,	562,	564,	566,	568,	570,	572,	574,	576, 
578,	580,	582, 584)]
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
write.table(L1HSall2, "ALLmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#exon 105312870
#read in exon 105312870
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon1<- read.table("105312870 exon.txt")
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon2<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon3<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon4<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon5<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon6<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon7<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon8<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon9<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon10<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon11<- read.table("105312870 exon2.txt")
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon12<- read.table("105312870 exon2.txt")
exonall<- rbind(exon1, exon2, exon3, exon4, exon5, exon6, exon7, exon8, exon9, exon10, exon11, exon12)
exon52betas<- exonall[,c(1,5)]
patients<- stri_sub(exonall[,1], 0, 16)
exon52betas2<- data.frame(cbind(patients, exon52betas[,2]))
outbetasexon52<- exon52betas2[which(exon52betas2[,1] %in% ssub[,1]),]
rownames(outbetasexon52)=outbetasexon52[,1]
L1HS52<- L1HSall2[c("L1HS_dup52:L1HS:L1:LINE:-"),]
tL1HS52<- data.frame(t(L1HS52))
library(stringr)
test<- str_replace_all(rownames(tL1HS52), "[.]", "-")
ttL1HS52<- cbind(test, tL1HS52)
rownames(ttL1HS52)=ttL1HS52[,1]
exon52<- merge(ttL1HS52, outbetasexon52, by="row.names", all.x=TRUE)
exon52a<- na.omit(exon52)
finalexon52<- exon52a[,c(1,3,5)]
colnames(finalexon52)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
write.table(finalexon52, "ALLexon52.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot
library("ggplot2")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
finalexon52<- read.table("ALLexon52.txt", header=TRUE)
legend_title= "ALL"
ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="ALL")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ALL L1HS_dup52")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon dup1134
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon1<- read.table("90699440 exon.txt")
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon2<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon3<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon4<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon5<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon6<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon7<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon8<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon9<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon10<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon11<- read.table("90699440 exon2.txt")
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon12<- read.table("90699440 exon2.txt")
exonall<- rbind(exon1, exon2, exon3, exon4, exon5, exon6, exon7, exon8, exon9, exon10, exon11, exon12)
exon52betas<- exonall[,c(1,5)]
patients<- stri_sub(exonall[,1], 0, 16)
exon1134betas<- exon2[,c(1,5)]
patients<- stri_sub(exon2[,1], 0, 16)
exon1134betas2<- data.frame(cbind(patients, exon1134betas[,2]))
outbetasexon1134<- exon1134betas2[which(exon1134betas2[,1] %in% ssub[,1]),]
outbetasexon1134<-na.omit(outbetasexon1134)
rownames(outbetasexon1134)=outbetasexon1134[,1]
L1HS1134<- L1HSall2[c("L1HS_dup1134:L1HS:L1:LINE:-"),]
tL1HS1134<- data.frame(t(L1HS1134))
library(stringr)
test<- str_replace_all(rownames(tL1HS1134), "[.]", "-")
ttL1HS1134<- cbind(test, tL1HS1134)
rownames(ttL1HS1134)=ttL1HS1134[,1]
exon1134<- merge(ttL1HS1134, outbetasexon1134, by="row.names", all.x=TRUE)
exon1134a<- na.omit(exon1134)
finalexon1134<- exon1134a[,c(1,3,5)]
colnames(finalexon1134)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
write.table(finalexon1134, "ALLexon1134.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
finalexon1134<- read.table("ALLexon1134.txt", header=TRUE)
library("ggplot2")
legend_title= "ALL"
ggplot(NULL) + 
  geom_point(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")], color="ALL")) + geom_smooth(aes(finalexon1134[,c("expression")], finalexon1134[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ALL L1HS_dup1134")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon795
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon1<- read.table("73787793 exon.txt")
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon2<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon3<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon4<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon5<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon6<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon7<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon8<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon9<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon10<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon11<- read.table("73787793 exon2.txt")
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon12<- read.table("73787793 exon2.txt")
exonall<- rbind(exon1, exon2, exon3, exon4, exon5, exon6, exon7, exon8, exon9, exon10, exon11, exon12)
exon795betas<- exonall[,c(1,5)]
patients<- stri_sub(exonall[,1], 0, 16)
exon795betas2<- data.frame(cbind(patients, exon795betas[,2]))
outbetasexon795<- exon795betas2[which(exon795betas2[,1] %in% ssub[,1]),]
outbetasexon795<-na.omit(outbetasexon795)
rownames(outbetasexon795)=outbetasexon795[,1]
L1HS795<- L1HSall2[c("L1HS_dup795:L1HS:L1:LINE:-"),]
tL1HS795<- data.frame(t(L1HS795))
library(stringr)
test<- str_replace_all(rownames(tL1HS795), "[.]", "-")
ttL1HS795<- cbind(test, tL1HS795)
rownames(ttL1HS795)=ttL1HS795[,1]
exon795<- merge(ttL1HS795, outbetasexon795, by="row.names", all.x=TRUE)
exon795a<- na.omit(exon795)
finalexon795<- exon795a[,c(1,3,5)]
colnames(finalexon795)= c("rownames", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
write.table(finalexon795, "ALLexon795.txt", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
finalexon795<- read.table("ALLexon795.txt", header=TRUE)
library("ggplot2")
legend_title= "ALL"
ggplot(NULL) + 
  geom_point(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")], color="ALL")) + geom_smooth(aes(finalexon795[,c("expression")], finalexon795[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ALL L1HS_dup795")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	
#read in exon355 (werid one)
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon1<- read.table("16943900 exon.txt")
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon2<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon3<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon4<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon5<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon6<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon7<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon8<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon9<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon10<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon11<- read.table("16943900 exon2.txt")
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon12<- read.table("16943900 exon2.txt")
exonall<- rbind(exon1, exon2, exon3, exon4, exon5, exon6, exon7, exon8, exon9, exon10, exon11, exon12)
exon355betas<- exonall[,c(1,5)]
patients<- stri_sub(exonall[,1], 0, 16)
exon355betas2<- data.frame(cbind(patients, exon355betas[,2]))
outbetasexon355<- exon355betas2[which(exon355betas2[,1] %in% ssub[,1]),]
outbetasexon355<-na.omit(outbetasexon355)
L1HS355<- L1HSall2[c("L1HS_dup355:L1HS:L1:LINE:+"),]
tL1HS355<- data.frame(t(L1HS355))
ttL1HS355<- cbind(rownames(tL1HS355), tL1HS355)
library(stringr)
test<- str_replace_all(rownames(tL1HS355), "[.]", "-")
ttL1HS355<- cbind(test, tL1HS355)
rownames(ttL1HS355)=ttL1HS355[,1]
exon355<- merge(ttL1HS355, outbetasexon355, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon355)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
write.table(exon355, "ALLexon355.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
finalexon355<- read.table("ALLexon355.txt", header=TRUE)
library("ggplot2")
legend_title= "ALL"
ggplot(NULL) + 
  geom_point(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")], color="ALL")) + geom_smooth(aes(finalexon355[,c("expression")], finalexon355[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ALL L1HS_dup355")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

#read in exon 1090 (same issue)
setwd("/home/quints1/Patient_Names/BLCA/edited_patients/byexon")
exon1<- read.table("13258876 exon.txt")
setwd("/home/quints1/Patient_Names/BRCA/edited_patients/byexon")
exon2<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/COAD/edited_patients/byexon")
exon3<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/ESCA/edited_patients/byexon")
exon4<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/HNSC/edited_patients/byexon")
exon5<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRC/edited_patients/byexon")
exon6<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/KIRP/edited_patients/byexon")
exon7<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/LIHC/edited_patients/byexon")
exon8<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/LUAD/edited_patients/byexon")
exon9<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/LUSC/edited_patients/byexon")
exon10<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/READ/edited_patients/byexon")
exon11<- read.table("13258876 exon2.txt")
setwd("/home/quints1/Patient_Names/THCA/edited_patients/byexon")
exon12<- read.table("13258876 exon2.txt")
exonall<- rbind(exon1, exon2, exon3, exon4, exon5, exon6, exon7, exon8, exon9, exon10, exon11, exon12)
exon1090betas<- exonall[,c(1,5)]
patients<- stri_sub(exonall[,1], 0, 16)
exon1090betas2<- data.frame(cbind(patients, exon1090betas[,2]))
outbetasexon1090<- exon1090betas2[which(exon1090betas2[,1] %in% ssub[,1]),]
outbetasexon1090<-na.omit(outbetasexon1090)
L1HS1090<- L1HSall2[c("L1HS_dup1090:L1HS:L1:LINE:+"),]
tL1HS1090<- data.frame(t(L1HS1090))
ttL1HS1090<- cbind(rownames(tL1HS1090), tL1HS1090)
library(stringr)
test<- str_replace_all(rownames(tL1HS1090), "[.]", "-")
ttL1HS1090<- cbind(test, tL1HS1090)
rownames(ttL1HS1090)=ttL1HS1090[,1]
exon1090<- merge(ttL1HS1090, outbetasexon1090, by.x="test", by.y="patients", all.x=TRUE)
colnames(exon1090)= c("patients", "expression", "betavalues")
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
write.table(exon1090, "ALLexon1090.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#plot in rstudio
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
finalexon1090<- read.table("ALLexon1090.txt", header=TRUE)
library("ggplot2")
legend_title= "ALL"
ggplot(NULL) + 
  geom_point(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")], color="ALL")) + geom_smooth(aes(finalexon1090[,c("expression")], finalexon1090[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title="ALL L1HS_dup1090")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))

##################Get P and R values
####BLCA
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
BLCA52<- read.table("BLCAexon52.txt")
BLCA795<- read.table("BLCAexon795.txt")
BLCA1134<- read.table("BLCAexon1134.txt")
#Function to calculate p-values. (lm(x ~ y))
BLCA52P<-anova(lm(BLCA52[,c("expression")] ~ BLCA52[,c("betavalues")]))$Pr[1]
BLCA795P<-anova(lm(BLCA795[,c("expression")] ~ BLCA795[,c("betavalues")]))$Pr[1]
BLCA1134P<-anova(lm(BLCA1134[,c("expression")] ~ BLCA1134[,c("betavalues")]))$Pr[1]
BLCAPvals<- rbind(BLCA52P, BLCA795P, BLCA1134P)
BLCA52R<- cor(BLCA52[,c("expression")], BLCA52[,c("betavalues")])
BLCA795R<- cor(BLCA795[,c("expression")], BLCA795[,c("betavalues")])
BLCA1134R<- cor(BLCA1134[,c("expression")], BLCA1134[,c("betavalues")])
BLCARvals<- rbind(BLCA52R, BLCA795R, BLCA1134R)
BLCAtotal<- cbind(BLCAPvals, BLCARvals)
colnames(BLCAtotal)= c("P-values", "R-values")
####BRCA
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
BRCA52<- read.table("BRCAexon52.txt", row.names=1, header=TRUE)
BRCA795<- read.table("BRCAexon795.txt", row.names=1, header=TRUE)
#Function to calculate p-values. (lm(x ~ y))
BRCA52P<-anova(lm(BRCA52[,c("expression")] ~ BRCA52[,c("betavalues")]))$Pr[1]
BRCA795P<-anova(lm(BRCA795[,c("expression")] ~ BRCA795[,c("betavalues")]))$Pr[1]
BRCAPvals<- rbind(BRCA52P, BRCA795P)
BRCA52R<- cor(BRCA52[,c("expression")], BRCA52[,c("betavalues")])
BRCA795R<- cor(BRCA795[,c("expression")], BRCA795[,c("betavalues")])
BRCARvals<- rbind(BRCA52R, BRCA795R)
BRCAtotal<- cbind(BRCAPvals, BRCARvals)
colnames(BRCAtotal)= c("P-values", "R-values")
###COAD
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
COAD795<- read.table("COADexon795.txt", header=TRUE, row.names=1)
COAD795P<-anova(lm(COAD795[,c("expression")] ~ COAD795[,c("betavalues")]))$Pr[1]
COAD795R<- cor(COAD795[,c("expression")], COAD795[,c("betavalues")])
COADtotal<- cbind(COAD795P, COAD795R)
colnames(COADtotal)= c("P-values", "R-values")
rownames(COADtotal)= c("COAD795")
###HNSC
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
HNSC52<- read.table("HNSCexon52.txt", header=TRUE, row.names=1)
HNSC795<- read.table("HNSCexon795.txt", header=TRUE, row.names=1)
HNSC1134<- read.table("HNSCexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
HNSC52P<-anova(lm(HNSC52[,c("expression")] ~ HNSC52[,c("betavalues")]))$Pr[1]
HNSC795P<-anova(lm(HNSC795[,c("expression")] ~ HNSC795[,c("betavalues")]))$Pr[1]
HNSC1134P<-anova(lm(HNSC1134[,c("expression")] ~ HNSC1134[,c("betavalues")]))$Pr[1]
HNSCPvals<- rbind(HNSC52P, HNSC795P, HNSC1134P)
HNSC52R<- cor(HNSC52[,c("expression")], HNSC52[,c("betavalues")])
HNSC795R<- cor(HNSC795[,c("expression")], HNSC795[,c("betavalues")])
HNSC1134R<- cor(HNSC1134[,c("expression")], HNSC1134[,c("betavalues")])
HNSCRvals<- rbind(HNSC52R, HNSC795R, HNSC1134R)
HNSCtotal<- cbind(HNSCPvals, HNSCRvals)
colnames(HNSCtotal)= c("P-values", "R-values")
###KIRC
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
KIRC795<- read.table("KIRCexon795.txt", header=TRUE, row.names=1)
KIRC1134<- read.table("KIRCexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
KIRC795P<-anova(lm(KIRC795[,c("expression")] ~ KIRC795[,c("betavalues")]))$Pr[1]
KIRC1134P<-anova(lm(KIRC1134[,c("expression")] ~ KIRC1134[,c("betavalues")]))$Pr[1]
KIRCPvals<- rbind(KIRC795P, KIRC1134P)
KIRC795R<- cor(KIRC795[,c("expression")], KIRC795[,c("betavalues")])
KIRC1134R<- cor(KIRC1134[,c("expression")], KIRC1134[,c("betavalues")])
KIRCRvals<- rbind(KIRC795R, KIRC1134R)
KIRCtotal<- cbind(KIRCPvals, KIRCRvals)
colnames(KIRCtotal)= c("P-values", "R-values")
###KIRP
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
KIRP1134<- read.table("KIRPexon1134.txt", header=TRUE, row.names=1)
KIRP1134P<-anova(lm(KIRP1134[,c("expression")] ~ KIRP1134[,c("betavalues")]))$Pr[1]
KIRP1134R<- cor(KIRP1134[,c("expression")], KIRP1134[,c("betavalues")])
KIRPtotal<- cbind(KIRP1134P, KIRP1134R)
colnames(KIRPtotal)= c("P-values", "R-values")
rownames(KIRPtotal)= c("KIRP1134")
###LUAD
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
LUAD795<- read.table("LUADexon795.txt", header=TRUE, row.names=1)
LUAD1134<- read.table("LUADexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
LUAD795P<-anova(lm(LUAD795[,c("expression")] ~ LUAD795[,c("betavalues")]))$Pr[1]
LUAD1134P<-anova(lm(LUAD1134[,c("expression")] ~ LUAD1134[,c("betavalues")]))$Pr[1]
LUADPvals<- rbind(LUAD795P, LUAD1134P)
LUAD795R<- cor(LUAD795[,c("expression")], LUAD795[,c("betavalues")])
LUAD1134R<- cor(LUAD1134[,c("expression")], LUAD1134[,c("betavalues")])
LUADRvals<- rbind(LUAD795R, LUAD1134R)
LUADtotal<- cbind(LUADPvals, LUADRvals)
colnames(LUADtotal)= c("P-values", "R-values")
###THCA
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
THCA1090<- read.table("THCAexon1090.txt", header=TRUE)
THCA1090P<-anova(lm(THCA1090[,c("expression")] ~ THCA1090[,c("betavalues")]))$Pr[1]
THCA1090R<- cor(THCA1090[,c("expression")], THCA1090[,c("betavalues")])
THCAtotal<- cbind(THCA1090P, THCA1090R)
colnames(THCAtotal)= c("P-values", "R-values")
rownames(THCAtotal)= c("THCA1090")
#ALL
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
ALL52<- read.table("ALLexon52.txt", header=TRUE, row.names=1)
ALL795<- read.table("ALLexon795.txt", header=TRUE, row.names=1)
ALL1134<- read.table("ALLexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
ALL52P<-anova(lm(ALL52[,c("expression")] ~ ALL52[,c("betavalues")]))$Pr[1]
ALL795P<-anova(lm(ALL795[,c("expression")] ~ ALL795[,c("betavalues")]))$Pr[1]
ALL1134P<-anova(lm(ALL1134[,c("expression")] ~ ALL1134[,c("betavalues")]))$Pr[1]
ALLPvals<- rbind(ALL52P, ALL795P, ALL1134P)
ALL52R<- cor(ALL52[,c("expression")], ALL52[,c("betavalues")])
ALL795R<- cor(ALL795[,c("expression")], ALL795[,c("betavalues")])
ALL1134R<- cor(ALL1134[,c("expression")], ALL1134[,c("betavalues")])
ALLRvals<- rbind(ALL52R, ALL795R, ALL1134R)
ALLtotal<- cbind(ALLPvals, ALLRvals)
colnames(ALLtotal)= c("P-values", "R-values")
###Combine them
combined<- rbind(ALLtotal, BLCAtotal, BRCAtotal, COADtotal, HNSCtotal, KIRCtotal, KIRPtotal, LUADtotal, THCAtotal)
sortedbypvalues<- combined[order(combined[,c("P-values")]),, drop=FALSE]
setwd("/home/quints1/Patient_Names/betavaluegraphs/")
write.table(sortedbypvalues, file="sorted by the negative slopes.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

###################################
############Do the same but get all the data for reference just in case
####BLCA
setwd("/home/quints1/Patient_Names/betavaluegraphs/BLCA")
BLCA52<- read.table("BLCAexon52.txt")
BLCA795<- read.table("BLCAexon795.txt")
BLCA1134<- read.table("BLCAexon1134.txt")
BLCA1090<-read.table("BLCAexon1090.txt")
#Function to calculate p-values. (lm(x ~ y))
BLCA52P<-anova(lm(BLCA52[,c("expression")] ~ BLCA52[,c("betavalues")]))$Pr[1]
BLCA795P<-anova(lm(BLCA795[,c("expression")] ~ BLCA795[,c("betavalues")]))$Pr[1]
BLCA1134P<-anova(lm(BLCA1134[,c("expression")] ~ BLCA1134[,c("betavalues")]))$Pr[1]
BLCA1090P<-anova(lm(BLCA1090[,c("expression")] ~ BLCA1090[,c("betavalues")]))$Pr[1]
BLCAPvals<- rbind(BLCA52P, BLCA795P, BLCA1134P, BLCA1090P)
BLCA52R<- cor(BLCA52[,c("expression")], BLCA52[,c("betavalues")])
BLCA795R<- cor(BLCA795[,c("expression")], BLCA795[,c("betavalues")])
BLCA1134R<- cor(BLCA1134[,c("expression")], BLCA1134[,c("betavalues")])
BLCA1090R<- cor(BLCA1090[,c("expression")], BLCA1090[,c("betavalues")])
BLCARvals<- rbind(BLCA52R, BLCA795R, BLCA1134R, BLCA1090R)
BLCAtotal<- cbind(BLCAPvals, BLCARvals)
colnames(BLCAtotal)= c("P-values", "R-values")
####BRCA
setwd("/home/quints1/Patient_Names/betavaluegraphs/BRCA")
BRCA52<- read.table("BRCAexon52.txt", row.names=1, header=TRUE)
BRCA795<- read.table("BRCAexon795.txt", row.names=1, header=TRUE)
BRCA1134<- read.table("BRCAexon1134.txt", row.names=1, header=TRUE)
BRCA1090<-read.table("BRCAexon1090.txt",header=TRUE)
#Function to calculate p-values. (lm(x ~ y))
BRCA52P<-anova(lm(BRCA52[,c("expression")] ~ BRCA52[,c("betavalues")]))$Pr[1]
BRCA795P<-anova(lm(BRCA795[,c("expression")] ~ BRCA795[,c("betavalues")]))$Pr[1]
BRCA1134P<-anova(lm(BRCA1134[,c("expression")] ~ BRCA1134[,c("betavalues")]))$Pr[1]
BRCA1090P<-anova(lm(BRCA1090[,c("expression")] ~ BRCA1090[,c("betavalues")]))$Pr[1]
BRCAPvals<- rbind(BRCA52P, BRCA795P, BRCA1134P, BRCA1090P)
BRCA52R<- cor(BRCA52[,c("expression")], BRCA52[,c("betavalues")])
BRCA795R<- cor(BRCA795[,c("expression")], BRCA795[,c("betavalues")])
BRCA1134R<- cor(BRCA1134[,c("expression")], BRCA1134[,c("betavalues")])
BRCA1090R<- cor(BRCA1090[,c("expression")], BRCA1090[,c("betavalues")])
BRCARvals<- rbind(BRCA52R, BRCA795R, BRCA1134R, BRCA1090R)
BRCAtotal<- cbind(BRCAPvals, BRCARvals)
colnames(BRCAtotal)= c("P-values", "R-values")
###COAD
setwd("/home/quints1/Patient_Names/betavaluegraphs/COAD")
COAD795<- read.table("COADexon795.txt", header=TRUE, row.names=1)
COAD795P<-anova(lm(COAD795[,c("expression")] ~ COAD795[,c("betavalues")]))$Pr[1]
COAD795R<- cor(COAD795[,c("expression")], COAD795[,c("betavalues")])
COADtotal<- cbind(COAD795P, COAD795R)
colnames(COADtotal)= c("P-values", "R-values")
rownames(COADtotal)= c("COAD795")
###HNSC
setwd("/home/quints1/Patient_Names/betavaluegraphs/HNSC")
HNSC52<- read.table("HNSCexon52.txt", header=TRUE, row.names=1)
HNSC795<- read.table("HNSCexon795.txt", header=TRUE, row.names=1)
HNSC1134<- read.table("HNSCexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
HNSC52P<-anova(lm(HNSC52[,c("expression")] ~ HNSC52[,c("betavalues")]))$Pr[1]
HNSC795P<-anova(lm(HNSC795[,c("expression")] ~ HNSC795[,c("betavalues")]))$Pr[1]
HNSC1134P<-anova(lm(HNSC1134[,c("expression")] ~ HNSC1134[,c("betavalues")]))$Pr[1]
HNSCPvals<- rbind(HNSC52P, HNSC795P, HNSC1134P)
HNSC52R<- cor(HNSC52[,c("expression")], HNSC52[,c("betavalues")])
HNSC795R<- cor(HNSC795[,c("expression")], HNSC795[,c("betavalues")])
HNSC1134R<- cor(HNSC1134[,c("expression")], HNSC1134[,c("betavalues")])
HNSCRvals<- rbind(HNSC52R, HNSC795R, HNSC1134R)
HNSCtotal<- cbind(HNSCPvals, HNSCRvals)
colnames(HNSCtotal)= c("P-values", "R-values")
###KIRC
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRC")
KIRC795<- read.table("KIRCexon795.txt", header=TRUE, row.names=1)
KIRC1134<- read.table("KIRCexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
KIRC795P<-anova(lm(KIRC795[,c("expression")] ~ KIRC795[,c("betavalues")]))$Pr[1]
KIRC1134P<-anova(lm(KIRC1134[,c("expression")] ~ KIRC1134[,c("betavalues")]))$Pr[1]
KIRCPvals<- rbind(KIRC795P, KIRC1134P)
KIRC795R<- cor(KIRC795[,c("expression")], KIRC795[,c("betavalues")])
KIRC1134R<- cor(KIRC1134[,c("expression")], KIRC1134[,c("betavalues")])
KIRCRvals<- rbind(KIRC795R, KIRC1134R)
KIRCtotal<- cbind(KIRCPvals, KIRCRvals)
colnames(KIRCtotal)= c("P-values", "R-values")
###KIRP
setwd("/home/quints1/Patient_Names/betavaluegraphs/KIRP")
KIRP1134<- read.table("KIRPexon1134.txt", header=TRUE, row.names=1)
KIRP1134P<-anova(lm(KIRP1134[,c("expression")] ~ KIRP1134[,c("betavalues")]))$Pr[1]
KIRP1134R<- cor(KIRP1134[,c("expression")], KIRP1134[,c("betavalues")])
KIRPtotal<- cbind(KIRP1134P, KIRP1134R)
colnames(KIRPtotal)= c("P-values", "R-values")
rownames(KIRPtotal)= c("KIRP1134")
###LUAD
setwd("/home/quints1/Patient_Names/betavaluegraphs/LUAD")
LUAD795<- read.table("LUADexon795.txt", header=TRUE, row.names=1)
LUAD1134<- read.table("LUADexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
LUAD795P<-anova(lm(LUAD795[,c("expression")] ~ LUAD795[,c("betavalues")]))$Pr[1]
LUAD1134P<-anova(lm(LUAD1134[,c("expression")] ~ LUAD1134[,c("betavalues")]))$Pr[1]
LUADPvals<- rbind(LUAD795P, LUAD1134P)
LUAD795R<- cor(LUAD795[,c("expression")], LUAD795[,c("betavalues")])
LUAD1134R<- cor(LUAD1134[,c("expression")], LUAD1134[,c("betavalues")])
LUADRvals<- rbind(LUAD795R, LUAD1134R)
LUADtotal<- cbind(LUADPvals, LUADRvals)
colnames(LUADtotal)= c("P-values", "R-values")
###THCA
setwd("/home/quints1/Patient_Names/betavaluegraphs/THCA")
THCA1090<- read.table("THCAexon1090.txt", header=TRUE)
THCA1090P<-anova(lm(THCA1090[,c("expression")] ~ THCA1090[,c("betavalues")]))$Pr[1]
THCA1090R<- cor(THCA1090[,c("expression")], THCA1090[,c("betavalues")])
THCAtotal<- cbind(THCA1090P, THCA1090R)
colnames(THCAtotal)= c("P-values", "R-values")
rownames(THCAtotal)= c("THCA1090")
#ALL
setwd("/home/quints1/Patient_Names/betavaluegraphs/ALL")
ALL52<- read.table("ALLexon52.txt", header=TRUE, row.names=1)
ALL795<- read.table("ALLexon795.txt", header=TRUE, row.names=1)
ALL1134<- read.table("ALLexon1134.txt", header=TRUE, row.names=1)
#Function to calculate p-values. (lm(x ~ y))
ALL52P<-anova(lm(ALL52[,c("expression")] ~ ALL52[,c("betavalues")]))$Pr[1]
ALL795P<-anova(lm(ALL795[,c("expression")] ~ ALL795[,c("betavalues")]))$Pr[1]
ALL1134P<-anova(lm(ALL1134[,c("expression")] ~ ALL1134[,c("betavalues")]))$Pr[1]
ALLPvals<- rbind(ALL52P, ALL795P, ALL1134P)
ALL52R<- cor(ALL52[,c("expression")], ALL52[,c("betavalues")])
ALL795R<- cor(ALL795[,c("expression")], ALL795[,c("betavalues")])
ALL1134R<- cor(ALL1134[,c("expression")], ALL1134[,c("betavalues")])
ALLRvals<- rbind(ALL52R, ALL795R, ALL1134R)
ALLtotal<- cbind(ALLPvals, ALLRvals)
colnames(ALLtotal)= c("P-values", "R-values")







