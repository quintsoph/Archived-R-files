#####hg19 with promoters will find exons
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/BRCAfiles/editedfiles/bedresults")
BRCAexons<- read.table("BRCAexons.txt", sep="", fill=TRUE)
BRCAexons2<- data.frame(BRCAexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% BRCAexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

###############DONT DO ONLY FOR RECORDS READ IN List_of_alldups.txt instead
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
comparison<- read.table("TCGA-ZU-A8S4-11A-60e3b4d8-081a-4247-8144-3cce35345dbf.discount.instance.cntTable")
L1HSlist=NULL
for (i in comparison[,1])
{
	k=NULL
	k<- data.frame(grep("^L1HS", i, value=TRUE))
	L1HSlist<- rbind(L1HSlist, k)
}
#####################

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/BRCAfiles/editedfiles/bedresults")
BRCAnames<- read.table("BRCAbedpatients2.txt")
newBRCA<- as.matrix(stri_sub(BRCAnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
BRCAsub<- as.matrix(combined[(combined[,1] %in% newBRCA[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSBRCA<- data.frame(cbind(L1HSBRCA, out))
	print(L1HSBRCA)
}

rownames(L1HSBRCA)=L1HSBRCA[,1]
L1HSBRCA2<- L1HSBRCA[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
62, 64, 66, 68, 70, 72, 74, 76, 78, 80,
82, 84, 86, 88, 90, 92, 94, 96, 98, 100,
102, 104, 106, 108, 110, 112, 114, 116, 118, 120,
122, 124, 126, 128, 130, 132, 134, 136, 138, 140,
142, 144, 146, 148, 150, 152, 154, 156, 158, 160,
162, 164, 166, 168)]
setwd("/home/quints1/Patient_Names/BRCAfiles/betaresults")
write.table(L1HSBRCA2, "BRCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/BRCAfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/BRCAfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/BRCAfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% BRCAsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSBRCA2[which(rownames(L1HSBRCA2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/BRCAfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/BRCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/BRCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/BRCAfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "BRCA"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="BRCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/BRCAfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/BRCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/BRCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
BRCAprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("BRCA", namer, sep=""))
	BRCAprvalues<- rbind(BRCAprvalues, exontotal)
}
sortedbypvalues<- BRCAprvalues[order(BRCAprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesBRCA.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

#############BLCA
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/BLCAfiles/editedfiles/bedresults")
BRCAexons<- read.table("BLCAexons.txt", sep="", fill=TRUE)
BRCAexons2<- data.frame(BRCAexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% BRCAexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/BLCAfiles/editedfiles/bedresults")
BLCAnames<- read.table("BLCAbedpatients.txt")
newBLCA<- as.matrix(stri_sub(BLCAnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
BLCAsub<- as.matrix(combined[(combined[,1] %in% newBLCA[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSBLCA<- data.frame(cbind(L1HSBLCA, out))
	print(L1HSBLCA)
}

rownames(L1HSBLCA)=L1HSBLCA[,1]
L1HSBLCA2<- L1HSBLCA[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34)]
setwd("/home/quints1/Patient_Names/BLCAfiles/betaresults")
write.table(L1HSBLCA2, "BLCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/BLCAfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/BLCAfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/BLCAfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% BLCAsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSBLCA2[which(rownames(L1HSBLCA2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/BLCAfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/BLCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/BLCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/BLCAfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "BLCA"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="BLCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/BLCAfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/BLCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/BLCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
BLCAprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("BLCA", namer, sep=""))
	BLCAprvalues<- rbind(BLCAprvalues, exontotal)
}
sortedbypvalues<- BLCAprvalues[order(BLCAprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesBLCA.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)


##############COAD
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/COADfiles/editedfiles/bedresults")
COADexons<- read.table("COADexons.txt", sep="", fill=TRUE)
COADexons2<- data.frame(COADexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% COADexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/COADfiles/editedfiles/bedresults")
COADnames<- read.table("COADbedpatients.txt")
newCOAD<- as.matrix(stri_sub(COADnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
COADsub<- as.matrix(combined[(combined[,1] %in% newCOAD[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSCOAD<- data.frame(cbind(L1HSCOAD, out))
	print(L1HSCOAD)
}

rownames(L1HSCOAD)=L1HSCOAD[,1]
L1HSCOAD2<- L1HSCOAD[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30)]
setwd("/home/quints1/Patient_Names/COADfiles/betaresults")
write.table(L1HSCOAD2, "COADmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/COADfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/COADfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/COADfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% COADsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSCOAD2[which(rownames(L1HSCOAD2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/COADfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/COADfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/COADfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/COADfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "COAD"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="COAD")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/COADfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/COADfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/COADfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
COADprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("COAD", namer, sep=""))
	COADprvalues<- rbind(COADprvalues, exontotal)
}
sortedbypvalues<- COADprvalues[order(COADprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesCOAD.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########ESCA
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/ESCAfiles/editedfiles/bedresults")
ESCAexons<- read.table("ESCAexons.txt", sep="", fill=TRUE)
ESCAexons2<- data.frame(ESCAexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% ESCAexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/ESCAfiles/editedfiles/bedresults")
ESCAnames<- read.table("ESCAbedpatients.txt")
newESCA<- as.matrix(stri_sub(ESCAnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
ESCAsub<- as.matrix(combined[(combined[,1] %in% newESCA[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSESCA<- data.frame(cbind(L1HSESCA, out))
	print(L1HSESCA)
}

rownames(L1HSESCA)=L1HSESCA[,1]
L1HSESCA2<- L1HSESCA[,c(2, 4, 6, 8, 10, 12, 14, 16, 18)]
setwd("/home/quints1/Patient_Names/ESCAfiles/betaresults")
write.table(L1HSESCA2, "ESCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/ESCAfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/ESCAfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/ESCAfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% ESCAsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSESCA2[which(rownames(L1HSESCA2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/ESCAfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/ESCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/ESCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/ESCAfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "ESCA"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="ESCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/ESCAfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/ESCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/ESCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
ESCAprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("ESCA", namer, sep=""))
	ESCAprvalues<- rbind(ESCAprvalues, exontotal)
}
sortedbypvalues<- ESCAprvalues[order(ESCAprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesESCA.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########HNSC
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/HNSCfiles/editedfiles/bedresults")
HNSCexons<- read.table("HNSCexons.txt", sep="", fill=TRUE)
HNSCexons2<- data.frame(HNSCexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% HNSCexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/HNSCfiles/editedfiles/bedresults")
HNSCnames<- read.table("HNSCbedpatients.txt")
newHNSC<- as.matrix(stri_sub(HNSCnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
HNSCsub<- as.matrix(combined[(combined[,1] %in% newHNSC[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSHNSC<- data.frame(cbind(L1HSHNSC, out))
	print(L1HSHNSC)
}

rownames(L1HSHNSC)=L1HSHNSC[,1]
L1HSHNSC2<- L1HSHNSC[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38)]
setwd("/home/quints1/Patient_Names/HNSCfiles/betaresults")
write.table(L1HSHNSC2, "HNSCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/HNSCfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/HNSCfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/HNSCfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% HNSCsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSHNSC2[which(rownames(L1HSHNSC2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/HNSCfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/HNSCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/HNSCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/HNSCfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "HNSC"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="HNSC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/HNSCfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/HNSCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/HNSCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
HNSCprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("HNSC", namer, sep=""))
	HNSCprvalues<- rbind(HNSCprvalues, exontotal)
}
sortedbypvalues<- HNSCprvalues[order(HNSCprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesHNSC.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########KIRC   
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/KIRCfiles/editedfiles/bedresults")
KIRCexons<- read.table("KIRCexons.txt", sep="", fill=TRUE)
KIRCexons2<- data.frame(KIRCexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% KIRCexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/KIRCfiles/editedfiles/bedresults")
KIRCnames<- read.table("KIRCbedpatients.txt")
newKIRC<- as.matrix(stri_sub(KIRCnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
KIRCsub<- as.matrix(combined[(combined[,1] %in% newKIRC[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSKIRC<- data.frame(cbind(L1HSKIRC, out))
	print(L1HSKIRC)
}

rownames(L1HSKIRC)=L1HSKIRC[,1]
L1HSKIRC2<- L1HSKIRC[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48)]
setwd("/home/quints1/Patient_Names/KIRCfiles/betaresults")
write.table(L1HSKIRC2, "KIRCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/KIRCfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/KIRCfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/KIRCfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% KIRCsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSKIRC2[which(rownames(L1HSKIRC2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/KIRCfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/KIRCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/KIRCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/KIRCfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "KIRC"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="KIRC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/KIRCfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/KIRCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/KIRCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
KIRCprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("KIRC", namer, sep=""))
	KIRCprvalues<- rbind(KIRCprvalues, exontotal)
}
sortedbypvalues<- KIRCprvalues[order(KIRCprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesKIRC.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########KIRP   
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/KIRPfiles/editedfiles/bedresults")
KIRPexons<- read.table("KIRPexons.txt", sep="", fill=TRUE)
KIRPexons2<- data.frame(KIRPexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% KIRPexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/KIRPfiles/editedfiles/bedresults")
KIRPnames<- read.table("KIRPbedpatients.txt")
newKIRP<- as.matrix(stri_sub(KIRPnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
KIRPsub<- as.matrix(combined[(combined[,1] %in% newKIRP[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSKIRP<- data.frame(cbind(L1HSKIRP, out))
	print(L1HSKIRP)
}

rownames(L1HSKIRP)=L1HSKIRP[,1]
L1HSKIRP2<- L1HSKIRP[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42)]
setwd("/home/quints1/Patient_Names/KIRPfiles/betaresults")
write.table(L1HSKIRP2, "KIRPmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/KIRPfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/KIRPfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/KIRPfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% KIRPsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSKIRP2[which(rownames(L1HSKIRP2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/KIRPfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/KIRPfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/KIRPfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/KIRPfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "KIRP"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="KIRP")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/KIRPfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/KIRPfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/KIRPfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
KIRPprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("KIRP", namer, sep=""))
	KIRPprvalues<- rbind(KIRPprvalues, exontotal)
}
sortedbypvalues<- KIRPprvalues[order(KIRPprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesKIRP.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########LIHC  
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/LIHCfiles/editedfiles/bedresults")
LIHCexons<- read.table("LIHCexons.txt", sep="", fill=TRUE)
LIHCexons2<- data.frame(LIHCexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% LIHCexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/LIHCfiles/editedfiles/bedresults")
LIHCnames<- read.table("LIHCbedpatients.txt")
newLIHC<- as.matrix(stri_sub(LIHCnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
LIHCsub<- as.matrix(combined[(combined[,1] %in% newLIHC[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSLIHC<- data.frame(cbind(L1HSLIHC, out))
	print(L1HSLIHC)
}

rownames(L1HSLIHC)=L1HSLIHC[,1]
L1HSLIHC2<- L1HSLIHC[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
62, 64, 66, 68, 70, 72, 74, 76, 78, 80,
82)]
setwd("/home/quints1/Patient_Names/LIHCfiles/betaresults")
write.table(L1HSLIHC2, "LIHCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/LIHCfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LIHCfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/LIHCfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% LIHCsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSLIHC2[which(rownames(L1HSLIHC2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/LIHCfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/LIHCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LIHCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/LIHCfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "LIHC"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="LIHC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/LIHCfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/LIHCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LIHCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
LIHCprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("LIHC", namer, sep=""))
	LIHCprvalues<- rbind(LIHCprvalues, exontotal)
}
sortedbypvalues<- LIHCprvalues[order(LIHCprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesLIHC.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########LUAD
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/LUADfiles/editedfiles/bedresults")
LUADexons<- read.table("LUADexons.txt", sep="", fill=TRUE)
LUADexons2<- data.frame(LUADexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% LUADexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/LUADfiles/editedfiles/bedresults")
LUADnames<- read.table("LUADbedpatients.txt")
newLUAD<- as.matrix(stri_sub(LUADnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
LUADsub<- as.matrix(combined[(combined[,1] %in% newLUAD[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSLUAD<- data.frame(cbind(L1HSLUAD, out))
	print(L1HSLUAD)
}

rownames(L1HSLUAD)=L1HSLUAD[,1]
L1HSLUAD2<- L1HSLUAD[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40)]
setwd("/home/quints1/Patient_Names/LUADfiles/betaresults")
write.table(L1HSLUAD2, "LUADmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/LUADfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LUADfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/LUADfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% LUADsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSLUAD2[which(rownames(L1HSLUAD2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/LUADfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/LUADfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LUADfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/LUADfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "LUAD"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="LUAD")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/LUADfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/LUADfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LUADfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
LUADprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("LUAD", namer, sep=""))
	LUADprvalues<- rbind(LUADprvalues, exontotal)
}
sortedbypvalues<- LUADprvalues[order(LUADprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesLUAD.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########LUSC
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/LUSCfiles/editedfiles/bedresults")
LUSCexons<- read.table("LUSCexons.txt", sep="", fill=TRUE)
LUSCexons2<- data.frame(LUSCexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% LUSCexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/LUSCfiles/editedfiles/bedresults")
LUSCnames<- read.table("LUSCbedpatients.txt")
newLUSC<- as.matrix(stri_sub(LUSCnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
LUSCsub<- as.matrix(combined[(combined[,1] %in% newLUSC[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSLUSC<- data.frame(cbind(L1HSLUSC, out))
	print(L1HSLUSC)
}

rownames(L1HSLUSC)=L1HSLUSC[,1]
L1HSLUSC2<- L1HSLUSC[,c(2, 4, 6, 8, 10, 12, 14, 16)]
setwd("/home/quints1/Patient_Names/LUSCfiles/betaresults")
write.table(L1HSLUSC2, "LUSCmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/LUSCfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LUSCfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/LUSCfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% LUSCsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSLUSC2[which(rownames(L1HSLUSC2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/LUSCfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/LUSCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LUSCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/LUSCfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "LUSC"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="LUSC")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/LUSCfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/LUSCfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/LUSCfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
LUSCprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("LUSC", namer, sep=""))
	LUSCprvalues<- rbind(LUSCprvalues, exontotal)
}
sortedbypvalues<- LUSCprvalues[order(LUSCprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesLUSC.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########READ
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/READfiles/editedfiles/bedresults")
READexons<- read.table("READexons.txt", sep="", fill=TRUE)
READexons2<- data.frame(READexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% READexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/READfiles/editedfiles/bedresults")
READnames<- read.table("READbedpatients.txt")
newREAD<- as.matrix(stri_sub(READnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
READsub<- as.matrix(combined[(combined[,1] %in% newREAD[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSREAD<- data.frame(cbind(L1HSREAD, out))
	print(L1HSREAD)
}

rownames(L1HSREAD)=L1HSREAD[,1]
L1HSREAD2<- L1HSREAD[,c(2, 4)]
setwd("/home/quints1/Patient_Names/READfiles/betaresults")
write.table(L1HSREAD2, "READmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/READfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/READfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/READfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% READsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSREAD2[which(rownames(L1HSREAD2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/READfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/READfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/READfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/READfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "READ"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="READ")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/READfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/READfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/READfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
READprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("READ", namer, sep=""))
	READprvalues<- rbind(READprvalues, exontotal)
}
sortedbypvalues<- READprvalues[order(READprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesREAD.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

###########THCA
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/THCAfiles/editedfiles/bedresults")
THCAexons<- read.table("THCAexons.txt", sep="", fill=TRUE)
THCAexons2<- data.frame(THCAexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% THCAexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/THCAfiles/editedfiles/bedresults")
THCAnames<- read.table("THCAbedpatients.txt")
newTHCA<- as.matrix(stri_sub(THCAnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
THCAsub<- as.matrix(combined[(combined[,1] %in% newTHCA[,1]),])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSTHCA<- data.frame(cbind(L1HSTHCA, out))
	print(L1HSTHCA)
}

rownames(L1HSTHCA)=L1HSTHCA[,1]
L1HSTHCA2<- L1HSTHCA[,c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
62, 64, 66, 68, 70, 72, 74, 76, 78, 80,
82, 84, 86, 88, 90, 92, 94, 96, 98, 100)]
setwd("/home/quints1/Patient_Names/THCAfiles/betaresults")
write.table(L1HSTHCA2, "THCAmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/THCAfiles/editedfiles/bedresults/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/THCAfiles/editedfiles/bedresults/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/THCAfiles/editedfiles/bedresults/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% THCAsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSTHCA2[which(rownames(L1HSTHCA2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/THCAfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/THCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/THCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/THCAfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "THCA"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="THCA")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/THCAfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/THCAfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/THCAfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
THCAprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("THCA", namer, sep=""))
	THCAprvalues<- rbind(THCAprvalues, exontotal)
}
sortedbypvalues<- THCAprvalues[order(THCAprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesTHCA.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)


###########ALL
library(stringi)
setwd("/home/quints1/Patient_Names/")
L1HShg<- read.table("hg19_edit_L1HS_promoters.gtf", header=FALSE, skip=1)
L1HShgcut<- data.frame(L1HShg[,c(4,5,12)])
setwd("/home/quints1/Patient_Names/ALLfiles/")
ALLexons<- read.table("ALLexons.txt", sep="", fill=TRUE)
ALLexons2<- data.frame(ALLexons[-65,])
exoncut<- L1HShgcut[which(L1HShgcut[,1] %in% ALLexons2[,1]),]
dupsout<- data.frame(exoncut[,c(1,3)])

setwd("/home/quints1/Patient_Names")
dupslist<-read.table("List_of_alldups.txt", header=FALSE)
dupcut<- stri_sub(dupslist[,1], 0, -16)
dupcut<- cbind(dupcut, dupslist)
newdups<- data.frame(dupslist[which(dupcut[,1] %in% dupsout[,2]),])
setwd("/home/quints1/Patient_Names/ALLfiles/")
ALLnames<- read.table("ALLbedpatients.txt")
newALL<- as.matrix(stri_sub(ALLnames[,1], 1, -27))
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

newfiles<- as.matrix(stri_sub(Instantfiles[,1], 0, 12))
combined<- cbind(newfiles, Instantfiles)
ALLsub<- as.matrix(combined[(combined[,1] %in% newALL[,1]),])
setwd("/DataDrives/dd7/pre-mrna/RawCnts")
L1HSALL=NULL
for (i in ALLsub[,2])
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
	out<- finalL1HS[which(finalL1HS[,1] %in% newdups[,1]),]
	colnames(out)=c("L1HS", name)
	out<- as.matrix(out)
	L1HSALL<- data.frame(cbind(L1HSALL, out))
	print(L1HSALL)
}

rownames(L1HSALL)=L1HSALL[,1]
L1HSALL2<- L1HSALL[,c(seq(2,620,2))]
setwd("/home/quints1/Patient_Names/ALLfiles/betaresults")
write.table(L1HSALL2, "ALLmethyl.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
###GRAPHS
setwd("/home/quints1/Patient_Names/ALLfiles/byexon")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/ALLfiles/byexon"))
rownames(dupsout)=dupsout[,2]
rownames(dupcut)= dupcut[,1]
totalmerge<- merge(dupsout, dupcut, by="row.names")
totalmerge<- totalmerge[,c(2,3,5)]
totalmerge2<- totalmerge[,c(1,3)]
x <- " exon.txt"
stri_sub(totalmerge2[,1], 13, -13)<- x
rownames(totalmerge2)=totalmerge2[,1]

finalexons=NULL
for (i in allfiles[,1])
{
	finalexons=NULL
	setwd("/home/quints1/Patient_Names/ALLfiles/byexon")
	exon1<- read.table(i, header=FALSE)
	exonbetas<- exon1[,c(1,5)]
	patients<- stri_sub(exon1[,1], 0, 12)
	exonbetas2<- data.frame(cbind(patients, exonbetas[,2]))
	outbetasexon<- data.frame(exonbetas2[which(exonbetas2[,1] %in% ALLsub[,1]),])
	L1HSdup<- totalmerge2[c(i),]
	filedup<- L1HSdup[,2]
	library(stringr)
	L1HSexoncut<- L1HSALL2[which(rownames(L1HSALL2) %in% filedup),]
	tL1HSexoncut<- data.frame(t(L1HSexoncut))
	test<- str_replace_all(rownames(tL1HSexoncut), "[.]", "-")
	test<- stri_sub(test, 0,12)
	ttL1HSexoncut<- data.frame(cbind(test, tL1HSexoncut))
	colnames(ttL1HSexoncut)=c("patients", "L1HS")
	exonout<- merge(ttL1HSexoncut, outbetasexon, by="patients", all.x=TRUE)
	exonouta<- na.omit(exonout)
	finalexons<- exonouta
	colnames(finalexons)= c("patients", "expression", "betavalues")
	name <- stri_sub(filedup, 0, -16)
	namer<- paste0(name, "file.txt", collapse="_")
	setwd("/home/quints1/Patient_Names/ALLfiles/betaresults")
	write.table(finalexons, file=namer, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#####for loop for graphs!!!
library("ggplot2")
setwd("/home/quints1/Patient_Names/ALLfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/ALLfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])

for (i in allfiles2[,1])
{
setwd("/home/quints1/Patient_Names/ALLfiles/betaresults")
library(stringi)
namer<- stri_sub(i, 0, -9)
finalexon52<- read.table(i, header=TRUE)
legend_title= "ALL"
p= ggplot(NULL) + 
  geom_point(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")], color="ALL")) + geom_smooth(aes(finalexon52[,c("expression")], finalexon52[,c("betavalues")]), method=lm, se=FALSE, color="blue", size=0.5) +
   labs(x = "expression levels", y = "beta values", title=namer)+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("black"))
	setwd("/home/quints1/Patient_Names/ALLfiles/betaresults/graphs")
	ggsave(p, filename=paste("graph", namer, ".png", sep=""))
}

#####overall p and r values
setwd("/home/quints1/Patient_Names/ALLfiles/betaresults")
allfiles<- as.matrix(dir(path="/home/quints1/Patient_Names/ALLfiles/betaresults"))
allfiles2<- as.matrix(allfiles[-c(1,2),])
ALLprvalues=NULL
for (i in allfiles2[,1])
{
	exontotal=NULL
	namer<- stri_sub(i, 0, -9)
	exonfile<- read.table(i, header=TRUE)
	exonfileP<-anova(lm(exonfile[,c("expression")] ~ exonfile[,c("betavalues")]))$Pr[1]
	exonfileR<- cor(exonfile[,c("expression")], exonfile[,c("betavalues")])
	exontotal<- cbind(exonfileP, exonfileR)
	colnames(exontotal)= c("P-values", "R-values")
	rownames(exontotal)= c(paste("ALL", namer, sep=""))
	ALLprvalues<- rbind(ALLprvalues, exontotal)
}
sortedbypvalues<- ALLprvalues[order(ALLprvalues[,c("P-values")]),, drop=FALSE]
write.table(sortedbypvalues, "PandRvaluesALL.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)