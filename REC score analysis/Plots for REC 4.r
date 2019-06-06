#######Graphs for REC Score 4
#edit the L1HS file
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4")
L1HS<- read.table("L1HS.VST.normalized.adjusted.txt", header=TRUE)
normal<- subset(L1HS, L1HS[,3]=="normal")
Finalized<- normal[,c(1, 2, 8)]
write.table(Finalized, "L1HS normal for Graphs.txt", quote=FALSE, row.names=FALSE, sep="\t")

####Now start graph
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4")
L1HS<- read.table("L1HS normal for Graphs.txt", header=TRUE)

#Negative Regular top 10
#CCT7
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
CCT<- read.table("CCT7.txt", header=TRUE)
CCT<- subset(CCT, CCT[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- CCT[which(CCT[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- CCT[which(CCT[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- CCT[which(CCT[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "CCT7 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and CCT7 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3"))
	
###COPS7A
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
COP<- read.table("COPS7A.txt", header=TRUE)
COP<- subset(COP, COP[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- COP[which(COP[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- COP[which(COP[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "COPS7A (rpm)", y = "L1HS Expression (rpm)", title="L1HS and COPS7A Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "lightsalmon3"))

###FAM195B
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
FAM<- read.table("FAM195B.txt", header=TRUE)
FAM<- subset(FAM, FAM[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- FAM[which(FAM[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- FAM[which(FAM[,1] %in% HNSC[,1]),]
HNSCwL1HS<- cbind(HNSCcut, HNSC)
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 6)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- FAM[which(FAM[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "FAM195B (rpm)", y = "L1HS Expression (rpm)", title="L1HS and FAM195B Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "lightsalmon3"))

###FLOT1
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
FLOT<- read.table("FLOT1.txt", header=TRUE)
FLOT<- subset(FLOT, FLOT[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- FLOT[which(FLOT[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- FLOT[which(FLOT[,1] %in% HNSC[,1]),]
HNSCwL1HS<- cbind(HNSCcut, HNSC)
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 6)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	labs(x = "FLOT1 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and FLOT1 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue"))
	
###MAF1
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
MAF<- read.table("MAF1.txt", header=TRUE)
MAF<- subset(MAF, MAF[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- MAF[which(MAF[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- MAF[which(MAF[,1] %in% HNSC[,1]),]
HNSCwL1HS<- cbind(HNSCcut, HNSC)
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 6)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- MAF[which(MAF[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "MAF1 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and MAF1 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "lightsalmon3"))

###PARK7
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
PAR<- read.table("PARK7.txt", header=TRUE)
PAR<- subset(PAR, PAR[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- PAR[which(PAR[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- PAR[which(PAR[,1] %in% HNSC[,1]),]
HNSCwL1HS<- cbind(HNSCcut, HNSC)
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 6)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- PAR[which(PAR[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "PARK7 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and PARK7 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "lightsalmon3"))

###PSMC3
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
PSM<- read.table("PSMC3.txt", header=TRUE)
PSM<- subset(PSM, PSM[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- PSM[which(PSM[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- PSM[which(PSM[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "PSMC3 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and PSMC3 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "lightsalmon3"))

###RAB1B
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
RAB<- read.table("RAB1B.txt", header=TRUE)
RAB<- subset(RAB, RAB[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- RAB[which(RAB[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- RAB[which(RAB[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- RAB[which(RAB[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "RAB1B (rpm)", y = "L1HS Expression (rpm)", title="L1HS and RAB1B Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "lightsalmon3", "brown"))

###SNX17
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Negative Genes")
SNX<- read.table("SNX17.txt", header=TRUE)
SNX<- subset(SNX, SNX[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- SNX[which(SNX[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- SNX[which(SNX[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- SNX[which(SNX[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- SNX[which(SNX[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "SNX17 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and SNX17 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3", "brown"))

###Positive Regular top 10
##ERVL-B4
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Positive Genes")
ERVLB4<- read.table("ERVL-B4-int%3AERVL%3ALTR.txt", header=TRUE)
ERVLB4<- subset(ERVLB4, ERVLB4[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- ERVLB4[which(ERVLB4[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#COAD
COAD<- subset(L1HS, L1HS[,2]=="COAD")
COADcut<- ERVLB4[which(ERVLB4[,1] %in% COAD[,1]),]
COADwL1HS<- cbind(COADcut, COAD)
COADaL1HS<- COADwL1HS[,c(1, 3, 6)]
colnames(COADaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- ERVLB4[which(ERVLB4[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- ERVLB4[which(ERVLB4[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- ERVLB4[which(ERVLB4[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "ERVL-B4-int:3AERVL:3ALTR (rpm)", y = "L1HS Expression (rpm)", title="L1HS and ERVL-B4-int:3AERVL:3ALTR Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "deeppink2", "darkkhaki", "lightsalmon3", "brown"))

####ERVL-E-int
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Positive Genes")
ERVLE<- read.table("ERVL-E-int%3AERVL%3ALTR.txt", header=TRUE)
ERVLE<- subset(ERVLE, ERVLE[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- ERVLE[which(ERVLE[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- ERVLE[which(ERVLE[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- ERVLE[which(ERVLE[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- ERVLE[which(ERVLE[,1] %in% LUSC[,1]),]
LUSCwL1HS<- cbind(LUSCcut, LUSC)
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 6)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- ERVLE[which(ERVLE[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "ERVL-E-int:3AERVL:3ALTR (rpm)", y = "L1HS Expression (rpm)", title="L1HS and ERVL-E-int:3AERVL:3ALTR Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###HERVH48
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Positive Genes")
HERV<- read.table("HERVH48-int%3AERV1%3ALTR.txt", header=TRUE)
HERV<- subset(HERV, HERV[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- HERV[which(HERV[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#COAD
COAD<- subset(L1HS, L1HS[,2]=="COAD")
COADcut<- HERV[which(HERV[,1] %in% COAD[,1]),]
COADwL1HS<- cbind(COADcut, COAD)
COADaL1HS<- COADwL1HS[,c(1, 3, 6)]
colnames(COADaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- HERV[which(HERV[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- HERV[which(HERV[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- HERV[which(HERV[,1] %in% LUSC[,1]),]
LUSCwL1HS<- cbind(LUSCcut, LUSC)
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 6)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- HERV[which(HERV[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "HERVH48-int:3AERV1:3ALTR (rpm)", y = "L1HS Expression (rpm)", title="L1HS and HERVH48-int:3AERV1:3ALTR Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "deeppink2", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###MLT1A
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Positive Genes")
MLT1A<- read.table("MLT1A%3AERVL-MaLR%3ALTR.txt", header=TRUE)
MLT1A<- subset(MLT1A, MLT1A[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- MLT1A[which(MLT1A[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#COAD
COAD<- subset(L1HS, L1HS[,2]=="COAD")
COADcut<- MLT1A[which(MLT1A[,1] %in% COAD[,1]),]
COADwL1HS<- cbind(COADcut, COAD)
COADaL1HS<- COADwL1HS[,c(1, 3, 6)]
colnames(COADaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- MLT1A[which(MLT1A[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- MLT1A[which(MLT1A[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- MLT1A[which(MLT1A[,1] %in% LUSC[,1]),]
LUSCwL1HS<- cbind(LUSCcut, LUSC)
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 6)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- MLT1A[which(MLT1A[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "MLT1A:3AERVL-MaLR:3ALTR (rpm)", y = "L1HS Expression (rpm)", title="L1HS and MLT1A:3AERVL-MaLR:3ALTR Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "deeppink2", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###MLT1E1A
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Positive Genes")
MLT1E<- read.table("MLT1E1A-int%3AERVL-MaLR%3ALTR.txt", header=TRUE)
MLT1E<- subset(MLT1E, MLT1E[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- MLT1E[which(MLT1E[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- MLT1E[which(MLT1E[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- MLT1E[which(MLT1E[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- MLT1E[which(MLT1E[,1] %in% LUSC[,1]),]
LUSCwL1HS<- cbind(LUSCcut, LUSC)
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 6)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- MLT1E[which(MLT1E[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "MLT1E1A-int:3AERVL-MaLR:3ALTR (rpm)", y = "L1HS Expression (rpm)", title="L1HS and MLT1E1A-int:3AERVL-MaLR:3ALTR Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

####Tigger1
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS/Positive Genes")
Tig<- read.table("Tigger1%3ATcMar-Tigger%3ADNA.txt", header=TRUE)
Tig- subset(Tig, Tig[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- Tig[which(Tig[,1] %in% BRCA[,1]),]
BRCAwL1HS<- cbind(BRCAcut, BRCA)
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 6)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- Tig[which(Tig[,1] %in% KIRC[,1]),]
KIRCwL1HS<- cbind(KIRCcut, KIRC)
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 6)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- Tig[which(Tig[,1] %in% LUAD[,1]),]
LUADwL1HS<- cbind(LUADcut, LUAD)
LUADaL1HS<- LUADwL1HS[,c(1, 3, 6)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- Tig[which(Tig[,1] %in% LUSC[,1]),]
LUSCwL1HS<- cbind(LUSCcut, LUSC)
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 6)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- Tig[which(Tig[,1] %in% THCA[,1]),]
THCAwL1HS<- cbind(THCAcut, THCA)
THCAaL1HS<- THCAwL1HS[,c(1, 3, 6)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "Tigger1:3ATcMar-Tigger:3ADNA (rpm)", y = "L1HS Expression (rpm)", title="L1HS and Tigger1:3ATcMar-Tigger:3ADNA Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

	
####miRNA Positives
### hsa-mir-141
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m141<- read.table("hsa-mir-141.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m141[which(m141[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- m141[which(m141[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 2, 4)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	labs(x = "hsa-mir-141 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-141 6Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue"))

###hsa-mir-200a
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m200a<- read.table("hsa-mir-200a.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m200a[which(m200a[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- m200a[which(m200a[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 2, 4)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	labs(x = "hsa-mir-200a (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-200a Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue"))

###hsa-mir-203
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m203<- read.table("hsa-mir-203.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m203[which(m203[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- m203[which(m203[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 2, 4)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	labs(x = "hsa-mir-203 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-203 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue"))

###hsa-mir-3065
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m3065<- read.table("hsa-mir-3065.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m3065[which(m3065[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- m3065[which(m3065[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 2, 4)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#PRAD
PRAD<- subset(L1HS, L1HS[,2]=="PRAD")
PRADcut<- m3065[which(m3065[,1] %in% PRAD[,1]),]
PRADwL1HS<- merge(PRADcut, PRAD, by.x="patient")
PRADaL1HS<- PRADwL1HS[,c(1, 2, 4)]
colnames(PRADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
    labs(x = "hsa-mir-3065 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-3065 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "deeppink2"))

###hsa-mir-33a
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m33a<- read.table("hsa-mir-33a.txt", header=TRUE)
#ESCA
ESCA<- subset(L1HS, L1HS[,2]=="ESCA")
ESCAcut<- m33a[which(m33a[,1] %in% ESCA[,1]),]
ESCAwL1HS<- merge(ESCAcut, ESCA, by.x="patient")
ESCAaL1HS<- ESCAwL1HS[,c(1, 2, 4)]
colnames(ESCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- m33a[which(m33a[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 2, 4)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#STAD
STAD<- subset(L1HS, L1HS[,2]=="STAD")
STADcut<- m33a[which(m33a[,1] %in% STAD[,1]),]
STADwL1HS<- merge(STADcut, STAD, by.x="patient")
STADaL1HS<- STADwL1HS[,c(1, 2, 4)]
colnames(STADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
    geom_jitter(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")], color="ESCA")) + geom_smooth(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="coral1") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")], color="STAD")) + geom_smooth(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkred") +
    labs(x = "hsa-mir-33a (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-33a Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("coral1", "darkblue", "darkred"))

###hsa-mir-375
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m375<- read.table("hsa-mir-375.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m375[which(m375[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#PRAD
PRAD<- subset(L1HS, L1HS[,2]=="PRAD")
PRADcut<- m375[which(m375[,1] %in% PRAD[,1]),]
PRADwL1HS<- merge(PRADcut, PRAD, by.x="patient")
PRADaL1HS<- PRADwL1HS[,c(1, 2, 4)]
colnames(PRADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
    labs(x = "hsa-mir-375 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-375 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "deeppink2"))

###hsa-mir-181c
#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m181c<- read.table("hsa-mir-181c.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m181c[which(m181c[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	labs(x = "hsa-mir-181c (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-181c Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen"))

###hsa-mir-200b
#BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m200b<- read.table("hsa-mir-200b.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m200b[which(m200b[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	labs(x = "hsa-mir-200b (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-200b Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen"))

###hsa-mir-642a
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m33a<- read.table("hsa-mir-33a.txt", header=TRUE)
#ESCA
ESCA<- subset(L1HS, L1HS[,2]=="ESCA")
ESCAcut<- m33a[which(m33a[,1] %in% ESCA[,1]),]
ESCAwL1HS<- merge(ESCAcut, ESCA, by.x="patient")
ESCAaL1HS<- ESCAwL1HS[,c(1, 2, 4)]
colnames(ESCAaL1HS)=c("patients", "gene", "L1HS")
#STAD
STAD<- subset(L1HS, L1HS[,2]=="STAD")
STADcut<- m33a[which(m33a[,1] %in% STAD[,1]),]
STADwL1HS<- merge(STADcut, STAD, by.x="patient")
STADaL1HS<- STADwL1HS[,c(1, 2, 4)]
colnames(STADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
    geom_jitter(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")], color="ESCA")) + geom_smooth(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="coral1") +
	geom_jitter(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")], color="STAD")) + geom_smooth(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkred") +
    labs(x = "hsa-mir-642a (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-642a Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("coral1", "darkblue", "darkred"))

###Negative miRNA
###hsa-mir-362
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Negative miRNA")
m362<- read.table("hsa-mir-362.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m362[which(m362[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- m362[which(m362[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 2, 4)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	labs(x = "hsa-mir-362 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-362 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue"))

###hsa-mir-652
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Negative miRNA")
m652<- read.table("hsa-mir-652.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m652[which(m652[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#PRAD
PRAD<- subset(L1HS, L1HS[,2]=="PRAD")
PRADcut<- m652[which(m652[,1] %in% PRAD[,1]),]
PRADwL1HS<- merge(PRADcut, PRAD, by.x="patient")
PRADaL1HS<- PRADwL1HS[,c(1, 2, 4)]
colnames(PRADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
    labs(x = "hsa-mir-652 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-652 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "deeppink2"))

###hsa-mir-215
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Negative miRNA")
m215<- read.table("hsa-mir-215.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m215[which(m215[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	labs(x = "hsa-mir-215 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-215 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen"))

###hsa-mir-486
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Negative miRNA")m486<- read.table("hsa-mir-486.txt", header=TRUE)
m486<- read.table("hsa-mir-486.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m486[which(m486[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	labs(x = "hsa-mir-486 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-486 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen"))

#hsa-mir-204
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Negative miRNA")
m204<- read.table("hsa-mir-204.txt", header=TRUE)
#PRAD
PRAD<- subset(L1HS, L1HS[,2]=="PRAD")
PRADcut<- m204[which(m204[,1] %in% PRAD[,1]),]
PRADwL1HS<- merge(PRADcut, PRAD, by.x="patient")
PRADaL1HS<- PRADwL1HS[,c(1, 2, 4)]
colnames(PRADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
    labs(x = "hsa-mir-204 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-204 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink2"))

###hsa-mir-363
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Negative miRNA")
m363<- read.table("hsa-mir-363.txt", header=TRUE)
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m363[which(m363[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 2, 4)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	labs(x = "hsa-mir-363 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and hsa-mir-363 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen"))

####GDC positive
####Now start graph
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC")
L1HS<- read.table("L1HS normal for Graphs.txt", header=TRUE)

###ARHGAP32
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
ARH<- read.table("ARHGAP32%7C9743.txt", header=TRUE)
ARH<- subset(ARH, ARH[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- ARH[which(ARH[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- ARH[which(ARH[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- ARH[which(ARH[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- ARH[which(ARH[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "ARHGAP32 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and ARHGAP32 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3"))

###BAT2L2
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
BAT<- read.table("BAT2L2%7C23215.txt", header=TRUE)
BAT<- subset(BAT, BAT[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- BAT[which(BAT[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- BAT[which(BAT[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- BAT[which(BAT[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- BAT[which(BAT[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- BAT[which(BAT[,1] %in% LUSC[,1]),]
LUSCwL1HS<- merge(LUSCcut, LUSC, by.x="patient")
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 5)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- BAT[which(BAT[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "BAT2L2 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and BAT2L2 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###FLJ45340
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
FLJ<- read.table("FLJ45340%7C402483.txt", header=TRUE)
FLJ<- subset(FLJ, FLJ[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- FLJ[which(FLJ[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- FLJ[which(FLJ[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- FLJ[which(FLJ[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- FLJ[which(FLJ[,1] %in% LUSC[,1]),]
LUSCwL1HS<- merge(LUSCcut, LUSC, by.x="patient")
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 5)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- FLJ[which(FLJ[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "FLJ45340 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and FLJ45340 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###LOC10019086
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
LOC<- read.table("LOC100190986%7C100190986.txt", header=TRUE)
LOC<- subset(LOC, LOC[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- LOC[which(LOC[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- LOC[which(LOC[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- LOC[which(LOC[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- LOC[which(LOC[,1] %in% LUSC[,1]),]
LUSCwL1HS<- merge(LUSCcut, LUSC, by.x="patient")
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 5)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- LOC[which(LOC[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "LOC100190986 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and LOC100190986 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###NKTR
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
NKT<- read.table("NKTR%7C4820.txt", header=TRUE)
NKT<- subset(NKT, NKT[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- NKT[which(NKT[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- NKT[which(NKT[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- NKT[which(NKT[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- NKT[which(NKT[,1] %in% LUSC[,1]),]
LUSCwL1HS<- merge(LUSCcut, LUSC, by.x="patient")
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 5)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- NKT[which(NKT[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "NKTR (rpm)", y = "L1HS Expression (rpm)", title="L1HS and NKTR Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###THOC2
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
THO<- read.table("THOC2%7C57187.txt", header=TRUE)
THO<- subset(THO, THO[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- THO[which(THO[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- THO[which(THO[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- THO[which(THO[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- THO[which(THO[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- THO[which(THO[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "THOC2 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and THOC2 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3", "brown"))

###UBN2
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
UBN<- read.table("UBN2%7C254048.txt", header=TRUE)
UBN<- subset(UBN, UBN[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- UBN[which(UBN[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- UBN[which(UBN[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- UBN[which(UBN[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- UBN[which(UBN[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- UBN[which(UBN[,1] %in% LUSC[,1]),]
LUSCwL1HS<- merge(LUSCcut, LUSC, by.x="patient")
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 5)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- UBN[which(UBN[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "UBN2 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and UBN2 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3","darkslategray", "brown"))

###WDR52
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
WDR<- read.table("WDR52%7C55779.txt", header=TRUE)
WDR<- subset(WDR, WDR[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- WDR[which(WDR[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- WDR[which(WDR[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- WDR[which(WDR[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- WDR[which(WDR[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "WDR52 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and WDR52 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3", "brown"))

###ZNF37B
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
ZNF<- read.table("ZNF37B%7C100129482.txt", header=TRUE)
ZNF<- subset(ZNF, ZNF[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- ZNF[which(ZNF[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- ZNF[which(ZNF[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- ZNF[which(ZNF[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- ZNF[which(ZNF[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "ZNF37B (rpm)", y = "L1HS Expression (rpm)", title="L1HS and ZNF37B Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3", "brown"))

###ZNF587
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Positive genes")
ZNF5<- read.table("ZNF587%7C84914.txt", header=TRUE)
ZNF5<- subset(ZNF5, ZNF5[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- ZNF5[which(ZNF5[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- ZNF5[which(ZNF5[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- ZNF5[which(ZNF5[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#LUSC
LUSC<- subset(L1HS, L1HS[,2]=="LUSC")
LUSCcut<- ZNF5[which(ZNF5[,1] %in% LUSC[,1]),]
LUSCwL1HS<- merge(LUSCcut, LUSC, by.x="patient")
LUSCaL1HS<- LUSCwL1HS[,c(1, 3, 5)]
colnames(LUSCaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
	labs(x = "ZNF587 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and ZNF587 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen",  "darkkhaki", "lightsalmon3","darkslategray"))

####Negative GDC genes
###AIP
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
AIP<- read.table("AIP%7C9049.txt", header=TRUE)
AIP<- subset(AIP, AIP[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- AIP[which(AIP[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- AIP[which(AIP[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- AIP[which(AIP[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- AIP[which(AIP[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "AIP (rpm)", y = "L1HS Expression (rpm)", title="L1HS and AIP Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3"))

	###ARL2
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
ARL<- read.table("ARL2%7C402.txt", header=TRUE)
ARL<- subset(ARL, ARL[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- ARL[which(ARL[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- ARL[which(ARL[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- ARL[which(ARL[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "ARL2 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and ARL2 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "lightsalmon3"))

###FAM195B
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
FAM<- read.table("FAM195B%7C348262.txt", header=TRUE)
FAM<- subset(FAM, FAM[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- FAM[which(FAM[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- FAM[which(FAM[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- FAM[which(FAM[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- FAM[which(FAM[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "FAM195B (rpm)", y = "L1HS Expression (rpm)", title="L1HS and FAM195B Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3"))

###FLOT1
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
FLOT<- read.table("FLOT1%7C10211.txt", header=TRUE)
FLOT<- subset(FLOT, FLOT[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- FLOT[which(FLOT[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- FLOT[which(FLOT[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- FLOT[which(FLOT[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- FLOT[which(FLOT[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "FLOT1 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and FLOT1 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3"))

###FUNDC2
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
FUND<- read.table("FUNDC2%7C65991.txt", header=TRUE)
FUND<- subset(FUND, FUND[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- FUND[which(FUND[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- FUND[which(FUND[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- FUND[which(FUND[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- FUND[which(FUND[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "FUNDC2 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and FUNDC2 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3"))

###PARK7
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
PARK<- read.table("PARK7%7C11315.txt", header=TRUE)
PARK<- subset(PARK, PARK[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- PARK[which(PARK[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- PARK[which(PARK[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- PARK[which(PARK[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- PARK[which(PARK[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "PARK7 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and PARK7 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "darkkhaki", "lightsalmon3"))

###RAB5C
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
RAB<- read.table("RAB5C%7C5878.txt", header=TRUE)
RAB<- subset(RAB, RAB[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- RAB[which(RAB[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- RAB[which(RAB[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "RAB5C (rpm)", y = "L1HS Expression (rpm)", title="L1HS and RAB5C Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "lightsalmon3"))

###SLC29A1
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
SLC<- read.table("SLC39A1%7C27173.txt", header=TRUE)
SLC<- subset(SLC, SLC[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- SLC[which(SLC[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- SLC[which(SLC[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- SLC[which(SLC[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- SLC[which(SLC[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "SLC39A1 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and SLC39A1 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3", "brown"))

###SNX17
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
SNX<- read.table("SNX17%7C9784.txt", header=TRUE)
SNX<- subset(SNX, SNX[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- SNX[which(SNX[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#KIRC
KIRC<- subset(L1HS, L1HS[,2]=="KIRC")
KIRCcut<- SNX[which(SNX[,1] %in% KIRC[,1]),]
KIRCwL1HS<- merge(KIRCcut, KIRC, by.x="patient")
KIRCaL1HS<- KIRCwL1HS[,c(1, 3, 5)]
colnames(KIRCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- SNX[which(SNX[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#THCA
THCA<- subset(L1HS, L1HS[,2]=="THCA")
THCAcut<- SNX[which(SNX[,1] %in% THCA[,1]),]
THCAwL1HS<- merge(THCAcut, THCA, by.x="patient")
THCAaL1HS<- THCAwL1HS[,c(1, 3, 5)]
colnames(THCAaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown") +
	labs(x = "SNX17 (rpm)", y = "L1HS Expression (rpm)", title="L1HS and SNX17 Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkkhaki", "lightsalmon3", "brown"))

###UROD
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC/Negative genes")
UROD<- read.table("UROD%7C7389.txt", header=TRUE)
UROD<- subset(UROD, UROD[,2]=="normal")
#BRCA
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- UROD[which(UROD[,1] %in% BRCA[,1]),]
BRCAwL1HS<- merge(BRCAcut, BRCA, by.x="patient")
BRCAaL1HS<- BRCAwL1HS[,c(1, 3, 5)]
colnames(BRCAaL1HS)=c("patients", "gene", "L1HS")
#HNSC
HNSC<- subset(L1HS, L1HS[,2]=="HNSC")
HNSCcut<- UROD[which(UROD[,1] %in% HNSC[,1]),]
HNSCwL1HS<- merge(HNSCcut, HNSC, by.x="patient")
HNSCaL1HS<- HNSCwL1HS[,c(1, 3, 5)]
colnames(HNSCaL1HS)=c("patients", "gene", "L1HS")
#LUAD
LUAD<- subset(L1HS, L1HS[,2]=="LUAD")
LUADcut<- UROD[which(UROD[,1] %in% LUAD[,1]),]
LUADwL1HS<- merge(LUADcut, LUAD, by.x="patient")
LUADaL1HS<- LUADwL1HS[,c(1, 3, 5)]
colnames(LUADaL1HS)=c("patients", "gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
	geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
	geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
	geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
	labs(x = "UROD (rpm)", y = "L1HS Expression (rpm)", title="L1HS and UROD Correlation")+
	theme_classic()+
	theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "darkblue", "lightsalmon3"))

######Figure out the patients missing from BRCA
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/miRNA/Positive miRNA")
m141<- read.table("hsa-mir-141.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 4/gene2L1HS.GDC")
L1HS<- read.table("L1HS normal for Graphs.txt", header=TRUE)
BRCA<- subset(L1HS, L1HS[,2]=="BRCA")
BRCAcut<- m141[which(m141[,1] %in% BRCA[,1]),]
inter<- intersect(BRCA[,1], BRCAcut[,1])
inter<- data.frame(inter)
colnames(inter)=c("patient")
merge(inter, BRCA, by.x="patient", all=TRUE)
library("dplyr")
bind_cols(inter, BRCA)