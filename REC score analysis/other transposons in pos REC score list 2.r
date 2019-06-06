###New transposon graphs with the top selected cancers
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs for REC scores top")
patients<- read.table("patient_type.txt", header=TRUE)

#####MLT transposon
MLT<- read.table("MLT1E1A-int%3AERVL-MaLR%3ALTR.txt", header=TRUE, row.names=1)
MLT<- cbind(MLT, rownames(MLT))
#BLCA
BLCA<- subset(patients, patients[,2]=="BLCA")
BLCAcut<- MLT[which(row.names(MLT) %in% BLCA[,1]),]
BLCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BLCAcut[,2]),]
rownames(BLCAwL1HS)=BLCAwL1HS[,1]
colnames(BLCAwL1HS)=c("patient", "L1HS")
BLCAaL1HS<-merge(BLCAcut, BLCAwL1HS, all=TRUE, by="row.names")
rownames(BLCAaL1HS)=BLCAaL1HS[,1]
BLCAaL1HS<- data.frame(BLCAaL1HS[,-c(1,3,4)])
colnames(BLCAaL1HS)=c("gene", "L1HS")
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- MLT[which(row.names(MLT) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#CHOL
CHOL<- subset(patients, patients[,2]=="CHOL")
CHOLcut<- MLT[which(row.names(MLT) %in% CHOL[,1]),]
CHOLwL1HS<- L1HSread[which(L1HSread[,1] %in% CHOLcut[,2]),]
rownames(CHOLwL1HS)=CHOLwL1HS[,1]
colnames(CHOLwL1HS)=c("patient", "L1HS")
CHOLaL1HS<-merge(CHOLcut, CHOLwL1HS, all=TRUE, by="row.names")
rownames(CHOLaL1HS)=CHOLaL1HS[,1]
CHOLaL1HS<- data.frame(CHOLaL1HS[,-c(1,3,4)])
colnames(CHOLaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- MLT[which(row.names(MLT) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#ESCA
ESCA<- subset(patients, patients[,2]=="ESCA")
ESCAcut<- MLT[which(row.names(MLT) %in% ESCA[,1]),]
ESCAwL1HS<- L1HSread[which(L1HSread[,1] %in% ESCAcut[,2]),]
rownames(ESCAwL1HS)=ESCAwL1HS[,1]
colnames(ESCAwL1HS)=c("patient", "L1HS")
ESCAaL1HS<-merge(ESCAcut, ESCAwL1HS, all=TRUE, by="row.names")
rownames(ESCAaL1HS)=ESCAaL1HS[,1]
ESCAaL1HS<- data.frame(ESCAaL1HS[,-c(1,3,4)])
colnames(ESCAaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- MLT[which(row.names(MLT) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KICH
KICH<- subset(patients, patients[,2]=="KICH")
KICHcut<- MLT[which(row.names(MLT) %in% KICH[,1]),]
KICHwL1HS<- L1HSread[which(L1HSread[,1] %in% KICHcut[,2]),]
rownames(KICHwL1HS)=KICHwL1HS[,1]
colnames(KICHwL1HS)=c("patient", "L1HS")
KICHaL1HS<-merge(KICHcut, KICHwL1HS, all=TRUE, by="row.names")
rownames(KICHaL1HS)=KICHaL1HS[,1]
KICHaL1HS<- data.frame(KICHaL1HS[,-c(1,3,4)])
colnames(KICHaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- MLT[which(row.names(MLT) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- MLT[which(row.names(MLT) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LIHC
LIHC<- subset(patients, patients[,2]=="LIHC")
LIHCcut<- MLT[which(row.names(MLT) %in% LIHC[,1]),]
LIHCwL1HS<- L1HSread[which(L1HSread[,1] %in% LIHCcut[,2]),]
rownames(LIHCwL1HS)=LIHCwL1HS[,1]
colnames(LIHCwL1HS)=c("patient", "L1HS")
LIHCaL1HS<-merge(LIHCcut, LIHCwL1HS, all=TRUE, by="row.names")
rownames(LIHCaL1HS)=LIHCaL1HS[,1]
LIHCaL1HS<- data.frame(LIHCaL1HS[,-c(1,3,4)])
colnames(LIHCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- MLT[which(row.names(MLT) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- MLT[which(row.names(MLT) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#PRAD
PRAD<- subset(patients, patients[,2]=="PRAD")
PRADcut<- MLT[which(row.names(MLT) %in% PRAD[,1]),]
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRADcut[,2]),]
rownames(PRADwL1HS)=PRADwL1HS[,1]
colnames(PRADwL1HS)=c("patient", "L1HS")
PRADaL1HS<-merge(PRADcut, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)=PRADaL1HS[,1]
PRADaL1HS<- data.frame(PRADaL1HS[,-c(1,3,4)])
colnames(PRADaL1HS)=c("gene", "L1HS")
#STAD
STAD<- subset(patients, patients[,2]=="STAD")
STADcut<- MLT[which(row.names(MLT) %in% STAD[,1]),]
STADwL1HS<- L1HSread[which(L1HSread[,1] %in% STADcut[,2]),]
rownames(STADwL1HS)=STADwL1HS[,1]
colnames(STADwL1HS)=c("patient", "L1HS")
STADaL1HS<-merge(STADcut, STADwL1HS, all=TRUE, by="row.names")
rownames(STADaL1HS)=STADaL1HS[,1]
STADaL1HS<- data.frame(STADaL1HS[,-c(1,3,4)])
colnames(STADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- MLT[which(row.names(MLT) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
#UCEC
UCEC<- subset(patients, patients[,2]=="UCEC")
UCECcut<- MLT[which(row.names(MLT) %in% UCEC[,1]),]
UCECwL1HS<- L1HSread[which(L1HSread[,1] %in% UCECcut[,2]),]
rownames(UCECwL1HS)=UCECwL1HS[,1]
colnames(UCECwL1HS)=c("patient", "L1HS")
UCECaL1HS<-merge(UCECcut, UCECwL1HS, all=TRUE, by="row.names")
rownames(UCECaL1HS)=UCECaL1HS[,1]
UCECaL1HS<- data.frame(UCECaL1HS[,-c(1,3,4)])
colnames(UCECaL1HS)=c("gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")], color="BLCA")) + geom_smooth(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="burlywood1") +
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="cadetblue1") +
  geom_jitter(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")], color="CHOL")) + geom_smooth(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chartreuse1") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")], color="ESCA")) + geom_smooth(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="coral1") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")], color="KICH")) + geom_smooth(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkgoldenrod3") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")], color="LIHC")) + geom_smooth(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkorange3") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
  geom_jitter(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")], color="STAD")) + geom_smooth(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkred") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  geom_jitter(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")], color="UCEC")) + geom_smooth(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  labs(x = "MLT1E1A-int:3AERVL-MaLR:3ALTR Expression (rpm)", y = "L1HS Expression (rpm)", title="L1HS and MLT1E1A-int:3AERVL-MaLR:3ALTR Correlation")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("burlywood1","cadetblue1", 
		"chartreuse1", "chocolate1", "coral1", "darkblue", "darkgoldenrod3", 
		"darkkhaki", "black", "darkorange3", "lightsalmon3", "darkslategray",
		"deeppink2", "darkred", "brown3", "springgreen3")) 

#####HERVH3
HERVH3<- read.table("HERVH48-int%3AERV1%3ALTR.txt", header=TRUE, row.names=1)
HERVH3<- cbind(HERVH3, rownames(HERVH3))
#BLCA
BLCA<- subset(patients, patients[,2]=="BLCA")
BLCAcut<- HERVH3[which(row.names(HERVH3) %in% BLCA[,1]),]
BLCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BLCAcut[,2]),]
rownames(BLCAwL1HS)=BLCAwL1HS[,1]
colnames(BLCAwL1HS)=c("patient", "L1HS")
BLCAaL1HS<-merge(BLCAcut, BLCAwL1HS, all=TRUE, by="row.names")
rownames(BLCAaL1HS)=BLCAaL1HS[,1]
BLCAaL1HS<- data.frame(BLCAaL1HS[,-c(1,3,4)])
colnames(BLCAaL1HS)=c("gene", "L1HS")
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- HERVH3[which(row.names(HERVH3) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#CHOL
CHOL<- subset(patients, patients[,2]=="CHOL")
CHOLcut<- HERVH3[which(row.names(HERVH3) %in% CHOL[,1]),]
CHOLwL1HS<- L1HSread[which(L1HSread[,1] %in% CHOLcut[,2]),]
rownames(CHOLwL1HS)=CHOLwL1HS[,1]
colnames(CHOLwL1HS)=c("patient", "L1HS")
CHOLaL1HS<-merge(CHOLcut, CHOLwL1HS, all=TRUE, by="row.names")
rownames(CHOLaL1HS)=CHOLaL1HS[,1]
CHOLaL1HS<- data.frame(CHOLaL1HS[,-c(1,3,4)])
colnames(CHOLaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- HERVH3[which(row.names(HERVH3) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#ESCA
ESCA<- subset(patients, patients[,2]=="ESCA")
ESCAcut<- HERVH3[which(row.names(HERVH3) %in% ESCA[,1]),]
ESCAwL1HS<- L1HSread[which(L1HSread[,1] %in% ESCAcut[,2]),]
rownames(ESCAwL1HS)=ESCAwL1HS[,1]
colnames(ESCAwL1HS)=c("patient", "L1HS")
ESCAaL1HS<-merge(ESCAcut, ESCAwL1HS, all=TRUE, by="row.names")
rownames(ESCAaL1HS)=ESCAaL1HS[,1]
ESCAaL1HS<- data.frame(ESCAaL1HS[,-c(1,3,4)])
colnames(ESCAaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- HERVH3[which(row.names(HERVH3) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KICH
KICH<- subset(patients, patients[,2]=="KICH")
KICHcut<- HERVH3[which(row.names(HERVH3) %in% KICH[,1]),]
KICHwL1HS<- L1HSread[which(L1HSread[,1] %in% KICHcut[,2]),]
rownames(KICHwL1HS)=KICHwL1HS[,1]
colnames(KICHwL1HS)=c("patient", "L1HS")
KICHaL1HS<-merge(KICHcut, KICHwL1HS, all=TRUE, by="row.names")
rownames(KICHaL1HS)=KICHaL1HS[,1]
KICHaL1HS<- data.frame(KICHaL1HS[,-c(1,3,4)])
colnames(KICHaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- HERVH3[which(row.names(HERVH3) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- HERVH3[which(row.names(HERVH3) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LIHC
LIHC<- subset(patients, patients[,2]=="LIHC")
LIHCcut<- HERVH3[which(row.names(HERVH3) %in% LIHC[,1]),]
LIHCwL1HS<- L1HSread[which(L1HSread[,1] %in% LIHCcut[,2]),]
rownames(LIHCwL1HS)=LIHCwL1HS[,1]
colnames(LIHCwL1HS)=c("patient", "L1HS")
LIHCaL1HS<-merge(LIHCcut, LIHCwL1HS, all=TRUE, by="row.names")
rownames(LIHCaL1HS)=LIHCaL1HS[,1]
LIHCaL1HS<- data.frame(LIHCaL1HS[,-c(1,3,4)])
colnames(LIHCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- HERVH3[which(row.names(HERVH3) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- HERVH3[which(row.names(HERVH3) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#PRAD
PRAD<- subset(patients, patients[,2]=="PRAD")
PRADcut<- HERVH3[which(row.names(HERVH3) %in% PRAD[,1]),]
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRADcut[,2]),]
rownames(PRADwL1HS)=PRADwL1HS[,1]
colnames(PRADwL1HS)=c("patient", "L1HS")
PRADaL1HS<-merge(PRADcut, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)=PRADaL1HS[,1]
PRADaL1HS<- data.frame(PRADaL1HS[,-c(1,3,4)])
colnames(PRADaL1HS)=c("gene", "L1HS")
#STAD
STAD<- subset(patients, patients[,2]=="STAD")
STADcut<- HERVH3[which(row.names(HERVH3) %in% STAD[,1]),]
STADwL1HS<- L1HSread[which(L1HSread[,1] %in% STADcut[,2]),]
rownames(STADwL1HS)=STADwL1HS[,1]
colnames(STADwL1HS)=c("patient", "L1HS")
STADaL1HS<-merge(STADcut, STADwL1HS, all=TRUE, by="row.names")
rownames(STADaL1HS)=STADaL1HS[,1]
STADaL1HS<- data.frame(STADaL1HS[,-c(1,3,4)])
colnames(STADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- HERVH3[which(row.names(HERVH3) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
#UCEC
UCEC<- subset(patients, patients[,2]=="UCEC")
UCECcut<- HERVH3[which(row.names(HERVH3) %in% UCEC[,1]),]
UCECwL1HS<- L1HSread[which(L1HSread[,1] %in% UCECcut[,2]),]
rownames(UCECwL1HS)=UCECwL1HS[,1]
colnames(UCECwL1HS)=c("patient", "L1HS")
UCECaL1HS<-merge(UCECcut, UCECwL1HS, all=TRUE, by="row.names")
rownames(UCECaL1HS)=UCECaL1HS[,1]
UCECaL1HS<- data.frame(UCECaL1HS[,-c(1,3,4)])
colnames(UCECaL1HS)=c("gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")], color="BLCA")) + geom_smooth(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="burlywood1") +
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="cadetblue1") +
  geom_jitter(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")], color="CHOL")) + geom_smooth(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chartreuse1") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")], color="ESCA")) + geom_smooth(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="coral1") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")], color="KICH")) + geom_smooth(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkgoldenrod3") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")], color="LIHC")) + geom_smooth(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkorange3") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
  geom_jitter(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")], color="STAD")) + geom_smooth(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkred") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  geom_jitter(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")], color="UCEC")) + geom_smooth(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  labs(x = "HERVH48-int:3AERV1:3ALTR Expression (rpm)", y = "L1HS Expression (rpm)", title="HERVH48-int:3AERV1:3ALTR Correlation")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("burlywood1","cadetblue1", 
		"chartreuse1", "chocolate1", "coral1", "darkblue", "darkgoldenrod3", 
		"darkkhaki", "black", "darkorange3", "lightsalmon3", "darkslategray",
		"deeppink2", "darkred", "brown3", "springgreen3")) 

####Tigger1
tig1<- read.table("Tigger1%3ATcMar-Tigger%3ADNA.txt", header=TRUE, row.names=1)
tig1<- cbind(tig1, rownames(tig1))
#BLCA
BLCA<- subset(patients, patients[,2]=="BLCA")
BLCAcut<- tig1[which(row.names(tig1) %in% BLCA[,1]),]
BLCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BLCAcut[,2]),]
rownames(BLCAwL1HS)=BLCAwL1HS[,1]
colnames(BLCAwL1HS)=c("patient", "L1HS")
BLCAaL1HS<-merge(BLCAcut, BLCAwL1HS, all=TRUE, by="row.names")
rownames(BLCAaL1HS)=BLCAaL1HS[,1]
BLCAaL1HS<- data.frame(BLCAaL1HS[,-c(1,3,4)])
colnames(BLCAaL1HS)=c("gene", "L1HS")
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- tig1[which(row.names(tig1) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#CHOL
CHOL<- subset(patients, patients[,2]=="CHOL")
CHOLcut<- tig1[which(row.names(tig1) %in% CHOL[,1]),]
CHOLwL1HS<- L1HSread[which(L1HSread[,1] %in% CHOLcut[,2]),]
rownames(CHOLwL1HS)=CHOLwL1HS[,1]
colnames(CHOLwL1HS)=c("patient", "L1HS")
CHOLaL1HS<-merge(CHOLcut, CHOLwL1HS, all=TRUE, by="row.names")
rownames(CHOLaL1HS)=CHOLaL1HS[,1]
CHOLaL1HS<- data.frame(CHOLaL1HS[,-c(1,3,4)])
colnames(CHOLaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- tig1[which(row.names(tig1) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#ESCA
ESCA<- subset(patients, patients[,2]=="ESCA")
ESCAcut<- tig1[which(row.names(tig1) %in% ESCA[,1]),]
ESCAwL1HS<- L1HSread[which(L1HSread[,1] %in% ESCAcut[,2]),]
rownames(ESCAwL1HS)=ESCAwL1HS[,1]
colnames(ESCAwL1HS)=c("patient", "L1HS")
ESCAaL1HS<-merge(ESCAcut, ESCAwL1HS, all=TRUE, by="row.names")
rownames(ESCAaL1HS)=ESCAaL1HS[,1]
ESCAaL1HS<- data.frame(ESCAaL1HS[,-c(1,3,4)])
colnames(ESCAaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- tig1[which(row.names(tig1) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KICH
KICH<- subset(patients, patients[,2]=="KICH")
KICHcut<- tig1[which(row.names(tig1) %in% KICH[,1]),]
KICHwL1HS<- L1HSread[which(L1HSread[,1] %in% KICHcut[,2]),]
rownames(KICHwL1HS)=KICHwL1HS[,1]
colnames(KICHwL1HS)=c("patient", "L1HS")
KICHaL1HS<-merge(KICHcut, KICHwL1HS, all=TRUE, by="row.names")
rownames(KICHaL1HS)=KICHaL1HS[,1]
KICHaL1HS<- data.frame(KICHaL1HS[,-c(1,3,4)])
colnames(KICHaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- tig1[which(row.names(tig1) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- tig1[which(row.names(tig1) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LIHC
LIHC<- subset(patients, patients[,2]=="LIHC")
LIHCcut<- tig1[which(row.names(tig1) %in% LIHC[,1]),]
LIHCwL1HS<- L1HSread[which(L1HSread[,1] %in% LIHCcut[,2]),]
rownames(LIHCwL1HS)=LIHCwL1HS[,1]
colnames(LIHCwL1HS)=c("patient", "L1HS")
LIHCaL1HS<-merge(LIHCcut, LIHCwL1HS, all=TRUE, by="row.names")
rownames(LIHCaL1HS)=LIHCaL1HS[,1]
LIHCaL1HS<- data.frame(LIHCaL1HS[,-c(1,3,4)])
colnames(LIHCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- tig1[which(row.names(tig1) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- tig1[which(row.names(tig1) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#PRAD
PRAD<- subset(patients, patients[,2]=="PRAD")
PRADcut<- tig1[which(row.names(tig1) %in% PRAD[,1]),]
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRADcut[,2]),]
rownames(PRADwL1HS)=PRADwL1HS[,1]
colnames(PRADwL1HS)=c("patient", "L1HS")
PRADaL1HS<-merge(PRADcut, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)=PRADaL1HS[,1]
PRADaL1HS<- data.frame(PRADaL1HS[,-c(1,3,4)])
colnames(PRADaL1HS)=c("gene", "L1HS")
#STAD
STAD<- subset(patients, patients[,2]=="STAD")
STADcut<- tig1[which(row.names(tig1) %in% STAD[,1]),]
STADwL1HS<- L1HSread[which(L1HSread[,1] %in% STADcut[,2]),]
rownames(STADwL1HS)=STADwL1HS[,1]
colnames(STADwL1HS)=c("patient", "L1HS")
STADaL1HS<-merge(STADcut, STADwL1HS, all=TRUE, by="row.names")
rownames(STADaL1HS)=STADaL1HS[,1]
STADaL1HS<- data.frame(STADaL1HS[,-c(1,3,4)])
colnames(STADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- tig1[which(row.names(tig1) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
#UCEC
UCEC<- subset(patients, patients[,2]=="UCEC")
UCECcut<- tig1[which(row.names(tig1) %in% UCEC[,1]),]
UCECwL1HS<- L1HSread[which(L1HSread[,1] %in% UCECcut[,2]),]
rownames(UCECwL1HS)=UCECwL1HS[,1]
colnames(UCECwL1HS)=c("patient", "L1HS")
UCECaL1HS<-merge(UCECcut, UCECwL1HS, all=TRUE, by="row.names")
rownames(UCECaL1HS)=UCECaL1HS[,1]
UCECaL1HS<- data.frame(UCECaL1HS[,-c(1,3,4)])
colnames(UCECaL1HS)=c("gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")], color="BLCA")) + geom_smooth(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="burlywood1") +
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="cadetblue1") +
  geom_jitter(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")], color="CHOL")) + geom_smooth(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chartreuse1") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")], color="ESCA")) + geom_smooth(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="coral1") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")], color="KICH")) + geom_smooth(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkgoldenrod3") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")], color="LIHC")) + geom_smooth(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkorange3") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
  geom_jitter(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")], color="STAD")) + geom_smooth(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkred") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  geom_jitter(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")], color="UCEC")) + geom_smooth(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  labs(x = "Tigger1:3ATcMar-Tigger:3ADNA Expression (rpm)", y = "L1HS Expression (rpm)", title="Tigger1:3ATcMar-Tigger:3ADNA Correlation")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("burlywood1","cadetblue1", 
		"chartreuse1", "chocolate1", "coral1", "darkblue", "darkgoldenrod3", 
		"darkkhaki", "black", "darkorange3", "lightsalmon3", "darkslategray",
		"deeppink2", "darkred", "brown3", "springgreen3")) 

###Tigger 2
tig2<- read.table("Tigger2%3ATcMar-Tigger%3ADNA.txt", header=TRUE, row.names=1)
tig2<- cbind(tig2, rownames(tig2))
#BLCA
BLCA<- subset(patients, patients[,2]=="BLCA")
BLCAcut<- tig2[which(row.names(tig2) %in% BLCA[,1]),]
BLCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BLCAcut[,2]),]
rownames(BLCAwL1HS)=BLCAwL1HS[,1]
colnames(BLCAwL1HS)=c("patient", "L1HS")
BLCAaL1HS<-merge(BLCAcut, BLCAwL1HS, all=TRUE, by="row.names")
rownames(BLCAaL1HS)=BLCAaL1HS[,1]
BLCAaL1HS<- data.frame(BLCAaL1HS[,-c(1,3,4)])
colnames(BLCAaL1HS)=c("gene", "L1HS")
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- tig2[which(row.names(tig2) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#CHOL
CHOL<- subset(patients, patients[,2]=="CHOL")
CHOLcut<- tig2[which(row.names(tig2) %in% CHOL[,1]),]
CHOLwL1HS<- L1HSread[which(L1HSread[,1] %in% CHOLcut[,2]),]
rownames(CHOLwL1HS)=CHOLwL1HS[,1]
colnames(CHOLwL1HS)=c("patient", "L1HS")
CHOLaL1HS<-merge(CHOLcut, CHOLwL1HS, all=TRUE, by="row.names")
rownames(CHOLaL1HS)=CHOLaL1HS[,1]
CHOLaL1HS<- data.frame(CHOLaL1HS[,-c(1,3,4)])
colnames(CHOLaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- tig2[which(row.names(tig2) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#ESCA
ESCA<- subset(patients, patients[,2]=="ESCA")
ESCAcut<- tig2[which(row.names(tig2) %in% ESCA[,1]),]
ESCAwL1HS<- L1HSread[which(L1HSread[,1] %in% ESCAcut[,2]),]
rownames(ESCAwL1HS)=ESCAwL1HS[,1]
colnames(ESCAwL1HS)=c("patient", "L1HS")
ESCAaL1HS<-merge(ESCAcut, ESCAwL1HS, all=TRUE, by="row.names")
rownames(ESCAaL1HS)=ESCAaL1HS[,1]
ESCAaL1HS<- data.frame(ESCAaL1HS[,-c(1,3,4)])
colnames(ESCAaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- tig2[which(row.names(tig2) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KICH
KICH<- subset(patients, patients[,2]=="KICH")
KICHcut<- tig2[which(row.names(tig2) %in% KICH[,1]),]
KICHwL1HS<- L1HSread[which(L1HSread[,1] %in% KICHcut[,2]),]
rownames(KICHwL1HS)=KICHwL1HS[,1]
colnames(KICHwL1HS)=c("patient", "L1HS")
KICHaL1HS<-merge(KICHcut, KICHwL1HS, all=TRUE, by="row.names")
rownames(KICHaL1HS)=KICHaL1HS[,1]
KICHaL1HS<- data.frame(KICHaL1HS[,-c(1,3,4)])
colnames(KICHaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- tig2[which(row.names(tig2) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- tig2[which(row.names(tig2) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LIHC
LIHC<- subset(patients, patients[,2]=="LIHC")
LIHCcut<- tig2[which(row.names(tig2) %in% LIHC[,1]),]
LIHCwL1HS<- L1HSread[which(L1HSread[,1] %in% LIHCcut[,2]),]
rownames(LIHCwL1HS)=LIHCwL1HS[,1]
colnames(LIHCwL1HS)=c("patient", "L1HS")
LIHCaL1HS<-merge(LIHCcut, LIHCwL1HS, all=TRUE, by="row.names")
rownames(LIHCaL1HS)=LIHCaL1HS[,1]
LIHCaL1HS<- data.frame(LIHCaL1HS[,-c(1,3,4)])
colnames(LIHCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- tig2[which(row.names(tig2) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- tig2[which(row.names(tig2) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#PRAD
PRAD<- subset(patients, patients[,2]=="PRAD")
PRADcut<- tig2[which(row.names(tig2) %in% PRAD[,1]),]
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRADcut[,2]),]
rownames(PRADwL1HS)=PRADwL1HS[,1]
colnames(PRADwL1HS)=c("patient", "L1HS")
PRADaL1HS<-merge(PRADcut, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)=PRADaL1HS[,1]
PRADaL1HS<- data.frame(PRADaL1HS[,-c(1,3,4)])
colnames(PRADaL1HS)=c("gene", "L1HS")
#STAD
STAD<- subset(patients, patients[,2]=="STAD")
STADcut<- tig2[which(row.names(tig2) %in% STAD[,1]),]
STADwL1HS<- L1HSread[which(L1HSread[,1] %in% STADcut[,2]),]
rownames(STADwL1HS)=STADwL1HS[,1]
colnames(STADwL1HS)=c("patient", "L1HS")
STADaL1HS<-merge(STADcut, STADwL1HS, all=TRUE, by="row.names")
rownames(STADaL1HS)=STADaL1HS[,1]
STADaL1HS<- data.frame(STADaL1HS[,-c(1,3,4)])
colnames(STADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- tig2[which(row.names(tig2) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
#UCEC
UCEC<- subset(patients, patients[,2]=="UCEC")
UCECcut<- tig2[which(row.names(tig2) %in% UCEC[,1]),]
UCECwL1HS<- L1HSread[which(L1HSread[,1] %in% UCECcut[,2]),]
rownames(UCECwL1HS)=UCECwL1HS[,1]
colnames(UCECwL1HS)=c("patient", "L1HS")
UCECaL1HS<-merge(UCECcut, UCECwL1HS, all=TRUE, by="row.names")
rownames(UCECaL1HS)=UCECaL1HS[,1]
UCECaL1HS<- data.frame(UCECaL1HS[,-c(1,3,4)])
colnames(UCECaL1HS)=c("gene", "L1HS")
#Graph!!
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")], color="BLCA")) + geom_smooth(aes(BLCAaL1HS[,c("gene")], BLCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="burlywood1") +
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="cadetblue1") +
  geom_jitter(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")], color="CHOL")) + geom_smooth(aes(CHOLaL1HS[,c("gene")], CHOLaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chartreuse1") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")], color="ESCA")) + geom_smooth(aes(ESCAaL1HS[,c("gene")], ESCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="coral1") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")], color="KICH")) + geom_smooth(aes(KICHaL1HS[,c("gene")], KICHaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkgoldenrod3") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")], color="LIHC")) + geom_smooth(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkorange3") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("gene")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
  geom_jitter(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")], color="STAD")) + geom_smooth(aes(STADaL1HS[,c("gene")], STADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkred") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  geom_jitter(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")], color="UCEC")) + geom_smooth(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  labs(x = "Tigger2:3ATcMar-Tigger:3ADNA Expression (rpm)", y = "L1HS Expression (rpm)", title="Tigger2:3ATcMar-Tigger:3ADNA Correlation")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("burlywood1","cadetblue1", 
		"chartreuse1", "chocolate1", "coral1", "darkblue", "darkgoldenrod3", 
		"darkkhaki", "black", "darkorange3", "lightsalmon3", "darkslategray",
		"deeppink2", "darkred", "brown3", "springgreen3")) 

