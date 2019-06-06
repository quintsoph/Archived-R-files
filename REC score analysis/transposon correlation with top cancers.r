####The transposon correlation with top cancers
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/Graph data for top genes")
L1HSread<- read.table("L1HS.VST.normal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 3/Graph data for top genes")
patients<- read.table("patient_type.txt", header=TRUE)

#####L1P
L1P<- read.table("L1P%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
L1P<- cbind(L1P, rownames(L1P))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- L1P[which(row.names(L1P) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="burlywood1") +
  labs(x = "L1P:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1HS and L1P:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("burlywood1"))

#####L1PA2
#BLCA
A2<- read.table("L1PA2%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A2<- cbind(A2, rownames(A2))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A2[which(row.names(A2) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A2[which(row.names(A2) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- A2[which(row.names(A2) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A2[which(row.names(A2) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LIHC
LIHC<- subset(patients, patients[,2]=="LIHC")
LIHCcut<- A2[which(row.names(A2) %in% LIHC[,1]),]
LIHCwL1HS<- L1HSread[which(L1HSread[,1] %in% LIHCcut[,2]),]
rownames(LIHCwL1HS)=LIHCwL1HS[,1]
colnames(LIHCwL1HS)=c("patient", "L1HS")
LIHCaL1HS<-merge(LIHCcut, LIHCwL1HS, all=TRUE, by="row.names")
rownames(LIHCaL1HS)=LIHCaL1HS[,1]
LIHCaL1HS<- data.frame(LIHCaL1HS[,-c(1,3,4)])
colnames(LIHCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A2[which(row.names(A2) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink2") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")], color="LIHC")) + geom_smooth(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkorange3") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA2:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA2:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"deeppink2", "darkblue","darkkhaki", "darkorange3",
		"black")) 

#####L1PA3
A3<- read.table("L1PA3%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A3<- cbind(A3, rownames(A3))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A3[which(row.names(A3) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A3[which(row.names(A3) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- A3[which(row.names(A3) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A3[which(row.names(A3) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- A3[which(row.names(A3) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A3[which(row.names(A3) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A3[which(row.names(A3) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A3[which(row.names(A3) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
#UCEC
UCEC<- subset(patients, patients[,2]=="UCEC")
UCECcut<- A3[which(row.names(A3) %in% UCEC[,1]),]
UCECwL1HS<- L1HSread[which(L1HSread[,1] %in% UCECcut[,2]),]
rownames(UCECwL1HS)=UCECwL1HS[,1]
colnames(UCECwL1HS)=c("patient", "L1HS")
UCECaL1HS<-merge(UCECcut, UCECwL1HS, all=TRUE, by="row.names")
rownames(UCECaL1HS)=UCECaL1HS[,1]
UCECaL1HS<- data.frame(UCECaL1HS[,-c(1,3,4)])
colnames(UCECaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  geom_jitter(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")], color="UCEC")) + geom_smooth(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  labs(x = "L1PA3:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA3:3AL1:ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkblue","darkkhaki", "black", "lightsalmon3", "darkslategray",
		"brown3", "springgreen3")) 

####L1PA4
A4<- read.table("L1PA4%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A4<- cbind(A4, rownames(A4))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A4[which(row.names(A4) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A4[which(row.names(A4) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A4[which(row.names(A4) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- A4[which(row.names(A4) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A4[which(row.names(A4) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A4[which(row.names(A4) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A4[which(row.names(A4) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
#UCEC
UCEC<- subset(patients, patients[,2]=="UCEC")
UCECcut<- A4[which(row.names(A4) %in% UCEC[,1]),]
UCECwL1HS<- L1HSread[which(L1HSread[,1] %in% UCECcut[,2]),]
rownames(UCECwL1HS)=UCECwL1HS[,1]
colnames(UCECwL1HS)=c("patient", "L1HS")
UCECaL1HS<-merge(UCECcut, UCECwL1HS, all=TRUE, by="row.names")
rownames(UCECaL1HS)=UCECaL1HS[,1]
UCECaL1HS<- data.frame(UCECaL1HS[,-c(1,3,4)])
colnames(UCECaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  geom_jitter(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")], color="UCEC")) + geom_smooth(aes(UCECaL1HS[,c("gene")], UCECaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  labs(x = "L1PA4:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA4:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "black", "lightsalmon3", "darkslategray",
		"brown3", "springgreen3")) 

####L1PA5
A5<- read.table("L1PA5%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A5<- cbind(A5, rownames(A5))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A5[which(row.names(A5) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A5[which(row.names(A5) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A5[which(row.names(A5) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A5[which(row.names(A5) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A5[which(row.names(A5) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A5[which(row.names(A5) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA5:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA5:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "lightsalmon3", "darkslategray",
		"brown3")) 

#######L1PA6
A6<- read.table("L1PA6%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A6<- cbind(A6, rownames(A6))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A6[which(row.names(A6) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A6[which(row.names(A6) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A6[which(row.names(A6) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A6[which(row.names(A6) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A6[which(row.names(A6) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A6[which(row.names(A6) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA6:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA6:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "lightsalmon3", "darkslategray",
		"brown3")) 

#######L1PA7
A7<- read.table("L1PA7%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A7<- cbind(A7, rownames(A7))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A7[which(row.names(A7) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A7[which(row.names(A7) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- A7[which(row.names(A7) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A7[which(row.names(A7) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A7[which(row.names(A7) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A7[which(row.names(A7) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA7:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA7:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "black", "lightsalmon3", "darkslategray",
		"brown3")) 

#####L1PA8
A8<- read.table("L1PA8%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A8<- cbind(A8, rownames(A8))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A8[which(row.names(A8) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A8[which(row.names(A8) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A8[which(row.names(A8) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A8[which(row.names(A8) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A8[which(row.names(A8) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  labs(x = "L1PA8:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA8:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "lightsalmon3", "darkslategray"
		)) 

#####L1PA8A
A8A<- read.table("L1PA8A%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A8A<- cbind(A8A, rownames(A8A))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A8A[which(row.names(A8A) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A8A[which(row.names(A8A) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A8A[which(row.names(A8A) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A8A[which(row.names(A8A) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA8A:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA8A:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "lightsalmon3",
		"brown3")) 

####L1PA10
A10<- read.table("L1PA10%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A10<- cbind(A10, rownames(A10))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A10[which(row.names(A10) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A10[which(row.names(A10) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- A10[which(row.names(A10) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LIHC
LIHC<- subset(patients, patients[,2]=="LIHC")
LIHCcut<- A10[which(row.names(A10) %in% LIHC[,1]),]
LIHCwL1HS<- L1HSread[which(L1HSread[,1] %in% LIHCcut[,2]),]
rownames(LIHCwL1HS)=LIHCwL1HS[,1]
colnames(LIHCwL1HS)=c("patient", "L1HS")
LIHCaL1HS<-merge(LIHCcut, LIHCwL1HS, all=TRUE, by="row.names")
rownames(LIHCaL1HS)=LIHCaL1HS[,1]
LIHCaL1HS<- data.frame(LIHCaL1HS[,-c(1,3,4)])
colnames(LIHCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A10[which(row.names(A10) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A10[which(row.names(A10) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A10[which(row.names(A10) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")], color="LIHC")) + geom_smooth(aes(LIHCaL1HS[,c("gene")], LIHCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA10:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA10:3AL1:ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "black", "springgreen3", "lightsalmon3", "darkslategray",
		"brown3")) 
 
######L1PA11
A11<- read.table("L1PA11%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A11<- cbind(A11, rownames(A11))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A11[which(row.names(A11) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A11[which(row.names(A11) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A11[which(row.names(A11) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A11[which(row.names(A11) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A11[which(row.names(A11) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA11:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA11:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "lightsalmon3", "brown3"
		)) 

#####L1PA12
A12<- read.table("L1PA12%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A12<- cbind(A12, rownames(A12))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A12[which(row.names(A12) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A12[which(row.names(A12) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A12[which(row.names(A12) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A12[which(row.names(A12) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A12[which(row.names(A12) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A12[which(row.names(A12) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA12:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA12:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "lightsalmon3", "darkslategray", "brown3"
		)) 

#####L1PA13
A13<- read.table("L1PA13%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A13<- cbind(A13, rownames(A13))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A13[which(row.names(A13) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A13[which(row.names(A13) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#KIRP
KIRP<- subset(patients, patients[,2]=="KIRP")
KIRPcut<- A13[which(row.names(A13) %in% KIRP[,1]),]
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRPcut[,2]),]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
colnames(KIRPwL1HS)=c("patient", "L1HS")
KIRPaL1HS<-merge(KIRPcut, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)=KIRPaL1HS[,1]
KIRPaL1HS<- data.frame(KIRPaL1HS[,-c(1,3,4)])
colnames(KIRPaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A13[which(row.names(A13) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A13[which(row.names(A13) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A13[which(row.names(A13) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")], color="KIRP")) + geom_smooth(aes(KIRPaL1HS[,c("gene")], KIRPaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA13:3AL1:3ALINE.normal Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA13:3AL1:3ALINE.normal Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "black", "lightsalmon3", "darkslategray",
		"brown3")) 

######L1PA14
A14<- read.table("L1PA14%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A14<- cbind(A14, rownames(A14))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A14[which(row.names(A14) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A14[which(row.names(A14) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A14[which(row.names(A14) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A14[which(row.names(A14) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A14[which(row.names(A14) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA14:3AL1:3ALINE.normal Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA14:3AL1:3ALINE.normal Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "lightsalmon3", "darkslategray",
		"brown3")) 

#####L1PA15
A15<- read.table("L1PA15%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A15<- cbind(A15, rownames(A15))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A15[which(row.names(A15) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A15[which(row.names(A15) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen" +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA15:3AL1:3ALINE.normal Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA15:3AL1:3ALINE.normal Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"brown3")) 

#####L1PA15-16
A1516<- read.table("L1PA15-16%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A1516<- cbind(A1516, rownames(A1516))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A1516[which(row.names(A1516) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A1516[which(row.names(A1516) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A1516[which(row.names(A1516) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A1516[which(row.names(A1516) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A1516[which(row.names(A1516) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA15-16:3AL1:3ALINE.normal Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA15-16:3AL1:3ALINE.normal Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "lightsalmon3", "darkslategray",
		"brown3")) 

#####L1PA16
A16<- read.table("L1PA16%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A16<- cbind(A16, rownames(A16))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A16[which(row.names(A16) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A16[which(row.names(A16) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- A16[which(row.names(A16) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A16[which(row.names(A16) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- A16[which(row.names(A16) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A16[which(row.names(A16) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightsalmon3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA16:3AL1:3ALINE.normal Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA16:3AL1:3ALINE.normal Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "lightsalmon3", "darkslategray",
		"brown3")) 

####L1PA17
A17<- read.table("L1PA17%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
A17<- cbind(A17, rownames(A17))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- A17[which(row.names(A17) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- A17[which(row.names(A17) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- A17[which(row.names(A17) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- A17[which(row.names(A17) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PA17:3AL1:3ALINE.normal Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PA17:3AL1:3ALINE.normal Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "springgreen3",
		"brown3")) 

#####L1PB
PB<- read.table("L1PB%3AL1%3ALINE.normal.txt", header=TRUE, row.names=1)
PB<- cbind(PB, rownames(PB))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- PB[which(row.names(PB) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- PB[which(row.names(PB) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- PB[which(row.names(PB) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- PB[which(row.names(PB) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "L1PB:3AL1:3ALINE Expression (rpm)", y = "L1HS Expression (rpm)", title="L1PB:3AL1:3ALINE Correlation")+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "springgreen3",
		"brown3")) 

######HERVH48
HERV<- read.table("HERVH48-int%3AERV1%3ALTR.normal.txt", header=TRUE, row.names=1)
HERV<- cbind(HERV, rownames(HERV))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- HERV[which(row.names(HERV) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- HERV[which(row.names(HERV) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- HERV[which(row.names(HERV) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- HERV[which(row.names(HERV) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- HERV[which(row.names(HERV) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- HERV[which(row.names(HERV) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "HERVH48-int:ERV1:LTR Expression (rpm)", y = "L1HS Expression (rpm)", title="HERVH48-int:ERV1:LTR Correlation")+
  xlim(5,11)+
  ylim(8,13)+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "springgreen3", "darkslategray",
		"brown3")) 

######MLT1E1A
MLT<- read.table("MLT1E1A%3AERVL-MaLR%3ALTR.normal.txt", header=TRUE, row.names=1)
MLT<- cbind(MLT, rownames(MLT))
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

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "MLT1E1A:ERVL-MaLR:LTR Expression (rpm)", y = "L1HS Expression (rpm)", title="MLT1E1A:ERVL-MaLR:LTR Correlation")+
  xlim(5,11)+
  ylim(8,13)+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "springgreen3", "darkslategray",
		"brown3")) 

######LTR9B
LT<- read.table("LTR9B%3AERV1%3ALTR.normal.txt", header=TRUE, row.names=1)
LT<- cbind(LT, rownames(LT))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- LT[which(row.names(LT) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- LT[which(row.names(LT) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- LT[which(row.names(LT) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- LT[which(row.names(LT) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- LT[which(row.names(LT) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- LT[which(row.names(LT) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- LT[which(row.names(LT) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "LTR9B:ERV1:LTR Expression (rpm)", y = "L1HS Expression (rpm)", title="LTR9B:ERV1:LTR Correlation")+
  xlim(5,11)+
  ylim(8,13)+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkblue","darkkhaki", "springgreen3", "darkslategray",
		"brown3")) 

#####MLT2A1
MLT2<- read.table("MLT2A1%3AERVL%3ALTR.normal.txt", header=TRUE, row.names=1)
MLT2<- cbind(MLT2, rownames(MLT2))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- MLT2[which(row.names(MLT2) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- MLT2[which(row.names(MLT2) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- MLT2[which(row.names(MLT2) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- MLT2[which(row.names(MLT2) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- MLT2[which(row.names(MLT2) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- MLT2[which(row.names(MLT2) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "MLT2A1:ERVL:LTR Expression (rpm)", y = "L1HS Expression (rpm)", title="MLT2A1:ERVL:LTR Correlation")+
  xlim(5,11)+
  ylim(8,13)+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"chocolate1", "darkkhaki", "springgreen3", "darkslategray",
		"brown3")) 

#######MSTB1
MST<- read.table("MSTB1-int%3AERVL-MaLR%3ALTR.normal.txt", header=TRUE, row.names=1)
MST<- cbind(MST, rownames(MST))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- MST[which(row.names(MST) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- MST[which(row.names(MST) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- MST[which(row.names(MST) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#LUSC
LUSC<- subset(patients, patients[,2]=="LUSC")
LUSCcut<- MST[which(row.names(MST) %in% LUSC[,1]),]
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSCcut[,2]),]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
colnames(LUSCwL1HS)=c("patient", "L1HS")
LUSCaL1HS<-merge(LUSCcut, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)=LUSCaL1HS[,1]
LUSCaL1HS<- data.frame(LUSCaL1HS[,-c(1,3,4)])
colnames(LUSCaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- MST[which(row.names(MST) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")], color="LUSC")) + geom_smooth(aes(LUSCaL1HS[,c("gene")], LUSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkslategray") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "MSTB1:ERVL-MaLR:LTR Expression (rpm)", y = "L1HS Expression (rpm)", title="MSTB1:ERVL-MaLR:LTR Correlation")+
  xlim(5,11)+
  ylim(8,13)+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", 
		"darkkhaki", "springgreen3", "darkslategray",
		"brown3")) 

#####ERVL
ERV<- read.table("ERVL-int%3AERVL%3ALTR.normal.txt", header=TRUE, row.names=1)
ERV<- cbind(ERV, rownames(ERV))
#BRCA
BRCA<- subset(patients, patients[,2]=="BRCA")
BRCAcut<- ERV[which(row.names(ERV) %in% BRCA[,1]),]
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCAcut[,2]),]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
colnames(BRCAwL1HS)=c("patient", "L1HS")
BRCAaL1HS<-merge(BRCAcut, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- data.frame(BRCAaL1HS[,-c(1,3,4)])
colnames(BRCAaL1HS)=c("gene", "L1HS")
#COAD
COAD<- subset(patients, patients[,2]=="COAD")
COADcut<- ERV[which(row.names(ERV) %in% COAD[,1]),]
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COADcut[,2]),]
rownames(COADwL1HS)=COADwL1HS[,1]
colnames(COADwL1HS)=c("patient", "L1HS")
COADaL1HS<-merge(COADcut, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- data.frame(COADaL1HS[,-c(1,3,4)])
colnames(COADaL1HS)=c("gene", "L1HS")
#HNSC
HNSC<- subset(patients, patients[,2]=="HNSC")
HNSCcut<- ERV[which(row.names(ERV) %in% HNSC[,1]),]
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSCcut[,2]),]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
colnames(HNSCwL1HS)=c("patient", "L1HS")
HNSCaL1HS<-merge(HNSCcut, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)=HNSCaL1HS[,1]
HNSCaL1HS<- data.frame(HNSCaL1HS[,-c(1,3,4)])
colnames(HNSCaL1HS)=c("gene", "L1HS")
#KIRC
KIRC<- subset(patients, patients[,2]=="KIRC")
KIRCcut<- ERV[which(row.names(ERV) %in% KIRC[,1]),]
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRCcut[,2]),]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
colnames(KIRCwL1HS)=c("patient", "L1HS")
KIRCaL1HS<-merge(KIRCcut, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- data.frame(KIRCaL1HS[,-c(1,3,4)])
colnames(KIRCaL1HS)=c("gene", "L1HS")
#LUAD
LUAD<- subset(patients, patients[,2]=="LUAD")
LUADcut<- ERV[which(row.names(ERV) %in% LUAD[,1]),]
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUADcut[,2]),]
rownames(LUADwL1HS)=LUADwL1HS[,1]
colnames(LUADwL1HS)=c("patient", "L1HS")
LUADaL1HS<-merge(LUADcut, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)=LUADaL1HS[,1]
LUADaL1HS<- data.frame(LUADaL1HS[,-c(1,3,4)])
colnames(LUADaL1HS)=c("gene", "L1HS")
#THCA
THCA<- subset(patients, patients[,2]=="THCA")
THCAcut<- ERV[which(row.names(ERV) %in% THCA[,1]),]
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCAcut[,2]),]
rownames(THCAwL1HS)=THCAwL1HS[,1]
colnames(THCAwL1HS)=c("patient", "L1HS")
THCAaL1HS<-merge(THCAcut, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- data.frame(THCAaL1HS[,-c(1,3,4)])
colnames(THCAaL1HS)=c("gene", "L1HS")

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("gene")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="lightseagreen") +
  geom_jitter(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")], color="COAD")) + geom_smooth(aes(COADaL1HS[,c("gene")], COADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="chocolate1") +
  geom_jitter(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("gene")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")], color="KIRC")) + geom_smooth(aes(KIRCaL1HS[,c("gene")], KIRCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="darkkhaki") +
  geom_jitter(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")], color="LUAD")) + geom_smooth(aes(LUADaL1HS[,c("gene")], LUADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")], color="THCA")) + geom_smooth(aes(THCAaL1HS[,c("gene")], THCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="brown3") +
  labs(x = "ERVL-B4-int:ERVL:LTR Expression (rpm)", y = "L1HS Expression (rpm)", title="ERVL-B4-int:ERVL:LTR Correlation")+
  xlim(5,11)+
  ylim(8,13)+
  theme_classic()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("lightseagreen", "chocolate1",
		"darkblue", "darkkhaki", "springgreen3", 
		"brown3")) 
