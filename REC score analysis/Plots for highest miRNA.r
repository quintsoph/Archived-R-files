###Plot highest miRNA
### First three negative 
# hsa-mir-452
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir452<- read.table("hsa-mir-452.normal.txt")
BLCA452<- subset(mir452, mir452[,2] == "BLCA", select=c("V1", "V3"))
BLCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BLCA452[,1]),]
rownames(BLCA452)=BLCA452[,1]
rownames(BLCAwL1HS)=BLCAwL1HS[,1]
BLCAaL1HS<-merge(BLCA452, BLCAwL1HS, all=TRUE, by="row.names")
rownames(BLCAaL1HS)=BLCAaL1HS[,1]
BLCAaL1HS<- BLCAaL1HS[,-c(1, 2, 4)]
BLCAaL1HS<- log(BLCAaL1HS, base=2)
colnames(BLCAaL1HS)=c("miRNA", "L1HS")

BRCA452<- subset(mir452, mir452[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA452[,1]),]
rownames(BRCA452)=BRCA452[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA452, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

CESC452<- subset(mir452, mir452[,2] == "CESC", select=c("V1", "V3"))
CESCwL1HS<- L1HSread[which(L1HSread[,1] %in% CESC452[,1]),]
rownames(CESC452)=CESC452[,1]
rownames(CESCwL1HS)=CESCwL1HS[,1]
CESCaL1HS<-merge(CESC452, CESCwL1HS, all=TRUE, by="row.names")
rownames(CESCaL1HS)=CESCaL1HS[,1]
CESCaL1HS<- CESCaL1HS[,-c(1, 2, 4)]
CESCaL1HS<- log(CESCaL1HS, base=2)
colnames(CESCaL1HS)=c("miRNA", "L1HS")

CHOL452<- subset(mir452, mir452[,2] == "CHOL", select=c("V1", "V3"))
CHOLwL1HS<- L1HSread[which(L1HSread[,1] %in% CHOL452[,1]),]
rownames(CHOL452)=CHOL452[,1]
rownames(CHOLwL1HS)=CHOLwL1HS[,1]
CHOLaL1HS<-merge(CHOL452, CHOLwL1HS, all=TRUE, by="row.names")
rownames(CHOLaL1HS)= CHOLaL1HS[,1]
CHOLaL1HS<- CHOLaL1HS[,-c(1, 2, 4)]
CHOLaL1HS<- log(CHOLaL1HS, base=2)
colnames(CHOLaL1HS)=c("miRNA", "L1HS")

COAD452<- subset(mir452, mir452[,2] == "COAD", select=c("V1", "V3"))
COADwL1HS<- L1HSread[which(L1HSread[,1] %in% COAD452[,1]),]
rownames(COAD452)=COAD452[,1]
rownames(COADwL1HS)=COADwL1HS[,1]
COADaL1HS<-merge(COAD452, COADwL1HS, all=TRUE, by="row.names")
rownames(COADaL1HS)=COADaL1HS[,1]
COADaL1HS<- COADaL1HS[,-c(1, 2, 4)]
COADaL1HS<- log(COADaL1HS, base=2)
colnames(COADaL1HS)=c("miRNA", "L1HS")

ESCA452<- subset(mir452, mir452[,2] == "ESCA", select=c("V1", "V3"))
ESCAwL1HS<- L1HSread[which(L1HSread[,1] %in% ESCA452[,1]),]
rownames(ESCA452)=ESCA452[,1]
rownames(ESCAwL1HS)=ESCAwL1HS[,1]
ESCAaL1HS<-merge(ESCA452, ESCAwL1HS, all=TRUE, by="row.names")
rownames(ESCAaL1HS)=ESCAaL1HS[,1]
ESCAaL1HS<- ESCAaL1HS[,-c(1, 2, 4)]
ESCAaL1HS<- log(ESCAaL1HS, base=2)
colnames(ESCAaL1HS)=c("miRNA", "L1HS")

HNSC452<- subset(mir452, mir452[,2] == "HNSC", select=c("V1", "V3"))
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSC452[,1]),]
rownames(HNSC452)=HNSC452[,1]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
HNSCaL1HS<-merge(HNSC452, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)= HNSCaL1HS[,1]
HNSCaL1HS<- HNSCaL1HS[,-c(1, 2, 4)]
HNSCaL1HS<- log(HNSCaL1HS, base=2)
colnames(HNSCaL1HS)=c("miRNA", "L1HS")

KICH452<- subset(mir452, mir452[,2] == "KICH", select=c("V1", "V3"))
KICHwL1HS<- L1HSread[which(L1HSread[,1] %in% KICH452[,1]),]
rownames(KICH452)=KICH452[,1]
rownames(KICHwL1HS)=KICHwL1HS[,1]
KICHaL1HS<-merge(KICH452, KICHwL1HS, all=TRUE, by="row.names")
rownames(KICHaL1HS)=KICHaL1HS[,1]
KICHaL1HS<- KICHaL1HS[,-c(1, 2, 4)]
KICHaL1HS<- log(KICHaL1HS, base=2)
colnames(KICHaL1HS)=c("miRNA", "L1HS")

KIRC452<- subset(mir452, mir452[,2] == "KIRC", select=c("V1", "V3"))
KIRCwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRC452[,1]),]
rownames(KIRC452)=KIRC452[,1]
rownames(KIRCwL1HS)=KIRCwL1HS[,1]
KIRCaL1HS<-merge(KIRC452, KIRCwL1HS, all=TRUE, by="row.names")
rownames(KIRCaL1HS)=KIRCaL1HS[,1]
KIRCaL1HS<- KIRCaL1HS[,-c(1, 2, 4)]
KIRCaL1HS<- log(KIRCaL1HS, base=2)
colnames(KIRCaL1HS)=c("miRNA", "L1HS")

KIRP452<- subset(mir452, mir452[,2] == "KIRP", select=c("V1", "V3"))
KIRPwL1HS<- L1HSread[which(L1HSread[,1] %in% KIRP452[,1]),]
rownames(KIRP452)=KIRP452[,1]
rownames(KIRPwL1HS)=KIRPwL1HS[,1]
KIRPaL1HS<-merge(KIRP452, KIRPwL1HS, all=TRUE, by="row.names")
rownames(KIRPaL1HS)= KIRPaL1HS[,1]
KIRPaL1HS<- KIRPaL1HS[,-c(1, 2, 4)]
KIRPaL1HS<- log(KIRPaL1HS, base=2)
colnames(KIRPaL1HS)=c("miRNA", "L1HS")

LIHC452<- subset(mir452, mir452[,2] == "LIHC", select=c("V1", "V3"))
LIHCwL1HS<- L1HSread[which(L1HSread[,1] %in% LIHC452[,1]),]
rownames(LIHC452)=LIHC452[,1]
rownames(LIHCwL1HS)=LIHCwL1HS[,1]
LIHCaL1HS<-merge(LIHC452, LIHCwL1HS, all=TRUE, by="row.names")
rownames(LIHCaL1HS)= LIHCaL1HS[,1]
LIHCaL1HS<- LIHCaL1HS[,-c(1, 2, 4)]
LIHCaL1HS<- log(LIHCaL1HS, base=2)
colnames(LIHCaL1HS)=c("miRNA", "L1HS")

LUAD452<- subset(mir452, mir452[,2] == "LUAD", select=c("V1", "V3"))
LUADwL1HS<- L1HSread[which(L1HSread[,1] %in% LUAD452[,1]),]
rownames(LUAD452)=LUAD452[,1]
rownames(LUADwL1HS)=LUADwL1HS[,1]
LUADaL1HS<-merge(LUAD452, LUADwL1HS, all=TRUE, by="row.names")
rownames(LUADaL1HS)= LUADaL1HS[,1]
LUADaL1HS<- LUADaL1HS[,-c(1, 2, 4)]
LUADaL1HS<- log(LUADaL1HS, base=2)
colnames(LUADaL1HS)=c("miRNA", "L1HS")

LUSC452<- subset(mir452, mir452[,2] == "LUSC", select=c("V1", "V3"))
LUSCwL1HS<- L1HSread[which(L1HSread[,1] %in% LUSC452[,1]),]
rownames(LUSC452)=LUSC452[,1]
rownames(LUSCwL1HS)=LUSCwL1HS[,1]
LUSCaL1HS<-merge(LUSC452, LUSCwL1HS, all=TRUE, by="row.names")
rownames(LUSCaL1HS)= LUSCaL1HS[,1]
LUSCaL1HS<- LUSCaL1HS[,-c(1, 2, 4)]
LUSCaL1HS<- log(LUSCaL1HS, base=2)
colnames(LUSCaL1HS)=c("miRNA", "L1HS")

PRAD452<- subset(mir452, mir452[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD452[,1]),]
rownames(PRAD452)=PRAD452[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD452, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

THCA452<- subset(mir452, mir452[,2] == "THCA", select=c("V1", "V3"))
THCAwL1HS<- L1HSread[which(L1HSread[,1] %in% THCA452[,1]),]
rownames(THCA452)=THCA452[,1]
rownames(THCAwL1HS)=THCAwL1HS[,1]
THCAaL1HS<-merge(THCA452, THCAwL1HS, all=TRUE, by="row.names")
rownames(THCAaL1HS)=THCAaL1HS[,1]
THCAaL1HS<- THCAaL1HS[,-c(1, 2, 4)]
THCAaL1HS<- log(THCAaL1HS, base=2)
colnames(THCAaL1HS)=c("miRNA", "L1HS")

UCEC452<- subset(mir452, mir452[,2] == "UCEC", select=c("V1", "V3"))
UCECwL1HS<- L1HSread[which(L1HSread[,1] %in% UCEC452[,1]),]
rownames(UCEC452)=UCEC452[,1]
rownames(UCECwL1HS)=UCECwL1HS[,1]
UCECaL1HS<-merge(UCEC452, UCECwL1HS, all=TRUE, by="row.names")
rownames(UCECaL1HS)=UCECaL1HS[,1]
UCECaL1HS<- UCECaL1HS[,-c(1, 2, 4)]
UCECaL1HS<- log(UCECaL1HS, base=2)
colnames(UCECaL1HS)=c("miRNA", "L1HS")

ggplot(NULL) + 
  geom_jitter(aes(BLCAaL1HS[,2], BLCAaL1HS[,3], color="BLCA")) + geom_smooth(aes(BLCAaL1HS[,2],BLCAaL1HS[,3]), method=lm, se=FALSE, color=#00FF99) +
  geom_jitter(aes(BRCAaL1HS[,2], BLCAaL1HS[,3], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,2], BRCAaL1HS[,3]), method=lm, se=FALSE, color=#FF00FF) +
  geom_jitter(aes(CESCaL1HS[,2], CESCaL1HS[,3], color="CESC")) + geom_smooth(aes(CESCaL1HS[,2], CESCaL1HS[,3]), method=lm, se=FALSE, color=#FF0066) +
  labs(x = "log2 normal of hsa-mir-452", y = "L1HS Expression(RPM)", title="hsa-mir-452")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c(#00FF99, #FF00FF, #FF0066)) 

###########START HERE
#I needed to just do the ones in the duplicate files
#minus files
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir3651<- read.table("hsa-mir-365-1.normal.txt")

BRCA3651<- subset(mir3651, mir3651[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA3651[,1]),]
rownames(BRCA3651)=BRCA3651[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA3651, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD3651<- subset(mir3651, mir3651[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD3651[,1]),]
rownames(PRAD3651)=PRAD3651[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD3651, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-365-1", y = "L1HS Expression(RPM)", title="hsa-mir-365-1")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-365-1 normal minus 2.png", width=9.27, height=10.2)

#mir-365-2
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir3652<- read.table("hsa-mir-365-2.normal.txt")

BRCA3652<- subset(mir3652, mir3652[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA3652[,1]),]
rownames(BRCA3652)=BRCA3652[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA3652, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD3652<- subset(mir3652, mir3652[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD3652[,1]),]
rownames(PRAD3652)=PRAD3652[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD3652, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-365-2", y = "L1HS Expression(RPM)", title="hsa-mir-365-2")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-365-2 normal minus 2.png", width=9.27, height=10.2)

#hsa-mir-378
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir378<- read.table("hsa-mir-378.normal.txt")

BRCA378<- subset(mir378, mir378[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA378[,1]),]
rownames(BRCA378)=BRCA378[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA378, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD378<- subset(mir378, mir378[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD378[,1]),]
rownames(PRAD378)=PRAD378[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD378, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-378", y = "L1HS Expression(RPM)", title="hsa-mir-378")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-378 normal minus 2.png", width=9.27, height=10.2)

#hsa-mir-484
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir484<- read.table("hsa-mir-484.normal.txt")

BRCA484<- subset(mir484, mir484[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA484[,1]),]
rownames(BRCA484)=BRCA484[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA484, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD484<- subset(mir484, mir484[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD484[,1]),]
rownames(PRAD484)=PRAD484[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD484, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-484", y = "L1HS Expression(RPM)", title="hsa-mir-484")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-484 normal minus 2.png", width=9.27, height=10.2)

#hsa-mir-486
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir486<- read.table("hsa-mir-486.normal.txt")

BRCA486<- subset(mir486, mir486[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA486[,1]),]
rownames(BRCA486)=BRCA486[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA486, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD486<- subset(mir486, mir486[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD486[,1]),]
rownames(PRAD486)=PRAD486[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD486, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-486", y = "L1HS Expression(RPM)", title="hsa-mir-486")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-486 normal minus 2.png", width=9.27, height=10.2)

#hsa-mir-652
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir652<- read.table("hsa-mir-652.normal.txt")

BRCA652<- subset(mir652, mir652[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA652[,1]),]
rownames(BRCA652)=BRCA652[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA652, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD652<- subset(mir652, mir652[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD652[,1]),]
rownames(PRAD652)=PRAD652[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD652, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-652", y = "L1HS Expression(RPM)", title="hsa-mir-652")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-652 normal minus 2.png", width=9.27, height=10.2)

##Plus graphs
#hsa-mir-141
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir141<- read.table("hsa-mir-141.normal.txt")

BRCA141<- subset(mir141, mir141[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA141[,1]),]
rownames(BRCA141)=BRCA141[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA141, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

HNSC141<- subset(mir141, mir141[,2] == "HNSC", select=c("V1", "V3"))
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSC141[,1]),]
rownames(HNSC141)=HNSC141[,1]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
HNSCaL1HS<-merge(HNSC141, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)= HNSCaL1HS[,1]
HNSCaL1HS<- HNSCaL1HS[,-c(1, 2, 4)]
HNSCaL1HS<- log(HNSCaL1HS, base=2)
colnames(HNSCaL1HS)=c("miRNA", "L1HS")

combined<- rbind(BRCAaL1HS, HNSCaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-141", y = "L1HS Expression(RPM)", title="hsa-mir-141")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "springgreen3")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-141 normal plus 2.png", width=9.27, height=10.2)

#hsa-mir-200a
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir200a<- read.table("hsa-mir-200a.normal.txt")

BRCA200a<- subset(mir200a, mir200a[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA200a[,1]),]
rownames(BRCA200a)=BRCA200a[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA200a, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

HNSC200a<- subset(mir200a, mir200a[,2] == "HNSC", select=c("V1", "V3"))
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSC200a[,1]),]
rownames(HNSC200a)=HNSC200a[,1]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
HNSCaL1HS<-merge(HNSC200a, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)= HNSCaL1HS[,1]
HNSCaL1HS<- HNSCaL1HS[,-c(1, 2, 4)]
HNSCaL1HS<- log(HNSCaL1HS, base=2)
colnames(HNSCaL1HS)=c("miRNA", "L1HS")

combined<- rbind(BRCAaL1HS, HNSCaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-200a", y = "L1HS Expression(RPM)", title="hsa-mir-200a")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "springgreen3")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-200a normal plus 2.png", width=9.27, height=10.2)

#hsa-mir-200c
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir200c<- read.table("hsa-mir-200c.normal.txt")

BRCA200c<- subset(mir200c, mir200c[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA200c[,1]),]
rownames(BRCA200c)=BRCA200c[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA200c, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

HNSC200c<- subset(mir200c, mir200c[,2] == "HNSC", select=c("V1", "V3"))
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSC200c[,1]),]
rownames(HNSC200c)=HNSC200c[,1]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
HNSCaL1HS<-merge(HNSC200c, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)= HNSCaL1HS[,1]
HNSCaL1HS<- HNSCaL1HS[,-c(1, 2, 4)]
HNSCaL1HS<- log(HNSCaL1HS, base=2)
colnames(HNSCaL1HS)=c("miRNA", "L1HS")

combined<- rbind(BRCAaL1HS, HNSCaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-200c", y = "L1HS Expression(RPM)", title="hsa-mir-200c")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "springgreen3")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-200c normal plus 2.png", width=9.27, height=10.2)

#hsa-mir-203
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir203<- read.table("hsa-mir-203.normal.txt")

BRCA203<- subset(mir203, mir203[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA203[,1]),]
rownames(BRCA203)=BRCA203[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA203, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

HNSC203<- subset(mir203, mir203[,2] == "HNSC", select=c("V1", "V3"))
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSC203[,1]),]
rownames(HNSC203)=HNSC203[,1]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
HNSCaL1HS<-merge(HNSC203, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)= HNSCaL1HS[,1]
HNSCaL1HS<- HNSCaL1HS[,-c(1, 2, 4)]
HNSCaL1HS<- log(HNSCaL1HS, base=2)
colnames(HNSCaL1HS)=c("miRNA", "L1HS")

combined<- rbind(BRCAaL1HS, HNSCaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-203", y = "L1HS Expression(RPM)", title="hsa-mir-203")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "springgreen3")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-203 normal plus 2.png", width=9.27, height=10.2)

#hsa-mir-205
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir205<- read.table("hsa-mir-205.normal.txt")

BRCA205<- subset(mir205, mir205[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA205[,1]),]
rownames(BRCA205)=BRCA205[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA205, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

HNSC205<- subset(mir205, mir205[,2] == "HNSC", select=c("V1", "V3"))
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSC205[,1]),]
rownames(HNSC205)=HNSC205[,1]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
HNSCaL1HS<-merge(HNSC205, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)= HNSCaL1HS[,1]
HNSCaL1HS<- HNSCaL1HS[,-c(1, 2, 4)]
HNSCaL1HS<- log(HNSCaL1HS, base=2)
colnames(HNSCaL1HS)=c("miRNA", "L1HS")

combined<- rbind(BRCAaL1HS, HNSCaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-205", y = "L1HS Expression(RPM)", title="hsa-mir-205")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "springgreen3")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-205 normal plus 2.png", width=9.27, height=10.2)

#hsa-mir-375
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir375<- read.table("hsa-mir-375.normal.txt")

BRCA375<- subset(mir375, mir375[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA375[,1]),]
rownames(BRCA375)=BRCA375[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA375, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD375<- subset(mir375, mir375[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD375[,1]),]
rownames(PRAD375)=PRAD375[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD375, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  labs(x = "log2 normal of hsa-mir-375", y = "L1HS Expression(RPM)", title="hsa-mir-375")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-375 normal plus 2.png", width=9.27, height=10.2)

#hsa-mir-3065
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs")
L1HSread<- read.table("L1HSnormal.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/normal miRNA")
mir3065<- read.table("hsa-mir-3065.normal.txt")

BRCA3065<- subset(mir3065, mir3065[,2] == "BRCA", select=c("V1", "V3"))
BRCAwL1HS<- L1HSread[which(L1HSread[,1] %in% BRCA3065[,1]),]
rownames(BRCA3065)=BRCA3065[,1]
rownames(BRCAwL1HS)=BRCAwL1HS[,1]
BRCAaL1HS<-merge(BRCA3065, BRCAwL1HS, all=TRUE, by="row.names")
rownames(BRCAaL1HS)=BRCAaL1HS[,1]
BRCAaL1HS<- BRCAaL1HS[,-c(1, 2, 4)]
BRCAaL1HS<- log(BRCAaL1HS, base=2)
colnames(BRCAaL1HS)=c("miRNA", "L1HS")

PRAD3065<- subset(mir3065, mir3065[,2] == "PRAD", select=c("V1", "V3"))
PRADwL1HS<- L1HSread[which(L1HSread[,1] %in% PRAD3065[,1]),]
rownames(PRAD3065)=PRAD3065[,1]
rownames(PRADwL1HS)=PRADwL1HS[,1]
PRADaL1HS<-merge(PRAD3065, PRADwL1HS, all=TRUE, by="row.names")
rownames(PRADaL1HS)= PRADaL1HS[,1]
PRADaL1HS<- PRADaL1HS[,-c(1, 2, 4)]
PRADaL1HS<- log(PRADaL1HS, base=2)
colnames(PRADaL1HS)=c("miRNA", "L1HS")

HNSC3065<- subset(mir3065, mir3065[,2] == "HNSC", select=c("V1", "V3"))
HNSCwL1HS<- L1HSread[which(L1HSread[,1] %in% HNSC3065[,1]),]
rownames(HNSC3065)=HNSC3065[,1]
rownames(HNSCwL1HS)=HNSCwL1HS[,1]
HNSCaL1HS<-merge(HNSC3065, HNSCwL1HS, all=TRUE, by="row.names")
rownames(HNSCaL1HS)= HNSCaL1HS[,1]
HNSCaL1HS<- HNSCaL1HS[,-c(1, 2, 4)]
HNSCaL1HS<- log(HNSCaL1HS, base=2)
colnames(HNSCaL1HS)=c("miRNA", "L1HS")

combined<- rbind(PRADaL1HS, BRCAaL1HS, HNSCaL1HS)

library("ggplot2")
legend_title<- "Cancer Types"
ggplot(NULL) + 
  geom_jitter(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")], color="PRAD")) + geom_smooth(aes(PRADaL1HS[,c("miRNA")], PRADaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="blueviolet") +
  geom_jitter(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")], color="BRCA")) + geom_smooth(aes(BRCAaL1HS[,c("miRNA")], BRCAaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="deeppink3") +
  geom_jitter(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")], color="HNSC")) + geom_smooth(aes(HNSCaL1HS[,c("miRNA")], HNSCaL1HS[,c("L1HS")]), method=lm, se=FALSE, color="springgreen3") +
  labs(x = "log2 normal of hsa-mir-3065", y = "L1HS Expression(RPM)", title="hsa-mir-3065")+
  theme_bw()+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("deeppink3", "springgreen3", "blueviolet")) 

setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/Graphs/Result Graphs")
ggsave("hsa-mir-3065 normal plus 2.png", width=9.27, height=10.2)
