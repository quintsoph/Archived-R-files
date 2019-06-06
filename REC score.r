Allraw <- subset(AllCancers, V1=="BRCA")
AllcancerL1HS<- L1HScancert[,c(Allraw[,2])]
AllcancerL1HS<- t(AllcancerL1HS)
AllnormalL1HS<- L1HSnormalt[,c(Allraw[,2])]
AllnormalL1HS<- t(AllnormalL1HS)

#BRCA MiRNA
Allcancermirna<- cancermirna[,c(Allraw[,2])]
Allcancermirna<- t(Allcancermirna)
Allnormalmirna<- framenormalmirna[,c(Allraw[,2])]
Allnormalmirna<- t(Allnormalmirna)

#All log2ratios
#find the dividen and log
DividensofL1HSAll<- AllcancerL1HS/AllnormalL1HS
logofL1HSAll <- log2(DividensofL1HSAll)
logofL1HSAll <- t(logofL1HSAll)
#find the dividensofmirna and log
DividensofmirnaAll<- Allcancermirna/Allnormalmirna
logofmirnaAll<- log2(DividensofmirnaAll)
#Change inf to NAs.
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.infinite))
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.nan))

 #L1HS and mirna
ratiopvaluesAll<- function(x)
{
	idea <- lm(logofL1HSAll[,c(1)] ~ logofmirnaAll[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaAll <- !apply (is.na(logofmirnaAll), 2, all)
framefalsenaAll<- data.frame(falsenaAll)
NONAAll <- subset(framefalsenaAll, falsenaAll==TRUE)
sublogmirnaAll<- logofmirnaAll[,c(rownames(NONAAll))]

#iterate
presultsAll=NULL
for (i in colnames(sublogmirnaAll))
{
	q = NULL
	
		q <-ratiopvaluesAll(i)

	presultsAll <- cbind(presultsAll, q) 
}
colnames(presultsAll)=colnames(sublogmirnaAll)
tpresults<- t(presultsAll)
tpresultsBRCA<- tpresults

#BLCA
Allraw <- subset(AllCancers, V1=="BLCA")
AllcancerL1HS<- L1HScancert[,c(Allraw[,2])]
AllcancerL1HS<- t(AllcancerL1HS)
AllnormalL1HS<- L1HSnormalt[,c(Allraw[,2])]
AllnormalL1HS<- t(AllnormalL1HS)

#BRCA MiRNA
Allcancermirna<- cancermirna[,c(Allraw[,2])]
Allcancermirna<- t(Allcancermirna)
Allnormalmirna<- framenormalmirna[,c(Allraw[,2])]
Allnormalmirna<- t(Allnormalmirna)

#All log2ratios
#find the dividen and log
DividensofL1HSAll<- AllcancerL1HS/AllnormalL1HS
logofL1HSAll <- log2(DividensofL1HSAll)
logofL1HSAll <- t(logofL1HSAll)
#find the dividensofmirna and log
DividensofmirnaAll<- Allcancermirna/Allnormalmirna
logofmirnaAll<- log2(DividensofmirnaAll)
#Change inf to NAs.
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.infinite))
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.nan))

 #L1HS and mirna
ratiopvaluesAll<- function(x)
{
	idea <- lm(logofL1HSAll[,c(1)] ~ logofmirnaAll[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaAll <- !apply (is.na(logofmirnaAll), 2, all)
framefalsenaAll<- data.frame(falsenaAll)
NONAAll <- subset(framefalsenaAll, falsenaAll==TRUE)
sublogmirnaAll<- logofmirnaAll[,c(rownames(NONAAll))]

#iterate
presultsAll=NULL
for (i in colnames(sublogmirnaAll))
{
	q = NULL
	
		q <-ratiopvaluesAll(i)

	presultsAll <- cbind(presultsAll, q) 
}
colnames(presultsAll)=colnames(sublogmirnaAll)
tpresults<- t(presultsAll)
tpresultsBLCA<- tpresults

Allraw <- subset(AllCancers, V1=="CESC")
AllcancerL1HS<- L1HScancert[,c(Allraw[,2])]
AllcancerL1HS<- t(AllcancerL1HS)
AllnormalL1HS<- L1HSnormalt[,c(Allraw[,2])]
AllnormalL1HS<- t(AllnormalL1HS)

#MiRNA
Allcancermirna<- cancermirna[,c(Allraw[,2])]
Allcancermirna<- t(Allcancermirna)
Allnormalmirna<- framenormalmirna[,c(Allraw[,2])]
Allnormalmirna<- t(Allnormalmirna)

#All log2ratios
#find the dividen and log
DividensofL1HSAll<- AllcancerL1HS/AllnormalL1HS
logofL1HSAll <- log2(DividensofL1HSAll)
logofL1HSAll <- t(logofL1HSAll)
#find the dividensofmirna and log
DividensofmirnaAll<- Allcancermirna/Allnormalmirna
logofmirnaAll<- log2(DividensofmirnaAll)
#Change inf to NAs.
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.infinite))
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.nan))

 #L1HS and mirna
ratiopvaluesAll<- function(x)
{
	idea <- lm(logofL1HSAll[,c(1)] ~ logofmirnaAll[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaAll <- !apply (is.na(logofmirnaAll), 2, all)
framefalsenaAll<- data.frame(falsenaAll)
NONAAll <- subset(framefalsenaAll, falsenaAll==TRUE)
sublogmirnaAll<- logofmirnaAll[,c(rownames(NONAAll))]

#iterate
presultsAll=NULL
for (i in colnames(sublogmirnaAll))
{
	q = NULL
	
		q <-ratiopvaluesAll(i)

	presultsAll <- cbind(presultsAll, q) 
}
colnames(presultsAll)=colnames(sublogmirnaAll)
tpresults<- t(presultsAll)
tpresultsCESC<- tpresults

#CHOL
Allraw <- subset(AllCancers, V1=="CHOL")
AllcancerL1HS<- L1HScancert[,c(Allraw[,2])]
AllcancerL1HS<- t(AllcancerL1HS)
AllnormalL1HS<- L1HSnormalt[,c(Allraw[,2])]
AllnormalL1HS<- t(AllnormalL1HS)

#BRCA MiRNA
Allcancermirna<- cancermirna[,c(Allraw[,2])]
Allcancermirna<- t(Allcancermirna)
Allnormalmirna<- framenormalmirna[,c(Allraw[,2])]
Allnormalmirna<- t(Allnormalmirna)

#All log2ratios
#find the dividen and log
DividensofL1HSAll<- AllcancerL1HS/AllnormalL1HS
logofL1HSAll <- log2(DividensofL1HSAll)
logofL1HSAll <- t(logofL1HSAll)
#find the dividensofmirna and log
DividensofmirnaAll<- Allcancermirna/Allnormalmirna
logofmirnaAll<- log2(DividensofmirnaAll)
#Change inf to NAs.
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.infinite))
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.nan))

 #L1HS and mirna
ratiopvaluesAll<- function(x)
{
	idea <- lm(logofL1HSAll[,c(1)] ~ logofmirnaAll[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaAll <- !apply (is.na(logofmirnaAll), 2, all)
framefalsenaAll<- data.frame(falsenaAll)
NONAAll <- subset(framefalsenaAll, falsenaAll==TRUE)
sublogmirnaAll<- logofmirnaAll[,c(rownames(NONAAll))]

#iterate
presultsAll=NULL
for (i in colnames(sublogmirnaAll))
{
	q = NULL
	
		q <-ratiopvaluesAll(i)

	presultsAll <- cbind(presultsAll, q) 
}
colnames(presultsAll)=colnames(sublogmirnaAll)
tpresults<- t(presultsAll)
tpresultsCHOL<- tpresults

#ESCA
Allraw <- subset(AllCancers, V1=="ESCA")
AllcancerL1HS<- L1HScancert[,c(Allraw[,2])]
AllcancerL1HS<- t(AllcancerL1HS)
AllnormalL1HS<- L1HSnormalt[,c(Allraw[,2])]
AllnormalL1HS<- t(AllnormalL1HS)

#BRCA MiRNA
Allcancermirna<- cancermirna[,c(Allraw[,2])]
Allcancermirna<- t(Allcancermirna)
Allnormalmirna<- framenormalmirna[,c(Allraw[,2])]
Allnormalmirna<- t(Allnormalmirna)

#All log2ratios
#find the dividen and log
DividensofL1HSAll<- AllcancerL1HS/AllnormalL1HS
logofL1HSAll <- log2(DividensofL1HSAll)
logofL1HSAll <- t(logofL1HSAll)
#find the dividensofmirna and log
DividensofmirnaAll<- Allcancermirna/Allnormalmirna
logofmirnaAll<- log2(DividensofmirnaAll)
#Change inf to NAs.
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.infinite))
 is.na(logofmirnaAll) <- do.call(cbind,lapply(logofmirnaAll, is.nan))

 #L1HS and mirna
ratiopvaluesAll<- function(x)
{
	idea <- lm(logofL1HSAll[,c(1)] ~ logofmirnaAll[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsenaAll <- !apply (is.na(logofmirnaAll), 2, all)
framefalsenaAll<- data.frame(falsenaAll)
NONAAll <- subset(framefalsenaAll, falsenaAll==TRUE)
sublogmirnaAll<- logofmirnaAll[,c(rownames(NONAAll))]

#iterate
presultsAll=NULL
for (i in colnames(sublogmirnaAll))
{
	q = NULL
	
		q <-ratiopvaluesAll(i)

	presultsAll <- cbind(presultsAll, q) 
}
colnames(presultsAll)=colnames(sublogmirnaAll)
tpresults<- t(presultsAll)
tpresultsESCA<- tpresults