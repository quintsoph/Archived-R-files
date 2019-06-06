setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
NewNames<- data.frame(read.table(file="SophiapatCanTypes"))
NewNames<- t(NewNames)

#format cancer
newL1HScancer<- data.frame(read.table(file="L1HS.Cancer.BaseMeansB", header=TRUE))
newL1HScancer<- data.matrix(newL1HScancer)
tnewL1HScancer<- t(newL1HScancer)
newmirnacancer<- read.csv(file="Mirna.Cancer.csv")
newmirnacancer<- data.matrix(newmirnacancer)
cancermirnanames<- read.csv(file="cancer mirna names.csv")
rownames(newmirnacancer)=cancermirnanames[,1]
tnewmirnacancer<- t(newmirnacancer)
tnewmirnacancer<- tnewmirnacancer[-1,]
tnewL1HScancer<- tnewL1HScancer[-1,]
cancer<- cbind(tnewL1HScancer, tnewmirnacancer)
cancermirna<- cancer[,-1]

#format normal
newL1HSnormal<- data.frame(read.table(file="L1HS.Norma.BaseMeansA", header=TRUE))
newL1HSnormal<- data.matrix(newL1HSnormal)
tnewL1HSnormal<- t(newL1HSnormal)
newmirnanormal<- read.csv(file="Mirna.Normal.csv")
newmirnanormal<- data.matrix(newmirnanormal)
normalmirnanames<- read.csv(file="normal mirna names.csv")
rownames(newmirnanormal)=normalmirnanames[,1]
tnewmirnanormal<- t(newmirnanormal)
tnewmirnanormal<- tnewmirnanormal[-1,]
tnewL1HSnormal<- tnewL1HSnormal[-1,]
normal<- cbind(tnewL1HSnormal, tnewmirnanormal)
normalmirna<- normal[,-1]

L1HScancer<- data.frame((cancer[,1]))
L1HSnormal<- data.frame(normal[,1])
mirnacancer<- data.frame(cancermirna)
mirnanormal<- data.frame(normalmirna)


#ALL CANCER FOLDCHANGE RATIO
Allraw<- subset(NewNames, NewNames[,2]=="LIHC")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))
 
#599
mir599mirna<- data.frame(logofmirnaall[,c("hsa.mir.599")])
rownames(mir599mirna)=rownames(logofmirnaall)
mir599L1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir599mirna))
mir599<- cbind(mir599L1HS, mir599mirna)
result599 <-subset(mir599, !is.na(mir599[,2]))
colnames(result599)=c("MiRNA", "L1HS")
reg=lm(result599[,2] ~ result599[,1])
qplot(data=result599, x=result599$MiRNA, y=result599$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-599 Correlation of MiRNA and L1HS in LIHC") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#3140
Allraw<- subset(NewNames, NewNames[,2]=="THCA")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))

mir3140mirna<- data.frame(logofmirnaall[,c("hsa.mir.3140")])
rownames(mir3140mirna)=rownames(logofmirnaall)
mir3140L1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir3140mirna))
mir3140<- cbind(mir3140L1HS, mir3140mirna)
result3140 <-subset(mir3140, !is.na(mir3140[,2]))
colnames(result3140)=c("MiRNA", "L1HS")
reg=lm(result3140[,2] ~ result3140[,1])
qplot(data=result3140, x=result3140$MiRNA, y=result3140$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-3140 Correlation of MiRNA and L1HS in THCA") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#656
Allraw<- subset(NewNames, NewNames[,2]=="STAD")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))

mir656mirna<- data.frame(logofmirnaall[,c("hsa.mir.656")])
rownames(mir656mirna)=rownames(logofmirnaall)
mir656L1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir656mirna))
mir656<- cbind(mir656L1HS, mir656mirna)
result656 <-subset(mir656, !is.na(mir656[,2]))
colnames(result656)=c("MiRNA", "L1HS")
reg=lm(result656[,2] ~ result656[,1])
qplot(data=result656, x=result656$MiRNA, y=result656$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-656 Correlation of MiRNA and L1HS in STAD") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#324
mir324mirna<- data.frame(logofmirnaall[,c("hsa.mir.324")])
rownames(mir324mirna)=rownames(logofmirnaall)
mir324L1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir324mirna))
mir324<- cbind(mir324L1HS, mir324mirna)
result324 <-subset(mir324, !is.na(mir324[,2]))
colnames(result324)=c("MiRNA", "L1HS")
reg=lm(result324[,2] ~ result324[,1])
qplot(data=result324, x=result324$MiRNA, y=result324$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-324 Correlation of MiRNA and L1HS in STAD") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#376c
Allraw<- subset(NewNames, NewNames[,2]=="CHOL")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))

mir376cmirna<- data.frame(logofmirnaall[,c("hsa.mir.376c")])
rownames(mir376cmirna)=rownames(logofmirnaall)
mir376cL1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir376cmirna))
mir376c<- cbind(mir376cL1HS, mir376cmirna)
result376c <-subset(mir376c, !is.na(mir376c[,2]))
colnames(result376c)=c("MiRNA", "L1HS")
reg=lm(result376c[,2] ~ result376c[,1])
qplot(data=result376c, x=result376c$MiRNA, y=result376c$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-376c Correlation of MiRNA and L1HS in CHOL") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#376b
mir376bmirna<- data.frame(logofmirnaall[,c("hsa.mir.376b")])
rownames(mir376bmirna)=rownames(logofmirnaall)
mir376bL1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir376bmirna))
mir376b<- cbind(mir376bL1HS, mir376bmirna)
result376b <-subset(mir376b, !is.na(mir376b[,2]))
colnames(result376b)=c("MiRNA", "L1HS")
reg=lm(result376b[,2] ~ result376b[,1])
qplot(data=result376b, x=result376b$MiRNA, y=result376b$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-376b Correlation of MiRNA and L1HS in CHOL") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#204
Allraw<- subset(NewNames, NewNames[,2]=="PRAD")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))

mir204mirna<- data.frame(logofmirnaall[,c("hsa.mir.204")])
rownames(mir204mirna)=rownames(logofmirnaall)
mir204L1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir204mirna))
mir204<- cbind(mir204L1HS, mir204mirna)
result204 <-subset(mir204, !is.na(mir204[,2]))
colnames(result204)=c("MiRNA", "L1HS")
reg=lm(result204[,2] ~ result204[,1])
qplot(data=result204, x=result204$MiRNA, y=result204$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-204 Correlation of MiRNA and L1HS in PRAD") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#374a
Allraw<- subset(NewNames, NewNames[,2]=="KICH")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))

mir374amirna<- data.frame(logofmirnaall[,c("hsa.mir.374a")])
rownames(mir374amirna)=rownames(logofmirnaall)
mir374aL1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir374amirna))
mir374a<- cbind(mir374aL1HS, mir374amirna)
result374a <-subset(mir374a, !is.na(mir374a[,2]))
colnames(result374a)=c("MiRNA", "L1HS")
reg=lm(result374a[,2] ~ result374a[,1])
qplot(data=result374a, x=result374a$MiRNA, y=result374a$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-374a Correlation of MiRNA and L1HS in KICH") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)

#3191
Allraw<- subset(NewNames, NewNames[,2]=="KIRP")
dfallraw<- data.frame(Allraw)

#cancer L1HS
matched<- match(dfallraw[,1], rownames(L1HScancer))
allL1HScancer <-L1HScancer[matched,]
allL1HScancer<- data.frame(allL1HScancer)
rownames(allL1HScancer)=dfallraw[,1]

#normal L1HS
matched<- match(dfallraw[,1], rownames(L1HSnormal))
allL1HSnormal <-L1HSnormal[matched,]
allL1HSnormal<- data.frame(allL1HSnormal)
rownames(allL1HSnormal)=dfallraw[,1]
dividenL1HS<- allL1HScancer/allL1HSnormal
logofL1HSall<- log2(dividenL1HS)

#cancermirna
matchedmirna<- match(dfallraw[,1], rownames(mirnacancer))
allmirnacancer<- mirnacancer[matched,]
allmirnacancer<- data.frame(allmirnacancer)
rownames(allmirnacancer)=dfallraw[,1]

#normal mirna
matchedmirna<- match(dfallraw[,1], rownames(mirnanormal))
allmirnanormal<- mirnanormal[matched,]
allmirnanormal<- data.frame(allmirnanormal)
rownames(allmirnacancer)=dfallraw[,1]
dividenmirna<- allmirnacancer/allmirnanormal
logofmirnaall<- log2(dividenmirna)
is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.infinite))
 is.na(logofmirnaall) <- do.call(cbind,lapply(logofmirnaall, is.nan))

mir3191mirna<- data.frame(logofmirnaall[,c("hsa.mir.3191")])
rownames(mir3191mirna)=rownames(logofmirnaall)
mir3191L1HS<- subset(logofL1HSall, rownames(logofL1HSall)==rownames(mir3191mirna))
mir3191<- cbind(mir3191L1HS, mir3191mirna)
result3191 <-subset(mir3191, !is.na(mir3191[,2]))
colnames(result3191)=c("MiRNA", "L1HS")
reg=lm(result3191[,2] ~ result3191[,1])
qplot(data=result3191, x=result3191$MiRNA, y=result3191$L1HS,
		xlab="log2 Fold Change MiRNA", ylab="log2 Fold Change L1HS",
		main="hsa-mir-3191 Correlation of MiRNA and L1HS in KIRP") +
geom_point(shape=1)+
geom_smooth(method=lm, se=FALSE)
