#Generate All cancer type log2 pvalues from new data
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results")
allcancer<- data.frame(read.table("log2 foldchange ratio only.txt", header=TRUE, row.names=1))
allcancercut<- data.frame(allcancer[,c("pvalues")])
rownames(allcancercut)=rownames(allcancer)
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
nullps<- read.table("NULLS p-values 1G NEW.txt", header=TRUE)
testm<- match(rownames(allcancercut), rownames(nullps))
testm<- testm[!is.na(testm)]
lists<-data.frame(rownames(nullps))
comb<- cbind(lists, nullps)
resultnames<- data.frame(comb[testm,])
match(rownames(allcancercut), rownames(resultnames))
together<- cbind(resultnames, allcancercut)
processed<- na.omit(together)

write.table(processed, file="correlation check.txt", sep="\t")
#graph!

plot(processed[,3], processed[,2], xlab="P-Values from Our MiRNA log2 Foldchang Ratio All Cancers", ylab="1G Nulls P-Values", main="Correlation", col="Green")
abline(lm(processed[,2] ~ processed[,3]))

#Try with the averages of the p-values

#read in the created files
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results/log2 foldchange ratio by each")
BLCA<- data.frame(read.table(file="BLCA pvalues and rvalues of log2 foldchange ratio.txt"))
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2 foldchange ratio.txt"))
CESC<- data.frame(read.table(file="CESC pvalues and rvalues of log2 foldchange ratio.txt"))
CHOL<- data.frame(read.table(file="CHOL pvalues and rvalues of log2 foldchange ratio.txt"))
ESCA<- data.frame(read.table(file="ESCA pvalues and rvalues of log2 foldchange ratio.txt"))
HNSC<- data.frame(read.table(file="HNSC pvalues and rvalues of log2 foldchange ratio.txt"))
KICH<- data.frame(read.table(file="KICH pvalues and rvalues of log2 foldchange ratio.txt"))
KIRC<- data.frame(read.table(file="KIRC pvalues and rvalues of log2 foldchange ratio.txt"))
KIRP<- data.frame(read.table(file="KIRP pvalues and rvalues of log2 foldchange ratio.txt"))
LIHC<- data.frame(read.table(file="LIHC pvalues and rvalues of log2 foldchange ratio.txt"))
LUAD<- data.frame(read.table(file="LUAD pvalues and rvalues of log2 foldchange ratio.txt"))
LUSC<- data.frame(read.table(file="LUSC pvalues and rvalues of log2 foldchange ratio.txt"))
PAAD<- data.frame(read.table(file="PAAD pvalues and rvalues of log2 foldchange ratio.txt"))
PCPG<- data.frame(read.table(file="PCPG pvalues and rvalues of log2 foldchange ratio.txt"))
PRAD<- data.frame(read.table(file="PRAD pvalues and rvalues of log2 foldchange ratio.txt"))
STAD<- data.frame(read.table(file="STAD pvalues and rvalues of log2 foldchange ratio.txt"))
THCA<- data.frame(read.table(file="THCA pvalues and rvalues of log2 foldchange ratio.txt"))
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2 foldchange ratio.txt"))
UCEC<- data.frame(read.table(file="UCEC pvalues and rvalues of log2 foldchange ratio.txt"))
READ<- data.frame(read.table(file="READ pvalues and rvalues of log2 foldchange ratio.txt"))
FPPP<- data.frame(read.table(file="FPPP pvalues and rvalues of log2 foldchange ratio.txt"))
COAD<- data.frame(read.table(file="COAD pvalues and rvalues of log2 foldchange ratio.txt"))

#merge all files
BLCA2<- data.frame(t(BLCA))
BLCA2<- BLCA2[-2,]
BRCA2<- data.frame(t(BRCA))
BRCA2<- BRCA2[-2,]
combine<- merge(BLCA2, BRCA2, all=TRUE)
CESC2<- data.frame(t(CESC))
CESC2<- CESC2[-2,]
combine2<- merge(combine, CESC2, all=TRUE)
CHOL2<- data.frame(t(CHOL))
CHOL2<- CHOL2[-2,]
combine3<- merge(combine2, CHOL2, all=TRUE)
ESCA2<- data.frame(t(ESCA))
ESCA2<- ESCA2[-2,]
combine4<- merge(combine3, ESCA2, all=TRUE)
HNSC2<- data.frame(t(HNSC))
HNSC2<- HNSC2[-2,]
combine5<- merge(combine4, HNSC2, all=TRUE)
KICH2<- data.frame(t(KICH))
KICH2<- KICH2[-2,]
combine6<- merge(combine5, KICH2, all=TRUE)
KIRC2<- data.frame(t(KIRC))
KIRC2<- KIRC2[-2,]
combine7<- merge(combine6, KIRC2, all=TRUE)
KIRP2<- data.frame(t(KIRP))
KIRP2<- KIRP2[-2,]
combine8<- merge(combine7, KIRP2, all=TRUE)
LIHC2<- data.frame(t(LIHC))
LIHC2<- LIHC2[-2,]
combine9<- merge(combine8, LIHC2, all=TRUE)
LUAD2<- data.frame(t(LUAD))
LUAD2<- LUAD2[-2,]
combine10<- merge(combine9, LUAD2, all=TRUE)
LUSC2<- data.frame(t(LUSC))
LUSC2<- LUSC2[-2,]
combine11<- merge(combine10, LUSC2, all=TRUE)
PAAD2<- data.frame(t(PAAD))
PAAD2<- PAAD2[-2,]
combine12<- merge(combine11, PAAD2, all=TRUE)
PCPG2<- data.frame(t(PCPG))
PCPG2<- PCPG2[-2,]
combine13<- merge(combine12, PCPG2, all=TRUE)
PRAD2<- data.frame(t(PRAD))
PRAD2<- PRAD2[-2,]
combine14<- merge(combine13, PRAD2, all=TRUE)
STAD2<- data.frame(t(STAD))
STAD2<- STAD2[-2,]
combine15<- merge(combine14, STAD2, all=TRUE)
THCA2<- data.frame(t(THCA))
THCA2<- THCA2[-2,]
combine16<- merge(combine15, THCA2, all=TRUE)
THYM2<- data.frame(t(THYM))
THYM2<- THYM2[-2,]
combine17<- merge(combine16, THYM2, all=TRUE)
UCEC2<- data.frame(t(UCEC))
UCEC2<- UCEC2[-2,]
combine18<- merge(combine17, UCEC2, all=TRUE)
COAD2<- data.frame(t(COAD))
COAD2<- COAD2[-2,]
combine19<- merge(combine18, COAD2, all=TRUE)
notna<- colSums(!is.na(combine19))
sum<- colSums(combine19, na.rm=TRUE)
avpval<- data.frame(sum/notna)
pvalav<- na.omit(avpval)

#create graph!
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
nullps<- read.table("results2.txt", header=TRUE)
testm<- match(rownames(pvalav), rownames(nullps))
testm<- testm[!is.na(testm)]
lists<-data.frame(rownames(nullps))
comb<- cbind(lists, nullps)
resultnames<- data.frame(comb[testm,])
match(rownames(pvalav), rownames(resultnames))
together<- cbind(resultnames, pvalav)
processed<- na.omit(together)

#graph!
plot(together[,3], together[,2], xlab="P-Values from Our MiRNA Cancer log2 foldchange total average", ylab="1G Nulls P-Values", main="Correlation", col="red")
abline(lm(together[,2] ~ together[,3]))

write.table(processed, file="correlation check 2.txt", sep="\t")

#Try with median
medpval<- data.frame(apply(combine19, 2, median, na.rm=TRUE))
colnames(medpval)=c("mut")
pvalmed<- na.omit(medpval)

#create graph!
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
nullps<- read.table("NULLS p-values 1G NEW.txt", header=TRUE)
testm<- match(rownames(pvalmed), rownames(nullps))
testm<- testm[!is.na(testm)]
lists<-data.frame(rownames(nullps))
comb<- cbind(lists, nullps)
resultnames<- data.frame(comb[testm,])
match(rownames(pvalmed), rownames(resultnames))
together<- cbind(resultnames, pvalmed)
processed<- na.omit(together)

#graph!
plot(processed[,3], processed[,2], xlab="P-Values from Our MiRNA Cancer log2 foldchange total median", ylab="1G Nulls P-Values", main="Correlation", col="orange")
abline(lm(processed[,2] ~ processed[,3]))
cor(processed[,3], processed[,2])

write.table(processed, file="correlation check 3.txt", sep="\t")


#cancer not fold change
#read in the created files
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results/log2 cancer by each")
BLCA<- data.frame(read.table(file="BLCA pvalues and rvalues of log2 cancer.txt"))
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2 cancer.txt"))
CESC<- data.frame(read.table(file="CESC pvalues and rvalues of log2 cancer.txt"))
CHOL<- data.frame(read.table(file="CHOL pvalues and rvalues of log2 cancer.txt"))
ESCA<- data.frame(read.table(file="ESCA pvalues and rvalues of log2 cancer.txt"))
HNSC<- data.frame(read.table(file="HNSC pvalues and rvalues of log2 cancer.txt"))
KICH<- data.frame(read.table(file="KICH pvalues and rvalues of log2 cancer.txt"))
KIRC<- data.frame(read.table(file="KIRC pvalues and rvalues of log2 cancer.txt"))
KIRP<- data.frame(read.table(file="KIRP pvalues and rvalues of log2 cancer.txt"))
LIHC<- data.frame(read.table(file="LIHC pvalues and rvalues of log2 cancer.txt"))
LUAD<- data.frame(read.table(file="LUAD pvalues and rvalues of log2 cancer.txt"))
LUSC<- data.frame(read.table(file="LUSC pvalues and rvalues of log2 cancer.txt"))
PAAD<- data.frame(read.table(file="PAAD pvalues and rvalues of log2 cancer.txt"))
PCPG<- data.frame(read.table(file="PCPG pvalues and rvalues of log2 cancer.txt"))
PRAD<- data.frame(read.table(file="PRAD pvalues and rvalues of log2 cancer.txt"))
STAD<- data.frame(read.table(file="STAD pvalues and rvalues of log2 cancer.txt"))
THCA<- data.frame(read.table(file="THCA pvalues and rvalues of log2 cancer.txt"))
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2 cancer.txt"))
UCEC<- data.frame(read.table(file="UCEC pvalues and rvalues of log2 cancer.txt"))

#merge all files
BLCA2<- data.frame(t(BLCA))
BLCA2<- BLCA2[-2,]
BRCA2<- data.frame(t(BRCA))
BRCA2<- BRCA2[-2,]
combine<- merge(BLCA2, BRCA2, all=TRUE)
CESC2<- data.frame(t(CESC))
CESC2<- CESC2[-2,]
combine2<- merge(combine, CESC2, all=TRUE)
CHOL2<- data.frame(t(CHOL))
CHOL2<- CHOL2[-2,]
combine3<- merge(combine2, CHOL2, all=TRUE)
ESCA2<- data.frame(t(ESCA))
ESCA2<- ESCA2[-2,]
combine4<- merge(combine3, ESCA2, all=TRUE)
HNSC2<- data.frame(t(HNSC))
HNSC2<- HNSC2[-2,]
combine5<- merge(combine4, HNSC2, all=TRUE)
KICH2<- data.frame(t(KICH))
KICH2<- KICH2[-2,]
combine6<- merge(combine5, KICH2, all=TRUE)
KIRC2<- data.frame(t(KIRC))
KIRC2<- KIRC2[-2,]
combine7<- merge(combine6, KIRC2, all=TRUE)
KIRP2<- data.frame(t(KIRP))
KIRP2<- KIRP2[-2,]
combine8<- merge(combine7, KIRP2, all=TRUE)
LIHC2<- data.frame(t(LIHC))
LIHC2<- LIHC2[-2,]
combine9<- merge(combine8, LIHC2, all=TRUE)
LUAD2<- data.frame(t(LUAD))
LUAD2<- LUAD2[-2,]
combine10<- merge(combine9, LUAD2, all=TRUE)
LUSC2<- data.frame(t(LUSC))
LUSC2<- LUSC2[-2,]
combine11<- merge(combine10, LUSC2, all=TRUE)
PAAD2<- data.frame(t(PAAD))
PAAD2<- PAAD2[-2,]
combine12<- merge(combine11, PAAD2, all=TRUE)
PCPG2<- data.frame(t(PCPG))
PCPG2<- PCPG2[-2,]
combine13<- merge(combine12, PCPG2, all=TRUE)
PRAD2<- data.frame(t(PRAD))
PRAD2<- PRAD2[-2,]
combine14<- merge(combine13, PRAD2, all=TRUE)
STAD2<- data.frame(t(STAD))
STAD2<- STAD2[-2,]
combine15<- merge(combine14, STAD2, all=TRUE)
THCA2<- data.frame(t(THCA))
THCA2<- THCA2[-2,]
combine16<- merge(combine15, THCA2, all=TRUE)
THYM2<- data.frame(t(THYM))
THYM2<- THYM2[-2,]
combine17<- merge(combine16, THYM2, all=TRUE)
UCEC2<- data.frame(t(UCEC))
UCEC2<- UCEC2[-2,]
combine18<- merge(combine17, UCEC2, all=TRUE)

#Try with median
medpval<- data.frame(apply(combine18, 2, median, na.rm=TRUE))
colnames(medpval)=c("mut")
pvalmed<- na.omit(medpval)

#create graph!
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
nullps<- read.table("NULLS p-values 1G.txt", header=TRUE)
testm<- match(rownames(pvalmed), rownames(nullps))
testm<- testm[!is.na(testm)]
lists<-data.frame(rownames(nullps))
comb<- cbind(lists, nullps)
resultnames<- data.frame(comb[testm,])
match(rownames(pvalmed), rownames(resultnames))
together<- cbind(resultnames, pvalmed)
processed<- na.omit(together)

#graph!
plot(together[,3], together[,2], xlab="P-Values from Our MiRNA Cancer log2 total median", ylab="1G Nulls P-Values", main="Correlation", col="pink")
abline(lm(together[,2] ~ together[,3]))
cor(together[,3], together[,2])

##
#
#Read in log2 foldchange by each from above

#merge all files
BLCA2<- data.frame(t(BLCA))
BLCA2<- BLCA2[-2,]
BRCA2<- data.frame(t(BRCA))
BRCA2<- BRCA2[-2,]
CESC2<- data.frame(t(CESC))
CESC2<- CESC2[-2,]
CHOL2<- data.frame(t(CHOL))
CHOL2<- CHOL2[-2,]
ESCA2<- data.frame(t(ESCA))
ESCA2<- ESCA2[-2,]
HNSC2<- data.frame(t(HNSC))
HNSC2<- HNSC2[-2,]
KICH2<- data.frame(t(KICH))
KICH2<- KICH2[-2,]
KIRC2<- data.frame(t(KIRC))
KIRC2<- KIRC2[-2,]
KIRP2<- data.frame(t(KIRP))
KIRP2<- KIRP2[-2,]
LIHC2<- data.frame(t(LIHC))
LIHC2<- LIHC2[-2,]
LUAD2<- data.frame(t(LUAD))
LUAD2<- LUAD2[-2,]
LUSC2<- data.frame(t(LUSC))
LUSC2<- LUSC2[-2,]
PAAD2<- data.frame(t(PAAD))
PAAD2<- PAAD2[-2,]
PCPG2<- data.frame(t(PCPG))
PCPG2<- PCPG2[-2,]
PRAD2<- data.frame(t(PRAD))
PRAD2<- PRAD2[-2,]
STAD2<- data.frame(t(STAD))
STAD2<- STAD2[-2,]
THCA2<- data.frame(t(THCA))
THCA2<- THCA2[-2,]
THYM2<- data.frame(t(THYM))
THYM2<- THYM2[-2,]
UCEC2<- data.frame(t(UCEC))
UCEC2<- UCEC2[-2,]
COAD2<- data.frame(t(COAD))
COAD2<- COAD2[-2,]

pvalmed=t(BRCA2)
#create graph!
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
nullps<- read.table("NULLS p-values 1G NEW.txt", header=TRUE)
pvalmed=t(COAD2)
testm<- match(rownames(pvalmed), rownames(nullps))
testm<- testm[!is.na(testm)]
lists<-data.frame(rownames(nullps))
comb<- cbind(lists, nullps)
resultnames<- data.frame(comb[testm,])
match(rownames(pvalmed), rownames(resultnames))
together<- cbind(resultnames, pvalmed)
processed<- na.omit(together)

plot(processed[,3], processed[,2], xlab="P-Values from Our MiRNA Cancer log2 foldchange COAD", ylab="10G Nulls P-Values", main="Correlation", col="orange")
abline(lm(processed[,2] ~ processed[,3]))
cor(processed[,3], processed[,2])


#create graph!
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
nullps<- read.table("NULLS p-values 1G NEW.txt", header=TRUE)
REC<- read.table("REC ONLY NEW remove.txt", row.names=1)
RECN<- -1 * REC
pREC<- 10 ^ REC
ppREC<- 1- pREC
testm<- match(rownames(pREC), rownames(nullps))
testm<- testm[!is.na(testm)]
lists<-data.frame(rownames(nullps))
comb<- cbind(lists, nullps)
resultnames<- data.frame(comb[testm,])
match(rownames(pREC), rownames(resultnames))
together<- cbind(resultnames, pREC)
processed<- na.omit(together)

#graph!
plot(processed[,3], processed[,2], xlab="1 minus of the Positive REC scores", ylab="1G Nulls P-Values", main="Correlation 9", col="green")
abline(lm(processed[,2] ~ processed[,3]))
cor(processed[,3], processed[,2])
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files/NEW REC correlation")
write.table(processed, file="1 minus of the Positve REC scores.txt", sep="\t")

setwd("F:/Bioinformatics Lab/Cancer Data/New REC correlation")
pos<- read.table("REC Power analysis Positive REC ONLY.txt", header=TRUE, row.names=1)
neg<- read.table("REC Power analysis Negative Rec ONLY.txt", header=TRUE, row.names=1)
together<- rbind(pos, neg)
processed<- na.omit(together)

#graph!
plot(processed[,3], processed[,2], xlab="Combined REC scores Power", ylab="1G Nulls P-Values", main="Correlation 8", col="purple")

#Test for doing the covariate
#combo of the graphs
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files/New REC correlation")
pos<- read.table("1 minus of the Positve REC scores.txt", header=TRUE, row.names=1)
neg<- read.table("REC Power analysis Negative Rec ONLY.txt", header=TRUE, row.names=1)
together<- rbind(pos, neg)
processed<- na.omit(together)
plot(processed[,3], processed[,2], xlab="Combined REC scores Power", ylab="1G Nulls P-Values", main="Correlation 8", col="purple")



