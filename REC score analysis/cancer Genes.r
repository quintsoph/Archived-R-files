####Cancer Genes
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis/GDG")
cancergenes<-read.table("cancergenes.txt")
rownames(cancergenes)=cancergenes[,1]
setwd("E:/Bioinformatics Lab/Cancer Data/REC Score Analysis 2/results.genes.VST")
GDG<- data.frame(read.csv("Ordered Gene REC Scores Negative 2.csv"))
GDGREC<- as.matrix(GDG[,c(1,2)])
rownames(GDGREC)=GDGREC[,1]

#match!!!
matchednum<- data.frame(match(cancergenes[,1], GDGREC[,1]))
rownames(matchednum)=rownames(cancergenes)
g1000<- ifelse(matchednum[,1] < 1000, TRUE, FALSE)
together<- cbind(g1000, matchednum)
togetherordered <- together[order(together[,2]),, drop=FALSE]
table(together[,1])["TRUE"]
test<- GDGREC[which(GDGREC[,1] %in% cancergenes[,1]),]
togetherordered<- na.omit(togetherordered)
results<- cbind(togetherordered, test)
colnames(results)=c("if less than 1000", "ranks in Neg", "Genes", "REC")
write.table(results, "Cancer Genes in the Neg REC score 2.txt", sep="\t", quote=FALSE)

#Regular analysis
setwd("F:/Bioinformatics Lab/Cancer Data/REC Score Analysis/Regular")
cancergenes<-read.table("cancergenes.txt")
rownames(cancergenes)=cancergenes[,1]
Regular<- data.frame(read.csv("Ordered Gene REC Scores Negative Regular.csv"))
RegularREC<- as.matrix(Regular[,c(1,2)])
rownames(RegularREC)=RegularREC[,1]

#match!!!
matchednum<- data.frame(match(cancergenes[,1], RegularREC[,1]))
rownames(matchednum)=rownames(cancergenes)
g1000<- ifelse(matchednum[,1] < 1000, TRUE, FALSE)
together<- cbind(g1000, matchednum)
togetherordered <- together[order(together[,2]),, drop=FALSE]
table(together[,1])["TRUE"]
test<- RegularREC[which(RegularREC[,1] %in% cancergenes[,1]),]
togetherordered<- na.omit(togetherordered)
results<- cbind(togetherordered, test)
colnames(results)=c("if less than 1000", "ranks in REC", "Genes", "REC")
write.table(results, "Cancer Genes in the Regular Neg REC score.txt", sep="\t", quote=FALSE)

