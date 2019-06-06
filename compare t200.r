#upload files
setwd("E:/Bioinformatics Lab/Cancer Data/Research")
genelist<- read.table(file="genes.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
results<- read.table("withmedian2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
matched<- matched[!is.na(matched)]
test<- tresult[,matched]
ttest<- t(test)
rownames(ttest)= ttest[,1]
numeric<- data.frame(as.numeric(ttest[,2]))
rownames(numeric)= ttest[,1]
colnames(numeric)= c("pvals")
newdata<- numeric[order(numeric[,1]) ,, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
write.table(newdata, file="withmedian2 t200 matches.txt", quote=FALSE, sep="\t")

#no median mutations
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
original<- read.table(file="cancertype.txt", header=TRUE, row.names=1)

#upload files
setwd("E:/Bioinformatics Lab/Cancer Data/Research")
genelist<- read.table(file="genes.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
results<- read.table("nomedian2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
matched<- matched[!is.na(matched)]
test<- tresult[,matched]
ttest<- t(test)
rownames(ttest)= ttest[,1]
numeric<- data.frame(as.numeric(ttest[,2]))
rownames(numeric)= ttest[,1]
colnames(numeric)= c("pvals")
newdata<- numeric[order(numeric[,1]) ,, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
write.table(newdata, file="t200 matches nomedian2.txt", quote=FALSE, sep="\t")

#rank only
#write in the files
setwd("E:/Bioinformatics Lab/Cancer Data/Research")
genelist<- read.table(file="genes.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
results<- read.table("nomedian2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
ranksnomed<- cbind(genelist, matched)
ranksnomed2<- data.frame(ranksnomed[,-c(1)])
rownames(ranksnomed2)= ranksnomed[,1]
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
results<- read.table("withmedian2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
rankswithmed<- cbind(genelist, matched)
rankswithmed2<- data.frame(rankswithmed[,-c(1)])
rownames(rankswithmed2)= rankswithmed[,1]
combined<- cbind(rankswithmed2, ranksnomed2)
colnames(combined)=c("withmedian", "nomedian")
write.table(combined, file="ranks of t200.txt", quote=FALSE, sep="\t")





#CANCER RPM
#pval results with median mutations
setwd("E:/Bioinformatics Lab/Cancer Data/Research")
genelist<- read.table(file="genes.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
results<- read.table("cancerwithmed2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
matched<- matched[!is.na(matched)]
test<- tresult[,matched]
ttest<- t(test)
rownames(ttest)= ttest[,1]
numeric<- data.frame(as.numeric(ttest[,2]))
rownames(numeric)= ttest[,1]
colnames(numeric)= c("pvals")
newdata<- numeric[order(numeric[,1]) ,, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
write.table(newdata, file="withmedian cancerrpm t200 matches.txt", quote=FALSE, sep="\t")


#pvals results with no median
setwd("E:/Bioinformatics Lab/Cancer Data/Research")
genelist<- read.table(file="genes.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
results<- read.table("cancernomed2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
matched<- matched[!is.na(matched)]
test<- tresult[,matched]
ttest<- t(test)
rownames(ttest)= ttest[,1]
numeric<- data.frame(as.numeric(ttest[,2]))
rownames(numeric)= ttest[,1]
colnames(numeric)= c("pvals")
newdata<- numeric[order(numeric[,1]) ,, drop=FALSE]
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
write.table(newdata, file="no median cancerrpm t200 matches.txt", quote=FALSE, sep="\t")

#create rpm ranks
setwd("E:/Bioinformatics Lab/Cancer Data/Research")
genelist<- read.table(file="genes.txt", header=TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
results<- read.table("cancernomed2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
ranksnomed<- cbind(genelist, matched)
ranksnomed2<- data.frame(ranksnomed[,-c(1)])
rownames(ranksnomed2)= ranksnomed[,1]
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
results<- read.table("cancerwithmed2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(genelist[,1], result[,1])
rankswithmed<- cbind(genelist, matched)
rankswithmed2<- data.frame(rankswithmed[,-c(1)])
rownames(rankswithmed2)= rankswithmed[,1]
combined<- cbind(rankswithmed2, ranksnomed2)
colnames(combined)=c("withmedian", "nomedian")
write.table(combined, file="ranks of t200 cancerrpm.txt", quote=FALSE, sep="\t")

#transposon matches
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations")
List<- read.table("transposon_genes 2.txt")
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
results<- read.table("cancernomed2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(List[,1], result[,1])
matched2<- matched[!is.na(matched)]
test<- tresult[,matched]
ttest<- t(test)
final<- cbind(ttest, matched)
write.table(final, "transposons no median cancerrpm all.txt", quote=FALSE, sep="\t")

setwd("E:/Bioinformatics Lab/Cancer Data/Mutations")
List<- read.table("transposon_genes 2.txt")
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
results<- read.table("nomedian2.txt", header=TRUE)
result<- results[-c(2)]
tresult<- t(result)
matched<- match(List[,1], result[,1])
matched2<- matched[!is.na(matched)]
test<- tresult[,matched]
ttest<- t(test)
final<- cbind(ttest, matched)
write.table(final, "transposons without median.txt", quote=FALSE, sep="\t")
