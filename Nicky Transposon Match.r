
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
