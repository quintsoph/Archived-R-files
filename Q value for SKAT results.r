#q-values for SKAT
#fold change
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results")
foldmed<- data.frame(read.table("withmedian2.txt", header=TRUE))
q<- data.frame(p.adjust(foldmed[,3], method="BH"))
result<- cbind(foldmed, q)
qfoldmed<- result[,c(1,4)]
colnames(qfoldmed)=c("qvalue", "Gene")
write.table(qfoldmed, "qvaluefoldmed.txt", quote=FALSE, sep="\t")

foldno<- data.frame(read.table("nomedian2.txt", header=TRUE))
q<- data.frame(p.adjust(foldno[,3], method="BH"))
result<- cbind(foldno, q)
qfoldno<- result[,c(1,4)]
colnames(qfoldno)=c("Gene", "qvalue")
write.table(qfoldno, "qvaluefoldnomed.txt", quote=FALSE, sep="\t")

#cancer rpm
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Final Results/SKAT Results/cancerrpm")
cancermed<-data.frame(read.table("cancerwithmed2.txt", header=TRUE))
q<- data.frame(p.adjust(cancermed[,3], method="BH"))
result<- cbind(cancermed, q)
qcancermed<- result[,c(1,4)]
colnames(qcancermed)=c("Gene", "qvalue")
write.table(qcancermed, "qvaluecancermed.txt", quote=FALSE, sep="\t")

cancerno<- data.frame(read.table("cancernomed2.txt", header=TRUE))
q<- data.frame(p.adjust(cancerno[,3], method="BH"))
result<- cbind(cancerno, q)
qcancerno<- result[,c(1,4)]
colnames(qcancerno)=c("Gene", "qvalue")
write.table(qcancerno, "qvaluecancernomed.txt", quote=FALSE, sep="\t")
