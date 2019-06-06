setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer")
meanstd<-read.csv(“Final2.csv”, header=TRUE, sep=”\t”)
 
x<-1
a<-meanstd[x,103]  #mean
b<-meanstd[x,104]  #standard deviation
c=a+2*b 
test<- meanstd[x,1:102]
sig<-replace(test,test<c,0)   

#Because I am only interested in those that are greater than 2 stds away (those that are >a+2*b), I replaced those that are less and '2' to '0'

repeat {
x<-1+x
a<-meanstd[x,103]
b<-meanstd[x,104]
c=a+2*b
test<- meanstd[x,1:102]
sig<-rbind(sig, replace (test, test <c, 0))
if (x >162) break
}

> write.table(sig,"significantones.txt",row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
