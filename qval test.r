#Cut out all Q-Values for the NEB gene 
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
raw_qs<- read.table("GeneQvals.txt", header=TRUE)
NEB_qs<- subset(raw_qs, raw_qs[,1] == "NEB", drop=FALSE)
write.table(NEB_qs, file="GENE NEB Pvals and Qvals.txt", sep="\t", quote=FALSE)