3rd column base mean A is cancer
4th column base mean B is normal
500 files
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
L1HS<- read.csv("L1HS.csv")
L1HScancer<- L1HS[,c(1,3)]
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/MiRNA files")
mirnafiles <- dir(pattern = "*.miRna", full.names = TRUE)
lst <- lapply(mirnafiles, read.table, header=TRUE, sep='')
mirna<- sapply(lst, '[[',3)
file_list <- list.files("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/MiRNA files")
colnames(mirna)=file_list
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
names<- read.table("miRNA names.txt")
tnames=t(names)
rownames(mirna)=tnames
tmirna=t(mirna)
file <-cbind(L1HScancer,tmirna)
files<-file[,-c(1)]
colnames(files)[1]<-"L1HS"
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
write.table(files, file="L1HS and miRNA", sep=",")


file2<- files[,-c(1)]
#make a function
testpvalues<- function(x)
{
value <- lm(files$L1HS ~ x)
anova(value)$Pr[1]
}

#make a for loop
pvalues=NULL
for (i in 1:ncol(file2))
{
p <-testpvalues(file2[,i])
pvalues <- cbind(pvalues, p) 
}

colnames(pvalues)=colnames(file2)
pvalues <0.05
cont=pvalues<0.05
cont2=t(cont)
colnames(cont2)= c("pvalues")
df<-data.frame(cont2)
df
sub<-subset(df, pvalues==TRUE)
pvaluesdf<- data.frame(pvalues)
pt<- t(pvaluesdf)
colnames(pt)= c("dataPV")
data<-subset(pt, pt[,1] < 0.05)
RealPV<- cbind(sub, data)
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
write.table(RealPV, file="PValues for all cancer", sep=",")
