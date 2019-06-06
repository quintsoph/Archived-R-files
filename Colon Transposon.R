setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells/Transposon Data Similar")
f1<- read.table("TCGA-A6-5662-01A.cntTable")
f1r<- f1c[25141,]
f1c <- f1r[,3]
f1c


setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells/Transposon Data Similar")
f2<- read.table("TCGA-F4-6704-01A.cntTable")
f2r<- f2[25141,]
f2c <- f2r[,3]
f2c

reads per million miRNA reads

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer")
file<- read.csv("TransposonVmiRNA Raw.csv")
plot(file$hsa.let.7a.1, file$L1HS, xlab="hsa.let.7a.1", ylab="L1HS", main="miRNA V. L1HS", col="green")
plot(file$hsa.let.7a.2, file$L1HS, xlab="hsa.let.7a.2", ylab="L1HS", main="miRNA V. L1HS", col="blue")
plot(file$hsa.let.7a.3, file$L1HS, xlab="hsa.let.7a.3", ylab="L1HS", main="miRNA V. L1HS", col="blue")
plot(file$hsa.mir.100, file$L1HS, xlab="hsa.mir.100", ylab="L1HS", main="miRNA V. L1HS", col="blue")

plot(file$hsa.mir.24.1, file$L1HS, xlab="hsa.mir.24.1", ylab="L1HS", main="miRNA V. L1HS", col="dark green")
abline(lm(file$L1HS~file$hsa.mir.24.1))
fit.2 <- lm(file$L1HS ~ file$hsa.mir.24.1)

significance of correlation r^2 or p
linear models

plot1= is the patients L1HS, and first MiRNA
plot1a= is without patient names

curve(coef(fit.2)[1] + coef(fit.2)[2]*x, add=TRUE)
In general, the higher the R-squared, the better the model fits your data.

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer")
files<- read.csv("TransposonVmiRNA Raw.csv")
file<- files[,-c(1)]
file2<- file[,-c(1)]
#make a function
testpvalues<- function(x)
{
value <- lm(file$L1HS ~ file$hsa.mir)
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
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer")
write.table(pvalues, file="p-values", sep=",")
write.table(cont, file="p-values conditions", sep=",")
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
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer")
write.table(RealPV, file="True P-Values", sep=",")


hsa.mir.1224     0.031875901
hsa.mir.124.2    0.032119996
hsa-mir-1251     0.024736957
hsa-mir-135a-1   0.036538536
hsa-mir-153-2    0.003903414
hsa-mir-181a-2   0.033005585
hsa-mir-2277     0.033548554
hsa-mir-24-1     0.000974501
hsa-mir-27b      0.026958708
hsa-mir-3065     0.021412869
hsa-mir-30b      0.038502373
hsa-mir-30c-2    0.040282478
hsa-mir-3128     0.007416874
hsa-mir-3138     0.007098966
hsa-mir-3158-1   0.040410516
hsa-mir-3164     0.007892258
hsa-mir-3179-2   0.007892258
hsa-mir-32       0.047286791
hsa-mir-320b-2   0.008238413
hsa-mir-324      0.037563392
hsa-mir-331      0.031989452
hsa-mir-33b      0.014703619
hsa-mir-345      0.04991828
hsa-mir-3607     0.001608531
hsa-mir-3647     0.033381932
hsa-mir-3651     0.027197092
hsa-mir-3653     0.013904727
hsa-mir-3654     0.019051408
hsa-mir-374a     0.01780868
hsa-mir-3910-2   0.004452029
hsa-mir-3920     0.014306502
hsa-mir-516a-1   0.04806702
hsa-mir-542      0.016379808
hsa-mir-574      0.007782052
hsa-mir-663b     0.007892258
hsa-mir-887      0.044926135
hsa-mir-934      0.011632546

name= c(hsa.mir.1224
,hsa.mir.124.2
,hsa-mir-1251
,hsa-mir-135a-1
,hsa-mir-153-2
,hsa-mir-181a-2
,hsa-mir-2277
,hsa-mir-24-1
,hsa-mir-27b
,hsa-mir-3065
,hsa-mir-30b
,hsa-mir-30c-2
,hsa-mir-3128
,hsa-mir-3138
,hsa-mir-3158-1
,hsa-mir-3164
,hsa-mir-3179-2
,hsa-mir-32
,hsa-mir-320b-2
,hsa-mir-324
,hsa-mir-331
,hsa-mir-33b
,hsa-mir-345
,hsa-mir-3607
,hsa-mir-3647
,hsa-mir-3651
,hsa-mir-3653
,hsa-mir-3654
,hsa-mir-374a
,hsa-mir-3910-2
,hsa-mir-3920
,hsa-mir-516a-1
,hsa-mir-542
,hsa-mir-574
,hsa-mir-663b
,hsa-mir-887
,hsa-mir-934

values= c(0.031875901
,0.032119996
,0.024736957
,0.03653853
,0.003903414
,0.033005585
,0.033548554
,0.000974501
,0.026958708
,0.021412869
,0.038502373
,0.040282478
,0.007416874
,0.007098966
,0.040410516
,0.007892258
,0.007892258
,0.047286791
,0.008238413
,0.037563392
,0.031989452
,0.014703619
,0.04991828
,0.001608531
,0.033381932
,0.027197092
,0.013904727
,0.019051408
,0.01780868
,0.004452029
,0.014306502
,0.04806702
,0.016379808
,0.007782052
,0.007892258
,0.044926135
,0.011632546


TrueP= matrix(data=name, byrow=FALSE)
TrueV= matrix(data=name, byrow=FALSE)
Real= cbind(TrueP, TrueV)

3rd column base mean A is cancer
4th column base mean B is normal
500 files



resultspvalues=NULL
for (i in colnames(tlogofmirna))
{
	a=NULL
	a <-ratiopvalues(tlogofmirna[,i])
	resultspvalues <- cbind(resultspvalues, q) 
}

resultspvalues

1.403156

-.0925141

