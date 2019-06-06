setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/formated L1HS and mirna")
L1HS<- read.csv("L1HS.csv")
L1HScancer<- L1HS[,c(1,4)]

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/MiRNA files")
mirnafiles <- dir(pattern = "*.miRna", full.names = TRUE)
lst <- lapply(mirnafiles, read.table, header=TRUE, sep='')
mirna<- sapply(lst, '[[',3)

#Adds names of patients.
file_list <- list.files("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/MiRNA files")
colnames(mirna)=file_list

#Adds miRna names to data.
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/formated L1HS and mirna")
names<- read.table("miRNA names.txt")
tnames=t(names)
rownames(mirna)=tnames
tmirna=t(mirna)

#Combine data into single dataframe.
combinedData <-cbind(L1HScancer,tmirna)

files <- combinedData[,-c(1)]
colnames(files)[1]<-"L1HS"

#files is both mirna and l1hs
files1<-log2(files)

#Change inf to NAs.
is.na(files1) <- sapply(files1, is.infinite)

#Files2 is our data without L1HS
files2 <- files1[,-c(1)]

#Function to calculate p-values.
testpvalues<- function(x)
{
	print(x)
	formula <- paste("files1$L1HS ~ ", files2[x], sep="")
	value <- lm(formula)
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(files2))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(files2[i]))){
		p <-testpvalues(i)
	}
	else{
		print(paste("Skipping", i))
	}
	pvalues <- cbind(pvalues, p) 
}

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/raw results")
names<- read.csv("names of log mirna.csv", header=FALSE)
tnames<- t(names)
colnames(pvalues)=tnames[1,]
pvalues<- t(pvalues)
qvalues<- p.adjust(pvalues, method="BH")
tqvalues<- t(qvalues)
colnames(tqvalues)=tnames[1,]
qvalues<- t(tqvalues)
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
write.table(RealPV, file="log True Pvalues.txt", sep=" ")


#set of first 10 MiRNA for comparison
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
ana<- read.table("PValues for all cancer raw", sep=",")

setoften<- files[,c(1:11)]

logoften<-log2(setoften)

is.na(logoften) <- sapply(logoften, is.infinite)
logoften

noL1HS<-logoften[,-c(1)]

#make a function
testpvalues<- function(x)
{
print(x)
formula <- paste("logoften$L1HS ~ ", noL1HS[x], sep="")
value <- lm(formula)
anova(value)$Pr[1]
}


#make a for loop
trypvalues=NULL
for (i in colnames(noL1HS))
{
try=NULL
if(! all(is.na(noL1HS[i]))){
	try <-testpvalues(i)
}
else{
	print(paste("Skipping", i))
	
}
trypvalues <- cbind(trypvalues, try) 
}

#Configures data for those files that are less than 0.05
tnames<-t(names)
subnames<-tnames[,c(1:10)]
tsubnames<- t(subnames)
colnames(trypvalues)=tsubnames[1,]
trypvalues <0.05
conttry=trypvalues<0.05
conttry2=t(conttry)
colnames(conttry2)= c("pvalues")
dft<-data.frame(conttry2)
dft
subtry<-subset(dft, conttry2==TRUE)
trypvaluesdf<- data.frame(trypvalues)
pttry<- t(trypvaluesdf)
colnames(pttry)= c("dataPV")
data<-subset(pttry, pttry[,1] < 0.05)
testRealPV<- cbind(subtry, data)
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
write.table(testRealPV, file="set of ten examples.txt", sep=" ")

write.table(trypvalues, file="set of pvalues for examples.txt", sep=" ")

#individual test for pvalues
anova(lm(files1$L1HS ~ files2[,c("hsa-mir-122")]))
anova(lm(files1$L1HS ~ files2[,c("hsa-mir-122")]))$Pr[1]

# of cancer/norm
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
L1HS<- read.csv("L1HS.csv")
L1HS
cL1HS<- L1HS[,4]
nL1HS<- L1HS[,3]
namesofpat<- L1HS[,1]
datacL1HS<- data.frame(cL1HS)
datanL1HS<- data.frame(nL1HS)
rownames(datacL1HS)=namesofpat
rownames(datanL1HS)= namesofpat

#find the dividen
DividensofL1HS<- datacL1HS/datanL1HS

#find the log2 of the divident
logofL1HS <- log2(DividensofL1HS)

#read in the normal mirna file
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
normalmirna<- read.table(file="normalMiRNA data ALL.ReadCountPerMillion", header=TRUE)
framenormalmirna<- data.matrix(normalmirna)
framecancermirna<- data.matrix(mirna)
framenormalmirna<- framenormalmirna[,-c(1)]

#find the dividensofmirna
Dividensofmirna<- framecancermirna/framenormalmirna

#find the log2 of mirna
logofmirna<- log2(Dividensofmirna)

#Adds miRna names to data.
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
names<- read.table("miRNA names.txt")
tnames=t(names)
rownames(logofmirna)=tnames
tlogofmirna<- t(logofmirna)

#Change inf to NAs.
 is.na(tlogofmirna) <- do.call(cbind,lapply(tlogofmirna, is.infinite))
 is.na(tlogofmirna) <- do.call(cbind,lapply(tlogofmirna, is.nan))

#make a function for ratio pvalues
ratiopvalues<- function(x)
{
	idea <- lm(logofL1HS[,c(1)] ~ tlogofmirna[,x])
	anova(idea)$Pr[1]
}

#remove columns with all NAs in the rows
falsena <- !apply (is.na(tlogofmirna), 2, all)
framefalsena<- data.frame(falsena)
NONA <- subset(framefalsena, falsena==TRUE)
sublogmirna<- tlogofmirna[,c(rownames(NONA))]

#iterate
presults=NULL
for (i in colnames(sublogmirna))
{
	q = NULL
	
		q <-ratiopvalues(i)

	presults <- cbind(presults, q) 
}
 #Organized results
colnames(presults)=colnames(sublogmirna)
tpresults<- t(presults)
goodpresults = presults< 0.05
tgoodpresults<- t(goodpresults)
frameps <-data.frame(tgoodpresults)
frameps
subcondit<-subset(frameps, tgoodpresults==TRUE)
subpvalues<- data.frame(tpresults)
subpvalues
colnames(subpvalues)= c("dataPV")
data<-subset(subpvalues, subpvalues[,1] < 0.05)
RealPV<- cbind(subcondit, data)
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons")
write.table(RealPV, file="All True Logs of AllCancer.txt", sep=" ")
write.table(tpresults, file="All Raw Logs of AllCancer.txt", sep=" ")

#Go from least to greatest
colnames(tpresults)=c("data")
pvaluesoflogsall<- tpresults[order(-tpresults),,drop=FALSE]
write.table(pvaluesoflogsall, file="All Logs of AllCancer Sorted.txt", sep=" ")


rvalues=NULL
for (i in colnames(sublogmirna))
{
	r = NULL
	
		r <-cor(logofL1HS[,c(1)], tlogofmirna[,i])

	rvalues <- cbind(rvalues, r) 
}

trvalues<- t(rvalues)
colnames(trvalues)=c("data")
rownames(trvalues)=colnames(sublogmirna)
rvaluesoflogsall<- trvalues[order(-trvalues),,drop=FALSE]
results<- cbind(tpresults, trvalues)
colnames(results)=c("pvalues", "rvalues")
sortedbypvalues<- results[order(-results[,c("pvalues")]),, drop=FALSE]
write.table(rvaluesoflogsall, file="All Logs of AllCancer Sorted rvalues.txt", sep=" ")
write.table(sortedbypvalues, file="All Logs of AllCancer with pvalues and rvalues.txt", sep=" ")
write.table(logofL1HS, file="log2 of L1HS.txt", sep=" ")
write.table(tlogofmirna, file="log2 of mirna.txt", sep=" ")