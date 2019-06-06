#read in patient names
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
names<- read.table("patientsofgenes.txt", header=TRUE)
library(stringi)
ednames<- data.frame(stri_sub(names[,1],9,-26))

#Read in L1HS
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
#format cancer
newL1HScancer<- data.frame(read.table(file="L1HS.Cancer.BaseMeansB", header=TRUE, row.names=1))
newL1HScancer<- data.matrix(newL1HScancer)
tnewL1HScancer<- t(newL1HScancer)
#format normal
newL1HSnormal<- data.frame(read.table(file="L1HS.Norma.BaseMeansA", header=TRUE, row.names=1))
newL1HSnormal<- data.matrix(newL1HSnormal)
tnewL1HSnormal<- t(newL1HSnormal)
dividenL1HS<- newL1HScancer/newL1HSnormal
foldedL1HS<- log2(dividenL1HS)
tfoldedL1HS<- data.frame(t(foldedL1HS))
edL1HS<- data.frame(stri_sub(rownames(tfoldedL1HS),9))
rownames(tfoldedL1HS)=edL1HS[,1]

#Get the patients to lineup
matchesL1HS<- intersect(ednames[,1], rownames(tfoldedL1HS))
setwd("F:/Bioinformatics Lab/Cancer Data/Mutations/Final Results")
MedMut<- data.frame(read.table("Mutationsnoout.txt", header=TRUE, row.names=1))
#L1HS matches
totalmatches<- intersect(matchesL1HS, rownames(MedMut))
testm<-match(totalmatches, rownames(tfoldedL1HS))
testm<- testm[!is.na(testm)]
comb<- cbind(rownames(tfoldedL1HS), tfoldedL1HS)
resultnamesL1HS<- data.frame(comb[testm,])
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
write.table(resultnamesL1HS, "L1HS_Gene_Matches.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
#Mutation matches
testm1<-match(totalmatches, rownames(MedMut))
testm1<- testm1[!is.na(testm1)]
comb<- cbind(rownames(MedMut), MedMut)
resultnamesMUT<- data.frame(comb[testm1,])
write.table(resultnamesMUT, "Mutation_Matches.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
#Names
Matchednames<- data.frame(totalmatches)
rownames(Matchednames)=Matchednames[,1]
filenames<- stri_join(Matchednames[,1], "-01A_gene_TE_analysis.txt", collapse=NULL)
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
#Combine L1HS and the mutations
testm2<- match(resultnamesMUT[,1], resultnamesL1HS[,1])
together<- cbind(resultnamesMUT, resultnamesL1HS)

#grep(*filenames, names)
names<- as.vector(names[,1])
filenames<- as.vector(filenames[,1])
total=NULL
for (i in filenames)
{
	y<- NULL
	y<-grep(paste0('*?',i), names, value=TRUE)
	total<- rbind(total, y)
}

setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
write.table(total, "Complete_Files_Match.txt", sep=" ", quote=FALSE, col.names=FALSE, row.names=FALSE)




#For the final program
#Get the log2 fold of the patients for all their genes by CANCER
library(stringi)
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis/Genes")
mirnafiles <- dir(pattern = "*.txt", full.names = TRUE)

agg<-data.frame(read.table("TCGA-AX-A05Y-01A_gene_TE_analysis.txt", header=TRUE))
agg<- agg[,c(1,6)]
for (file in mirnafiles)
{
	z<- NULL
	y<- NULL
	z<- data.frame(read.table(file, header=TRUE))
	y<- z[,c(1,6)]
	colnames(y)=c("id", file)
	agg<- merge(agg, y, all=TRUE)
}

#read in L1HS matches
setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
L1HSfold<- read.table("L1HS_Gene_Matches.txt", header=FALSE, row.names=1)
tL1HS<- data.frame(t(L1HSfold))
colnames(tL1HS)=rownames(L1HSfold)
#read in the Mutation file
Mutations<- data.frame(read.table("Mutation_Matches.txt", header=FALSE, row.names=1))
tMutations<- data.frame(t(Mutations))
colnames(tMutations)=rownames(Mutations)
rownames(agg)=agg[,1]
foldgroup<- data.frame(agg[,-c(1,2)])
colnames(foldgroup)=stri_sub(colnames(foldgroup),11,-26)
tfoldgroup<- t(foldgroup)
if(!require("ppcor")){
  install.packages("ppcor", repos='http://cran.us.r-project.org')
  library(ppcor)
}
library(dplyr)
comp<- bind_rows(tL1HS, foldgroup)
genes<- rownames(agg)
rownames(comp)=c("L1HS", genes)
comp2<- bind_rows(tMutations, comp)
rownames(comp2)=c("Mutations", "L1HS", genes)
tcomp<- data.frame(t(comp2))
Genes<- tcomp[,-c(1,2)]
#read in the cancertype file
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/organized by cancer types")
AllCancers <-data.frame(read.table(file="Cancer Types Sorted"))

setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
NewNames<- data.frame(read.table(file="SophiapatCanTypes"))
NewNames<- t(NewNames)
edits<- cbind(NewNames, stri_sub(NewNames[,1], 9))
CancerPatients<- edits[,-c(1)]
test<-intersect(CancerPatients[,2], rownames(Genes))
match(test, rownames(Genes))
testm<- match(test, CancerPatients[,2])
testm<- testm[!is.na(testm)]
comb<- CancerPatients
rownames(comb)=CancerPatients[,2]
CPatients<- data.frame(comb[testm,])

practice<- tcomp[c("A05Y", "5768", "5987", "AA2U"),]
practice2<- lm(practice[,c("L1HS")] ~ practice[,c("AADAT")] + practice[,c("Mutations")], sep="")
practice2$coefficients[3]

library(ppcor)
is.na(practice)<- sapply(practice, is.infinite)
practice1<- data.matrix(practice)
practice2<- data.matrix(practice[,-c(1,2)])
partialpearresults=NULL
for (i in colnames(practice1[3:ncol]))
{
	r = NULL
	if(! all(is.na(practice1[i]))){
		r <-pcor.test(practice1[,c("L1HS")], practice2[,c(i)], practice1[,c("Mutations")])
}
	else{
		print(paste("Skipping", i))
	}
	partialpearresults <- cbind(partialpearresults, r) 
}
tpartialpearresults<- t(partialpearresults)
skip=NULL
for (i in colnames(practice))
{
	skipped=NULL
	if(!all(is.na(practice[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}
skip<- t(skip)
rownames(tpartialpear)= skip[,1]



library(ppcor)
library(matrixcalc)
is.na(practice)<- sapply(practice, is.infinite)
practice1<- data.frame(practice)
practice2<- data.frame(practice[,-c(1,2)])
partialpear<- function(x)
{
	if(((nrow(data.frame(na.omit(practice2[,x])))) > 1 ) == TRUE) {
	m<- NULL
	k<- NULL
	m<- cbind(practice1[,c("L1HS")], practice2[,x], practice1[,c("Mutations")])
	k<- m[complete.cases(m),]
			r <-pcor.test(k[,1], k[,2], k[,3], method="p", na.rm=TRUE)
			rcors<- r$estimate[1]
			return(rcors)}
	else{
	return(NA)
		}
}

partialpear("ZG16")
partialpear("MYH15")
partialpear("ZNFX1")
m<- cbind(practice1[,c("L1HS")], practice1[,c("ZNFX1")], practice1[,c("Mutations")])
	k<- m[complete.cases(m),]

partialpearresults=NULL
for (i in colnames(practice2))
{
	r = NULL
	if(! all(is.na(practice2[i]))){
	r<- partialpear(i)
	}
	else{
		print(paste("Skipping", i))
	}
	partialpearresults <- cbind(partialpearresults, r) 
}
tpartialpearresults<- t(partialpearresults)
skip=NULL
for (i in colnames(practice2))
{
	skipped=NULL
	if(!all(is.na(practice2[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}
skip<- t(skip)
rownames(tpartialpearresults)= skip[,1]



#read in the partial coefficient function script
library(matrixcalc)
is.na(practice)<- sapply(practice, is.infinite)
practice1<- data.frame(practice)
practice2<- data.frame(practice[,-c(1,2)])
partialpear<- function(x)
{
	if(((nrow(data.frame(na.omit(practice2[,x])))) > 1 ) == TRUE) {
	m<- NULL
	k<- NULL
	m<- cbind(practice1[,c("L1HS")], practice2[,x], practice1[,c("Mutations")])
	k<- as.matrix(m[complete.cases(m),])
		if(is.square.matrix(k) == TRUE){
			if(ifelse(det(k) == 0, FALSE, TRUE) == TRUE){
				
				r <-pcor.test(k[,1], k[,2], k[,3], method="p", na.rm=TRUE)
				rcors<- r$estimate[1]
				return(rcors)}
				
			else{return(NA)}}
		else{
			r <-pcor.test(k[,1], k[,2], k[,3], method="p", na.rm=TRUE)
			rcors<- r$estimate[1]
			return(rcors)}}
	else{
	return(NA)
		}
}
partialpear("ZG16")
partialpear("MYH15")

#save just in case
partialpearresults=NULL
for (i in colnames(ALLGENES))
{ tryCatch({
	r = NULL
	if(! all(is.na(ALLGENES[i]))){
	r<- partialpear(i)
	}
	else{
		print(paste("Skipping", i))
	}
	partialpearresults <- cbind(partialpearresults, r) 
	}, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
}



test1k<- tpartialpearresults
test2<- tpartialpearresults

setwd("F:/Bioinformatics Lab/Cancer Data/Gene L1HS analysis")
write.table(example, file="match example.txt", quote=FALSE, sep="\t", row.names=FALSE)