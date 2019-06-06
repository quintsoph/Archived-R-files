#READ
library(methods)
library(stringi)
library(dplyr)
library(matrixcalc)
setwd("/home/quints1/logfoldgenes/matchedgenes")
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
setwd("/home/quints1/logfoldgenes")
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
library(dplyr)
comp<- bind_rows(tL1HS, foldgroup)
genes<- rownames(agg)
rownames(comp)=c("L1HS", genes)
comp2<- bind_rows(tMutations, comp)
rownames(comp2)=c("Mutations", "L1HS", genes)
tcomp<- data.frame(t(comp2))
Genes<- tcomp[,-c(1,2)]
#read in the cancertype file
detach("package:dplyr")
NewNames<- data.frame(read.table(file="SophiaPatCanTypes.txt"))
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

##READ
READ<- subset(CPatients, CPatients[,1]=="READ")
test<- intersect(READ[,2], rownames(Genes))
READcut<- tcomp[which(row.names(tcomp) %in% READ[,2]),]
is.na(READcut)<- sapply(READcut, is.infinite)
ALLL1HSandMUT<- READcut
ALLGENES<- data.frame(READcut[,-c(1,2)])

#pvalues
testpvalues<- function(x)
{
	print(x)
	value <-lm(READcut[,c("L1HS")] ~ ALLGENES[,x] + READcut[,c("Mutations")], sep="")
	anova(value)$Pr[1]
}

#Iterates over miRnas and computes p-values.
pvalues=NULL
for (i in colnames(ALLGENES))
{
	p = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(ALLGENES[i]))){
		p <-testpvalues(i)
	}
	else{
		print(paste("Skipping", i))
	}

	pvalues <- cbind(pvalues, p) 
}
pvalues<- data.frame(pvalues)
tpvalues<- t(pvalues)
skip=NULL
for (i in colnames(ALLGENES))
{
	skipped=NULL
	if(!all(is.na(ALLGENES[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}
skip<- t(skip)
rownames(tpvalues)= skip[,1]

#rvalues
#Pearson
pear=NULL
for (i in colnames(ALLGENES))
{
	r = NULL
		#cor(x,y)
		#x is the genes and y is L1HS
		r <-cor(ALLGENES[,i], READcut[,c(2)], method="pearson")

	pear <- cbind(pear, r) 
}
tpear<- t(pear)
rownames(tpear)=colnames(ALLGENES)

#Spearman
spear=NULL
for (i in colnames(ALLGENES))
{
	r = NULL
		#cor(x,y)
		r <-cor(ALLGENES[,i], READcut[,c(2)], method="spearman")


	spear <- cbind(spear, r) 
}
tspear<- t(spear)
rownames(tspear)=colnames(ALLGENES)

#partial for pearson
pcor.test <- function(x,y,z,use="mat",method="p",na.rm=T){
	# The partial correlation coefficient between x and y given z
	#
	# pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
	#
	# x and y should be vectors
	#
	# z can be either a vector or a matrix
	#
	# use: There are two methods to calculate the partial correlation coefficient.
	#	 One is by using variance-covariance matrix ("mat") and the other is by using recursive formula ("rec").
	#	 Default is "mat".
	#
	# method: There are three ways to calculate the correlation coefficient, 
	#	    which are Pearson's ("p"), Spearman's ("s"), and Kendall's ("k") methods.
	# 	    The last two methods which are Spearman's and Kendall's coefficient are based on the non-parametric analysis.
	#	    Default is "p".
	#
	# na.rm: If na.rm is T, then all the missing samples are deleted from the whole dataset, which is (x,y,z).
	#        If not, the missing samples will be removed just when the correlation coefficient is calculated.
	#	   However, the number of samples for the p-value is the number of samples after removing 
	#	   all the missing samples from the whole dataset.
	#	   Default is "T".

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(use == "mat"){
		p.use <- "Var-Cov matrix"
		pcor = pcor.mat(x,y,z,method=method,na.rm=na.rm)
	}else if(use == "rec"){
		p.use <- "Recursive formula"
		pcor = pcor.rec(x,y,z,method=method,na.rm=na.rm)
	}else{
		stop("\'use\' should be either \"rec\" or \"mat\"!\n")
	}

	# print the method
	if(gregexpr("p",method)[[1]][1] == 1){
		p.method <- "Pearson"
	}else if(gregexpr("s",method)[[1]][1] == 1){
		p.method <- "Spearman"
	}else if(gregexpr("k",method)[[1]][1] == 1){
		p.method <- "Kendall"
	}else{
		stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
	}

	# sample number
	n <- dim(na.omit(data.frame(x,y,z)))[1]
	
	# given variables' number
	gn <- dim(z)[2]

	# p-value
	if(p.method == "Kendall"){
		statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
		p.value <- 2*pnorm(-abs(statistic))

	}else{
		statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  		p.value <- 2*pnorm(-abs(statistic))
	}

	data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use)
}			

# By using var-cov matrix
pcor.mat <- function(x,y,z,method="p",na.rm=T){

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z)

	if(na.rm == T){
		data = na.omit(data)
	}

	xdata <- na.omit(data.frame(data[,c(1,2)]))
	Sxx <- cov(xdata,xdata,m=method)

	xzdata <- na.omit(data)
	xdata <- data.frame(xzdata[,c(1,2)])
	zdata <- data.frame(xzdata[,-c(1,2)])
	Sxz <- cov(xdata,zdata,m=method)

	zdata <- na.omit(data.frame(data[,-c(1,2)]))
	Szz <- cov(zdata,zdata,m=method)

	# is Szz positive definite?
	zz.ev <- eigen(Szz)$values
	if(min(zz.ev)[1]<0){
		stop("\'Szz\' is not positive definite!\n")
	}

	# partial correlation
	Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
	
	rxx.z <- cov2cor(Sxx.z)[1,2]

	rxx.z
}

# By using recursive formula
pcor.rec <- function(x,y,z,method="p",na.rm=T){
	# 

	x <- c(x)
	y <- c(y)
	z <- as.data.frame(z)

	if(dim(z)[2] == 0){
		stop("There should be given data\n")
	}

	data <- data.frame(x,y,z)

	if(na.rm == T){
		data = na.omit(data)
	}

	# recursive formula
	if(dim(z)[2] == 1){
		tdata <- na.omit(data.frame(data[,1],data[,2]))
		rxy <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,1],data[,-c(1,2)]))
		rxz <- cor(tdata[,1],tdata[,2],m=method)

		tdata <- na.omit(data.frame(data[,2],data[,-c(1,2)]))
		ryz <- cor(tdata[,1],tdata[,2],m=method)

		rxy.z <- (rxy - rxz*ryz)/( sqrt(1-rxz^2)*sqrt(1-ryz^2) )
		
		return(rxy.z)
	}else{
		x <- c(data[,1])
		y <- c(data[,2])
		z0 <- c(data[,3])
		zc <- as.data.frame(data[,-c(1,2,3)])

		rxy.zc <- pcor.rec(x,y,zc,method=method,na.rm=na.rm)
		rxz0.zc <- pcor.rec(x,z0,zc,method=method,na.rm=na.rm)
		ryz0.zc <- pcor.rec(y,z0,zc,method=method,na.rm=na.rm)
		
		rxy.z <- (rxy.zc - rxz0.zc*ryz0.zc)/( sqrt(1-rxz0.zc^2)*sqrt(1-ryz0.zc^2) )
		return(rxy.z)
	}			
}	

partialpear<- function(x)
{
	if(((nrow(data.frame(na.omit(ALLGENES[,x])))) > 1 ) == TRUE) {
	m<- NULL
	k<- NULL
	m<- cbind( ALLGENES[,x], READcut[,c("L1HS")], READcut[,c("Mutations")])
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

trycatcher<- function(x)
{
	tryCatch({x},
	error=function(e){cat("ERROR :", conditionMessage(e), "\n")
	cnd<- conditionMessage(e)
	return(cnd)})
}
partialpearresults=NULL
for (i in colnames(ALLGENES))
{
	r = NULL
	if(! all(is.na(ALLGENES[i]))){
	r<- trycatcher(partialpear(i))
	}
	else{
		print(paste("Skipping", i))
	}
	partialpearresults <- cbind(partialpearresults, r) 
	}


tpartialpearresults<- t(partialpearresults)
skip=NULL
for (i in colnames(ALLGENES))
{
	skipped=NULL
	if(!all(is.na(ALLGENES[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}
skip<- t(skip)
rownames(tpartialpearresults)= skip[,1]



#partial for spearman
partialspear<- function(x)
{
	if(((nrow(data.frame(na.omit(ALLGENES[,x])))) > 1 ) == TRUE) {
	m<- NULL
	k<- NULL
	m<- cbind( ALLGENES[,x], READcut[,c("L1HS")], READcut[,c("Mutations")])
	k<- as.matrix(m[complete.cases(m),])
		if(is.square.matrix(k) == TRUE){
			if(ifelse(det(k) == 0, FALSE, TRUE) == TRUE){
				r <-pcor.test(k[,1], k[,2], k[,3], method="s", na.rm=TRUE)
				rcors<- r$estimate[1]
				return(rcors)}
			else{return(NA)}}
		else{
			r <-pcor.test(k[,1], k[,2], k[,3], method="s", na.rm=TRUE)
			rcors<- r$estimate[1]
			return(rcors)}}
	else{
	return(NA)
		}
}

partialspearresults=NULL
for (i in colnames(ALLGENES))
{
	r = NULL
	if(! all(is.na(ALLGENES[i]))){
	r<- trycatcher(partialspear(i))
	}
	else{
		print(paste("Skipping", i))
	}
	partialspearresults <- cbind(partialspearresults, r) 
}

tpartialspearresults<- t(partialspearresults)
skip=NULL
for (i in colnames(ALLGENES))
{
	skipped=NULL
	if(!all(is.na(ALLGENES[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}
skip<- t(skip)
rownames(tpartialspearresults)= skip[,1]


#Coefficient of the genes
testcoefgene<- function(x)
{
	print(x)
	value <-lm(READcut[,c("L1HS")] ~ ALLGENES[,x] + READcut[,c("Mutations")], sep="")
	value$coefficients[2]
}

coefgene=NULL
for (i in colnames(ALLGENES))
{
	rg = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(ALLGENES[i]))){
		rg <-testcoefgene(i)
	}
	else{
		print(paste("Skipping", i))
	}

	coefgene <- cbind(coefgene, rg) 
}
coefgene<- data.frame(coefgene)
tcoefgene<- t(coefgene)
skip=NULL
for (i in colnames(ALLGENES))
{
	skipped=NULL
	if(!all(is.na(ALLGENES[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}
skip<- t(skip)
rownames(tcoefgene)= skip[,1]

#coefficient of the Mutations
testcoefmut<- function(x)
{
	print(x)
	value <-lm(READcut[,c("L1HS")] ~ ALLGENES[,x] + READcut[,c("Mutations")], sep="")
	value$coefficients[3]
}

coefmut=NULL
for (i in colnames(ALLGENES))
{
	rmm = NULL
	#If not entire row is NAs, run.
	if(! all(is.na(ALLGENES[i]))){
		rmm <-testcoefmut(i)
	}
	else{
		print(paste("Skipping", i))
	}

	coefmut <- cbind(coefmut, rmm) 
}
coefmut<- data.frame(coefmut)
tcoefmut<- t(coefmut)
skip=NULL
for (i in colnames(ALLGENES))
{
	skipped=NULL
	if(!all(is.na(ALLGENES[i]))){
	skipped<- print(paste(i))
	}
	skip<- cbind(skip, skipped)
	}
skip<- t(skip)
rownames(tcoefmut)= skip[,1]

tpvalues<- data.frame(tpvalues)
tpear<- data.frame(tpear)
tspear<- data.frame(tspear)
tpartialpearresults<- data.frame(tpartialpearresults)
tpartialspearresults<- data.frame(tpartialspearresults)
tcoefgene<- data.frame(tcoefgene)
tcoefmut<- data.frame(tcoefmut)
READtogether<- merge(tpvalues, tpear, by=0, all=TRUE)
READtogether<- merge(READtogether, tspear, by.y=0, by.x="Row.names", all=TRUE)
READtogether<- merge(READtogether, tpartialpearresults, by.y=0, by.x="Row.names", all=TRUE)
READtogether<- merge(READtogether, tpartialspearresults, by.y=0, by.x="Row.names", all=TRUE)
READtogether<- merge(READtogether, tcoefgene, by.y=0, by.x="Row.names", all=TRUE)
READtogether<- merge(READtogether, tcoefmut, by.y=0, by.x="Row.names", all=TRUE)

colnames(READtogether)=c("Rownames", "Pvalues", "cor_pearson", "cor_spearman", "partial_pearson", "partial_spearman", "Summary_Gene_Coeff", "Summary_Mutation_Coeff")
READtogether <- READtogether[order(READtogether[,c("Pvalues")]),, drop=FALSE]
setwd("/home/quints1/logfoldgenes/matchedgenes/CancerResults")
write.table(READtogether, file="READresults.txt", sep="\t", quote=FALSE, row.names=FALSE)
