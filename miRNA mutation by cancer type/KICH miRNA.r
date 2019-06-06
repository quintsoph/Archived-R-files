###KICH miRNA
####read in partial coeficient function
####read in tcompfiles from "By cancer log2 fold with median mutations.r"


##KICH
is.na(tcomp)<- sapply(tcomp, is.infinite)
is.na(tcomp)<- sapply(tcomp, is.nan)
mirnaset<- tcomp[,-c(1,2)]
edits<- cbind(NewNames, stri_sub(NewNames[,1], 9))
CancerPatients<- data.frame(edits[,-c(1)])
test<-intersect(CancerPatients[,2], rownames(mirnaset))
Listmatch<- data.frame(test)
match(test, rownames(mirnaset))
testm<- match(test, CancerPatients[,2])
testm<- testm[!is.na(testm)]
comb<- CancerPatients
rownames(comb)=CancerPatients[,2]
CPatients<- data.frame(comb[testm,])
KICH<- subset(CPatients, CPatients[,1]=="KICH")
KICHcut<- tcomp[which(row.names(tcomp) %in% KICH[,2]),]
ALLL1HSandMUT<- KICHcut
ALLGENES<- data.frame(KICHcut[,-c(1,2)])


#pvalues
testpvalues<- function(x)
{
	print(x)
	value <-lm(KICHcut[,c("L1HS")] ~ ALLGENES[,x] + KICHcut[,c("Mutations")], sep="")
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
		r <-cor(ALLGENES[,i], KICHcut[,c(2)], method="pearson")

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
		r <-cor(ALLGENES[,i], KICHcut[,c(2)], method="spearman")


	spear <- cbind(spear, r) 
}
tspear<- t(spear)
rownames(tspear)=colnames(ALLGENES)


partialpear<- function(x)
{
	if(((nrow(data.frame(na.omit(ALLGENES[,x])))) > 1 ) == TRUE) {
	m<- NULL
	k<- NULL
	m<- cbind( ALLGENES[,x], KICHcut[,c("L1HS")], KICHcut[,c("Mutations")])
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
	m<- cbind( ALLGENES[,x], KICHcut[,c("L1HS")], KICHcut[,c("Mutations")])
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
	value <-lm(KICHcut[,c("L1HS")] ~ ALLGENES[,x] + KICHcut[,c("Mutations")], sep="")
	value$coefficients[2]
}

coefgene=NULL
for (i in colnames(ALLGENES))
{
	rg = NULL
	#If not entire column is NAs, run.
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
	value <-lm(KICHcut[,c("L1HS")] ~ ALLGENES[,x] + KICHcut[,c("Mutations")], sep="")
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
KICHtogether<- merge(tpvalues, tpear, by=0, all=TRUE)
KICHtogether<- merge(KICHtogether, tspear, by.y=0, by.x="Row.names", all=TRUE)
KICHtogether<- merge(KICHtogether, tpartialpearresults, by.y=0, by.x="Row.names", all=TRUE)
KICHtogether<- merge(KICHtogether, tpartialspearresults, by.y=0, by.x="Row.names", all=TRUE)
KICHtogether<- merge(KICHtogether, tcoefgene, by.y=0, by.x="Row.names", all=TRUE)
KICHtogether<- merge(KICHtogether, tcoefmut, by.y=0, by.x="Row.names", all=TRUE)

colnames(KICHtogether)=c("Rownames", "Pvalues", "cor_pearson", "cor_spearman", "partial_pearson", "partial_spearman", "Summary_Gene_Coeff", "Summary_Mutation_Coeff")
KICHtogether <- KICHtogether[order(KICHtogether[,c("Pvalues")]),, drop=FALSE]
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Added Mutation Cov")
write.table(KICHtogether, file="KICHresults.txt", sep="\t", quote=FALSE, row.names=FALSE)
