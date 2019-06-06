#BRCA
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2ratios.txt"))
BRCAorder <- BRCA[order(-BRCA[,c("pvalues")]),, drop=FALSE]

#rid NA rows
BRCAorderp<- BRCAorder[,1]
BRCAorderNA<- BRCAorderp[1:554]
#rank
BRCAorderrank<- data.frame(rank(BRCAorderNA, na.last = TRUE, ties.method=c("min")))
rownames(BRCAorderrank)= rownames(BRCAorder[c(1:554),])

#find rr
#554 is the number of mirna that is not NA 
#function for rr in BRCA
rrfunction<- function(x)
{
	(x/554)-(1/(2*554))
}

#find Hknot
Hknot<- function(x)
{
	-2*(sum(log(x)))
}

BRCAH<- Hknot(rrBRCA)

#create an Hknot data frame
Hknotresults=NULL
Hknotresults <-cbind(Hknotresults, BLCAH)

#rrfunction
rrBRCA <- rrfunction(BRCAorderrank)
BRCAH<- Hknot(rrBRCA)
Hknotresults<- cbind(Hknotresults, BRCAH)

