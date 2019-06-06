#read in files
setwd("/home/quints1/workDrive/sophia/1Gpvaltest")
Nullfiles <- dir(pattern = "*_stripped.txt", full.names = TRUE)
lst <- lapply(Nullfiles, read.table, header=TRUE, sep="\t")
null <- do.call(rbind, lst)
names<- read.table("dot mirna names.txt", row.names=1)

#Format the null rows
microcut<- function(x)
{
	subs=NULL
	subs<- subset(null, null[,1] == x)
	final<- data.frame(subs[,2])
	colnames(final)= c(x)
	return(final)
}


final=data.frame(matrix(NA, nrow = 4, ncol = 0))
for (i in rownames(names))
{
	y=NULL
	y<-microcut(i)
	final <- cbind(final, y)

}

tfinal<- data.frame(t(final))


#find pvalue
#number
less<- function(x)
{
	Nulls<-as.numeric(tfinal[c(x),])
	k<- as.numeric(names[c(x),])
	y<- sum(Nulls < k)
	return(y)
}

lessthans=NULL
for (i in rownames(names))
{
	st=NULL
	st<-less(i)
	lessthans<- cbind(lessthans, st)
}

colnames(lessthans)=rownames(names)
tlessthans<- t(lessthans)


pvalue=NULL
for (i in rownames(tlessthans))
{
	y=NULL
	y<- (tlessthans[c(i),])/100000
	pvalue<- cbind(pvalue, y)
}


tpvalue<- t(pvalue)
rownames(tpvalue)=rownames(tlessthans)
write.table(tpvalue, sep="\t", quote=FALSE)