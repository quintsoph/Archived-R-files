#read in files
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
Nullfiles <- dir(pattern = "*_stripped.txt", full.names = TRUE)
lst <- lapply(Nullfiles, read.table, header=TRUE, sep="\t")
null <- do.call(cbind, lst)
names<- read.table("REC ONLY NEW.txt", row.names=1)

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
tfinals<- tfinal * -1

#find pvalue
#number
less<- function(x)
{
	Nulls<-as.numeric(tfinals[c(x),])
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
	y<- (tlessthans[c(i),])/1000
	pvalue<- cbind(pvalue, y)
}


tpvalue<- t(pvalue)
rownames(tpvalue)=rownames(tlessthans)
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
write.table(tpvalue, file="NULLS p-values 1G NEW.txt", sep="\t")
