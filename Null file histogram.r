#congregate the files
Nullfiles <- dir(pattern = "*_stripped.txt", full.names = TRUE)
lst <- lapply(Nullfiles, read.table, header=TRUE, sep="\t")
null <- do.call(rbind, lst)
names<- read.table("dot mirna names.txt", row.names=1)

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

k=names
histo<- function(x)
{
png(paste(rownames(x), ".png"))
r<- as.numeric(x)
hist(r, xlim=range(-3:4), main=rownames(x))
abline(v=names[c(i),], lwd=2, col="purple")
abline(v=rowMeans(tfinal[c(i),]), lwd=2, col="magenta")
dev.off()
}

for (i in rownames(names))
{
	y<- histo(tfinal[c(i),])
}




