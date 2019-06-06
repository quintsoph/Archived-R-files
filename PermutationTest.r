#premutation test

alt.diff<- data.frame(read.table(file.choose(), header=TRUE, row.names=1))
nulldata<- read.table("Null1new.txt", row.names=1)

args = commandArgs(trailingOnly = TRUE)
alt.diff = args[1]
nulldata = args[2]
outputfile = args[3]

many.falsenull<- replicate(1000, nulldata[,1])
test<- hist(many.falsenull, xlim=range(-3:4), main="histogram of Nul1")
abline(v=mean(alt.diff[,1]), lwd=2, col="magenta")
abline(v=max(alt.diff[,1]), lwd=2, col="purple")
abline(v=min(alt.diff[,1]), lwd=2, col="darkgreen")



for i in Sys.glob(paste(nulldata, "/*.txt", sep=""))
{
hist$i
many.falsenull<- replicate(1000, i[,1])
hist(many.falsenull, xlim=range(-3:3), main="histogram of i")
abline(v=mean(alt.diff[,1]), lwd=2, col="magenta")
abline(v=max(alt.diff[,1]), lwd=2, col="purple")
abline(v=min(alt.diff[,1]), lwd=2, col="darkgreen")
}


nulldata<- setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/test")
for (i in Sys.glob(paste(nulldata, "/*.txt", sep="")))
{
hist[i]=NULL
many.falsenull<- replicate(1000, i[,1])
result<-hist(many.falsenull, xlim=range(-3:3), main="histogram of i")
abline(v=mean(alt.diff[,1]), lwd=2, col="magenta")
abline(v=max(alt.diff[,1]), lwd=2, col="purple")
abline(v=min(alt.diff[,1]), lwd=2, col="darkgreen")
}
