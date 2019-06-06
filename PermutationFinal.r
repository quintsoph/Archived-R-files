#final permutation
#have the one program lead into the other but save the output

args = commandArgs(trailingOnly = TRUE)
alt = args[1]
alt.diff<- read.table(alt, header=TRUE, row.names=1)
nulldata = args[2]

png("")
many.falsenull<- replicate(1000, nulldata[,1])
result<-hist(many.falsenull, xlim=range(-3:4), main="histogram of i")
abline(v=mean(alt.diff[,1]), lwd=2, col="magenta")
abline(v=max(alt.diff[,1]), lwd=2, col="purple")
abline(v=min(alt.diff[,1]), lwd=2, col="darkgreen")
dev.off()
