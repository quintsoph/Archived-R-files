#filenametxt is the folder where they are contained
#cbind keeps adding them to the file
probabilities <- c()
#filenametxt is a file with all the headers


testcancer <- function(fnames)
{
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer")
G<-read.csv("Final2.csv", stringsAsFactors=FALSE)
Averages <- G[,c("Averages")]
Standard.Deviation<- G[,c("Standard.Deviation")]
#Shows the normal distribution for all the rows
#dnorm(what you are comparing, Average, with the standard deviation set)
#below is for calling just one column
probs <- dnorm(fnames, Averages, Standard.Deviation)
return(probs)
}

#i=1 for testing loop before actually doing loop

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells/txtfolder")
files <- dir(pattern = "*.txt", full.names = TRUE)
lst <- lapply(files, read.table, header=TRUE, sep='')
fnames<- sapply(lst, '[[',3)
d=NULL
for (i in 1:ncol(fnames))
{
p <-testcancer(fnames[i,])
d <- cbind(d, p) 
}

setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells/BCGSC__IlluminaHiSeq_miRNASeq")
filesS2 <- dir(pattern = "*.txt", full.names = TRUE)
lstS2 <- lapply(filesS2, read.table, header=TRUE, sep='')
fnamesS2<- sapply(lstS2, '[[',3)
e=NULL
for (i in 1:ncol(fnamesS2))
{
q <-testcancer(fnames[i,])
e <- cbind(e, q) 
}
 
sum(d<.00005)
number/1000=
with the first sample location 137.623 killomiRNA
with the second sample location 130.606 killomiRNA
with both sample locations 260.061 killomiRNA



#junk data
colnames(fnames)<- c("2670", "2674", "2677", "2678")
files <- dir(pattern = "*.txt", full.names = TRUE)
files <- dir(pattern = "TCGA-A6-2670-01A-02T-0822-13.hg19.mirna.quantification.txt", full.names = TRUE)
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells")
write(p, "normaldtest.txt", append=TRUE, sep=",")
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells")
d <-write(p, file="reald.txt", append=TRUE, sep=",")
nextcol <-  data.frame(p)
colnames(nextcol) <- c(paste("col", i, sep="")) # rename the comlum
df <- cbind(df, nextcol)

#
#adding the column names seperately
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells")
pnames<- read.table(file.choose())
normalfile<- read.table("normald.txt")
final <- rbind(pnames, d)
#


setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells")
write.table(e, file="realdistS2.txt", sep=",")

1st one is cancer
the 2nd one is normal
normal deviation for the TE data?
