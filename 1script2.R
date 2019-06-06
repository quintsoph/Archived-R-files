#filenametxt is the folder where they are contained
#cbind keeps adding them to the file
probabilities <- c()
#filenametxt is a file with all the headers
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells")
fnames <- read.table(filenametxt)
for (i in 1:nrow(fnames))
{
p <- testcancer(fnames[i])
cbind(probs, p)
}
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells/txtfolder")
testcancer<- list.files(pattern="*.txt")
{
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/CSV file")
G<-read.csv("averagesofmirnaHG19.csv", stringsAsFactors=FALSE)
G<-G[-1]
Averages <- rowMeans(G)
Standard.Deviation<- apply(G, 1, sd)
#Shows the normal distribution for all the rows
#dnorm(what you are comparing, Average, with the standard deviation set)
#below is for calling just one column
setwd("E:/Bioinformatics Lab/Cancer Data/miRNASeq/Colon Cancer/Cancer Cells/txtfolder")
H<- read.table(list.file(pattern="*.txt"))
cancer.set <- H[,3]
probs <- dnorm(cancer.set, Averages, Standard.Deviation)
probs
}
sum(probabilities<.00005)
