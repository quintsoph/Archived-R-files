#Read in Nicky Files
setwd("E:/Bioinformatics Lab/Cancer Data/Mutations/Processed Data Issue")
compare <- read.table("comparemutation.txt", header=TRUE)
counts<- read.table("mutationcounts.txt", header=TRUE)
totaldap<- read.table("totalmutationsdaphnie.txt", header=TRUE)
