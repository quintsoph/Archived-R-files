args = commandArgs(trailingOnly = TRUE)
#intersectFilename = args[1]
#phenotypeFilename = args[2]
#mutationsPath = args[3]
#outputfile = args[4]

TEDcounts = as.matrix(read.table("L1HSData2.txt", sep="\t"))
TEDcounts = TEDcounts[2:nrow(TEDcounts),]
commonIDs = data.frame(read.table("newintersect.txt"))

TEDhash = new.env()
for (i in seq(nrow(TEDcounts))) {
  TEDhash[[ TEDcounts[i,1] ]] = as.numeric(TEDcounts[i,2])
}

TEs = rbind()
for (i in seq(nrow(commonIDs))) {
  TEs = rbind(TEs, TEDhash[[ toString(commonIDs[i,]) ]])
}

mutations = t(as.matrix(read.table("CCDC162P.txt", sep="\t")))
mutations = mutations[2:nrow(mutations),]
mutHash = new.env()
for (i in seq(nrow(mutations))) {
  mutHash[[ mutations[i,1] ]] = as.numeric(mutations[i,2:ncol(mutations)])
}

muts = rbind()
for (i in seq(nrow(commonIDs))) {
  muts = rbind(muts, mutHash[[ toString(commonIDs[i,]) ]])
}

mutstogether = as.matrix(rowSums(muts[,1:ncol(muts)]))
data = cbind(TEs, mutstogether)
colnames(data) = c("L1HS ratio", "Mutations")
data
write.table(data, file = "Result Files/C10orf99.txt", sep = "\t", quote = FALSE, row.names = FALSE)
boxplot(TEs~mutstogether, main = "CCDC162P", xlab = "Number of Mutations", ylab = "log Foldchange L1HS")
