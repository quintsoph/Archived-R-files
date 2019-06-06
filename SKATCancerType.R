library(SKAT)

allZero = function(x) length(x[which(x == 0)]) == nrow(x)

args = commandArgs(trailingOnly = TRUE)
intersectFilename = args[1]
phenotypeFilename = args[2]
mutationsPath = args[3]
outputfile = args[4]

write(paste("Gene", "Chromosome", "pvalue", sep="\t"), file=outputfile)

#Read in the TE file as a matrix. File should be formatted beginning with 3 columns: patient ID, cancer type, and TE counts.
phenotypes = as.matrix(read.table(phenotypeFilename, sep="\t"))
TEDcounts = phenotypes[2:nrow(phenotypes),c(1,3)]
commonIDs = data.frame(read.table(intersectFilename))

TEDhash = new.env()
for (i in seq(nrow(TEDcounts))) {
  TEDhash[[ TEDcounts[i,1] ]] = as.numeric(TEDcounts[i,2])
}
#TEs will be the dependent variable
TEs = rbind()
for (i in seq(nrow(commonIDs))) {
  TEs = rbind(TEs, TEDhash[[ toString(commonIDs[i,]) ]])
}
#Create the covariate matrix with cancer types
cancertypes = phenotypes[2:nrow(phenotypes),2]
columnnames = sort(unique(cancertypes))

covariate = as.data.frame(matrix(0,ncol=length(columnnames), nrow=(nrow(commonIDs))))
colnames(covariate) = columnnames
rownames(covariate) = commonIDs[,1]

for (i in seq(2, nrow(phenotypes))) {
    id = phenotypes[i,1]
    canc = phenotypes[i,2]
    if (id %in% row.names(covariate)) {
	#Can change 1 
        covariate[id,canc] = 1
    }
}
covariate = as.matrix(covariate)

#Read in the mutations file by file (by gene)
for (mutFilename in Sys.glob(paste(mutationsPath, "/*/*.txt", sep=""))) {
  ID = strsplit(strsplit(mutFilename, "/")[[1]][3], '[.]')[[1]][1]
  chromosome = strsplit(mutFilename, "/")[[1]][2]
  mutations = t(as.matrix(read.table(mutFilename, sep="\t", colClasses='character')))
  mutations = mutations[2:nrow(mutations),]
  mutHash = new.env()
  for (i in seq(nrow(mutations))) {
    mutHash[[ mutations[i,1] ]] = as.numeric(mutations[i,2:ncol(mutations)])
  }
#muts will be the independent variable
  muts = rbind()
  for (i in seq(nrow(commonIDs))) {
    muts = rbind(muts, mutHash[[ toString(commonIDs[i,]) ]])
  }
#Build the null model for SKAT, run SKAT, and write the pvalues to a file
  if (!allZero(muts)) {
    nullModel = SKAT_Null_Model(TEs~covariate, out_type="C")
    write(paste(ID, chromosome, SKAT(muts, nullModel)$p.value, sep="\t"), file=outputfile, append=TRUE)
  }
  else {
    write(paste(ID, chromosome, "NA", sep="\t"), file=outputfile, append=TRUE)
  }
}
