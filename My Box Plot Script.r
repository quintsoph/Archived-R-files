#my way to create box plots
#set directory
L1HS<- read.table("L1HSnoout.txt", header=TRUE, row.names=1)
tL1HS<- t(L1HS)
commonIDs<- data.frame(read.table("nooutintersect.txt"))
matched<- match(rownames(L1HS), commonIDs[,1])
matched2<- match(commonIDs[,1], rownames(L1HS))
#Both L1HS and Intersect Patients match!

#read in mutation file
#When there are more than 1 row
#read in mutation file
mutation<- read.table("TDRD1.txt", check.names=FALSE, header=TRUE)
mutation<- mutation[,-1]
mutation<- data.frame(colSums(mutation))
tmutation<- t(mutation)
matched<- match(rownames(L1HS), rownames(mutation))
matched<- matched[!is.na(matched)]
mutmatched<- data.frame(tmutation[,matched])
tmutmatched<- t(mutmatched)
combineddata<- merge(tL1HS, tmutmatched, all=TRUE)
data<- t(combineddata)
colnames(data)= c("Mutations", "L1HS")

#Create BoxPlot
write.table(data, file = "TDRD1.txt", sep = "\t", quote = FALSE)
boxplot(data[,2]~data[,1], main = "TDRD1.txt", xlab = "Number of Mutations", ylab = "log Foldchange L1HS")
