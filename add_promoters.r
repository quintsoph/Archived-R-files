####add promoters to the genome range
##Subtract 500 from + on the left hand
##Add 500 to - on the right hand
setwd("/home/quints1/Patient_Names")
files<- read.table("hg19_rmsk_L1HS.gtf")

Addpro<- function(x)
{
	if(files[x, 7] == "+"){
	y<- files[x, 4] - 500 
	result<- cbind(files[x ,c(1,2,3)], y, files[x ,c(5:18)])
	colnames(result)=c("")}
	else {y<- files[x, 5] + 500
		  result<- cbind(files[x ,c(1:4)], y, files[x ,c(6:18)])
		  colnames(result)=c("")}
	return(result)
	}

New=NULL
for (row in 1:1544)
{
	k<- Addpro(row)
	New<- rbind(New, k)
}

setwd("/home/quints1/Patient_Names")
write.table(New, file="hg19_edit_L1HS_promoters.gtf", sep="\t", quote=FALSE, row.names=FALSE)
