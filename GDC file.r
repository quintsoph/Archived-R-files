setwd("/home/quints1/Patient_Names")
trycatcher<- function(x)
{
	tryCatch({x},
	error=function(e){cat("ERROR :", conditionMessage(e), "\n")
	cnd<- conditionMessage(e)
	return(cnd)})
}

BLCA<- read.table("BLCA patients.txt")
rownames(BLCA)=BLCA[,1]
for (i in rownames(BLCA))
{	print(i)
	querytest<- GDCquery(project="TCGA-BLCA", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }

BRCA<- read.table("BRCA patients.txt")
rownames(BRCA)=BRCA[,1]
for (i in rownames(BRCA))
{	print(i)
	querytest<- GDCquery(project="TCGA-BRCA", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
	
COAD<- read.table("COAD patients.txt")
rownames(COAD)=COAD[,1]
for (i in rownames(COAD))
{	print(i)
	querytest<- GDCquery(project="TCGA-COAD", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
	
ESCA<- read.table("ESCA patients.txt")
rownames(ESCA)=ESCA[,1]
for (i in rownames(ESCA))
{	print(i)
	querytest<- GDCquery(project="TCGA-ESCA", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
	
HNSC<- read.table("HNSC patients.txt")
rownames(HNSC)=HNSC[,1]
for (i in rownames(HNSC))
{	print(i)
	querytest<- GDCquery(project="TCGA-HNSC", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
	
KIRC<- read.table("KIRC patients.txt")
rownames(KIRC)=KIRC[,1]
for (i in rownames(KIRC))
{	print(i)
	querytest<- GDCquery(project="TCGA-KIRC", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }	
	
KIRP<- read.table("KIRP patients.txt")
rownames(KIRP)=KIRP[,1]
for (i in rownames(KIRP))
{	print(i)
	querytest<- GDCquery(project="TCGA-KIRP", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }	

LIHC<- read.table("LIHC patients.txt")
rownames(LIHC)=LIHC[,1]
for (i in rownames(LIHC))
{	print(i)
	querytest<- GDCquery(project="TCGA-LIHC", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
	
LUAD<- read.table("LUAD patients.txt")
rownames(LUAD)=LUAD[,1]
for (i in rownames(LUAD))
{	print(i)
	querytest<- GDCquery(project="TCGA-LUAD", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }	
	
LUSC<- read.table("LUSC patients.txt")
rownames(LUSC)=LUSC[,1]
for (i in rownames(LUSC))
{	print(i)
	querytest<- GDCquery(project="TCGA-LUSC", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
	
READ<- read.table("READ patients.txt")
rownames(READ)=READ[,1]
for (i in rownames(READ))
{	print(i)
	querytest<- GDCquery(project="TCGA-READ", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
	
STAD<- read.table("STAD patients.txt")
rownames(STAD)=STAD[,1]
for (i in rownames(STAD))
{	print(i)
	querytest<- GDCquery(project="TCGA-STAD", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }	

THCA<- read.table("THCA patients.txt")
rownames(THCA)=THCA[,1]
for (i in rownames(THCA))
{	print(i)
	querytest<- GDCquery(project="TCGA-THCA", legacy=TRUE, data.category="DNA methylation", barcode=c(i))
	trycatcher(GDCdownload(querytest, method="api")) }
