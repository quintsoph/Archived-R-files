#Checking high lows
#hsa.mir.708
#hsa.mir.621
#read in the created files
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)/Final analysis results/log2 foldchange ratio by each")
BLCA<- data.frame(read.table(file="BLCA pvalues and rvalues of log2 foldchange ratio.txt"))
BRCA<- data.frame(read.table(file="BRCA pvalues and rvalues of log2 foldchange ratio.txt"))
CESC<- data.frame(read.table(file="CESC pvalues and rvalues of log2 foldchange ratio.txt"))
CHOL<- data.frame(read.table(file="CHOL pvalues and rvalues of log2 foldchange ratio.txt"))
ESCA<- data.frame(read.table(file="ESCA pvalues and rvalues of log2 foldchange ratio.txt"))
HNSC<- data.frame(read.table(file="HNSC pvalues and rvalues of log2 foldchange ratio.txt"))
KICH<- data.frame(read.table(file="KICH pvalues and rvalues of log2 foldchange ratio.txt"))
KIRC<- data.frame(read.table(file="KIRC pvalues and rvalues of log2 foldchange ratio.txt"))
KIRP<- data.frame(read.table(file="KIRP pvalues and rvalues of log2 foldchange ratio.txt"))
LIHC<- data.frame(read.table(file="LIHC pvalues and rvalues of log2 foldchange ratio.txt"))
LUAD<- data.frame(read.table(file="LUAD pvalues and rvalues of log2 foldchange ratio.txt"))
LUSC<- data.frame(read.table(file="LUSC pvalues and rvalues of log2 foldchange ratio.txt"))
PAAD<- data.frame(read.table(file="PAAD pvalues and rvalues of log2 foldchange ratio.txt"))
PCPG<- data.frame(read.table(file="PCPG pvalues and rvalues of log2 foldchange ratio.txt"))
PRAD<- data.frame(read.table(file="PRAD pvalues and rvalues of log2 foldchange ratio.txt"))
STAD<- data.frame(read.table(file="STAD pvalues and rvalues of log2 foldchange ratio.txt"))
THCA<- data.frame(read.table(file="THCA pvalues and rvalues of log2 foldchange ratio.txt"))
THYM<- data.frame(read.table(file="THYM pvalues and rvalues of log2 foldchange ratio.txt"))
UCEC<- data.frame(read.table(file="UCEC pvalues and rvalues of log2 foldchange ratio.txt"))
READ<- data.frame(read.table(file="READ pvalues and rvalues of log2 foldchange ratio.txt"))
FPPP<- data.frame(read.table(file="FPPP pvalues and rvalues of log2 foldchange ratio.txt"))
COAD<- data.frame(read.table(file="COAD pvalues and rvalues of log2 foldchange ratio.txt"))

#merge all files
BLCA2<- data.frame(t(BLCA))
BLCA2<- BLCA2[-2,]
BRCA2<- data.frame(t(BRCA))
BRCA2<- BRCA2[-2,]
CESC2<- data.frame(t(CESC))
CESC2<- CESC2[-2,]
CHOL2<- data.frame(t(CHOL))
CHOL2<- CHOL2[-2,]
ESCA2<- data.frame(t(ESCA))
ESCA2<- ESCA2[-2,]
HNSC2<- data.frame(t(HNSC))
HNSC2<- HNSC2[-2,]
KICH2<- data.frame(t(KICH))
KICH2<- KICH2[-2,]
KIRC2<- data.frame(t(KIRC))
KIRC2<- KIRC2[-2,]
KIRP2<- data.frame(t(KIRP))
KIRP2<- KIRP2[-2,]
LIHC2<- data.frame(t(LIHC))
LIHC2<- LIHC2[-2,]
LUAD2<- data.frame(t(LUAD))
LUAD2<- LUAD2[-2,]
LUSC2<- data.frame(t(LUSC))
LUSC2<- LUSC2[-2,]
PAAD2<- data.frame(t(PAAD))
PAAD2<- PAAD2[-2,]
PCPG2<- data.frame(t(PCPG))
PCPG2<- PCPG2[-2,]
PRAD2<- data.frame(t(PRAD))
PRAD2<- PRAD2[-2,]
STAD2<- data.frame(t(STAD))
STAD2<- STAD2[-2,]
THCA2<- data.frame(t(THCA))
THCA2<- THCA2[-2,]
THYM2<- data.frame(t(THYM))
THYM2<- THYM2[-2,]
UCEC2<- data.frame(t(UCEC))
UCEC2<- UCEC2[-2,]
COAD2<- data.frame(t(COAD))
COAD2<- COAD2[-2,]

#708
BLCA708<- BLCA2[,c("hsa.mir.708")]
BRCA708<- BRCA2[,c("hsa.mir.708")]
CESC708<- CESC2[,c("hsa.mir.708")]
CHOL708<- CHOL2[,c("hsa.mir.708")]
ESCA708<- ESCA2[,c("hsa.mir.708")]
HNSC708<- HNSC2[,c("hsa.mir.708")]
KICH708<- KICH2[,c("hsa.mir.708")]
KIRC708<- KIRC2[,c("hsa.mir.708")]
KIRP708<- KIRP2[,c("hsa.mir.708")]
LIHC708<- LIHC2[,c("hsa.mir.708")]
LUAD708<- LUAD2[,c("hsa.mir.708")]
LUSC708<- LUSC2[,c("hsa.mir.708")]
PAAD708<- PAAD2[,c("hsa.mir.708")]
PCPG708<- PCPG2[,c("hsa.mir.708")]
PRAD708<- PRAD2[,c("hsa.mir.708")]
STAD708<- STAD2[,c("hsa.mir.708")]
THCA708<- THCA2[,c("hsa.mir.708")]
THYM708<- THYM2[,c("hsa.mir.708")]
UCEC708<- UCEC2[,c("hsa.mir.708")]
COAD708<- COAD2[,c("hsa.mir.708")]
all708<- rbind(BLCA708, BRCA708, CESC708, CHOL708, ESCA708, HNSC708, KICH708,
KIRC708, KIRP708, LIHC708, LUAD708, LUSC708, PAAD708, PCPG708, PRAD708, STAD708,
THCA708, THYM708, UCEC708, COAD708)
colnames(all708)=c("P-Values")
setwd("F:/Bioinformatics Lab/Cancer Data/1G Files")
write.table(all708, file="hsa.mir.708 check.txt", sep="\t", quote=FALSE)

#621
BLCA219.2<- BLCA2[,c("hsa.mir.219.2")]
BRCA219.2<- BRCA2[,c("hsa.mir.219.2")]
CESC219.2<- CESC2[,c("hsa.mir.219.2")]
CHOL219.2<- CHOL2[,c("hsa.mir.219.2")]
ESCA219.2<- ESCA2[,c("hsa.mir.219.2")]
HNSC219.2<- HNSC2[,c("hsa.mir.219.2")]
KICH219.2<- KICH2[,c("hsa.mir.219.2")]
KIRC219.2<- KIRC2[,c("hsa.mir.219.2")]
KIRP219.2<- KIRP2[,c("hsa.mir.219.2")]
LIHC219.2<- LIHC2[,c("hsa.mir.219.2")]
LUAD219.2<- LUAD2[,c("hsa.mir.219.2")]
LUSC219.2<- LUSC2[,c("hsa.mir.219.2")]
PAAD219.2<- PAAD2[,c("hsa.mir.219.2")]
PCPG219.2<- PCPG2[,c("hsa.mir.219.2")]
PRAD219.2<- PRAD2[,c("hsa.mir.219.2")]
STAD219.2<- STAD2[,c("hsa.mir.219.2")]
THCA219.2<- THCA2[,c("hsa.mir.219.2")]
THYM219.2<- THYM2[,c("hsa.mir.219.2")]
UCEC219.2<- UCEC2[,c("hsa.mir.219.2")]
COAD219.2<- COAD2[,c("hsa.mir.219.2")]
all219.2<- rbind(BLCA219.2, BRCA219.2, CESC219.2, CHOL219.2, ESCA219.2, HNSC219.2,
KIRC219.2, KIRP219.2, LIHC219.2, LUAD219.2, LUSC219.2, PAAD219.2, PRAD219.2, STAD219.2,
THCA219.2, THYM219.2, UCEC219.2, COAD219.2)
colnames(all219.2)=c("P-Values")
write.table(all219.2, file="hsa.mir.219.2 check.txt", sep="\t", quote=FALSE)


#get the rankings in each cancer
#read in REC score stuff
rankBLCA<- BLCANA[c("hsa.mir.708"),]
rankBRCA<- BRCANA[c("hsa.mir.708"),]
rankCESC<- CESCNA[c("hsa.mir.708"),]
rankCHOL<- CHOLNA[c("hsa.mir.708"),]
rankESCA<- ESCANA[c("hsa.mir.708"),]
rankHNSC<- HNSCNA[c("hsa.mir.708"),]
rankKICH<- KICHNA[c("hsa.mir.708"),]
rankKIRC<- KIRCNA[c("hsa.mir.708"),]
rankKIRP<- KIRPNA[c("hsa.mir.708"),]
rankLIHC<- LIHCNA[c("hsa.mir.708"),]
rankLUAD<- LUADNA[c("hsa.mir.708"),]
rankLUSC<- LUSCNA[c("hsa.mir.708"),]
rankPAAD<- PAADNA[c("hsa.mir.708"),]
rankPCPG<- PCPGNA[c("hsa.mir.708"),]
rankPRAD<- PRADNA[c("hsa.mir.708"),]
rankSTAD<- STADNA[c("hsa.mir.708"),]
rankTHCA<- THCANA[c("hsa.mir.708"),]
rankUCEC<- UCECNA[c("hsa.mir.708"),]
rankCOAD<- COADNA[c("hsa.mir.708"),]
ranks708<- rbind(rankBLCA, rankBRCA, rankCESC, rankCHOL, rankESCA,
rankHNSC, rankKICH, rankKIRC, rankKIRP, rankLIHC, rankLUAD, rankLIHC,
rankLUAD, rankLUSC, rankPAAD, rankPCPG, rankPRAD, rankSTAD, rankTHCA, rankUCEC, rankCOAD)

rankBLCA<- BLCANA[c("hsa.mir.219.2"),]
rankBRCA<- BRCANA[c("hsa.mir.219.2"),]
rankCESC<- CESCNA[c("hsa.mir.219.2"),]
rankCHOL<- CHOLNA[c("hsa.mir.219.2"),]
rankESCA<- ESCANA[c("hsa.mir.219.2"),]
rankHNSC<- HNSCNA[c("hsa.mir.219.2"),]
rankKICH<- KICHNA[c("hsa.mir.219.2"),]
rankKIRC<- KIRCNA[c("hsa.mir.219.2"),]
rankKIRP<- KIRPNA[c("hsa.mir.219.2"),]
rankLIHC<- LIHCNA[c("hsa.mir.621"),]
rankLUAD<- LUADNA[c("hsa.mir.219.2"),]
rankLUSC<- LUSCNA[c("hsa.mir.219.2"),]
rankPAAD<- PAADNA[c("hsa.mir.219.2"),]
rankPCPG<- PCPGNA[c("hsa.mir.219.2"),]
rankPRAD<- PRADNA[c("hsa.mir.219.2"),]
rankSTAD<- STADNA[c("hsa.mir.219.2"),]
rankTHCA<- THCANA[c("hsa.mir.219.2"),]
rankUCEC<- UCECNA[c("hsa.mir.219.2"),]
rankCOAD<- COADNA[c("hsa.mir.219.2"),]
ranks219.2<- rbind(rankBLCA, rankBRCA, rankCESC, rankCHOL, rankESCA,
rankHNSC, rankKICH, rankKIRC, rankKIRP, rankLIHC, rankLUAD, rankLIHC,
rankLUAD, rankLUSC, rankPAAD, rankPCPG, rankPRAD, rankSTAD, rankTHCA, rankUCEC, rankCOAD)








