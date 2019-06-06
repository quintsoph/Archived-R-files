#for creating a combined graph of the top miRNA
#read in library ggplot2
#read in the files for creating the graphs
#Find patients per cancer patients
setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
Names<- read.table("SophiaPatCanTypes", header=FALSE)
tNames<- t(Names)
FPPP<- subset(tNames, tNames[,2] == "FPPP", drop=FALSE)


setwd("F:/Bioinformatics Lab/Cancer Data/miRNASeq/All Cancer Transposons/New Files (need org)")
#format cancer
newL1HScancer<- data.frame(read.table(file="L1HS.Cancer.BaseMeansB", header=TRUE))
newL1HScancer<- data.matrix(newL1HScancer)
tnewL1HScancer<- t(newL1HScancer)
newmirnacancer<- read.csv(file="Mirna.Cancer.csv")
newmirnacancer<- data.matrix(newmirnacancer)
cancermirnanames<- read.csv(file="cancer mirna names.csv")
rownames(newmirnacancer)=cancermirnanames[,1]
tnewmirnacancer<- t(newmirnacancer)
tnewmirnacancer<- tnewmirnacancer[-1,]
tnewL1HScancer<- tnewL1HScancer[-1,]
cancer<- cbind(tnewL1HScancer, tnewmirnacancer)
cancermirna<- cancer[,-1]

#format normal
newL1HSnormal<- data.frame(read.table(file="L1HS.Norma.BaseMeansA", header=TRUE))
newL1HSnormal<- data.matrix(newL1HSnormal)
tnewL1HSnormal<- t(newL1HSnormal)
newmirnanormal<- read.csv(file="Mirna.Normal.csv")
newmirnanormal<- data.matrix(newmirnanormal)
normalmirnanames<- read.csv(file="normal mirna names.csv")
rownames(newmirnanormal)=normalmirnanames[,1]
tnewmirnanormal<- t(newmirnanormal)
tnewmirnanormal<- tnewmirnanormal[-1,]
tnewL1HSnormal<- tnewL1HSnormal[-1,]
normal<- cbind(tnewL1HSnormal, tnewmirnanormal)
normalmirna<- normal[,-1]
dividenL1HS<- newL1HScancer/newL1HSnormal
foldedL1HS<- log2(dividenL1HS)
dividenmirna<- newmirnacancer/newmirnanormal
foldedmirna<- t(log2(dividenmirna))
foldedmirna<- foldedmirna[-1,]
#ratio<- cancer/normal

#PRAD
PRAD<- subset(tNames, tNames[,2] == "PRAD", drop=FALSE)
matched<- match(PRAD[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(PRAD[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.204")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.378c")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.487b")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.204")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.378c")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.487b")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in PRAD", y = "L1HS Expression(RPM)", title="PRAD")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#BLCA
BLCA<- subset(tNames, tNames[,2] == "BLCA", drop=FALSE)
matched<- match(BLCA[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(BLCA[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.195")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.143")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.497")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.195")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.143")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.497")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in BLCA", y = "L1HS Expression(RPM)", title="BLCA")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 


#BRCA
BRCA<- subset(tNames, tNames[,2] == "BRCA", drop=FALSE)
matched<- match(BRCA[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(BRCA[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.22")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.486")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.378")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.22")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.486")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.378")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in BRCA", y = "L1HS Expression(RPM)", title="BRCA")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#HNSC
HNSC<- subset(tNames, tNames[,2] == "HNSC", drop=FALSE)
matched<- match(HNSC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(HNSC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.29a")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.424")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.450a.2")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.29a")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.424")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.450a.2")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in HNSC", y = "L1HS Expression(RPM)", title="HNSC")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#KICH
KICH<- subset(tNames, tNames[,2] == "KICH", drop=FALSE)
matched<- match(KICH[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(KICH[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.374a")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.503")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.505")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(NULL) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.374a")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.503")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.505")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in KICH", y = "L1HS Expression(RPM)", title="KICH")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#KIRC
KIRC<- subset(tNames, tNames[,2] == "KIRC", drop=FALSE)
matched<- match(KIRC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(KIRC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.379")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.127")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.21")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.379")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.127")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.21")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in KIRC", y = "L1HS Expression(RPM)", title="KIRC")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#KIRP
KIRP<- subset(tNames, tNames[,2] == "KIRP", drop=FALSE)
matched<- match(KIRP[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(KIRP[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.let.7d")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.181b.2")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.148a")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.let.7d")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.181b.2")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.148a")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in KIRP", y = "L1HS Expression(RPM)", title="KIRP")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#LIHC
LIHC<- subset(tNames, tNames[,2] == "LIHC", drop=FALSE)
matched<- match(LIHC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(LIHC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.29a")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.29b.2")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.29b.1")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.29a")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.29b.2")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.29b.1")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in LIHC", y = "L1HS Expression(RPM)", title="LIHC")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#LUAD
LUAD<- subset(tNames, tNames[,2] == "LUAD", drop=FALSE)
matched<- match(LUAD[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(LUAD[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.542")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.10b")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.101.1")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.542")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.10b")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.101.1")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in LUAD", y = "L1HS Expression(RPM)", title="LUAD")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 


#LUSC
LUSC<- subset(tNames, tNames[,2] == "LUSC", drop=FALSE)
matched<- match(LUSC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(LUSC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.let.7e")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.181b.1")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.153.2")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.let.7e")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.181b.1")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.153.2")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in LUSC", y = "L1HS Expression(RPM)", title="LUSC")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 


#THCA
THCA<- subset(tNames, tNames[,2] == "THCA", drop=FALSE)
matched<- match(THCA[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(THCA[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.514.3")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.508")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.514.1")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.514.3")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.508")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.514.1")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in THCA", y = "L1HS Expression(RPM)", title="THCA")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 

#UCEC
UCEC<- subset(tNames, tNames[,2] == "UCEC", drop=FALSE)
matched<- match(UCEC[,1], colnames(foldedL1HS))
L1HS <-foldedL1HS[,matched]
L1HS<- data.frame(L1HS)
matched<- match(UCEC[,1], rownames(foldedmirna))
mirna<- foldedmirna[matched,]
mirna<- data.frame(mirna)
mirna1<- data.frame(mirna[,c("hsa.mir.133a.1")])
first<- cbind(L1HS, mirna1)
mirna2<- data.frame(mirna[,c("hsa.mir.1.2")])
second<- cbind(L1HS, mirna2)
mirna3<- data.frame(mirna[,c("hsa.mir.145")])
third<- cbind(L1HS, mirna3)
together<- cbind(first, second[,2], third[,2])
colnames(together)=c("L1HS", "first", "second", "third")
legend_title<- "MIRNA Names"
ggplot(together) + 
  geom_jitter(aes(first,L1HS, color="hsa.mir.133a.1")) + geom_smooth(aes(first,L1HS), method=lm, se=FALSE, color="darkblue") +
  geom_jitter(aes(second,L1HS, color="hsa.mir.1.2")) + geom_smooth(aes(second,L1HS), method=lm, se=FALSE, color="black") +
  geom_jitter(aes(third,L1HS, color="hsa.mir.145")) + geom_smooth(aes(third,L1HS), method=lm, se=FALSE, color="darkred") +
  labs(x = "Top Three MiRNA in UCEC", y = "L1HS Expression(RPM)", title="UCEC")+
  theme(legend.position="right")+ 
	scale_color_manual(legend_title, values = c("darkblue","black","darkred")) 




  
  