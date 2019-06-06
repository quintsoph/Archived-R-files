mir651<- microcut(x)
mir651<- na.omit(mir651)

#rank miRNA 
rr651 <- rrfunction(mir651)
invrr651<- invrrfunction(mir651)
invrr651<- rrfunction(invrr651)
Hknot651 <- Hknot(rr651)
Hknotinv651<- Hknot(invrr651)
Hknots<- cbind(Hknot651, Hknots)
Hknotsinv<- cbind(Hknotinv651, Hknotsinv)

#microRNA let 7a 1
x="hsa-let-7a-1"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(w)
invrrw<- rrfunction(invrrfunction(w))
Hknotlet7a1 <- Hknot(rrw)
Hknots<- cbind(Hknotlet7a1, Hknots)
Hknotinvlet7a1<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7a1, Hknotsinv)
Hknots
Hknotsinv

#microRNA let 7a 2
x="hsa-let-7a-2"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7a2 <- Hknot(rrw)
Hknots<- cbind(Hknotlet7a1, Hknots)
Hknotinvlet7a2<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7a2, Hknotsinv)
Hknots
Hknotsinv

#microRNA let 7a 3
x="hsa-let-7a-3"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7a3 <- Hknot(rrw)
Hknots<- cbind(Hknotlet7a3, Hknots)
Hknotinvlet7a3<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7a3, Hknotsinv)
Hknots
Hknotsinv

#mirna let 7b
x="hsa-let-7b"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7b <- Hknot(rrw)
Hknots<- cbind(Hknotlet7b, Hknots)
Hknotinvlet7b<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7b, Hknotsinv)
Hknots
Hknotsinv

#hsa-let-7c
x="hsa-let-7c"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7c <- Hknot(rrw)
Hknots<- cbind(Hknotlet7c, Hknots)
Hknotinvlet7c<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7c, Hknotsinv)
Hknots
Hknotsinv

#hsa-let-7d
x="hsa-let-7d"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7d <- Hknot(rrw)
Hknots<- cbind(Hknotlet7c, Hknots)
Hknotinvlet7d <- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7d, Hknotsinv)
Hknots
Hknotsinv

#hsa-let-7e
x="hsa-let-7e"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7e <- Hknot(rrw)
Hknots<- cbind(Hknotlet7e, Hknots)
Hknotinvlet7e <- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7e, Hknotsinv)
Hknots
Hknotsinv

#hsa-let-7f-1
x="hsa-let-7f-1"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7f1 <- Hknot(rrw)
Hknots<- cbind(Hknotlet7f1, Hknots)
Hknotinvlet7f1 <- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7f1, Hknotsinv)
Hknots
Hknotsinv

#hsa-let-7f-2
x="hsa-let-7f-2"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7f2 <- Hknot(rrw)
Hknots<- cbind(Hknotlet7f2, Hknots)
Hknotinvlet7f2 <- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7f2, Hknotsinv)
Hknots
Hknotsinv

#hsa-let-7g
x="hsa-let-7g"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7g <- Hknot(rrw)
Hknots<- cbind(Hknotlet7g, Hknots)
Hknotinvlet7g <- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7g, Hknotsinv)
Hknots
Hknotsinv

#hsa-let-7i
x="hsa-let-7i"
w <- microcut(x)
w<- na.omit(w)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotlet7i <- Hknot(rrw)
Hknots<- cbind(Hknotlet7i, Hknots)
Hknotinvlet7i <- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvlet7i, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 1 1
x= "hsa-mir-1-1"
w <- microcut(x)
#rank miRNA 
w<-t(w)
w<- na.omit(w)
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 3)
invrrw<- invrrfunction(rankw, 3)
Hknotmir1.1<- Hknot(rrw)
Hknots<- cbind(Hknotmir1.1, Hknots)
Hknotinvmir1.1<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir1.1, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 1 2
x= "hsa-mir-1-2"
w <- microcut(x)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotmir1.2<- Hknot(rrw)
Hknots<- cbind(Hknotmir1.2, Hknots)
Hknotinvmir1.2<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir1.2, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 100
x= "hsa-mir-100"
w <- microcut(x)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotmir100<- Hknot(rrw)
Hknots<- cbind(Hknotmir100, Hknots)
Hknotinvmir100<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir100, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 101-1
x= "hsa-mir-101-1"
w <- microcut(x)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotmir101.1<- Hknot(rrw)
Hknots<- cbind(Hknotmir101.1, Hknots)
Hknotinvmir101.1<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir101.1, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 101-2
x= "hsa-mir-101-2"
w <- microcut(x)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotmir101.2<- Hknot(rrw)
Hknots<- cbind(Hknotmir101.2, Hknots)
Hknotinvmir101.2<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir101.2, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 103-1
x= "hsa-mir-103-1"
w <- microcut(x)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotmir103.1<- Hknot(rrw)
Hknots<- cbind(Hknotmir103.1, Hknots)
Hknotinvmir103.1<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir103.1, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 103-2
x= "hsa-mir-103-2"
w <- microcut(x)
#rank miRNA 
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min")))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotmir103.2<- Hknot(rrw)
Hknots<- cbind(Hknotmir103.2, Hknots)
Hknotinvmir103.2<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir103.2, Hknotsinv)
Hknots
Hknotsinv

#mirna mir 105-1
x= "hsa-mir-105-1"
w <- microcut(x)
#rank miRNA 
w<-t(w)
w<- na.omit(w)
rankw <- data.frame(rank(w, na.last = TRUE, ties.method=c("min"), header= TRUE))
rrw <- rrfunction(rankw, 18)
invrrw<- invrrfunction(rankw, 18)
Hknotmir105.1<- Hknot(rrw)
Hknots<- cbind(Hknotmir105.1, Hknots)
Hknotinvmir105.1<- Hknot(invrrw)
Hknotsinv<- cbind(Hknotinvmir105.1, Hknotsinv)
Hknots
Hknotsinv

 x="hsa-mir-9-3"
 test1<- microcut(x)
 colnames(test1)<- c("UCEC", "THCA", "STAD", "PRAD", "PCPG", "PAAD", "LUSC", "LUAD", "LIHC", "KIRP", "KIRC", "KICH", "HNSC", "ESCA", "CHOL", "CESC", "BRCA", "BLCA")
 test2<-invrrfunction(test1)
 test2.1<- rrfunction(test2)
 Hknot(test2.1)