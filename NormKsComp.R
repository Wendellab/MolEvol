library(seqinr)

#use kaks function of seqinr or use KaKs_Calculator files?
ad1.norm.factor = 0.6 #This is the normalization factor for proper comparison of
#At and Dt sequences.  A to At rate/D to Dt rate for AD1.
ad2.norm.factor = 5/9 #normalization factor for AD2, same calculation as above

files <- list.files(pattern = "Gorai.*fasta$")
key <- read.table("key_anthocyanin.orthoMCL.txt", header = TRUE)

antho.files <- files[unlist(lapply(key[,1],grep,x=files))]

ad1 <- c("ARK2402","CL7","Coker315","CRB252","FM958","LKT511","Maxxa",
         "PM145","TAMCOT","TM1","TX1009","TX1037","TX1055","TX1107","TX1110",
         "TX1120","TX1182","TX1226","TX1228","TX1236","TX1748","TX1982",
         "TX1988","TX1996","TX2002","TX2089","TX2090","TX2091","TX2092",
         "TX2094","TX2095","TX44","TX480","TX665","TX672","TX786")
ad2 <- c("GB0262","GB0287","GB0303","GB0369","GPS52","K101","Phy76","PS6")
ad1.a <- lapply(ad1,paste,"A", sep = ".")
ad1.d <- lapply(ad1,paste,"D", sep = ".")
ad2.a <- lapply(ad2,paste,"A", sep = ".")
ad2.d <- lapply(ad2,paste,"D", sep = ".")

for(i in 1:length(antho.files)){
  fas.aln <- read.alignment(file = antho.files[i],format = "fasta")
  At.index.ad1 <- unlist(lapply(ad1.a,grep,x = fas.aln$nam))
  Dt.index.ad1 <- unlist(lapply(ad1.d,grep,x = fas.aln$nam))
  At.index.ad2 <- unlist(lapply(ad2.a,grep,x = fas.aln$nam))
  Dt.index.ad2 <- unlist(lapply(ad2.d,grep,x = fas.aln$nam))  
  D5.index <- grep("^D5",fas.aln$nam)
  A2.index <- grep("^A2",fas.aln$nam)
  K.table <- kaks(fas.aln)
  A.Ks.ad1 <- as.matrix(K.table$ks)[A2.index,c(At.index.ad1)]
  A.Ks.ad1.avg <- mean(A.Ks.ad1)
  D.Ks.ad1 <- as.matrix(K.table$ks)[D5.index,c(Dt.index.ad1)]
  D.Ks.ad1.avg <- mean(D.Ks.ad1)
  A.Ks.ad2 <- as.matrix(K.table$ks)[A2.index,c(At.index.ad2)]
  A.Ks.ad2.avg <- mean(A.Ks.ad2)
  D.Ks.ad2 <- as.matrix(K.table$ks)[D5.index,c(Dt.index.ad2)]
  D.Ks.ad2.avg <- mean(D.Ks.ad2)
  D.Ks.ad1.adj <- D.Ks.ad1*ad1.norm.factor
  D.Ks.ad1.adj.avg <- mean(D.Ks.ad1.adj)
  D.Ks.ad2.adj <- D.Ks.ad2*ad2.norm.factor
  D.Ks.ad2.adj.avg <- mean(D.Ks.ad2.adj)
  ad1.result <- try(wilcox.test(x = A.Ks.ad1, y = D.Ks.ad1.adj))
  if(inherits(ad1.result,"try-error")){
    ad1.result$statistic <- NA
    ad1.result$p.value <- NA
    ad1.result$method <- NA
  }
  ad2.result <- try(wilcox.test(x = A.Ks.ad2, y = D.Ks.ad2.adj))
  if(inherits(ad2.result,"try-error")){
    ad2.result$statistic <- NA
    ad2.result$p.value <- NA
    ad2.result$method <- NA
  }
  write(c(antho.files[i],A.Ks.ad1.avg,D.Ks.ad1.avg,D.Ks.ad1.adj.avg,ad1.result$statistic,ad1.result$p.value,ad1.result$method),file = "NormKsCompResultsAD1.txt",ncolumns = 7,sep = "\t",append =TRUE)
  write(c(antho.files[i],A.Ks.ad2.avg,D.Ks.ad2.avg,D.Ks.ad2.adj.avg,ad2.result$statistic,ad2.result$p.value,ad2.result$method),file = "NormKsCompResultsAD2.txt",ncolumns = 7,sep = "\t",append =TRUE)
}

ad1.norm.factor <- 2/3
ad2.norm.factor <- 2/3

for(i in 1:length(antho.files)){
  fas.aln <- read.alignment(file = antho.files[i],format = "fasta")
  At.index.ad1 <- unlist(lapply(ad1.a,grep,x = fas.aln$nam))
  Dt.index.ad1 <- unlist(lapply(ad1.d,grep,x = fas.aln$nam))
  At.index.ad2 <- unlist(lapply(ad2.a,grep,x = fas.aln$nam))
  Dt.index.ad2 <- unlist(lapply(ad2.d,grep,x = fas.aln$nam))  
  D5.index <- grep("^D5",fas.aln$nam)
  A2.index <- grep("^A2",fas.aln$nam)
  K.table <- kaks(fas.aln)
  A.Ka.ad1 <- as.matrix(K.table$ka)[A2.index,c(At.index.ad1)]
  A.Ka.ad1.avg <- mean(A.Ka.ad1)
  D.Ka.ad1 <- as.matrix(K.table$ka)[D5.index,c(Dt.index.ad1)]
  D.Ka.ad1.avg <- mean(D.Ka.ad1)
  A.Ka.ad2 <- as.matrix(K.table$ka)[A2.index,c(At.index.ad2)]
  A.Ka.ad2.avg <- mean(A.Ka.ad2)
  D.Ka.ad2 <- as.matrix(K.table$ka)[D5.index,c(Dt.index.ad2)]
  D.Ka.ad2.avg <- mean(D.Ka.ad2)
  D.Ka.ad1.adj <- D.Ka.ad1*ad1.norm.factor
  D.Ka.ad1.adj.avg <- mean(D.Ka.ad1.adj)
  D.Ka.ad2.adj <- D.Ka.ad2*ad2.norm.factor
  D.Ka.ad2.adj.avg <- mean(D.Ka.ad2.adj)
  ad1.result <- try(wilcox.test(x = A.Ka.ad1, y = D.Ka.ad1.adj))
  if(inherits(ad1.result,"try-error")){
    ad1.result$statistic <- NA
    ad1.result$p.value <- NA
    ad1.result$method <- NA
  }
  ad2.result <- try(wilcox.test(x = A.Ka.ad2, y = D.Ka.ad2.adj))
  if(inherits(ad2.result,"try-error")){
    ad2.result$statistic <- NA
    ad2.result$p.value <- NA
    ad2.result$method <- NA
  }
  write(c(antho.files[i],A.Ka.ad1.avg,D.Ka.ad1.avg,D.Ka.ad1.adj.avg,ad1.result$statistic,ad1.result$p.value,ad1.result$method),file = "NormKaCompResultsAD1.txt",ncolumns = 7,sep = "\t",append =TRUE)
  write(c(antho.files[i],A.Ka.ad2.avg,D.Ka.ad2.avg,D.Ka.ad2.adj.avg,ad2.result$statistic,ad2.result$p.value,ad2.result$method),file = "NormKaCompResultsAD2.txt",ncolumns = 7,sep = "\t",append =TRUE)
}
