#Q2: How do the subgenomes differ in their evolutionary rates in the ABP?
#H0: The homoeologous genes of both subgenomes evolve at similar rates due to 
#similar selective forces.
#HA: The homoeologous genes from each subgenome evolve at different rates due to
#selective forces unique to each subgenome.
#How: Using prior rates of substitutions to normalize for differential distances 
#between diploids (in Flagel et al 2012), determine if the distance is greater 
#between D-DT or A-AT.  

library(seqinr)

#use kaks function of seqinr or use KaKs_Calculator files?
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
  At.index <- unlist(lapply(ad1.a,grep,x = fas.aln$nam))
  Dt.index <- unlist(lapply(ad1.d,grep,x = fas.aln$nam))
  D5.index <- grep("^D5",fas.aln$nam)
  A2.index <- grep("^A2",fas.aln$nam)
  K.table <- kaks(fas.aln)
  A.Ks <- as.matrix(K.table$ks)[A2.index,c(At.index)]
  D.Ks <- as.matrix(K.table$ks)[D5.index,c(Dt.index)]
  D.Ks.adj <- D.Ks*0.6
  result <- wilcox.test(x = A.Ks, y = D.Ks.adj,)
  print(c(antho.files[i],result$statistic,result$p.value,result$method))
}

