#Q1: How does polyploidy affect the evolutionary rate of the genes of the ABP?
#H0: The genes continue to evolve at the same rate as prior to polyploidy.  
#Ha: The genes evolve at rate different from those prior to polyploidy due to 
#changes in selective pressures.
#How: relative rate test, absolute number of substitutions between diploids and
#their subgenomes

library(seqinr)
library(plyr)

#read in alignment files
files <- list.files(pattern = "Gorai.*fasta$")
key <- read.table("key_anthocyanin.orthoMCL.txt", header = TRUE)

antho.files <- files[unlist(lapply(key[,1],grep,x=files))]

ad1 <- c("ARK2402","CL7","Coker315","CRB252","FM958","LKT511","Maxxa",
         "PM145","TAMCOT","TM1","TX1009","TX1037","TX1055","TX1107","TX1110",
         "TX1120","TX1182","TX1226","TX1228","TX1236","TX1748","TX1982",
         "TX1988","TX1996","TX2002","TX2089","TX2090","TX2091","TX2092",
         "TX2094","TX2095","TX44","TX480","TX665","TX672","TX786")
ad2 <- c("GB0262","GB0287","GB0303","GB0369","GPS52","K101","Phy76","PS6")

monomorphic <- function(pos,ambi="N") {
  #This function returns TRUE if a site is monomorphic or FALSE if a site is
  #polymorphic
  #position = a vector of sites that correspond with all the residues in a 
  #position
  #p = proportion of sites necessary to call as polymorphic
  #ambi = ambiguous site designation
  #Checks site to see if all the characters are the same
  require(plyr)
  pos.res.table <- count(pos)
  ambi.index <- grep(ambi,pos.res.table[,1],ignore.case=TRUE)
  if(length(ambi.index) >= 1){
    pos.res.table <- pos.res.table[-ambi.index,]
  }    
  if(nrow(pos.res.table) == 1){
    result <- TRUE
  } else{
      result <- FALSE
  }
  result
}

#Tajima's 1D Relative Rate Test
for(i in 1:length(antho.files)) {
  fas.aln <- read.alignment(file = antho.files[i],format = "fasta")
  At.index <- grep("A$",fas.aln$nam)
  Dt.index <- grep("D$",fas.aln$nam)
  D5.index <- grep("^D5",fas.aln$nam)
  A2.index <- grep("^A2",fas.aln$nam)
  #monomorphic sites between A2 and D5
  A2.to.D5 <- grep(pattern = TRUE,x = apply(rbind(s2c(unlist(fas.aln$seq[D5.index])),s2c(unlist(fas.aln$seq[A2.index]))),2,monomorphic))
  K.A2.to.D5 <- length(s2c(unlist(fas.aln$seq[D5.index]))) - length(A2.to.D5)
  for(j in 1:length(At.index)){
    #find the indices of their sequences
    og.index <- D5.index
    ig.index <- A2.index
    seq.index <- At.index[j]
    #find the monomorphic sites between the og and seq of interest
    seq.to.og <- grep(pattern = TRUE,x = apply(rbind(s2c(unlist(fas.aln$seq[og.index])),s2c(unlist(fas.aln$seq[seq.index]))),2,monomorphic))
    K.seq.to.og <- length(s2c(unlist(fas.aln$seq[og.index]))) - length(seq.to.og)
    seq.to.ig <- grep(pattern = TRUE,x = apply(rbind(s2c(unlist(fas.aln$seq[ig.index])),s2c(unlist(fas.aln$seq[seq.index]))),2,monomorphic))
    K.seq.to.ig <- length(s2c(unlist(fas.aln$seq[ig.index]))) - length(seq.to.ig)
    K.seq <- (K.seq.to.og + K.seq.to.ig - K.A2.to.D5)/2
    K.og <- (K.A2.to.D5 + K.seq.to.og - K.seq.to.ig)/2
    K.ig <- (K.seq.to.ig + K.A2.to.D5 - K.seq.to.og)/2
    d <- K.seq - K.ig
    #find sites that are monomorphic in one group, but not in the other
    dip.unique <- setdiff(seq.to.og,A2.to.D5)
    pp.unique <- setdiff(A2.to.D5,seq.to.og)
    m1 <- length(dip.unique)
    m2 <- length(pp.unique)
    chisq.val <- (m1 - m2)^2/(m1 + m2)
    p.val <- 1-pchisq(chisq.val,df = 1)
    write(x = c(antho.files[i],"D5","A2",fas.aln$nam[seq.index],m1,m2,chisq.val,p.val),file = "Taj1DRelRateTest.txt",ncol = 8, append = TRUE,sep = "\t")
  }
  for(j in 1:length(Dt.index)){
    #find the indices of their sequences
    og.index <- A2.index
    seq.index <- Dt.index[j]
    #find the pm sites between the og and seq of interest
    seq.to.og <- grep(pattern = TRUE,x = apply(rbind(s2c(unlist(fas.aln$seq[og.index])),s2c(unlist(fas.aln$seq[seq.index]))),2,monomorphic))
    #find sites that are monomorphic in one group, but not in the other
    dip.unique <- setdiff(seq.to.og,A2.to.D5)
    pp.unique <- setdiff(A2.to.D5,seq.to.og)
    m1 <- length(dip.unique)
    m2 <- length(pp.unique)
    chisq.val <- (m1 - m2)^2/(m1 + m2)
    p.val <- 1-pchisq(chisq.val,df = 1)
    write(x = c(antho.files[i],"A2","D5",fas.aln$nam[seq.index],m1,m2,chisq.val,p.val),file = "Taj1DRelRateTest.txt",ncol = 8, append = TRUE,sep = "\t")
  }
}

#Tajima's 2D Relative Rate Test
#?

