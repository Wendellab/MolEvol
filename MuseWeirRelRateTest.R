#This just for filtering HyPhy's Muse and Weir relative rate test results
library(seqinr)
library(plyr)

files <- list.files(pattern = "Gorai.*csv$")
key <- read.table("key_anthocyanin.orthoMCL.txt", header = TRUE)

antho.files <- files[unlist(lapply(key[,1],grep,x=files))]

name <- vector()
gen <- vector()
num <- vector()
for(i in 1:length(antho.files)){
  LRT.table <- read.csv(file = antho.files[i],header = TRUE)
  if(grepl("A2toAgenome",x = antho.files[i])){
    index <- intersect(grep("^A2",LRT.table[,2]),grep(".A$",LRT.table[,3]))
    genome <- "A"
  }
  else{
    index <- union(intersect(grep("^D5",LRT.table[,2]),grep(".D$",LRT.table[,3])),
                   intersect(grep("^D5",LRT.table[,3]),grep(".D$",LRT.table[,2])))
    genome <- "D"
  }
  LRT.table.trunc <- LRT.table[index,]
  count <- 0
  for(j in 1:length(LRT.table.trunc$P.Value)){
    if(LRT.table.trunc[j,15] < 0.05){
      count <- count + 1
    }
  }
  name <- append(name,antho.files[i])
  gen <- append(gen,genome)
  num <- append(num,count)  
}
d <- data.frame(name,gen,num)
write.table(x = d,file = "MuseWeir.LRT.table.txt",sep = "\t")