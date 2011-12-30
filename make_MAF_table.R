#!/usr/bin/Rscript
args <- commandArgs(TRUE)
tpedfile <- args[1]

files <- list.files(path = ".", pattern ="frq")

master <- read.table(tpedfile,header=FALSE)[c(1,2,4)]
colnames(master) <- c("chr","SNP","pos")

alleles <- read.table(files[1],header=TRUE)[c(2,3,4)]
master <- merge(master, alleles)

for (i in c(1:length(files))) {
  pop <- unlist(strsplit(files[i],"\\."))
  pop <- (pop[[1]])
  tmp <- read.table(files[i],header=TRUE)[c(2,5,6)]
  colnames(tmp) <- c("SNP",paste(pop,"_MAF",sep=""),paste(pop,"_NChr",sep=""))
  master <- merge(master,tmp)
}

outfile <- gsub(".tped",".popMAFs",tpedfile)

write.table(master, file = outfile,
    row.names = FALSE,
    col.names  = TRUE,
    sep = "\t",
    quote = FALSE)
