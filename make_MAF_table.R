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

master <- subset(master, 
          select=c("euro_MAF","CEU_MAF","GBR_MAF","TSI_MAF","IBS_MAF","FIN_MAF",
                   "asian_MAF","CHB_MAF","CHS_MAF","JPT_MAF",
                   "african_MAF","YRI_MAF","LWK_MAF","ASW_MAF",
                   "latino_MAF","MXL_MAF","CLM_MAF","PUR_MAF",
                   "euro_NChr","CEU_NChr","GBR_NChr","TSI_NChr","IBS_NChr","FIN_NChr",
                   "asian_NChr","CHB_NChr","CHS_NChr","JPT_NChr",
                   "african_NChr","YRI_NChr","LWK_NChr","ASW_NChr",
                   "latino_NChr","MXL_NChr","CLM_NChr","PUR_NChr"))

outfile <- gsub(".tped",".popMAFs",tpedfile)

write.table(master, file = outfile,
    row.names = FALSE,
    col.names  = TRUE,
    sep = "\t",
    quote = FALSE)
