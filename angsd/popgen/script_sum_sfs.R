#input <- "data/thetas/AI_cat.sfs"
#df.sfs <- read.table(input)

#sfs concatanated across all chromosomes:
## Rscript script_sum_sfs.R AI_cat.sfs
#args <- commandArgs()
args <- commandArgs(trailingOnly = TRUE) #apparently trailing is critial here

input<-args[1]
df.sfs <- read.table(input)
out.name <- gsub(".sfs","_rsum.sfs",basename(input))
sfs <- as.numeric(colSums(df.sfs))
write.table(t(sfs), out.name, sep = "\t", quote=F, row.names = F, col.names = F)
