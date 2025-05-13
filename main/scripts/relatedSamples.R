rm(list=ls())
gc()

if(!require(plinkQC)){
    install.packages("plinkQC")
}
library("plinkQC")

args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
name <- args[2]
fail <- args[3]
plink <- args[4]
genome_build <- args[5]
ibd <- as.numeric(args[6])
sub=strsplit(name,"\\.")[[1]][1]
output <- paste0(indir,"/",sub,".",fail)

res = check_relatedness(indir = indir,name,filter_high_ldregion=F,genomebuild=genome_build, imissTh = 1, path2plink = plink, run.check_relatedness = T, highIBDTh = ibd)

write.table(res$failIDs, output, quote = FALSE, row.names = FALSE, col.names = FALSE )
