rm(list=ls())
gc()

library("plinkQC")

args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
name <- args[2]
fail <- args[3]
plink <- args[4]
het <- as.numeric(args[5])
sub=strsplit(name,"\\.")[[1]][1]
output <- paste0(indir,"/",sub,".",fail)

res = check_het_and_miss(indir = indir,name, imissTh = 1, hetTh=het, path2plink = plink, run.check_het_and_miss = T)

dt = rbind(data.frame(FID=res$fail_imiss$FID,IID=res$fail_imiss$IID),data.frame(FID=res$fail_het$FID,IID=res$fail_het$IID))

write.table(dt, output, quote=F, row.names=F,col.name=F, sep="\t") # print FID and IID for valid samples
