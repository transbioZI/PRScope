rm(list=ls())
gc()

library("data.table")
args = commandArgs(trailingOnly=TRUE)
maf_file = args[1]
ss_path = args[2]
suffix = args[3]

maf_file = "/data/references/1000Genome/reference/EUR/eur_hg38.phase3.frq"
ss_path = "/data/sumstats/transformed"
suffix = "transformed.h.tsv"

print("read maf file")
maf = as.data.frame(fread(maf_file))
ss = list.files(ss_path,suffix)
lenght_of_na = c()
print("start")
foreach(xy = c(1:length(ss)) ) {
  s = ss[xy]
  print(s)
  ss_file_path = paste0(ss_path,"/",s)
  current_s = as.data.frame(fread(ss_file_path, nThread = 10))
  matches = match(current_s$VARID,maf$SNP)
  matches_na = which(is.na(matches))
  matches_not_na = which(!is.na(matches))
  current_s$BASEMAF = NA
  current_s[matches_not_na,]$BASEMAF = maf$MAF[matches[!is.na(matches)]]
  current_s[matches_na,]$BASEMAF = 0.000000001
  
  write.table(current_s, paste0(ss_file_path,".new"), quote = F, row.names =F, sep = " ")
  system(paste0("rm ",ss_file_path))
  system(paste0("mv ",paste0(ss_file_path,".new"), " ",ss_file_path))
  lenght_of_na = c(lenght_of_na,length(matches_na))
  print("done")
}
print(lenght_of_na)

jpeg(paste0(ss_path,"/lenght_of_na.jpg"))
hist(lenght_of_na)
dev.off()
