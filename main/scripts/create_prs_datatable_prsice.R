rm(list=ls())
gc()
library(stringr)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

prs_path = args[1]
threshold = as.numeric(args[2])
up_threshold = as.numeric(args[3])
study_list = args[4]
name = args[5]

study_metadata = fread(study_list)
files = study_metadata$study_id

studies_to_remove = c()
studies_empty = c()
for(f in files) {
  ff= sprintf("%s/%s.prsice", prs_path, f)
  if(file.exists(ff) & file.size(ff) != 0L) {
    nsnp <- read.table(ff, header = T, colClasses = c("character", "character","character", "numeric"))
    nsnp$use = nsnp[,"Num_SNP"] > up_threshold
    nsnp = nsnp[nsnp$use == T,]
    if(dim(nsnp)[1] == 0) {
      studies_to_remove = c(studies_to_remove,f)
    }
  } else {
    if(file.exists(ff) & file.size(ff) == 0L){
      studies_empty = c(studies_empty,f)
    } else {
      studies_to_remove = c(studies_to_remove,f)
    }
  }
}

files_exist = c()
for(f in files) {
  if(!(f %in% studies_to_remove) & !(f %in% studies_empty)) {
    files_exist = c(files_exist, f)
  }
}

df <- read.table(sprintf("%s/%s.all_score", prs_path, files_exist[1]), header = T)
colnames(df) = gsub("e.", "e-", colnames(df))
number_of_snps = read.table(sprintf("%s/%s.prsice", prs_path, files_exist[1]), header = T, colClasses = c("character", "character","character", "numeric"))
number_of_snps$use = number_of_snps[,"Num_SNP"] > threshold
number_of_snps = number_of_snps[number_of_snps$use == T,]
th_columns = c("FID","IID",paste0("Pt_",as.character(number_of_snps$Threshold)))
df = df[,th_columns]
colnames(df) <- c("FID","IID",paste0(th_columns[3:ncol(df)], sprintf("_%s",files_exist[1])))

if(length(files_exist) > 1) {
  for(f in files_exist[2:length(files_exist)]) {
    df2 <- read.table(sprintf("%s/%s.all_score", prs_path, f), header = T)
    colnames(df2) = gsub("e.", "e-", colnames(df2))
    n_snps = read.table(sprintf("%s/%s.prsice", prs_path, f), header = T, colClasses = c("character", "character","character", "numeric"))
    n_snps$use = n_snps[,"Num_SNP"] > threshold
    n_snps = n_snps[n_snps$use == T,]
    tc = c("FID","IID",paste0("Pt_",as.character(n_snps$Threshold)))
    df2 = df2[,tc]
    colnames(df2) <- c("FID","IID",paste0(colnames(df2)[3:ncol(df2)], sprintf("_%s",f)))
    df <- merge(df,df2, by = c("FID","IID"))
  }
}

prs_calculated = rep(TRUE,length(files))
prs_calculated[match(studies_to_remove,files)] = FALSE
prs_calculated[match(studies_empty,files)] = FALSE

study_metadata$prs_calculated = prs_calculated

write.table(apply(study_metadata,2,as.character),paste0(study_list,".prs_calculated"), quote = F,col.names = T, row.names = F, sep = "\t")

write.table(df[,-1], paste0(prs_path,"/",name,"_",threshold,".tsv"), row.names = F, col.names = T, quote = F,sep="\t")
