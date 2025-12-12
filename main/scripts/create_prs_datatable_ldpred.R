rm(list=ls())
gc()

library(stringr)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

prs_path = args[1]
mode = args[2]
study_list = args[3]

studi_df = fread(study_list)

files = studi_df$study_id

studies_to_remove = c()
studies_empty = c()
for(f in files) {
  ff= sprintf("%s/%s.%s", prs_path,f, mode)
  if(!file.exists(ff) || file.size(ff) == 0L || is.na(file.size(ff))) {
    if(file.size(ff) == 0L || is.na(file.size(ff))) {
      studies_empty = c(studies_empty,f)
    } else {
      studies_to_remove = c(studies_to_remove,f)
    }
  } else {
    df <- read.table(ff, header = T)
    if(dim(df)[1] == 0) {
      studies_empty = c(studies_empty,f)
    } else {
      if(sum(is.na(df[,3])) == dim(df)[1]) {
        studies_empty = c(studies_empty,f)
      }
    }
  }
}


files_exist = c()
for(f in files) {
  if(!(f %in% studies_to_remove) & !(f %in% studies_empty)) {
    files_exist = c(files_exist, f)
  }
}

if(length(files_exist) != 0) {
  df <- read.table(sprintf("%s/%s.%s", prs_path,files_exist[1],mode), header = T)
  if(length(files_exist) >1) {
    for(f in files_exist[2:length(files_exist)]) {
      df2 <- read.table(sprintf("%s/%s.%s", prs_path,f,mode), header = T)
      df <- merge(df,df2, by = c("FID","IID"))
    }
  }
  write.table(df[,-1], paste0(prs_path,"/results.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
} else {
  write.table(data.frame(x="empty"), paste0(prs_path,"/results.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
}

