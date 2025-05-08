rm(list=ls())
gc()

library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

prs_path = args[1]
mode = args[2]
study_list = args[3]
intercept_upper_thr = 100
intercept_lower_thr = 0
heritability_thr = 0

files = unique(readLines(study_list))

get_filename = function(x) {
  return(str_split(x,"\t")[[1]][1])
}

files = unlist(lapply(files, get_filename))

studies_to_remove = c()
studies_empty = c()
for(f in files) {
  ff= sprintf("%s/%s.%s", prs_path, f,mode)
  if(!file.exists(ff) || file.size(ff) == 0L || is.na(file.size(ff))) {
    if(file.size(ff) == 0L || is.na(file.size(ff))) {
      studies_empty = c(studies_empty,f)
    } else {
      studies_to_remove = c(studies_to_remove,f)
    }
  } else {
    df <- read.table(sprintf("%s/%s.%s", prs_path, f,mode), header = T)
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

remove_low_her_and_ind = c()
for(f in files_exist) {
  her = readLines(sprintf("%s/%s.her", prs_path, f))
  intercept = as.numeric(strsplit(her[2],":")[[1]][2])
  heritability = as.numeric(strsplit(her[4],":")[[1]][2])
  heritability_se = as.numeric(strsplit(her[5],":")[[1]][2])

  if(heritability/heritability_se < heritability_thr) {
    remove_low_her_and_ind = c(remove_low_her_and_ind, f)
  } else if(intercept > intercept_upper_thr || intercept < intercept_lower_thr ) {
    remove_low_her_and_ind = c(remove_low_her_and_ind, f)
  }
}

files_exist_2 = c()
for(f in files_exist) {
  if(!(f %in% remove_low_her_and_ind)) {
    files_exist_2 = c(files_exist_2, f)
  }
}

df <- read.table(sprintf("%s/%s.%s", prs_path, files_exist_2[1],mode), header = T)
for(f in files_exist_2[2:length(files_exist_2)]) {
  df2 <- read.table(sprintf("%s/%s.%s", prs_path, f, mode), header = T)
  df <- merge(df,df2, by = c("FID","IID"))
}

write.table(df[,-1], paste0(prs_path,"/resultPRSs.tsv"), row.names = F, col.names = T, quote = F)

PRS_data = read.table(paste0(prs_path,"/resultPRSs.tsv"), sep=" ", header=TRUE)
print(dim(PRS_data))
#select(PRS_data, starts_with("Pt_0.0001"))

pca = prcomp(PRS_data[,-1])$x[,1:2]
pdf(paste0(prs_path,"/PCA_PRS.pdf"))
plot(pca[,1], pca[,2])
dev.off()


if(length(studies_to_remove) !=0) {
  writeLines(studies_to_remove, paste0(prs_path,"/studies_not_included.txt"))
}

if(length(studies_empty) !=0) {
  writeLines(studies_empty, paste0(prs_path,"/studies_empty.txt"))
}

if(length(remove_low_her_and_ind) !=0) {
  writeLines(remove_low_her_and_ind, paste0(prs_path,"/studies_threshold.txt"))
}

writeLines(files_exist_2, paste0(prs_path,"/studies_included.txt"))
print("Number of studies included:")
print(length(files_exist_2))
print("Number of studies not included:")
print(length(studies_to_remove))
print("Number of studies empty:")
print(length(studies_empty))
print("Number of studies thresholding:")
print(length(remove_low_her_and_ind))
print("Samples X PRSs:")
print(dim(df[,-1]))
