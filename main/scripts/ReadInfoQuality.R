rm(list=ls())
gc()

library(tidyverse)
library(data.table)
library(optparse)
library(dplyr)
library(lattice)
library(doParallel)
args = commandArgs(trailingOnly=TRUE)

infoFilesPath = "/data/data/TRR/imputation/QCed/qc/qc/test_impute/qc/imputation/dasuqc1_deleteme.hg19.ch.fl/info"
outputPath = "/data/data/TRR/imputation/QCed/qc/qc/test_impute/qc/imputation/cobg_dir_genome_wide"
outputName = "infos"

list.filenames <- list.files(path = infoFilesPath,pattern="\\.info$")

cl <- makeCluster(30)
registerDoParallel(cl)
trialsdata = foreach(i=1:length(list.filenames)) %dopar% {
  read.table(paste0(infoFilesPath,"/",list.filenames[i]),sep="\t",header = TRUE)
}
stopCluster(cl)
infos <- do.call(rbind , trialsdata)

saveRDS(infos, paste0(outputPath,"/",outputName))

infos_filtered_08 = infos[infos$info > 0.8,]
infos_filtered_06 = infos[infos$info > 0.6,]
infos_filtered_04 = infos[infos$info > 0.4,]

writeLines(infos_filtered_08$SNP,paste0(outputPath,"/infos_08.extract"))
writeLines(infos_filtered_06$SNP,paste0(outputPath,"/infos_06.extract"))
writeLines(infos_filtered_04$SNP,paste0(outputPath,"/infos_04.extract"))

