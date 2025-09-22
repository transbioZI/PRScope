rm(list=ls())
gc()

library(tidyverse)
library(data.table)
library(optparse)
library(dplyr)
library(lattice)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

imputation_directory = args[1]

outputPath = paste0(imputation_directory,"/cobg_dir_genome_wide")
outputName = "infos.RDS"

list.filenames <- list.files(path = imputation_directory,pattern="dos.*\\.out\\.dosage\\.info$", recursive = T, full.names = T)

read_all = list()
for(i in c(1:length(list.filenames))) {
  read_all[[i]] = read.table(list.filenames[i],sep="\t",header = TRUE)
}

infos <- do.call(rbind , read_all)


system(paste0("mkdir -p ",outputPath,"/infos"))

saveRDS(infos, paste0(outputPath,"/infos/",outputName))

infos_filtered_08 = infos[infos$info > 0.8,]
infos_filtered_06 = infos[infos$info > 0.6,]
infos_filtered_04 = infos[infos$info > 0.4,]

writeLines(infos_filtered_08$SNP,paste0(outputPath,"/infos/infos_08.extract"))
writeLines(infos_filtered_06$SNP,paste0(outputPath,"/infos/infos_06.extract"))
writeLines(infos_filtered_04$SNP,paste0(outputPath,"/infos/infos_04.extract"))

plink_file = list.files(path = outputPath,pattern="\\.bg.bim$", recursive = T, full.names = T)
filename = basename(plink_file)
filename = str_split(filename,"\\.bim")[[1]][1]
system(paste0("/zi/home/ersoy.kocak/Desktop/Tools/plink --bfile ",paste0(outputPath,"/",filename)," --make-bed --extract ",paste0(outputPath,"/infos/infos_06.extract")," --out ",paste0(outputPath,"/infos/",filename,"_info_06")))
