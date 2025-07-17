rm(list=ls())
gc()
library(gwasrapidd)
library(data.table)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

x = fread("/zi/home/ersoy.kocak/Desktop/Projects/CoviDrug/study_lists/GWAS_Sample_Sizes_ALL.tsv")

trait.names = x$V1

ld <- "/zi/home/ersoy.kocak/Desktop/Projects/Gen_Cor/references/eur_w_ld_chr_hg38/"
wld <- ld

list_of_munged_files = c()
for(f in trait.names) {
  list_of_munged_files = c(list_of_munged_files,paste0("/zi/home/ersoy.kocak/Desktop/Data/munged_new/", f,".sumstats.gz") )
}

lines = c()
for(i in c(1:length(list_of_munged_files))) {

  for(j in c(i:length(list_of_munged_files))) {
    munged_file_1 = list_of_munged_files[i]
    munged_file_2 = list_of_munged_files[j]
    name_1 = trait.names[i]
    name_2 = trait.names[j]
    if(name_1 != name_2) {
      lines = c(lines,paste0("/zi/home/ersoy.kocak/Desktop/Tools/ldsc/ldsc.py --rg ", paste0(munged_file_1,",",munged_file_2), " --ref-ld-chr ", ld, " --w-ld-chr ", wld, " --out  /zi/home/ersoy.kocak/Desktop/Projects/Gen_Cor/results/",name_1,"_",name_2," &"))
      if(length(lines) %% 150 == 0) {
        lines = c(lines,"wait")
      }
    }
  }
}

writeLines(lines,"/zi/home/ersoy.kocak/Desktop/Projects/Gen_Cor/gen_cor_command_new.sh")


