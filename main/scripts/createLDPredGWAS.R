rm(list=ls())
gc()

library(dplyr)
library(data.table)
library(doParallel)
library(bigsnpr)
library(bigreadr)

cl <- makeCluster(5)
registerDoParallel(cl)

ldpred_ref = "/data/projects/on_going/PRScope/input/reference/ldpred2_ref/map_hm3_plus.rds"
gwas_table = fread("/data/projects/on_going/PRScope/output/gwas_list/gwas_search.qced_all.heritability.tsv")
gwas_table = gwas_table[gwas_table$heritability_passed == T,]
hg38 = TRUE

setwd("/data/projects/on_going/PRScope")

map_ldref = readRDS(ldpred_ref)

if(hg38==TRUE) {
  map_ldref$pos_hg19 = map_ldref$pos
  map_ldref$pos = map_ldref$pos_hg38
}

foreach(i=c(1:nrow(gwas_table))) %dopar% {
  require(data.table)
  require(bigsnpr)
  require(bigreadr)
  path = gwas_table[i,]$path
  study_id = gwas_table[i,]$study_id
  gwas = paste0(path,"/",study_id,".qced.h.tsv.gz")
  ss = fread(gwas,data.table = FALSE)
  ss$NEFF = gwas_table[i,]$neff
  colnames(ss) = tolower(colnames(ss))
  colnames(ss)[colnames(ss) == "a2"] = "a0"
  colnames(ss)[colnames(ss) == "bp"] = "pos"
  
  ss <- snp_match(ss, map_ldref)
  
  sd_ldref <- with(ss, sqrt(2*maf*(1 - maf)))
  sd_ss <- with(ss, 2 / sqrt(neff * se^2 + beta^2))
  
  is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05
  ss <- ss[!is_bad, ]
  ss <- ss[ , c("varid","snp","chr","pos","a1","a0","p","maf","beta","or","se","n","neff")] 
  colnames(ss)[colnames(ss) == "a0"] = "a2"
  colnames(ss)[colnames(ss) == "pos"] = "bp"
  colnames(ss) = toupper(colnames(ss))
  ldpred_ss = paste0(dirname(path),"/ldpred")
  
  if(!file.exists(ldpred_ss)) {
    dir.create(ldpred_ss)
  }
  
  fwrite(ss, file = paste0(ldpred_ss,"/",study_id,".ldpred.h.tsv.gz"), quote = F, row.names = F, sep = "\t", compress = "gzip")
}

stopCluster(cl)
