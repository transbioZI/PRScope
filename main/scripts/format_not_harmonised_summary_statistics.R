rm(list=ls())

library(tidyverse)
library(data.table)
library(optparse)
library(dplyr)


library("stringr")
library("readxl")

format_numeric <- function(df, column) {
  df[,column] = gsub("[^0-9| ^e | ^E | ^.| ^-]", "", df[,column])
  df = df[df[,column] != "[e | E]",] 
  df = df[df[,column] != ".",] 
  df = df[df[,column] != "-",]
  df = df[df[,column] != "",]
  df = df[!is.na(df[,column]),]
  df = df[!is.na(as.numeric(df[,column])),]
  return(df)
}
coll = c()
download_path="/data/GWAS_Catalog_Not_Harmonized"
formatted = "/data/GWAS_Catalog_Not_Harmonized/formatted"
setwd(formatted)
files = list.files(formatted, pattern = "*.tsv")
c = c()

for(f in files) {
  if(!grepl("meta",f,fixed = T)) {
    c = c(c, paste0("nextflow run  EBISPOT/gwas-sumstats-harmoniser --ref /data/harmonizer_reference --harm --file ",paste0(formatted,"/",f)," -profile standard,docker -r v1.1.1"))
    
  }
  #print()
 
}
c = c[-c(1:28)]
writeLines(c,"/data/GWAS_Catalog_Not_Harmonized/harmonizing_command_5.sh")

  df = as_tibble(fread(f))

  colnames(df)[colnames(df) == "variant_id"] = "rsid"
  colnames(df)[colnames(df) == "rs_id_all"] = "rsid"
  colnames(df)[colnames(df) == "rs_id"] = "rsid"
  colnames(df)[colnames(df) == "A1"] = "effect_allele"
  colnames(df)[colnames(df) == "EA"] = "effect_allele"
  colnames(df)[colnames(df) == "A2"] = "other_allele"
  colnames(df)[colnames(df) == "OA"] = "other_allele"
  colnames(df)[colnames(df) == "non_effect_allele"] = "other_allele"
  colnames(df)[colnames(df) == "SE"] = "standard_error"
  if("BETA" %in% colnames(df)) {
    colnames(df)[colnames(df) == "BETA"] = "beta"
  } else if("Z" %in% colnames(df) | "z" %in% colnames(df)) {
    colnames(df)[colnames(df) == "Z"] = "beta"
    colnames(df)[colnames(df) == "z"] = "beta"
  }
  
  colnames(df)[colnames(df) == "OR"] = "odds_ratio"
  colnames(df)[colnames(df) == "chr"] = "chromosome"
  colnames(df)[colnames(df) == "CHR"] = "chromosome"
  colnames(df)[colnames(df) == "CHROM"] = "chromosome"
  colnames(df)[colnames(df) == "chrom"] = "chromosome"
  colnames(df)[colnames(df) == "CHROMOSOME"] = "chromosome"
  colnames(df)[colnames(df) == "BP"] = "base_pair_location"
  colnames(df)[colnames(df) == "POS"] = "base_pair_location"
  colnames(df)[colnames(df) == "bp"] = "base_pair_location"
  colnames(df)[colnames(df) == "pos"] = "base_pair_location"
  colnames(df)[colnames(df) == "pvalue"] = "p_value"
  colnames(df)[colnames(df) == "EAF"] = "effect_allele_frequency"
  colnames(df)[colnames(df) == "info_all"] = "info"
  mandatory_columns = c('chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'p_value','rsid')
  
  columns = c()
  if(all(mandatory_columns %in% names(df))) {
    columns = c(mandatory_columns, columns)
  } else {
    #file.remove(f)
    #file.remove(paste0(f,"-meta.yaml"))
    coll = c(coll,colnames(df))
    print(f)
    print(colnames(df))
    next
  }
  
  if("beta" %in% names(df)) {
    columns = c("beta", columns)
    #df = format_numeric(df,"beta")
  }
  
  if("odds_ratio" %in% names(df)) {
    columns = c("odds_ratio", columns)
    #df = format_numeric(df,"odds_ratio")
  }
  if("standard_error" %in% names(df)) {
    columns = c("standard_error", columns)
    #df = format_numeric(df,"standard_error")
  }
  if("effect_allele_frequency" %in% names(df)) {
    columns = c("effect_allele_frequency", columns)
    #df = format_numeric(df,"effect_allele_frequency")
  }
  if("info" %in% names(df)) {
    columns = c("info", columns)
    #df = format_numeric(df,"info")
    #df[df$info > 1.0,]$info = 1.0
  }
  
  df = df[nchar(as.character(df$effect_allele))==1,]
  df = df[nchar(as.character(df$other_allele))==1,]
  
  df = df[grep('^rs', df$rsid),]
  
  
  
  df = df[with(df, order(as.numeric(chromosome), as.numeric(base_pair_location))), ]
  write.table(df, paste0(formatted,"/",f), quote = F, row.names =F, sep = "\t")
  #file.copy(paste0(f,"-meta.yaml"), formatted)
  file.remove(f)
  #file.remove(paste0(f,"-meta.yaml"))
}
