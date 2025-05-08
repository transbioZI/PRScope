rm(list=ls())
gc()

library(MASS)
library(igraph)
library(dplyr)
library(liver)
library(stringr)
library(data.table)
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))
suppressMessages(library(rtracklayer))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))
input = "/data/projects/on_going/multiprs_patient_stratification/study_lists/pipeline_studies_final_list_without_other_sources.txt"
output = paste0(input,".tsv")
pgc_list = "/data/sumstats/PGC_GWASs/pgc_list.csv"
studies = unique(readLines(input))

pgc_studies = fread(pgc_list, header = T, data.table = F)

complete_results = data.frame(matrix(ncol = 10, nrow = length(studies)))
colnames(complete_results) <- c('study_id', 'efo_trait', 'reported_trait', "initial_sample", "replication_sample","snp_count",
                                "publication_date","pubmed_id","title","publication")
complete_results$study_id = studies
complete_results$efo_trait = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE) | grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    pgc_studies[pgc_studies$id == x,]$study
  }  else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    trait = get_traits(study_id = st)@traits$trait
    paste(trait, collapse = '/')
  }  
  else{
    trait = get_traits(study_id = x)@traits$trait
    paste(trait, collapse = '/')
  }
  
})
complete_results$reported_trait = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    pgc_studies[pgc_studies$id == x,]$link
  }  else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    get_studies(study_id = st)@studies$reported_trait
  }  
  else{
    get_studies(study_id = x)@studies$reported_trait
  }
})
complete_results$initial_sample = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    x
  }  else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    get_studies(study_id = st)@studies$initial_sample_size
  } 
  else{
    get_studies(study_id = x)@studies$initial_sample_size
  }
})
complete_results$replication_sample = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    x
  }  else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    get_studies(study_id = st)@studies$replication_sample_size
  }  
  else{
    get_studies(study_id = x)@studies$replication_sample_size
  }
})
complete_results$snp_count = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    x
  }  else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    get_studies(study_id = st)@studies$snp_count
  }  
  else{
    get_studies(study_id = x)@studies$snp_count
  }
})
complete_results$publication_date = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    x
  }  else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    as.character(get_studies(study_id = st)@publications$publication_date)
  }  
  else{
    as.character(get_studies(study_id = x)@publications$publication_date)
  }
})
complete_results$pubmed_id = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    x
  } else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    get_studies(study_id = st)@publications$pubmed_id
  }  
  else{
    get_studies(study_id = x)@publications$pubmed_id
  }
})
complete_results$title = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    x
  } else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    get_studies(study_id = st)@publications$title
  }  
  else{
    get_studies(study_id = x)@publications$title
  }
})
complete_results$publication = sapply(studies, function(x) {
  if(grepl("transbio", x, fixed = TRUE) | grepl("PGC", x, fixed = TRUE) | grepl("CTG", x, fixed = TRUE)| grepl("BIG", x, fixed = TRUE) | grepl("ALT", x, fixed = TRUE)) {
    x
  } else if(grepl("build", x, fixed = TRUE)) {
    st = str_split(x,"_")[[1]][1]
    get_studies(study_id = st)@publications$publication
  } 
  else{
    get_studies(study_id = x)@publications$publication
  }
})

write.table(apply(complete_results,2,as.character),output, quote = F,col.names = T, row.names = F, sep = "\t")
