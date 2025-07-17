rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)
path_to_sumstats <- args[1]
studies_list <- args[2]
number_of_snps <- as.numeric(args[3])
complete_results = fread(studies_list)
studies = complete_results$study_id

if(length(studies) != 0) {
  number_of_lines = c()
  z_score_transformed = c()
  hm_readed = c()
  problem_beta = c()
  for(f in studies) {
    number_of_lines = c(number_of_lines, as.numeric(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.snpcount"))))
    z_score_transformed = c(z_score_transformed, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.z_score_converted"))))
    hm_readed = c(hm_readed, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.hm_readed"))))
    problem_beta = c(problem_beta, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.problematic_beta"))))
  }

  complete_results$include_subjects_from_UK = sapply(complete_results$countries_of_recruitment, function(x) {
    grepl("U.K.",x, fixed = TRUE)
  })

  complete_results$problematic_beta = problem_beta
  complete_results$z_score_transformed = z_score_transformed
  complete_results$hm_readed = hm_readed
  qc_p = rep(FALSE,length(studies))
  qc_p[which((number_of_lines > number_of_snps) & (z_score_transformed == FALSE) & (problem_beta == FALSE) & !(is.na(complete_results$sample_size)) & (complete_results$sample_size > 0))] = TRUE
  complete_results$qc_passed = qc_p
  complete_results$snp_count_current = number_of_lines
  write.table(complete_results,paste0(studies_list,".qced"), quote = F,col.names = T, row.names = F, sep = "\t")
} else {
  write.table(complete_results,paste0(studies_list,".qced"), quote = F,col.names = T, row.names = F, sep = "\t")
}
