rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))
suppressMessages(library(stringr))
suppressMessages(library(data.table))


apply_which_false = function(condition, threshold, criteria) {
    return(sapply(condition,function(x) {
        ifelse(x, "", paste0(criteria,": FAILED threshold: ",as.character(threshold)))
    }))
}

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
  problem_p_value = c()
  problematic_N = c()
  for(f in studies) {
    problematic_N = c(problematic_N, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.problematic_N"))))
    problem_p_value = c(problem_p_value, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.problematic_p_value"))))
    number_of_lines = c(number_of_lines, as.numeric(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.snpcount"))))
    z_score_transformed = c(z_score_transformed, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.z_score_converted"))))
    hm_readed = c(hm_readed, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.hm_readed"))))
    problem_beta = c(problem_beta, as.character(readLines(paste0(path_to_sumstats,"/",f,".qced.h.tsv.gz.problematic_beta"))))
  }

  complete_results$include_subjects_from_UK = sapply(complete_results$countries_of_recruitment, function(x) {
    grepl("U.K.",x, fixed = TRUE)
  })

  complete_results$problematic_p_value = problem_p_value
  complete_results$problematic_beta = problem_beta
  complete_results$z_score_transformed = z_score_transformed
  complete_results$hm_readed = hm_readed
  complete_results$problematic_N = problematic_N
  qc_p = rep(FALSE,length(studies))
  qc_p[which((number_of_lines > number_of_snps) & (z_score_transformed == FALSE) & (problem_p_value == FALSE) & (problem_beta == FALSE) & !(is.na(complete_results$sample_size)) & (complete_results$sample_size > 0))] = TRUE

  a = apply_which_false(number_of_lines > number_of_snps, number_of_snps, "number_of_snps")
  b = apply_which_false(z_score_transformed == FALSE, "TRUE", "z_score_transformed")
  c = apply_which_false(problem_p_value == FALSE, "TRUE", "problem_p_value")
  d = apply_which_false(problem_beta == FALSE, "TRUE", "problem_beta")
  e = apply_which_false(!(is.na(complete_results$sample_size)), "NA", "sample_size")
  f = apply_which_false(complete_results$sample_size > 0, 0, "sample_size")
  g = apply_which_false(problematic_N == FALSE, "TRUE", "problematic_N")

  all_criteria = data.frame(a= a,b =b, c = c, d=d,e=e,f=f, g = g)

  comment_qc = unlist(apply(all_criteria, 1, function(x) {
      str = as.character(x[ x!=""])
      xtr = paste0(str, collapse = " and ")
      if(!is.na(xtr) & xtr !="") {
        paste0(xtr," (first QC) ")
      } else {
         xtr
      }

  }))

  complete_results$qc_passed_comment = comment_qc

  complete_results$qc_passed = qc_p
  complete_results$snp_count_current = number_of_lines
  write.table(complete_results,paste0(studies_list,".qced"), quote = F,col.names = T, row.names = F, sep = "\t")
} else {
  write.table(complete_results,paste0(studies_list,".qced"), quote = F,col.names = T, row.names = F, sep = "\t")
}
