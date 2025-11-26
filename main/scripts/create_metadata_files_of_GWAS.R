rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

apply_which_false = function(condition, threshold, criteria) {
    return(sapply(condition,function(x) {
        if(x == TRUE) {
          ""
        } else {
          if(!is.na(threshold)) {
            paste0(criteria,": FAIL, threshold: ",as.character(threshold))
          } else {
            paste0(criteria,": FAIL")
          }
        }
    }))
}

args = commandArgs(trailingOnly=TRUE)

path_to_sumstats <- args[1]
studies_list <- args[2]
number_of_snps <- as.numeric(args[3])
output_path = args[4]
complete_results = fread(studies_list,data.table = F)
studies = complete_results$study_id
paths = complete_results$path
if(length(studies) != 0) {
  metadata_list = list()
  for(f in studies) {
    ai = which(f == studies)
    if(is.na(paths[ai]) | length(paths[ai]) == 0 ) {
      path_ss = paste0(path_to_sumstats,"/",f)
    } else {
      path_ss = paths[ai]
    }
    res = fread(paste0(path_ss,"/qced/",f,".metadata.tsv"))

    colnames(res) = c("N_failed", "P_failed", "BETA_failed","BETA_failed_reason", "MAF_matching_failed", "HM_prefix_readed","SNP_count_after_QC",
                                "BETA_is_Z_SCORE", "SE_failed", "CHR_failed","SNP_count_before_QC","duplicate_columns_dropped","columns_before_QC","columns_after_QC","ld_pred_failed")

    md1 = readLines(paste0(path_ss,"/raw/",f,"_md5sum_download.txt"))
    md1 = str_split(md1[1],"\\s+")[[1]][1]
    if(file.exists(paste0(path_ss,"/raw/",f,"_md5sum_real.txt"))) {
    	md2 = readLines(paste0(path_ss,"/raw/",f,"_md5sum_real.txt"))
    	md2 = md2[grepl("*.h.tsv.gz$",md2)]
    	md2 = str_split(md2[1],"\\s+")[[1]][1]
    } else {
    	md2 = "unknown"
    }

    res$md5sum_download = md1
    res$md5sum_original = md2
    res$identical_md5sum = identical(md1,md2)
    res$study_id = f
    res$path = path_ss
    metadata_list[[ai]] = res
  }

  qc_report_gwas = do.call(rbind,metadata_list)
  complete_results$path = NULL
  complete_results$include_subjects_from_UK = sapply(complete_results$countries_of_recruitment, function(x) {
    grepl("U.K.",x, fixed = TRUE)
  })

  merge_with_study_file = left_join(complete_results, qc_report_gwas, by="study_id")

  SNP_count_after_QC = as.numeric(merge_with_study_file$SNP_count_after_QC)
  P_failed = as.logical(merge_with_study_file$P_failed)
  BETA_failed = as.logical(merge_with_study_file$BETA_failed)
  CHR_failed = as.logical(merge_with_study_file$CHR_failed)
  SE_failed = as.logical(merge_with_study_file$SE_failed)
  SAMPLE_SIZE = as.numeric(merge_with_study_file$sample_size)
  N_failed = as.logical(merge_with_study_file$N_failed)
  MAF_matching_failed = as.logical(merge_with_study_file$MAF_matching_failed)
  BETA_is_Z_SCORE = as.logical(merge_with_study_file$BETA_is_Z_SCORE)
  ld_pred_failed = as.logical(merge_with_study_file$ld_pred_failed)
  qc_p = rep(FALSE,length(studies))
  qc_p[which((SNP_count_after_QC > number_of_snps) &
             (P_failed == FALSE) &
             (BETA_failed == FALSE) &
             (CHR_failed == FALSE) &
             !(is.na(SAMPLE_SIZE)) & (MAF_matching_failed == FALSE) &
             (SAMPLE_SIZE > 0) )] = TRUE
  qc_p_ldpred = rep(FALSE,length(studies))
  qc_p_ldpred[which((SNP_count_after_QC > number_of_snps) &
             (P_failed == FALSE) &
             (BETA_failed == FALSE) &
             (CHR_failed == FALSE) &
             (SE_failed == FALSE) &
             (N_failed == FALSE) &
             !(is.na(SAMPLE_SIZE)) & (MAF_matching_failed == FALSE) & (ld_pred_failed == FALSE) &
             (SAMPLE_SIZE > 0))] = TRUE
  a = apply_which_false(SNP_count_after_QC > number_of_snps, number_of_snps, "Number of SNPs")
  b = apply_which_false(CHR_failed == FALSE, NA, "CHR")
  c = apply_which_false(P_failed == FALSE, NA, "P")
  d = apply_which_false(BETA_failed == FALSE, NA, "BETA")
  e = apply_which_false(!(is.na(SAMPLE_SIZE)), "missingness", "SAMPLE SIZE")
  f = apply_which_false(SAMPLE_SIZE > 0, 0, "SAMPLE SIZE")
  g = apply_which_false(N_failed == FALSE, NA, "N")
  k = apply_which_false(SE_failed == FALSE, NA, "SE")
  kk = apply_which_false(ld_pred_failed == FALSE, NA, "ld_pred_failed")
  kka = apply_which_false(MAF_matching_failed == FALSE, NA, "MAF_matching_failed")
  all_criteria = data.frame(a = a, b = b, c = c, d = d, e = e, f = f, g = g, k = k, kk = kk, kka = kka)

  comment_qc = unlist(apply(all_criteria, 1, function(x) {
      str = as.character(x[ x!=""])
      xtr = paste0(str, collapse = " and ")
      if(!is.na(xtr) & xtr !="") {
        paste0(xtr," (first QC) ")
      } else {
         xtr
      }
  }))

  merge_with_study_file$qc_passed_comment = comment_qc
  merge_with_study_file$qc_passed = qc_p
  merge_with_study_file$qc_passed_ldpred = qc_p_ldpred
  write.table(merge_with_study_file,paste0(output_path,".qced.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")
} else {
  write.table(merge_with_study_file,paste0(output_path,".qced.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")
}
