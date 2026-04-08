rm(list=ls())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(gwasrapidd))
suppressMessages(library(bigsnpr))
suppressMessages(library(stats))
args = commandArgs(trailingOnly=TRUE)
    
input <- args[1]
output <- args[2]
maf_file <- args[3]
N = as.numeric(args[4])
output_other_files = args[5]
NEFF = as.numeric(args[6])

base <- fread(input, showProgress = FALSE, data.table = F)

write.table(summary(base),paste0(output_other_files,".summary.before_qc.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

problematic_beta = FALSE
beta_reason = NA
columns_dropped = FALSE
columns_before = colnames(base)
if(any(table(colnames(base)) > 1)) {
  base[,names(table(colnames(base))[table(colnames(base)) > 1])] = NULL # drop duplicate columns
  columns_dropped = TRUE
}

drop_columns = c("BETA","Z_SCORE") # it should be lowercase beta or z_score

df_colnames = colnames(base)
if(any(is.element(df_colnames, drop_columns))) {
  base[,which(is.element(df_colnames, drop_columns) == TRUE)] = NULL
  columns_dropped = TRUE
}

colnames(base) = toupper(colnames(base))

if(any(table(colnames(base)) > 1)) {
  base[,names(table(colnames(base))[table(colnames(base)) > 1])] = NULL # drop duplicate columns
  columns_dropped = TRUE
}

drop_columns = c("SNP","CHR","BP","A1","A2","OR","P","MAF","VARID")
df_colnames = colnames(base)
if(any(is.element(df_colnames, drop_columns))) {
  base[,which(is.element(df_colnames, drop_columns) == TRUE)] = NULL
  columns_dropped = TRUE
}

columns_after = colnames(base)
snpcount_before = as.character(dim(base)[1])

hm_readed = FALSE
if("HM_BETA" %in% colnames(base)) {
  colnames(base)[which(names(base) == "HM_RSID")] <- "SNP"
  colnames(base)[which(names(base) == "HM_CHROM")] <- "CHR"
  colnames(base)[which(names(base) == "HM_POS")] <- "BP"
  colnames(base)[which(names(base) == "HM_OTHER_ALLELE")] <- "A2"
  colnames(base)[which(names(base) == "HM_EFFECT_ALLELE")] <- "A1"
  if("BETA" %in% colnames(base)) {
    if(all(is.na(base$HM_BETA))) {
      problematic_beta = TRUE
      beta_reason = "hm_beta is NA"
    }
    base$BETA = NULL
  }
  colnames(base)[which(names(base) == "HM_BETA")] <- "BETA"
  colnames(base)[which(names(base) == "HM_ODDS_RATIO")] <- "OR"
  hm_readed = TRUE
} else {

  if(!("BETA" %in% colnames(base))) {
    base$BETA = NA
  }
  if(!("ODDS_RATIO" %in% colnames(base))) {
    base$ODDS_RATIO = NA
  }

  colnames(base)[which(names(base) == "RSID")] <- "SNP"
  colnames(base)[which(names(base) == "CHROMOSOME")] <- "CHR"
  colnames(base)[which(names(base) == "BASE_PAIR_LOCATION")] <- "BP"
  colnames(base)[which(names(base) == "OTHER_ALLELE")] <- "A2"
  colnames(base)[which(names(base) == "EFFECT_ALLELE")] <- "A1"
  colnames(base)[which(names(base) == "BETA")] <- "BETA"
  colnames(base)[which(names(base) == "ODDS_RATIO")] <- "OR"
}

problematic_CHR = FALSE
if(length(table(base$CHR)) < 22) {
  problematic_CHR = TRUE
}

if(problematic_beta == FALSE & problematic_CHR == FALSE ) {

  if(("STANDARD_ERROR" %in% colnames(base)) & ("SE" %in% colnames(base))) {
    if(all(is.na(base$STANDARD_ERROR)) & all(is.na(base$SE))) {
      base$SE = NULL
    }

    if(!all(is.na(base$STANDARD_ERROR)) & !all(is.na(base$SE))) {
      base$SE = NULL
    }

    if(!all(is.na(base$STANDARD_ERROR)) & all(is.na(base$SE))) {
      base$SE = NULL
    }

    if(all(is.na(base$STANDARD_ERROR)) & !all(is.na(base$SE))) {
      base$STANDARD_ERROR = base$SE
      base$SE = NULL
    }
  }

  if(!("STANDARD_ERROR" %in% colnames(base))) {
    if("SE" %in% colnames(base)) {
      colnames(base)[which(names(base) == "SE")] <- "STANDARD_ERROR"
    } else {
      base$STANDARD_ERROR = NA
    }
  }

  colnames(base)[which(names(base) == "P_VALUE")] <- "P"
  colnames(base)[which(names(base) == "STANDARD_ERROR")] <- "SE"

  base$MAF = 0
  base$VARID <- str_c(base$CHR, ":", base$BP,"_",base$A1,"_",base$A2)

  if(dim(base)[1] != 0 ) {
    maf = as.data.frame(fread(maf_file, select = c("SNP","MAF_hg38"), showProgress = FALSE))
    colnames(maf)[colnames(maf) == "MAF_hg38"] = "MAF"
    matches = match(base$SNP,maf$SNP)
    matches_na = which(is.na(matches))
    matches_not_na = which(!is.na(matches))
    if(length(matches_not_na) > 0) {
      base[matches_not_na,]$MAF = maf$MAF[matches[!is.na(matches)]]
    }
  }
  rm(maf)
  tmp_l = dim(base)[1]
  base <- dplyr::filter(base, MAF > 0.001)
  tmp_l2 = dim(base)[1]

  problematic_MAF_match = FALSE
  if(tmp_l2 == 0 & tmp_l > 0) {
    problematic_MAF_match = TRUE
  }

  if ("N" %in% colnames(base)) {
    if(all(is.na(base$N))) {
      base$N = N
    }
  } else {
    base$N = N
  }

  base$N = suppressWarnings(as.numeric(base$N))
  tmp_l = dim(base)[1]
  base = dplyr::filter(base, !(is.na(N)))
  base = dplyr::filter(base, N>0)
  tmp_l2 = dim(base)[1]
  problematic_N = FALSE

  if(tmp_l2 == 0 & tmp_l > 0) {
    problematic_N = TRUE
  }

  z_score_studies = FALSE
  if(all(is.na(base$BETA)) & all(is.na(base$OR))) {
     if(("Z_SCORE" %in% colnames(base)) & problematic_N == FALSE & problematic_MAF_match == FALSE) {
       base$Z_SCORE = suppressWarnings(as.numeric(base$Z_SCORE))
       base = dplyr::filter(base, !(is.na(Z_SCORE)))
       base = dplyr::filter(base, (is.finite(Z_SCORE) == TRUE))
       if( dim(base)[1] > 0 ) {
        base$BETA = base$Z_SCORE / sqrt( 2 * base$MAF * ( 1 - base$MAF) * ( base$N + base$Z_SCORE^2 ) )
        base$SE = 1 / sqrt( 2 * base$MAF * ( 1 - base$MAF ) * ( base$N + base$Z_SCORE^2 ) )
        z_score_studies = TRUE
       }
     }
  }

  if ("NEFF" %in% colnames(base)) {
    if(all(is.na(base$N))) {
      base$NEFF = NEFF
    }
  } else {
    base$NEFF = NEFF
  }

  base$NEFF = suppressWarnings(as.numeric(base$NEFF))
  tmp_l = dim(base)[1]
  base = dplyr::filter(base, !(is.na(NEFF)))
  base = dplyr::filter(base, NEFF>0)
  tmp_l2 = dim(base)[1]
  problematic_N = FALSE

  if(tmp_l2 == 0 & tmp_l > 0) {
    problematic_N = TRUE
  }

  base = base %>% select(SNP,CHR,BP,A1,A2,BETA,OR,P,SE,N,NEFF,MAF,VARID)
  
  #### Remove SNPs with no beta or OR - these cannot be used by PRSice ####
  base <- dplyr::filter(base, !(is.na(BETA) == TRUE & is.na(OR) == TRUE))

  if(dim(base)[1] == 0 ) {
      problematic_beta = TRUE
      beta_reason = "BETA is NA"
  }

  base <- subset(base, nchar(as.character(A1)) == 1)
  base <- subset(base, nchar(as.character(A2)) == 1)

  base <- dplyr::filter(base, !(A1 == "A" & A2 == "T"))
  base <- dplyr::filter(base, !(A1 == "T" & A2 == "A"))
  base <- dplyr::filter(base, !(A1 == "G" & A2 == "C"))
  base <- dplyr::filter(base, !(A1 == "C" & A2 == "G"))

  base$BETA <- suppressWarnings(as.numeric(base$BETA))
  base$OR <- suppressWarnings(as.numeric(base$OR))
  base$CHR = suppressWarnings(as.numeric(base$CHR))
  base = dplyr::filter(base, !(is.na(CHR)))
  base = dplyr::filter(base, (CHR %in% c(1:22)))

  base <- base %>%
    mutate(BETA = if_else(is.na(BETA), log(OR), BETA),
           OR = if_else(is.na(OR), exp(BETA), OR))

  base <- dplyr::filter(base, (is.finite(BETA) == TRUE))
  base <- dplyr::filter(base, (is.finite(OR) == TRUE))

  if(dim(base)[1] == 0 ) {
      problematic_beta = TRUE
      beta_reason = "BETA is infinite"
  }

  #### For SNPs with no p-value, replace the NA by a 1 ####
  nan_values = which(base$P == "NaN")
  if(length(nan_values) > 0) {
    base[nan_values,]$P = NA
  }
  base$P = gsub("e-$","e-01",base$P)
  base$P = gsub("e$","e-01",base$P)
  base <- base %>% mutate(P = if_else(is.na(P), 1, as.numeric(P)))
  #### Remove duplicate SNPs - these cannot be used by PRSice (an error will be thrown) ####
  base <- distinct(base, SNP, .keep_all = TRUE)

  if(dim(base)[1] != 0 ) {
    if(abs(median(base$BETA)) > 0.5 ) {
      problematic_beta = TRUE
      beta_reason = "BETA valut does not have expected range."
    }
  }

  problematic_p_value = FALSE
  if(dim(base)[1] != 0 ) {
    p_1 = c(which(base$P > 1), which(base$P < 0))
    if(length(p_1) != 0) {
      base = base[-p_1,]
    }
    
    if(dim(base)[1] == 0) {
      problematic_p_value = TRUE
    }
  }

  problematic_SE = FALSE
  if(all(is.na(base$SE))) {
      problematic_SE = TRUE
  }

  if(length(table(base$CHR)) !=22) {
    problematic_CHR = TRUE
  }

  dfx = data.frame(matrix(ncol = 15, nrow = 1))
  colnames(dfx) <- c("problematic_N", "problematic_p_value", "problematic_beta","problematic_beta_reason", "problematic_MAF_match", "hm_readed","snpcount",
                                  "z_score_converted", "problematic_SE", "problematic_CHR","snpcount_in_raw_file","columns_dropped","columns_before","columns_after","ld_pred_failed")

  dfx[1,]= c(problematic_N,problematic_p_value, problematic_beta,beta_reason,problematic_MAF_match,hm_readed ,dim(base)[1], z_score_studies,problematic_SE,problematic_CHR,snpcount_before,columns_dropped,paste0(columns_before,collapse = ","),paste0(columns_after,collapse = ","),FALSE)

  write.table(dfx,paste0(output_other_files,".metadata.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

  base = base[,c("VARID","SNP","CHR","BP","A1","A2","P","MAF","BETA","OR","SE","N","NEFF")]

  write.table(summary(base),paste0(output_other_files,".summary.qced.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

  fwrite(base, file = output, quote = F, row.names = F, sep = "\t", compress = "gzip")

  writeLines(as.character(warnings()),paste0(output_other_files,".warnings"))
} else {

  dfx = data.frame(matrix(ncol = 15, nrow = 1))
  colnames(dfx) <- c("problematic_N", "problematic_p_value", "problematic_beta","problematic_beta_reason", "problematic_MAF_match", "hm_readed","snpcount",
                                "z_score_converted", "problematic_SE", "problematic_CHR","snpcount_in_raw_file","columns_dropped","columns_before","columns_after","ld_pred_failed")

  dfx[1,]= c(FALSE,FALSE, problematic_beta,beta_reason,FALSE,hm_readed ,dim(base)[1], FALSE,FALSE,problematic_CHR,snpcount_before,columns_dropped,paste0(columns_before,collapse = ","),paste0(columns_after,collapse = ","),FALSE)

  write.table(dfx,paste0(output_other_files,".metadata.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

  write.table(summary(base[1:5,]),paste0(output_other_files,".summary.qced.tsv"), quote = F, sep = "\t", col.names = T, row.names = F)

  fwrite(base, file = output, quote = F, row.names = F, sep = "\t", compress = "gzip")

  writeLines(as.character(warnings()),paste0(output_other_files,".warnings"))
}
