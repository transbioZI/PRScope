rm(list=ls())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(gwasrapidd))

args = commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]
maf_file <- args[3]

base <- fread(input, showProgress = FALSE, data.table = F)
hm_readed = FALSE
if("hm_beta" %in% colnames(base)) {
  
  if(!("hm_effect_allele_frequency" %in% colnames(base))) {
    base$hm_effect_allele_frequency = NA
  }
  
  colnames(base)[which(names(base) == "hm_rsid")] <- "SNP"
  colnames(base)[which(names(base) == "hm_chrom")] <- "CHR"
  colnames(base)[which(names(base) == "hm_pos")] <- "BP"
  colnames(base)[which(names(base) == "hm_other_allele")] <- "A2"
  colnames(base)[which(names(base) == "hm_effect_allele")] <- "A1"
  colnames(base)[which(names(base) == "hm_beta")] <- "BETA"
  colnames(base)[which(names(base) == "hm_odds_ratio")] <- "OR"
  colnames(base)[which(names(base) == "hm_effect_allele_frequency")] <- "EAF"
  hm_readed = TRUE
} else {
  if(!("beta" %in% colnames(base))) {
    base$beta = NA
  }
  if(!("odds_ratio" %in% colnames(base))) {
    base$odds_ratio = NA
  }
  if(!("effect_allele_frequency" %in% colnames(base))) {
    base$effect_allele_frequency = NA
  }
  colnames(base)[which(names(base) == "rsid")] <- "SNP"
  colnames(base)[which(names(base) == "chromosome")] <- "CHR"
  colnames(base)[which(names(base) == "base_pair_location")] <- "BP"
  colnames(base)[which(names(base) == "other_allele")] <- "A2"
  colnames(base)[which(names(base) == "effect_allele")] <- "A1"
  colnames(base)[which(names(base) == "beta")] <- "BETA"
  colnames(base)[which(names(base) == "odds_ratio")] <- "OR"
  colnames(base)[which(names(base) == "effect_allele_frequency")] <- "EAF"
}

colnames(base)[which(names(base) == "p_value")] <- "P"
colnames(base)[which(names(base) == "standard_error")] <- "SE"

z_score_studies = FALSE
if((sum(is.na(base$BETA)) == length(base$BETA)) & (sum(is.na(base$OR)) == length(base$OR))) {
  if(("z_score" %in% colnames(base))) {
    if((sum(is.na(base$SE)) != length(base$SE)) & (sum(is.na(base$z_score)) != length(base$z_score))) {
      base = base[!is.na(base$SE),]
      base = base[!is.na(base$z_score),]
      base$BETA = base$SE*base$z_score
      z_score_studies = TRUE
    }
  }
}

base = base %>% select(SNP,CHR,BP,A1,A2,BETA,OR,EAF,P,SE)

base$VARID <- str_c(base$CHR, ":", base$BP)
base$MAF = NA

#### Remove SNPs with no beta or OR - these cannot be used by PRSice ####
base <- dplyr::filter(base, !(is.na(BETA) == TRUE & is.na(OR) == TRUE))

problematic_beta = FALSE
if(dim(base)[1] == 0 ) {
    problematic_beta = TRUE
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

#### For SNPs with no p-value, replace the NA by a 1 ####
base <- base %>% mutate(P = if_else(is.na(P), 1, as.numeric(P)))
#### Remove duplicate SNPs - these cannot be used by PRSice (an error will be thrown) ####
base <- distinct(base, SNP, .keep_all = TRUE)
base <- distinct(base, VARID, .keep_all = TRUE)

if(dim(base)[1] != 0 ) {
  if(abs(median(base$BETA)) > 0.5 ) {
    problematic_beta = TRUE
  }
}

problematic_p_value = FALSE
if(dim(base)[1] != 0 ) {
  if((sum(base$P > 1) > 0) | (sum(base$P < 0) > 0)) {
    problematic_p_value = TRUE
  }
}

if(dim(base)[1] != 0 ) {
  maf = as.data.frame(fread(maf_file, select = c("SNP","MAF"), showProgress = FALSE))
  matches = match(base$VARID,maf$SNP)
  matches_na = which(is.na(matches))
  matches_not_na = which(!is.na(matches))
  base[matches_not_na,]$MAF = maf$MAF[matches[!is.na(matches)]]
  base[matches_na,]$MAF = 0.000000001
}

writeLines(as.character(problematic_p_value), paste0(output,"problematic_p_value"))
writeLines(as.character(problematic_beta), paste0(output,".problematic_beta"))
writeLines(as.character(hm_readed), paste0(output,".hm_readed"))
writeLines(as.character(dim(base)[1]), paste0(output,".snpcount"))
writeLines(as.character(z_score_studies), paste0(output,".z_score_converted"))
fwrite(base[,c("VARID","SNP","CHR","BP","A1","A2","P","MAF","BETA","OR","SE","EAF")], file = output, quote = F, row.names = F, sep = " ", compress = "gzip")
