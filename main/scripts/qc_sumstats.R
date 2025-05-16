#!/usr/bin/env Rscript
rm(list=ls())
gc()

###################
# Import packages #
###################

library(tidyverse)
library(data.table)
library(optparse)
library(dplyr)
suppressMessages(library(gwasrapidd))

calculate_maf <- function(heaf,eaf) {
  if(is.na(eaf) & is.na(heaf)) {
    return(1) # or NA
  }
  if(is.na(heaf)) {
    return(0)
  }
  return(heaf)
}

#####################
# Parsing arguments #
#####################

args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
gwas_id <- args[3]
number_of_individuals <- args[4]
maf_file <- args[5]

gwas_catalogue_file <- as_tibble(fread(input))

#### Keep harmonized data only ####

base <- gwas_catalogue_file

#### Remove SNPs with no beta or OR - these cannot be used by PRSice ####
if("hm_beta" %in% colnames(base)) {
  colnames(base)[which(names(base) == "hm_rsid")] <- "SNPID"
  colnames(base)[which(names(base) == "hm_chrom")] <- "CHR"
  colnames(base)[which(names(base) == "hm_pos")] <- "POS"
  colnames(base)[which(names(base) == "hm_other_allele")] <- "Allele2"
  colnames(base)[which(names(base) == "hm_effect_allele")] <- "Allele1"
  colnames(base)[which(names(base) == "hm_beta")] <- "BETA"
  colnames(base)[which(names(base) == "hm_odds_ratio")] <- "OR"

} else {
  if(!("beta" %in% colnames(base))) {
    base$beta = NA
  }
  if(!("odds_ratio" %in% colnames(base))) {
    base$odds_ratio = NA
  }
  
  colnames(base)[which(names(base) == "rsid")] <- "SNPID"
  colnames(base)[which(names(base) == "chromosome")] <- "CHR"
  colnames(base)[which(names(base) == "base_pair_location")] <- "POS"
  colnames(base)[which(names(base) == "other_allele")] <- "Allele2"
  colnames(base)[which(names(base) == "effect_allele")] <- "Allele1"
  colnames(base)[which(names(base) == "beta")] <- "BETA"
  colnames(base)[which(names(base) == "odds_ratio")] <- "OR"
}

base <- dplyr::filter(base, !(is.na(BETA) == TRUE & is.na(OR) == TRUE))
base <- subset(base, nchar(as.character(Allele1)) == 1)
base <- subset(base, nchar(as.character(Allele2)) == 1)
if("info" %in% colnames(base)) {
  base <- dplyr::filter(base, is.na(info) | info > 0.8)
}
base <- dplyr::filter(base, !(Allele1 == "A" & Allele2 == "T"))
base <- dplyr::filter(base, !(Allele1 == "T" & Allele2 == "A"))
base <- dplyr::filter(base, !(Allele1 == "G" & Allele2 == "C"))
base <- dplyr::filter(base, !(Allele1 == "C" & Allele2 == "G"))
#### For SNPs with at least a beta or OR, alternatively use the beta or OR to calculate the other ####

base$VARID <- str_c(base$CHR, ":", base$POS)

base$BETA <- suppressWarnings(as.numeric(base$BETA))
base$OR <- suppressWarnings(as.numeric(base$OR))

base <- base %>%
  mutate(BETA = if_else(is.na(BETA), log(OR), BETA),
         OR = if_else(is.na(OR), exp(BETA), OR))
#### For SNPs with no p-value, replace the NA by a 1 ####

base <- base %>% mutate(p_value = if_else(is.na(p_value), 1, as.numeric(p_value)))
#### Remove duplicate SNPs - these cannot be used by PRSice (an error will be thrown) ####

base <- distinct(base, SNPID, .keep_all = TRUE)
base <- distinct(base, VARID, .keep_all = TRUE)
if("n" %in% colnames(base)) {
  colnames(base)[which(names(base) == "n")] <- "N"
}

if("N" %in% colnames(base)) {
  base <- dplyr::filter(base, !(is.na(N) == TRUE ))
}

if(!is.na(number_of_individuals) & number_of_individuals == "NA") {
  number_of_individuals = NA
}

if(!("N" %in% colnames(base))) {
  if(is.na(number_of_individuals)) {
    s = get_studies(study_id = strsplit(gwas_id,"_",fixed = T)[[1]][1])
    base$N = sum(s@ancestries$number_of_individuals, na.rm = T)
  }
  else {
    base$N = as.numeric(number_of_individuals)
  }
}
#### Calculate MAF using effect allele frequency ####

if("beta" %in% trimws(colnames(base))) {
  base$beta = NULL
}

maf = as.data.frame(fread(maf_file))

current_s = as.data.frame(base)
matches = match(current_s$VARID,maf$SNP)
matches_na = which(is.na(matches))
matches_not_na = which(!is.na(matches))
current_s$BASEMAF = NA
current_s[matches_not_na,]$BASEMAF = maf$MAF[matches[!is.na(matches)]]
current_s[matches_na,]$BASEMAF = 0.000000001
#### Save data ####

write.table(current_s[,c("VARID","SNPID","CHR","POS","Allele1","Allele2","p_value","BASEMAF","BETA","OR","standard_error","N")], output, quote = F, row.names =F, sep = " ")
