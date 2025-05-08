rm(list=ls())

###################
# Import packages #
###################
library(tidyverse)
library(data.table)
library(optparse)
library(dplyr)

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

######################################
# Importing data and transforming it #
######################################

gwas_catalogue_file <- as_tibble(fread(input))

if(sum(is.element(c("other_allele","effect_allele", "effect_weight", "hm_chr", "hm_pos"),colnames(gwas_catalogue_file))) == 5) {
  
  #### Keep harmonized data only ####
  
  base <- gwas_catalogue_file %>% dplyr::select(starts_with("hm_"), "effect_allele","other_allele", "effect_weight")
  
  #### Remove SNPs with no beta or OR - these cannot be used by PRSice ####
  
  base <- dplyr::filter(base, is.na(hm_chr) == FALSE)
  base <- dplyr::filter(base, is.na(hm_pos) == FALSE)
  base <- dplyr::filter(base, is.na(effect_allele) == FALSE)
  base <- dplyr::filter(base, is.na(other_allele) == FALSE)
  base <- dplyr::filter(base, is.na(effect_weight) == FALSE)
  base <- subset(base, nchar(as.character(effect_allele)) == 1)
  base <- subset(base, nchar(as.character(other_allele)) == 1)
  base <- dplyr::filter(base, !(effect_allele == "A" & other_allele == "T"))
  base <- dplyr::filter(base, !(effect_allele == "T" & other_allele == "A"))
  base <- dplyr::filter(base, !(effect_allele == "G" & other_allele == "C"))
  base <- dplyr::filter(base, !(effect_allele == "C" & other_allele == "G"))
  #base <- base %>% group_by(hm_rsid) %>% slice(which.min(p_value))
  
  #### For SNPs with at least a beta or OR, alternatively use the beta or OR to calculate the other ####
  
  base$VARID <- str_c(base$hm_chr, ":", base$hm_pos)
  base$p_value <- rep(1,nrow(base))
  
  #### Remove duplicate SNPs - these cannot be used by PRSice (an error will be thrown) ####
  
  base <- distinct(base, VARID, .keep_all = TRUE)
  
  #### Calculate MAF using effect allele frequency ####
  
  #if ("effect_allele_frequency" %in% colnames(base)) {
  #  base$MAF <- mapply(calculate_maf, base$hm_effect_allele_frequency, base$effect_allele_frequency)
  #} else {
  #  base$MAF <- mapply(calculate_maf, base$hm_effect_allele_frequency, NA)
  #}
  
  #### Change column names to match SAIGE output ####
  base$hm_source = NULL
  base$hm_inferOtherAllele = NULL
  base$hm_rsID = NULL
  colnames(base)[which(names(base) == "hm_chr")] <- "CHR"
  colnames(base)[which(names(base) == "hm_pos")] <- "POS"
  colnames(base)[which(names(base) == "effect_weight")] <- "BETA"
  colnames(base)[which(names(base) == "other_allele")] <- "Allele2"
  colnames(base)[which(names(base) == "effect_allele")] <- "Allele1"
  colnames(base)[which(names(base) == "p_value")] <- "p.value"
  
  #### Save data ####
  
  write.table(base, output, quote = F, row.names =F, sep = " ")
} else {
  header = c('CHR', 'POS', 'BETA', 'Allele2', 'Allele1', 'p.value', 'VARID')
  write.table(t(data.frame(header)), output, quote = F, row.names =F, sep = " ", col.names = F)
}
