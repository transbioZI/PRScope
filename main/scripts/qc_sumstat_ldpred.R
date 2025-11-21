rm(list=ls())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(gwasrapidd))
suppressMessages(library(bigsnpr))

args = commandArgs(trailingOnly=TRUE)

input <- args[1]
ldpred_output = args[2]
ldpred_ref = args[3]
genome_build = args[4]

map_ldref = readRDS(ldpred_ref)

if(genome_build=="hg38") {
  map_ldref$pos_hg19 = map_ldref$pos
  map_ldref$pos = map_ldref$pos_hg38
}

base <- fread(input, showProgress = FALSE, data.table = F)

if(dim(base)[1] > 1) {

    colnames(base) = tolower(colnames(base))
    colnames(base)[colnames(base) == "a2"] = "a0"
    colnames(base)[colnames(base) == "bp"] = "pos"

    base <- snp_match(base, map_ldref)

    sd_ldref <- with(base, sqrt(2*maf*(1 - maf)))
    sd_ss <- with(base, 2 / sqrt(neff * se^2 + beta^2))

    is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05
    base <- base[!is_bad, ]
    base <- base[ , c("varid","snp","chr","pos","a1","a0","p","maf","beta","or","se","n","neff")]
    colnames(base)[colnames(base) == "a0"] = "a2"
    colnames(base)[colnames(base) == "pos"] = "bp"
    colnames(base) = toupper(colnames(base))

    fwrite(base, file = ldpred_output, quote = F, row.names = F, sep = "\t", compress = "gzip")
} else {
    fwrite(base, file = ldpred_output, quote = F, row.names = F, sep = "\t", compress = "gzip")
}

