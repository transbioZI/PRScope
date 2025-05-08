rm(list=ls())
gc()

library(tidyverse)
library(data.table)
library(optparse)
library(dplyr)
library(lattice)

args = commandArgs(trailingOnly=TRUE)

their_imputation = paste0("/data/Psycourse/genotype/imputed/vcf_maf0.01_rsID_nodup_snps/psycourse_maf0.01_rsID_nodup_hg19_gsa_id.dose",".bim")
their_maf = "/data/Psycourse/genotype/imputed/vcf_maf0.01_rsID_nodup_snps/psycourse_maf0.01_rsID_nodup_hg19_gsa_id.dose.frq"
my_imputation = paste0("/data/Psycourse/genotype/not_imputed/qc/imputation/cobg_dir_genome_wide/psy_psy1_eur_ek-qc1.hg19.ch.fl.bgn",".bim")
my_maf = "/data/Psycourse/genotype/not_imputed/qc/imputation/cobg_dir_genome_wide/psy_psy1_eur_ek-qc1.hg19.ch.fl.bgn.frq"
infos = readRDS("/data/Psycourse/genotype/not_imputed/qc/imputation/infos")
output = "/data/Psycourse/genotype/comparison/maf.pdf"

imputed_by_them = as_tibble(read.table(their_imputation))
imputed_by_me = as_tibble(read.table(my_imputation))

imputed_by_them$converted_variant_id = paste0(imputed_by_them$V1,":",imputed_by_them$V4,"_",imputed_by_them$V6,"_",imputed_by_them$V5)

inter = intersect( imputed_by_them$converted_variant_id, imputed_by_me$V2)

imputed_by_them_int = imputed_by_them[imputed_by_them$converted_variant_id %in% inter,]

imputed_by_me_int = imputed_by_me[imputed_by_me$V2 %in% inter,]

stopifnot(identical(imputed_by_them_int$converted_variant_id,imputed_by_me_int$V2))

imputed_by_them_maf = read.table(their_maf, header = T)
imputed_by_me_maf = read.table(my_maf, header = T)

imputed_by_them_maf = imputed_by_them_maf[match(imputed_by_them_int$V2,imputed_by_them_maf$SNP),]
imputed_by_me_maf = imputed_by_me_maf[match(imputed_by_me_int$V2,imputed_by_me_maf$SNP),]

imputed_by_them_maf_filter = imputed_by_them_maf[imputed_by_them_maf$MAF > 0.1,]

imputed_by_me_maf_filter = imputed_by_me_maf[imputed_by_me_maf$MAF > 0.1,]

x = imputed_by_them_int[match(imputed_by_them_maf_filter$SNP,imputed_by_them_int$V2),]
imputed_by_them_maf_filter$SNP_converted = x$converted_variant_id
inter_maf = intersect(imputed_by_me_maf_filter$SNP, imputed_by_them_maf_filter$SNP_converted)

infos_filtered = infos[infos$info > 0.8,]

inter_maf_filtered = intersect(inter_maf,infos_filtered$SNP)

imputed_by_them_maf_inter = imputed_by_them_maf_filter[match(inter_maf_filtered,imputed_by_them_maf_filter$SNP_converted) ,]
imputed_by_me_maf_inter = imputed_by_me_maf_filter[match(inter_maf_filtered,imputed_by_me_maf_filter$SNP),]
imputed_by_them_maf_inter_sorted = imputed_by_them_maf_inter[order(imputed_by_them_maf_inter$SNP_converted),]
imputed_by_me_maf_inter_sorted = imputed_by_me_maf_inter[order(imputed_by_me_maf_inter$SNP),]
dim(imputed_by_me_maf_inter_sorted)
stopifnot(identical(imputed_by_me_maf_inter_sorted$SNP, imputed_by_them_maf_inter_sorted$SNP_converted))
pdf(output)
plot(imputed_by_me_maf_inter_sorted$MAF, imputed_by_them_maf_inter_sorted$MAF)
dev.off()

