rm(list=ls())
gc()

library(tidyverse)
library(data.table)
library(optparse)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
convert <- as.logical(args[3])
plink <- args[4]
plink2 <- args[5]
resources <- args[6]
target_data = as_tibble(read.table(paste0(input,".bim")))

target_data_filtered =  subset(target_data, nchar(as.character(V5)) == 1)
target_data_filtered =  subset(target_data_filtered, nchar(as.character(V6)) == 1)
target_data_filtered =  dplyr::filter(target_data_filtered, !(V5 == "." | V6 == "."))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V5 == "I" & V6 == "D"))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V5 == "D" & V6 == "I"))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V5 == "A" & V6 == "T"))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V5 == "T" & V6 == "A"))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V5 == "G" & V6 == "C"))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V5 == "C" & V6 == "G"))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V1 == 6 & V4 > 25000000 & V4 < 35000000))
target_data_filtered =  dplyr::filter(target_data_filtered, !(V1 == 8 & V4 > 7000000 & V4 < 13000000))
target_data_filtered = distinct(target_data_filtered, V2, .keep_all = TRUE)

writeLines(target_data_filtered$V2, paste0(input,".extract1"))

fam = as_tibble(read.table(paste0(input,".fam"), header = F, colClasses = c("character","character","character","character","character","character")))
total_sample = dim(fam)[1]
selected_samples = sample(c(1:total_sample), sum(runif(n=total_sample,min = 0, max = 1) < min((1000/total_sample),1.1)))
a = fam[selected_samples,1]
b = fam[selected_samples,2]
writeLines(paste0(a$V1," ",b$V2), paste0(input,".keep_samples"))

system(paste0(plink, " --bfile ", input, " --make-bed --allow-no-sex --autosome --extract ", paste0(input,".extract1")," --out ", paste0(input,".temp1")))
if(convert == TRUE) {
system(paste0(plink2, " --bfile ", paste0(input,".temp1"), " --make-bed --allow-no-sex --set-all-var-ids @:# --out ", paste0(input,".temp2")))
} else {
system(paste0(plink2, " --bfile ", paste0(input,".temp1"), " --make-bed --out ", paste0(input,".temp2")))
}

system( paste0("rm ", input,".temp1*") )

system(paste0(plink2, " --bfile ", paste0(input,".temp2"), " --make-bed --allow-no-sex --rm-dup exclude-all list --freq --het --missing --out ", output))

system( paste0("rm ", input,".temp2*") )

target_data = as_tibble(read.table(paste0(output,".bim")))
total_snp = length(target_data$V2)
selected_snps = sample(target_data$V2, sum(runif(n=total_snp,min = 0, max = 1) < min((1000000/total_snp),1.1)))
writeLines(selected_snps, paste0(input,".extract2"))

system(paste0(plink, " --bfile ", output, " --make-bed --allow-no-sex --out ", paste0(input,".temp4")," --extract ",paste0(input,".extract2")," --keep ", paste0(input,".keep_samples")))
system(paste0(plink, " --bfile ", paste0(input,".temp4")," --hwe 0.001 --geno 0.02 --make-founders  --make-bed --out ", paste0(input,".filtered2")))
system( paste0("rm ", input,".temp4*") )
system(paste0(plink, " --bfile ", paste0(input,".filtered2"), " --indep-pairwise 200 100 0.2 --maf 0.05 --out ", paste0(input,".filtered3")))
system(paste0(plink, " --bfile ", paste0(input,".filtered2"), " --make-bed --extract ",paste0(input,".filtered3.prune.in"), " --out ", paste0(input,".filtered4")))
system( paste0("rm ", input,".filtered2*") )
system(paste0(plink, " --bfile ", paste0(input,".filtered4"), " --indep-pairwise 200 100 0.2 --out ", paste0(input,".filtered5")))
system( paste0("rm ", input,".filtered4*") )

prune_in = readLines(paste0(input,".filtered5.prune.in"))

if(length(prune_in) > 100000) {
    prune_in = sample(prune_in, 100000, replace = FALSE, prob = NULL)
    writeLines(prune_in,paste0(input,".filtered5.prune.in"))
}
writeLines(prune_in,paste0(output,".POP_STRATIFICATION.SNPS"))
system(paste0(plink, " --bfile ", output, " --extract ",paste0(input,".filtered5.prune.in"), " --pca --autosome --out ", paste0(output,".POP_STRATIFICATION")))

system( paste0("rm ", input,".filtered5*") )
pop_stratification_plink = read.table(paste0(output,".POP_STRATIFICATION.eigenvec"))
colnames(pop_stratification_plink) = c("FID","IID",paste0("PC",c(1:20)))
jpeg(paste0(output,".POP_STRATIFICATION.jpeg"))
pairs(pop_stratification_plink[,3:8])
dev.off()

system( paste0("rm ", input,".filtered*") )
system( paste0("rm ", input,".temp*") )
system( paste0("rm ", input,".keep_s*") )
system( paste0("rm ", input,".extract*") )
system( paste0("rm ", input,"*") )
