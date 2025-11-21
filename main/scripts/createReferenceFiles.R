library(data.table)
library(dplyr)

setwd("/data/projects/on_going/PRScope/input/reference")

hg37 = fread("/data/references/1000Genome/reference/EUR_RS_HG37/eur_hg37.phase3.frq", data.table = F)
hg38 = fread("/data/references/1000Genome/reference/EUR_NEW_RS/eur_hg38.phase3.frq", data.table = F)
hg38_bim = fread("/data/references/1000Genome/reference/EUR_NEW_RS/eur_hg38.phase3.bim",data.table = F)
hg37_bim = fread("/data/references/1000Genome/reference/EUR_RS_HG37/eur_hg37.phase3.bim",data.table = F)

colnames(hg38_bim)[colnames(hg38_bim) == "V2"] ="SNP"
colnames(hg37_bim)[colnames(hg37_bim) == "V2"] ="SNP"

hg37 = hg37[!duplicated(hg37$SNP),]
hg37_bim = hg37_bim[!duplicated(hg37_bim$SNP),]

stopifnot(identical(hg37$SNP,hg37_bim$SNP))
stopifnot(identical(hg38$SNP,hg38_bim$SNP))

hg38 = inner_join(hg38,hg38_bim,by ="SNP", suffix = c("_frq", "_bim"))
hg37 = inner_join(hg37,hg37_bim,by ="SNP", suffix = c("_frq", "_bim"))
  
merged = inner_join(hg38,hg37,by =c("SNP"), suffix = c("_hg38", "_hg37"))
merged = merged[,c("CHR_hg38","SNP","A1_hg38","A2_hg38","A1_hg37","A2_hg37","MAF_hg38","MAF_hg37")]
colnames(merged)[colnames(merged) =="CHR_hg38"] = "CHR"
write.table(merged,"eur_hg38_hg37.frq",row.names = F, col.names = T, quote = F, sep = "\t")
