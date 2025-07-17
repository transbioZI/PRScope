rm(list=ls())
gc()

library(gwasrapidd)
library(data.table)
library(stringr)
library(doParallel)

x = fread("/zi/home/ersoy.kocak/Desktop/Projects/CoviDrug/study_lists/GWAS_Sample_Sizes_ALL.tsv")

trait.names = x$V1

results_path = "/zi/home/ersoy.kocak/Desktop/Projects/Gen_Cor/results/"

gen_cor = list.files(results_path)
id = c()
not_munged = c()
correlation_res = c()
heritibality = rep(NA,length(trait.names))
intercepts = rep(NA,length(trait.names))
se = rep(NA,length(trait.names))
pvalue = rep(NA,length(trait.names))
zvalue = rep(NA,length(trait.names))
inter_se = rep(NA,length(trait.names))

for(f in gen_cor) {
  st1 = str_split(f,"_")[[1]][1]
  st2 = str_split(f,"_")[[1]][2]
  st2 = str_split(st2,"\\.")[[1]][1]
  id = c(id, paste(sort(c(st1,st2)),collapse = "_"))
  
  current_file_read = readLines(paste0(results_path,st1,"_",st2,".log"))
  ind = grep('^Summary of Genetic Correlation Results', current_file_read)
  her = grep('^Heritability of phenotype 1', current_file_read)
  
  if(length(her) != 0) {
    heri = current_file_read[her+2]
    interi = current_file_read[her+5]
    heritibality[which(trait.names == st1)] = as.numeric(str_split(heri," ")[[1]][5])
    intercepts[which(trait.names == st1)] = as.numeric(str_split(interi," ")[[1]][2])
    str_p = str_split(heri, " ")[[1]][6]
    se[which(trait.names == st1)] = as.numeric(gsub(")", "", gsub("\\(", "", str_p)) )
    zvalue[which(trait.names == st1)] = as.numeric(str_split(heri," ")[[1]][5]) / as.numeric(gsub(")", "", gsub("\\(", "", str_p)))
    pvalue[which(trait.names == st1)] = exp(-0.717*zvalue[which(trait.names == st1)] - 0.416*zvalue[which(trait.names == st1)]^2)
    inter_se[which(trait.names == st1)] = as.numeric(gsub("\\)","",gsub("\\(","",str_split(interi," ")[[1]][3])))
  } 
  
  if(length(ind) == 0) {
    print("######1")
    print(f)
    not_munged = c(not_munged,f)
    next
  }
  
  rg = str_split(current_file_read[ind+2],"\\s+")[[1]][3]
  
  correlation_res = c(correlation_res,rg)
}

heri = "Total Observed scale h2: 0.0072 (0.0005)"
interi = "Intercept: 1.0412 (0.0073)"
heritibality[which(trait.names == "PGC045")] = as.numeric(str_split(heri," ")[[1]][5])
intercepts[which(trait.names == "PGC045")] = as.numeric(str_split(interi," ")[[1]][2])
str_p = str_split(heri, " ")[[1]][6]
se[which(trait.names == "PGC045")] = as.numeric(gsub(")", "", gsub("\\(", "", str_p)) )
zvalue[which(trait.names == "PGC045")] = as.numeric(str_split(heri," ")[[1]][5]) / as.numeric(gsub(")", "", gsub("\\(", "", str_p)))
pvalue[which(trait.names == "PGC045")] = exp(-0.717*zvalue[which(trait.names == "PGC045")] - 0.416*zvalue[which(trait.names == "PGC045")]^2)
inter_se[which(trait.names == "PGC045")] = as.numeric(gsub("\\)","",gsub("\\(","",str_split(interi," ")[[1]][3])))

df = data.frame(gwas_1 = rep(NA,length(id)), gwas_2 = rep(NA,length(id)), 
           intercept_1 = rep(NA,length(id)) , intercept_2 = rep(NA,length(id)),
           intercept_se_1 = rep(NA,length(id)),intercept_se_2 = rep(NA,length(id)),
           heritibality_1 =rep(NA,length(id)),heritibality_2 =rep(NA,length(id)),
           heritibality_se_1 =rep(NA,length(id)),heritibality_se_2 =rep(NA,length(id)),
           zvalue_1 =rep(NA,length(id)),zvalue_2 =rep(NA,length(id)),
           pvalue_1 =rep(NA,length(id)),pvalue_2 =rep(NA,length(id)), 
           correlation = rep(NA,length(id)))

for(i in id) {
  st1 = str_split(i,"_")[[1]][1]
  st2 = str_split(i,"_")[[1]][2]
  df[which(i == id),] = c(st1,st2,
                          intercepts[which(trait.names == st1)],intercepts[which(trait.names == st2)],
                          inter_se[which(trait.names == st1)],inter_se[which(trait.names == st2)],
                          heritibality[which(trait.names == st1)],heritibality[which(trait.names == st2)],
                          se[which(trait.names == st1)],se[which(trait.names == st2)],
                          zvalue[which(trait.names == st1)],zvalue[which(trait.names == st2)],
                          pvalue[which(trait.names == st1)],pvalue[which(trait.names == st2)],
                          correlation_res[which(i == id)])
}

df[which(df$correlation == "NA"),]$correlation = NA
df$correlation = as.numeric(df$correlation)
write.table(df, "/zi/home/ersoy.kocak/Desktop/Projects/CoviDrug/heritability_table_all.tsv", row.names = F, col.names = T, quote = F, sep = "\t")


