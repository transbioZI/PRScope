rm(list=ls())
gc()
library(gwasrapidd)
studies = readLines("/zi/flstorage/HITKIP/share/Ersoy/sumstats/pipeline_studies_final_list.txt")
studies = studies[grepl("GCST",studies)]

efo_study_matching <- get_studies(study_id =studies, warnings = F)
efo_study_ids <-unique(efo_study_matching@publications$study_id)

sample_sizes = c()
for (x in efo_study_ids) {
  sample_sizes = c(sample_sizes,sum(efo_study_matching[x,]@ancestries$number_of_individuals, na.rm = T))
  
}

df = data.frame(id = studies, sizes = sample_sizes)

studies_2 = readLines("/zi/flstorage/HITKIP/share/Ersoy/sumstats/pipeline_studies_final_list.txt")
studies_2 = studies_2[!grepl("GCST",studies_2)]

pgc = fread("/data/sumstats/PGC_GWASs/pgc_list.csv",data.table = T)
pgc = pgc[,c("New_Tag","N")]
pgc = pgc[match(studies_2,pgc$New_Tag),]

df_1 = data.frame(id = pgc$New_Tag, sizes = pgc$N)

x = rbind(df,df_1)

write.table(x,"/zi/flstorage/HITKIP/share/Ersoy/MultiPRS_SCZ_Clustering_Project/sumstats/GWAS_Sample_Sizes.tsv", col.names = F,row.names = F, sep = "\t", quote = F)


