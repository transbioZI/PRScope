rm(list=ls())
gc()

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))

args = commandArgs(trailingOnly=TRUE)

input = args[1]
output = paste0(input,".tsv")
studies = unique(readLines(input))

complete_results = data.frame(matrix(ncol = 10, nrow = length(studies)))
colnames(complete_results) <- c('study_id', 'efo_trait', 'reported_trait', "initial_sample", "replication_sample","snp_count",
                                "publication_date","pubmed_id","title","publication")
complete_results$study_id = studies
complete_results$efo_trait = sapply(studies, function(x) {
  trait = get_traits(study_id = x)@traits$trait
  paste(trait, collapse = '/')

})
complete_results$reported_trait = sapply(studies, function(x) {
  get_studies(study_id = x)@studies$reported_trait
})
complete_results$initial_sample = sapply(studies, function(x) {
  get_studies(study_id = x)@studies$initial_sample_size
})
complete_results$replication_sample = sapply(studies, function(x) {
  get_studies(study_id = x)@studies$replication_sample_size
})
complete_results$snp_count = sapply(studies, function(x) {
  get_studies(study_id = x)@studies$snp_count
})
complete_results$publication_date = sapply(studies, function(x) {
  as.character(get_studies(study_id = x)@publications$publication_date)
})
complete_results$pubmed_id = sapply(studies, function(x) {
  get_studies(study_id = x)@publications$pubmed_id
})
complete_results$title = sapply(studies, function(x) {
  get_studies(study_id = x)@publications$title
})
complete_results$publication = sapply(studies, function(x) {
  get_studies(study_id = x)@publications$publication
})

write.table(apply(complete_results,2,as.character),output, quote = F,col.names = T, row.names = F, sep = "\t")
