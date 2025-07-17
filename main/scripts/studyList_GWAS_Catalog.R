rm(list=ls())
gc()

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))

args = commandArgs(trailingOnly=TRUE)

input = args[1]
output = paste0(input,".tsv")
studies = unique(readLines(input))

if(length(studies) != 0) {
studies_read = get_studies(study_id = studies)

complete_results = data.frame(matrix(ncol = 14, nrow = length(studies)))

colnames(complete_results) <- c('study_id', 'efo_trait', 'reported_trait', "initial_sample", "replication_sample","snp_count",
                                "publication_date","pubmed_id","title","publication", "sample_size", "number_of_cases","number_of_controls", "countries_of_recruitment")

complete_results$study_id = studies

complete_results$efo_trait = sapply(studies, function(x) {
  trait = get_traits(study_id = x)@traits$trait
  paste(trait, collapse = '/')

})

complete_results$reported_trait = studies_read@studies$reported_trait
complete_results$initial_sample = studies_read@studies$initial_sample_size
complete_results$replication_sample = studies_read@studies$replication_sample_size
complete_results$snp_count = studies_read@studies$snp_count
complete_results$publication_date = as.character(studies_read@publications$publication_date)
complete_results$pubmed_id = studies_read@publications$pubmed_id
complete_results$title = studies_read@publications$title
complete_results$publication = studies_read@publications$publication
complete_results$sample_size = sapply(c(1:length(studies)), function(x) {
  sum(studies_read[x]@ancestries$number_of_individuals, na.rm = T)
})

complete_results$countries_of_recruitment = sapply(c(1:length(studies)), function(x) {
  paste(studies_read[x]@countries_of_recruitment$country_name, collapse = ",")
})

initial_sample_size_info = studies_read@studies$initial_sample_size

complete_results$number_of_cases = sapply(initial_sample_size_info, function(sample_size_info) {
  ts = sample_size_info
  rs = as.character(unlist(strsplit(ts,", ")))
  rs = rs[str_detect(rs, regex("cases"))]
  total_cases = sum(as.numeric(sapply(rs,function(x) {
    x = gsub(",","",x)
    a = unlist(str_extract_all(x, "\\d+"))
    if(length(a) > 1) {
      a[1]
    } else if(length(a) == 0) {
      0
    } else {
      a
    }
  })), na.rm = T)
  total_cases
})

complete_results$number_of_controls = sapply(initial_sample_size_info, function(sample_size_info) {
  ts = sample_size_info
  rs = as.character(unlist(strsplit(ts,", ")))
  rs = rs[str_detect(rs, regex("controls"))]
  total_controls = sum(as.numeric(sapply(rs,function(x) {
    x = gsub(",","",x)
    a = unlist(str_extract_all(x, "\\d+"))
    if(length(a) > 1) {
      a[1]
    } else if(length(a) == 0) {
      0
    } else {
      a
    }
  })), na.rm = T)
  total_controls
})

write.table(apply(complete_results,2,as.character),output, quote = F,col.names = T, row.names = F, sep = "\t")
} else {

complete_results = data.frame(matrix(ncol = 14, nrow = 0))

colnames(complete_results) <- c('study_id', 'efo_trait', 'reported_trait', "initial_sample", "replication_sample","snp_count",
                                "publication_date","pubmed_id","title","publication", "sample_size", "number_of_cases","number_of_controls", "countries_of_recruitment")
write.table(apply(complete_results,2,as.character),output, quote = F,col.names = T, row.names = F, sep = "\t")

}
