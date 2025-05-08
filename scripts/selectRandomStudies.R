rm(list=ls())
library(stringr)
library(gwasrapidd)

tab = read.csv("/zi/home/ersoy.kocak/Downloads/gwas-efo-trait-mappings.tsv", sep = "\t")
existing = unique(readLines("/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/studies_to_calculate_new.txt.childEfoIds"))

efo_ids = lapply(tab$EFO.URI, function(x) {
  res = str_split(x,pattern = "/")
  res[[1]][length(res[[1]])]
})

efo_ids = unique(unlist(efo_ids))

efo_ids = efo_ids[!(efo_ids %in%existing)]
study = get_studies(efo_id = efo_ids, warnings = F)

studies = lapply(efo_ids, function(x) {
  study = get_studies(efo_id = x, warnings = F)
  if(isS4(study) == T) {
    study = study@studies
    if(dim(study)[1] != 0) {
      study$study_id
    }
  }
})

studies_1 = unique(unlist(studies))

harm = unique(readLines("/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/harmonised_list.txt"))

harm = unlist(lapply(harm, function(x) {
  str_split(x,"/")[[1]][3]
}))

studies_1 = studies_1[studies_1 %in% harm]

res = c()
for (x in studies_1) {
  st = get_studies(study_id = x)
  if(dim(st@ancestries)[1] == 1) {
    if(!is.na(st@ancestries$number_of_individuals) & st@ancestries$number_of_individuals >= 10000) {
      if(st@publications$publication_date >= "2018-01-01") {
        if(is.na(st@studies$snp_count) | st@studies$snp_count > 100000) {
          res = c(x,res)
        }
      }
    }
  }
}

res_1 = c()
for (x in res) {
  st = get_studies(study_id = x)
  if(length(st@ancestral_groups$ancestral_group) == 1) {
    if(!is.na(st@ancestral_groups$ancestral_group) & st@ancestral_groups$ancestral_group == "European") {
      res_1 = c(res_1,x)
    } else {
      print(x)
    }
  }
}

writeLines(res_1,"/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/randomStudies_all.txt")

