rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))

collectChildEfos <- function(efoIds) {
  efos <- get_child_efo(efoIds, warnings = F)
  child_efos <- unlist(efos, use.names = FALSE)
  all_child_efos <- get_child_efo(unique(child_efos), warnings = F)
  #filtered_child_efos <- sapply(names(all_child_efos), function(x) length(all_child_efos[[x]]) == 0)
  return(unique(c(names(all_child_efos), efoIds)))
}

args = commandArgs(trailingOnly=TRUE)
harmonised_path <- args[1]
efoIdsPath <- args[2]
efoStudies <- paste0(args[3],"/",args[4])
number_of_controls <- as.numeric(args[5])
number_of_snps <- as.numeric(args[6])
population <- as.character(args[7])
pubmed_ids <- as.numeric(gsub(" ", "", unlist(strsplit(args[8],",")) , fixed = TRUE))
number_of_cases <- as.numeric(args[9])

efo <- readLines(efoIdsPath)
harmonised = readLines(harmonised_path)
harm = c()
for(h in harmonised) {
  harm = c(str_split(h, pattern = "/")[[1]][3],harm)
}

child_efos <- collectChildEfos(efo)

efo_study_matching <- get_studies(efo_id=child_efos, warnings = F)
efo_study_ids <-unique(efo_study_matching@publications$study_id)
res = c()
no_harm = c()
less_cases = c()
less_control = c()
less_cases_non_harm = c()
less_control_non_harm = c()
for (x in efo_study_ids) {
    if(sum(efo_study_matching[x,]@ancestries$number_of_individuals, na.rm = T) >= number_of_cases + number_of_controls) {

        if((is.na(efo_study_matching[x,]@studies$snp_count) | efo_study_matching[x,]@studies$snp_count > number_of_snps) & !(as.numeric(efo_study_matching[x,]@publications$pubmed_id) %in% pubmed_ids)) {
          if(x %in% harm) {
            st = efo_study_matching[x,]
            ts = st@studies$initial_sample_size
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
              }})), na.rm = T)

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

            if(is.na(total_cases)) {
              total_cases = 0
            }
            if(is.na(total_controls)) {
              total_controls = 0
            }

            if(total_cases == 0 | total_cases > number_of_cases) {
              if(total_controls == 0 | total_controls > number_of_controls) {
                res = c(x,res)
              } else {
                less_control = c(x, less_control)
              }
            } else {
              less_cases = c(x, less_cases)
            }
          } else {
            st = efo_study_matching[x,]
            ts = st@studies$initial_sample_size
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
              }})), na.rm = T)

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

            if(is.na(total_cases)) {
              total_cases = 0
            }
            if(is.na(total_controls)) {
              total_controls = 0
            }

            if(total_cases == 0 | total_cases > number_of_cases) {
              if(total_controls == 0 | total_controls > number_of_controls) {
                no_harm = c(x,no_harm)
              } else {
                less_control_non_harm = c(x, less_control_non_harm)
              }
            } else {
              less_cases_non_harm = c(x, less_cases_non_harm)
            }

          }
        }
    }
}

res_1 = c()
other_ancestry = c()
for (x in res) {
    pops = unique(efo_study_matching[x,]@ancestral_groups$ancestral_group)
    pops = pops[!is.na(pops)]
    if(length(pops) == 1 & population %in% pops ) {
      res_1 = c(res_1,x)
    } else {
      other_ancestry = c(other_ancestry,x)
    }
}

non_harm_1 = c()
other_ancestry = c()
for (x in no_harm) {
    pops = unique(efo_study_matching[x,]@ancestral_groups$ancestral_group)
    pops = pops[!is.na(pops)]
    if(length(pops) == 1 & population %in% pops ) {
      non_harm_1 = c(non_harm_1,x)
    } else {
      other_ancestry = c(other_ancestry,x)
    }
}

less_control_1 = c()
for (x in less_control) {
    pops = unique(efo_study_matching[x,]@ancestral_groups$ancestral_group)
    pops = pops[!is.na(pops)]
    if(length(pops) == 1 & population %in% pops ) {
      less_control_1 = c(less_control_1,x)
    }
}

less_cases_1 = c()
for (x in less_cases) {
    pops = unique(efo_study_matching[x,]@ancestral_groups$ancestral_group)
    pops = pops[!is.na(pops)]
    if(length(pops) == 1 & population %in% pops ) {
      less_cases_1 = c(less_cases_1,x)
    }
}


less_control_non_harm_1 = c()
for (x in less_control_non_harm) {
    pops = unique(efo_study_matching[x,]@ancestral_groups$ancestral_group)
    pops = pops[!is.na(pops)]
    if(length(pops) == 1 & population %in% pops ) {
      less_control_non_harm_1 = c(less_control_non_harm_1,x)
    }
}

less_cases_non_harm_1 = c()
for (x in less_cases_non_harm) {
    pops = unique(efo_study_matching[x,]@ancestral_groups$ancestral_group)
    pops = pops[!is.na(pops)]
    if(length(pops) == 1 & population %in% pops ) {
      less_cases_non_harm_1 = c(less_cases_non_harm_1,x)
    }
}

writeLines(as.character(less_cases_non_harm_1), paste0(efoStudies,"_less_cases_no_harm.txt"))
writeLines(as.character(less_control_non_harm_1), paste0(efoStudies,"_less_control_no_harm.txt"))
writeLines(as.character(less_control_1), paste0(efoStudies,"_less_control.txt"))
writeLines(as.character(less_cases_1), paste0(efoStudies,"_less_cases.txt"))
writeLines(as.character(non_harm_1), paste0(efoStudies,"_not_harmonized.txt"))
writeLines(as.character(res_1), paste0(efoStudies,".txt"))
writeLines(as.character(child_efos), paste0(efoStudies,"_all_child_efo_ids.txt"))
writeLines(as.character(other_ancestry), paste0(efoStudies,"_other_ancestry_gwas.txt"))
