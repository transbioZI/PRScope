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
sample_size <- as.numeric(args[5])
number_of_snps <- as.numeric(args[6])
publication_date <- as.character(args[7])
population <- as.character(args[8])

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
for (x in efo_study_ids) {
    if(sum(efo_study_matching[x,]@ancestries$number_of_individuals, na.rm = T) >= sample_size) {
      if(efo_study_matching[x,]@publications$publication_date >= publication_date) {
        if(is.na(efo_study_matching[x,]@studies$snp_count) | efo_study_matching[x,]@studies$snp_count > number_of_snps) {
          if(x %in% harm) {
            res = c(x,res)
          } else {
            no_harm = c(x, no_harm)
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

write.table(res_1, file=paste0(efoStudies,".txt"), sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
write.table(no_harm, file=paste0(efoStudies,"_not_harmonized.txt"), sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
write.table(child_efos, file=paste0(efoStudies,"_all_child_efo_ids.txt"), sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
write.table(other_ancestry, file=paste0(efoStudies,"_other_ancestry_gwas.txt"), sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
