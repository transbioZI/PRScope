rm(list=ls())

# imports

if(!require(gwasrapidd)){
    install.packages("gwasrapidd")
}
if(!require(dplyr)){
    install.packages("dplyr")
}
if(!require(tidyverse)){
    install.packages("tidyverse")
}

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

efoIdsPath <- args[1]
efoStudies <- args[2]
harmonised_path <- args[3]

efo <- readLines(efoIdsPath)
harmonised = readLines(harmonised_path)
harm = c()
for(h in harmonised) {
  harm = c(str_split(h, pattern = "/")[[1]][3],harm)
}
message("collectChildEfos")
child_efos <- collectChildEfos(efo)
message("get_studies")
efo_study_matching <- get_studies(efo_id=child_efos, warnings = F)
efo_study_ids <-unique(efo_study_matching@publications$study_id)
message("filter")
res = c()
no_harm = c()
for (x in efo_study_ids) {
    if(sum(efo_study_matching[x,]@ancestries$number_of_individuals, na.rm = T) >= 10000) {
      if(efo_study_matching[x,]@publications$publication_date >= "2018-01-01") {
        if(is.na(efo_study_matching[x,]@studies$snp_count) | efo_study_matching[x,]@studies$snp_count > 250000) {
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
    if(length(pops) == 1 & "European" %in% pops ) {
      res_1 = c(res_1,x)
    } else {
      other_ancestry = c(other_ancestry,x)
    }
}

write.table(res_1, file=efoStudies, sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
write.table(no_harm, file=paste0(efoStudies,".no_harm"), sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
write.table(child_efos, file=paste0(efoStudies,".childEfoIds"), sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
write.table(other_ancestry, file=paste0(efoStudies,".otherAncestry"), sep="\t", row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8", col.names = FALSE)
