
rm(list=ls())
gc()

library(dplyr)
library(stringr)
path_to_sumstats = "/data/sumstats/transformed"
save_ = "/data/projects/on_going/ld_pred_psycourse/less_than_250000.txt"
# Count your lines of R code
files = list.files(path = path_to_sumstats , recursive = T, full.names = T)
number_of_lines = c()
for(f in files) {
  number_of_lines = c(number_of_lines, as.numeric(str_split(system(paste0("wc -l ",f),intern = TRUE), " ")[[1]][1]))
}

files[which(number_of_lines<250000)]
less_than = sapply(files[which(number_of_lines<250000)], function(x) {
  a = str_split(x,"/")[[1]][5]
  str_split(a,"\\.")[[1]][1]
})

writeLines(less_than, save_)
