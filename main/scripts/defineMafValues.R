rm(list=ls())

# imports
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))

args = commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
default_maf <- args[3]
study_ids = readLines(input)

trmap = c("transbio001"=0.01,"transbio002"=0.01,"transbio003"=0.01,"transbio004"=0.01,
  "transbio005"=0.01,"transbio006"=0.01,"transbio007"=0.01,"transbio008"=0.1,
  "transbio009"=0.1,"transbio010"=0.1,"transbio011"=0.1,"transbio012"=0.1,
  "transbio013"=0.01,"transbio014"=0.01,"transbio015"=0.1,"transbio016"=0.05,
  "transbio017"=0.01,"transbio018"=0.01,"transbio019"=0.01,"transbio020"=0.01,"transbio022"=0.01,
  "transbio023"=0.01,"transbio024"=0.01,"transbio025"=0.01,"transbio026"=0.01,"transbio027"=0.01,
  "transbio028"=0.01,"transbio029"=0.01,"transbio030"=0.01,"transbio031"=0.01,"transbio032"=0.01,"transbio033"=0.01)

st = str_split(study_ids, "_", simplify = FALSE)

st_l = sapply(c(1:length(st)), function(x) {
  st[[x]][1]
})
print("starting")
mafs = rep(0.01,length(st_l))
#for(s in st_l) {
#  if(grepl("transbio",s, fixed = T)) {
#    mafs = c(mafs, trmap[s])
#  } else {
#    study = get_studies(study_id = s)
#
#    if(is.null(study) || is.null(study@ancestries) || is.na(sum(study@ancestries$number_of_individuals,rm.na = T))) {
#      print(s)
#      mafs = c(mafs, default_maf)
#    } else if(sum(study@ancestries$number_of_individuals, rm.na = T) >= 10000 & sum(study@ancestries$number_of_individuals, rm.na = T) <= 50000) {
#      mafs = c(mafs, 0.1)
#    } else if(sum(study@ancestries$number_of_individuals, rm.na = T) > 50000 & sum(study@ancestries$number_of_individuals, rm.na = T) <= 100000) {
#      mafs = c(mafs, 0.05)
#    } else if(sum(study@ancestries$number_of_individuals, rm.na = T) > 100000 ) {
#      mafs = c(mafs, 0.01)
#    }
#
#  }
#}


write.table(data.frame(x =study_ids,y=mafs ),output, quote = F, sep = "\t", row.names = F, col.names = F)
