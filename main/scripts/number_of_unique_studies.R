suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(gwasrapidd))
suppressMessages(library(rtracklayer))
suppressMessages(library(magrittr))
suppressMessages(library(purrr))

list_of_ids = "/zi/flstorage/group_transbio/ersoy.kocak/Projects/PRSCalculator/resources/studies_172.txt"

res = get_studies(study_id = readLines(list_of_ids))

length(unique(res@publications$pubmed_id))
