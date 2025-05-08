rm(list=ls())
library("RCurl")

files = readLines("/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/studies_no_harm.txt.noharm")
ftp_link = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

get_number = function(x) {
  as.numeric(gsub("\\D", "", x))
}
res = c()
for(file in files) {
  number = get_number(file)
  mod = number%%1000
  lower = number-mod+1
  upper = number - mod +1000
  path1 = paste0("GCST",lower,"-","GCST",upper)
  path2 = paste0(file,"_buildGRCh38.txt")
  path3 = paste0(file,"_buildGRCh38.txt-meta.yaml")
  path4 = paste0(file,"_buildGRCh38.tsv.gz")
  path5 = paste0(file,"_buildGRCh38.tsv.gz-meta.yaml")
  path6 = paste0(file,"_buildGRCh37.tsv.gz")
  path7 = paste0(file,"_buildGRCh37.tsv.gz-meta.yaml")
  path8 = paste0(file,"_buildGRCh37.txt")
  path9 = paste0(file,"_buildGRCh37.txt.gz-meta.yaml")
  path10 = paste0(file,"_buildGRCh38.tsv")
  path11 = paste0(file,"_buildGRCh37.tsv")
  path12 = paste0(file,"_buildGRCh38.txt.gz")
  path13 = paste0(file,"_buildGRCh37.txt.gz")
  base = paste0(ftp_link,"/",path1,"/",file)
  if(url.exists(base)) {
    url1 = paste0(base,"/",path2)
    url2 = paste0(base,"/",path4)
    url3 = paste0(base,"/",path6)
    url4 = paste0(base,"/",path8)
    url5 = paste0(base,"/",path10)
    url6 = paste0(base,"/",path11)
    url7 = paste0(base,"/",path12)
    url8 = paste0(base,"/",path13)
    if(url.exists(url1)) {
      res = c(res, url1)
    } else if(url.exists(url2)) {
      res = c(res, url2)
    } else if(url.exists(url3)) {
      res = c(res, url3)
    } else if(url.exists(url4)) {
      res = c(res, url4)
    } else if(url.exists(url5)) {
      res = c(res, url5)
    } else if(url.exists(url6)) {
      res = c(res, url6)
    } else if(url.exists(url7)) {
      res = c(res, url7)
    } else if(url.exists(url8)) {
      res = c(res, url8)
    } else {
      res = c(res, file)
    }
  }
}
res = readRDS("/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/rtt.rds")
ftp_link = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

get_number = function(x) {
  as.numeric(gsub("\\D", "", x))
}
rr = c()
for(r in res) {
  number = get_number(r)
  mod = number%%1000
  lower = number-mod+1
  upper = number - mod +1000
  path0 = paste0("GCST",lower,"-","GCST",upper)
  
  base = paste0(ftp_link,"/",path0,"/",r)
  
  path1 = paste0(r,".txt")
  path2 = paste0(r,".tsv")
  path3 = paste0(r,".txt.gz")
  path4 = paste0(r,".tsv.gz")
  
  if(!grepl("http",r, fixed = T)) {
    url1 = paste0(base,"/",path1)
    url2 = paste0(base,"/",path2)
    url3 = paste0(base,"/",path3)
    url4 = paste0(base,"/",path4)
    if(url.exists(url1)) {
      rr = c(rr, url1)
    } else if(url.exists(url2)) {
      rr = c(rr, url2)
    } else if(url.exists(url3)) {
      rr = c(rr, url3)
    } else if(url.exists(url4)) {
      rr = c(rr, url4)
    } else {
      print(r)
    }
  } else {
    rr = c(rr,r)
  }
}

writeLines(rr,"/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/not_harm.txt")
require(doParallel)
registerDoParallel(cores=25)
formatted = "/data/GWAS_Catalog_Not_Harmonized"
library("stringr")
setwd(formatted)
foreach(i=1:length(rr)) %dopar% try(download.file(paste0(rr[i],"-meta.yaml"),destfile=paste0(str_split(rr[i],"/")[[1]][10],"-meta.yaml")))

saveRDS(res,"/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator/resources/rtt.rds")
