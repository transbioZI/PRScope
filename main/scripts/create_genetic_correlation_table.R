rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(plyr))

args = commandArgs(trailingOnly=TRUE)
genetic_cor_study_list = args[1]
rg_max = as.numeric(args[2])
output_genetic_correlation_gwas_list = args[3]

st_df = fread(genetic_cor_study_list)
st_df = st_df[st_df$heritability_passed == TRUE, ]
st_list = st_df$study_id

for(studi in st_list) {
  genetic_cor_results = st_df[which(st_df$study_id == studi),]$path
  current_file_read = readLines(paste0(genetic_cor_results,"/genetic_correlation/",studi,".log"))
  if(length(current_file_read) != 0 ) {
    ind_start = grep('^Summary of Genetic Correlation Results', current_file_read)+1
    ind_end = grep('^Analysis finished', current_file_read)-2

    rf_to_df = current_file_read[ind_start:ind_end]
    writeLines(rf_to_df, paste0(genetic_cor_results,"/genetic_correlation/",studi,".result.log"))
    df = fread(paste0(genetic_cor_results,"/genetic_correlation/",studi,".result.log"))

    df$p1 = unlist(sapply(sapply(df$p1,basename),function(x) {
      str_split(x,"\\.")[[1]][1]
    }))

    df$p2 = unlist(sapply(sapply(df$p2,basename),function(x) {
      str_split(x,"\\.")[[1]][1]
    }))

    df$id = sapply(c(1:nrow(df)),function(x) {
      paste(sort(c(df$p1[x],df$p2[x])),collapse = "_")
    })

    dire = paste0(genetic_cor_results,"/genetic_correlation/pairs")

    dir.create(dire, showWarnings = FALSE)

    for(i in 1:nrow(df)) {
      write.table(df[i,], paste0(dire,"/",df[i,]$p1,"_",df[i,]$p2,".cor"), row.names = F, col.names = T, quote = F, sep = "\t")
    }
  }
}

df_list = list()
counter = 1
for(i in 1:length(st_list)) {
  stu1 = st_list[i]
  genetic_cor_results = st_df[which(st_df$study_id == stu1),]$path
  for(j in i:length(st_list)) {
    stu2 = st_list[j]
    if(stu1 != stu2) {
      genetic_cor_results2 = st_df[which(st_df$study_id == stu2),]$path
      if(file.exists(paste0(genetic_cor_results,"/genetic_correlation/pairs/",stu1,"_",stu2,".cor"))) {
        df_list[[counter]] = fread(paste0(genetic_cor_results,"/genetic_correlation/pairs/",stu1,"_",stu2,".cor"), data.table = F)
      } else {
        df_list[[counter]] = fread(paste0(genetic_cor_results2,"/genetic_correlation/pairs/",stu2,"_",stu1,".cor"), data.table = F)
      }
      counter = counter + 1
    }
  }
}

dtf = do.call(rbind,df_list)
dtf = dtf[which(dtf$p1 %in% st_list),]
dtf = dtf[which(dtf$p2 %in% st_list),]
dtf = dtf[!duplicated(dtf$id),]

write.table(dtf,paste0(output_genetic_correlation_gwas_list,".all_genetic_correlations.tsv"), sep="\t", row.names = F, col.names = T, quote=F)

high_rg <- which(abs(dtf$rg) >= rg_max)

df_high <- dtf[high_rg,]

if(dim(df_high)[1]!=0) {

# Filtering scheme
iter <- 1000
set.seed(123293)
seeds <- sample(1:100000, iter)

rm_list <- list()
for(i in 1:iter) {
  
  set.seed(seeds[i])
  shuf <- df_high[sample(1:nrow(df_high), nrow(df_high)),]
  removed <- c()

  for(j in 1:nrow(df_high)) {
      g1 = shuf$p1[j]
      g2 = shuf$p2[j]
    
      if(length(intersect(c(g1,g2), removed)) == 0) {
        st1 = st_df[st_df$study_id == g1,]
        st2 = st_df[st_df$study_id == g2,]
        z1_score = st1$z_score
        z2_score = st2$z_score
        
        removed[length(removed)+1] <- ifelse( z1_score < z2_score , g1, g2)
        if (z1_score == z2_score ) {
          removed[length(removed)] <- sample(c(g1,g2), 1)
        }
      }
  }

  rm_list[[i]] <- removed
}

rm_lengths <- sapply(rm_list, length)

shortest <- which(rm_lengths == min(rm_lengths))

# All are the same

final_list <- rm_list[[shortest[1]]]
final_df <- dtf[-which(dtf$p1 %in% final_list | dtf$p2 %in% final_list),]
df_high = as.data.frame(df_high)
new_df = data.frame(study_id = st_list, correlated_with = NA)
for(i in c(1:length(st_list))) {
  studi = st_list[i]
  dxc = df_high[grepl(studi,df_high$id),]
  corre_with = unique(c(dxc$p1,dxc$p2))
  corre_with = corre_with[-which(corre_with == studi)]
  new_df[i,]$correlated_with = paste0(corre_with,collapse = ",")
}

stopifnot(identical(st_df$study_id,new_df$study_id))

df = join(st_df,new_df,by = "study_id")

stopifnot(identical(df$study_id,st_list))

res = unique(c(final_df$p1, final_df$p2))

passed = rep(FALSE,nrow(df))
passed[match(res,st_list)] = TRUE
df$genetic_correlation_passed = passed

} else {
  df = st_df
  df$genetic_correlation_passed = TRUE
  df$correlated_with = ""
}

apply_which_false = function(condition, criteria) {
  return(sapply(condition,function(x) {
    ifelse(x, "", paste0(criteria,": FAIL"))
  }))
}

a = apply_which_false(df$genetic_correlation_passed == TRUE, "GC")

cm = sapply(c(1:length(a)), function(x) {
  res=""
  if(!is.na(df$qc_passed_comment[x]) &  df$qc_passed_comment[x] != "") {
    res = df$qc_passed_comment[x]
  }
  
  if(!is.na(a[x]) &  a[x] != "") {
    if(res != "") {
      res = paste0(res, "#", a[x], " (genetic correlation) ")
    } else {
      res = paste0(a[x], " (genetic correlation) ")
    }
  }
  res
})

df$qc_passed_comment = cm
df = df[,c("study_id","qc_passed_comment","genetic_correlation_passed","correlated_with")]

st_df_2 = fread(genetic_cor_study_list, data.table = F)
st_df_2$genetic_correlation_passed = FALSE
st_df_2$correlated_with = ""
st_df_2[match(df$study_id,st_df_2$study_id),]$qc_passed_comment = as.character(df$qc_passed_comment)
st_df_2[match(df$study_id,st_df_2$study_id),]$genetic_correlation_passed = df$genetic_correlation_passed
st_df_2[match(df$study_id,st_df_2$study_id),]$correlated_with = df$correlated_with

write.table(st_df_2,paste0(output_genetic_correlation_gwas_list,".genetic_correlation.tsv"), quote = F,col.names = T, row.names = F, sep = "\t")

for(studi in st_list) {
  genetic_cor_results = st_df[which(st_df$study_id == studi),]$path
  system(paste0("rm ",genetic_cor_results,"/genetic_correlation/",studi,".log"))
  system(paste0("rm ",genetic_cor_results,"/genetic_correlation/",studi,".result.log"))
}
