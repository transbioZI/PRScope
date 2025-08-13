rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(plyr))

args = commandArgs(trailingOnly=TRUE)
genetic_cor_results =  args[1]
genetic_cor_study_list = args[2]
rg_max = as.numeric(args[3])

st_df = fread(genetic_cor_study_list)
st_df = st_df[st_df$heritability_passed == TRUE, ]
st_list = st_df$study_id
df_list = list()

for(studi in st_list[-length(st_list)]) {
  current_file_read = readLines(paste0(genetic_cor_results,"/",studi,".log"))
  ind_start = grep('^Summary of Genetic Correlation Results', current_file_read)+1
  ind_end = grep('^Analysis finished', current_file_read)-2

  rf_to_df = current_file_read[ind_start:ind_end]
  writeLines(rf_to_df, paste0(genetic_cor_results,"/",studi,".result.log"))
  df = fread(paste0(genetic_cor_results,"/",studi,".result.log"))

  df$p1 = unlist(sapply(sapply(df$p1,basename),function(x) {
    str_split(x,"\\.")[[1]][1]
  }))

  df$p2 = unlist(sapply(sapply(df$p2,basename),function(x) {
    str_split(x,"\\.")[[1]][1]
  }))

  
  df$id = sapply(c(1:nrow(df)),function(x) {
    paste(sort(c(df$p1[x],df$p2[x])),collapse = "_")
    
  })
  
  df_list[[which(studi == st_list)]] = df
  
}

dtf = do.call(rbind,df_list)
dtf = dtf[which(dtf$p1 %in% st_list),]
dtf = dtf[which(dtf$p2 %in% st_list),]
dtf = dtf[!duplicated(dtf$id),]

write.table(dtf,paste0(genetic_cor_study_list,".all_genetic_correlations"), sep="\t", row.names = F, col.names = T, quote=F)

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

apply_which_false = function(condition, threshold, criteria) {
  return(sapply(condition,function(x) {
    ifelse(x, "", paste0(criteria,": FALSE threshold: ",as.character(threshold)))
  }))
}

a = apply_which_false(df$genetic_correlation_passed == TRUE, "FALSE", "genetic_correlation")

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

st_df_2 = fread(genetic_cor_study_list)
st_df_2$genetic_correlation_passed = FALSE
st_df_2$correlated_with = ""
st_df_2[match(df$study_id,st_df_2$study_id),]$qc_passed_comment = df$qc_passed_comment
st_df_2[match(df$study_id,st_df_2$study_id),]$genetic_correlation_passed = df$genetic_correlation_passed
st_df_2[match(df$study_id,st_df_2$study_id),]$correlated_with = df$correlated_with

write.table(st_df_2,paste0(genetic_cor_study_list,".genetic_correlation"), quote = F,col.names = T, row.names = F, sep = "\t")
