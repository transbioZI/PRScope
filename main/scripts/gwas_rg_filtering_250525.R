# Filter out GWAS based on genetic correlations and heritabilities
# 20.05.2025 Joonas Naamanka

library(data.table)
library(tidyverse)

df <- fread("/zi/flstorage/HITKIP/share/Ersoy/heritability_table_all.tsv", data.table = F)

# Low z score
min_z <- 4

low_z <- unique(c(df$gwas_1[which(df$zvalue_1 < min_z)], df$gwas_2[which(df$zvalue_2 < min_z)]))

length(low_z)

df <- df[-which(df$gwas_1 %in% low_z | df$gwas_2 %in% low_z),]

# High intercept
intercept_max <- 1.05

high_int <- unique(c(df$gwas_1[which(df$intercept_1 > intercept_max)], df$gwas_2[which(df$intercept_2 > intercept_max)]))

# High rg
rg_max <- 0.99

high_rg <- which(abs(df$correlation)> rg_max)


View(df[high_rg,])

df_high <- df[high_rg,]

# Filtering scheme
iter <- 1000
set.seed(123293)
seeds <- sample(1:100000, iter)

rm_list <- list()
for(i in 1:iter) {
  print(i)
  set.seed(seeds[i])
  shuf <- df_high[sample(1:nrow(df_high), nrow(df_high)),]
  removed <- c()
  for(j in 1:nrow(df_high)) {
      if(length(intersect(c(unlist(shuf[j,1:2])), removed)) == 0) {
      removed[length(removed)+1] <- ifelse(shuf$zvalue_1[j] < shuf$zvalue_2[j], shuf$gwas_1[j], shuf$gwas_2[j])
      if (shuf$zvalue_1[j] == shuf$zvalue_2[j] ) {
        removed[length(removed)] <- sample(c(shuf$gwas_1[j],shuf$gwas_2[j]), 1)
      }
      }
  }
  rm_list[[i]] <- removed
}

rm_lengths <- sapply(rm_list, length)

shortest <- which(rm_lengths == min(rm_lengths))

# All are the same

final_list <- rm_list[[shortest[1]]]

final_df <- df[-which(df$gwas_1 %in% final_list | df$gwas_2 %in% final_list),]

# Sanity check
min(final_df$zvalue_1)
min(final_df$zvalue_2)
max(abs(final_df$correlation))

n_distinct(c(final_df$gwas_1, final_df$gwas_2))

# 229 remain





############################################################################
# The code below is only need if there was no unique selection
#
#
#




# No unique optimum, let's use another criterion
z_means <- c()
for(i in rm_list[shortest]) {
  rems <- c(unlist(i))
  myrows <- which(df$gwas_1 %in% rems)
  zs <- aggregate(zvalue_1 ~ gwas_1, data = df[myrows,], FUN = mean)
  if(nrow(zs) != length(rems)) {
    print("WARNING: lengths do not match")
  }
  z_means[length(z_means)+1] <- mean(zs[,"zvalue_1"])
}

final_removal <- shortest[which(z_means == min(z_means))]

final_list <- rm_list[final_removal]

final_df <- df[-which(df$gwas_1 %in% final_list | df$gwas_2 %in% final_list),]

# Sanity check
min(final_df$zvalue_1)
min(final_df$zvalue_2)
max(final_df$correlation)

n_distinct(c(final_df$gwas_1, final_df$gwas_2))

# 225 remain