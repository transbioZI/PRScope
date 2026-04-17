rm(list=ls())

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(bigsnpr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

input <- args[1]
ldpred_output = args[2]
ldpred_ref = args[3]
genome_build = args[4]
metadata = args[5]
snps_list = args[6]

map_ldref = readRDS(ldpred_ref)

if(genome_build=="hg38") {
  map_ldref$pos_hg19 = map_ldref$pos
  map_ldref$pos = map_ldref$pos_hg38
}

base <- fread(input, showProgress = FALSE, data.table = F)
res = fread(metadata, showProgress = FALSE, data.table = F)
snps = readLines(snps_list)

problematic_p_value = !as.logical(res$problematic_p_value)
problematic_beta = !as.logical(res$problematic_beta)
problematic_MAF_match = !as.logical(res$problematic_MAF_match)
problematic_SE = !as.logical(res$problematic_SE)
problematic_CHR = !as.logical(res$problematic_CHR)

if(all(c(problematic_p_value, problematic_beta, problematic_MAF_match, problematic_SE, problematic_CHR))) {

    colnames(base) = tolower(colnames(base))
    colnames(base)[colnames(base) == "a2"] = "a0"
    colnames(base)[colnames(base) == "bp"] = "pos"
    base_other = base[which(base$snp %in% snps),]
    base = tryCatch({
      snp_match(base, map_ldref)
    }, error = function(err) {
      print(err)
      data.frame(x=c(),y=c())
    })

    if(dim(base)[1] == 0) {
      res$ld_pred_failed = TRUE
    } else {
      res$ld_pred_failed = FALSE
      columns = c("varid","snp","chr","pos","a1","a0","p","maf","beta","or","se","neff")
      base <- base[ , columns]

      if(dim(base_other)[1] != 0) {
        base_other = base_other[,columns]
      }
      
      cat("Test ",dim(base), "\n")
      base = rbind(base,base_other)
      colnames(base)[colnames(base) == "a0"] = "a2"
      colnames(base)[colnames(base) == "pos"] = "bp"
      
      if(length(unique(base$neff)) == 1) {
        cat("Start ",dim(base), "\n")
        TotalNeff = base$neff[1]
        base$neff <- 4/((2*base$maf*(1-base$maf))*base$se^2)
        base$neff <-ifelse(base$neff  > 1.1*TotalNeff, 1.1*TotalNeff, base$neff)
        base$neff<-ifelse(base$neff < 0.5*TotalNeff, 0.5*TotalNeff, base$neff)
        
        sd_val <- with(base, sqrt(2 * maf * (1 - maf)))
        sd_y_est = median(sd_val * base$se * sqrt(base$neff))
        sd_ss = with(base, sd_y_est / sqrt(neff * se^2))
        is_bad <-sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.10) | sd_ss < 0.1 | sd_val < 0.05
        base = base[!is_bad, ]
        
        png(paste0(ldpred_output,".png"), res=300, unit='px',height=2000, width=2000)
          plot_obj <- qplot(sd_val, sd_ss, color = is_bad) +
            theme_bigstatsr() +
            coord_equal() +
            scale_color_viridis_d(direction = -1) +
            geom_abline(linetype = 2, color = "red") +
            labs(x = "Standard deviations in the validation set",
            y = "Standard deviations derived from the summary statistics",
            color = "Removed?")
          print(plot_obj)
        dev.off()
        
        cat("End ",dim(base), "\n")
      }
        
      colnames(base) = toupper(colnames(base))
    }
    
    write.table(res, metadata, quote = F, sep = "\t", col.names = T, row.names = F)
    fwrite(base, file = ldpred_output, quote = F, row.names = F, sep = "\t", compress = "gzip")
    #sd_ldref <- with(base, sqrt(2*maf*(1 - maf)))
    #sd_ss <- with(base, 2 / sqrt(neff * se^2 + beta^2))

    #is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05
    #base <- base[!is_bad, ]
    #print(paste0("3:",dim(base)))

} else {
    res$ld_pred_failed = TRUE
    write.table(res, metadata, quote = F, sep = "\t", col.names = T, row.names = F)
    base = data.frame(x=c("empty"),y=c("empty"))
    fwrite(base, file = ldpred_output, quote = F, row.names = F, sep = "\t", compress = "gzip")
}

