rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(data.table))

apply_which_false = function(condition, threshold, criteria) {
    return(sapply(condition,function(x) {
        ifelse(x, "", paste0(criteria,": FALSE threshold: ",as.character(threshold)))
    }))
}

args = commandArgs(trailingOnly=TRUE)
heritability_results_path = args[1]
studies_list_path = args[2]
heritability_threshold_zscore = as.numeric(args[3])
output_path = args[4]
munged_file_path = args[5]

snp_used_thr = 200000
mean_chi_munged_thr = 1.02
complete_results = fread(studies_list_path)

studies = complete_results$study_id

complete_results$mean_chi_munged = sapply(studies, function(study) {
    file_path = paste0(munged_file_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^^Mean chi\\^2 = ', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            as.numeric(str_split(heritability_result[ind]," = ")[[1]][2])
        }
    } else {
        NA
    }
})

complete_results$observed_h2 = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Total Observed scale h2', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            as.numeric(str_split(heritability_result[ind]," ")[[1]][5])
        }
    } else {
        NA
    }
})

complete_results$observed_h2_se = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Total Observed scale h2', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            se_temp = str_split(heritability_result[ind]," ")[[1]][6]
            as.numeric(gsub(")", "", gsub("\\(", "", se_temp)))
        }
    } else {
        NA
    }
})

complete_results$lambda_gc = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Lambda GC', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            as.numeric(str_split(heritability_result[ind]," ")[[1]][3])
        }

    } else {
        NA
    }

})

complete_results$mean_chi2 = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Mean Chi', heritability_result)

        if(length(ind) == 0) {
            NA
        } else {
            as.numeric(str_split(heritability_result[ind]," ")[[1]][3])
        }
    } else {
        NA
    }
})

complete_results$intercept = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Intercept', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            as.numeric(str_split(heritability_result[ind]," ")[[1]][2])
        }

    } else {
        NA
    }
})

complete_results$intercept_se = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Intercept', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            intercept_temp = str_split(heritability_result[ind]," ")[[1]][3]
            as.numeric(gsub(")", "", gsub("\\(", "", intercept_temp)))
        }

    } else {
        NA
    }
})

complete_results$ratio = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Ratio', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {

            suppressWarnings(as.numeric(str_split(heritability_result[ind]," ")[[1]][2]))

        }

    } else {
        NA
    }
})

complete_results$ratio_se = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^Ratio', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            intercept_temp = str_split(heritability_result[ind]," ")[[1]][3]
            suppressWarnings(as.numeric(gsub(")", "", gsub("\\(", "", intercept_temp))))
        }
    } else {
        NA
    }
})

complete_results$snps_used = sapply(studies, function(study) {
    file_path = paste0(heritability_results_path,"/",study,".log")
    if (file.exists(file_path)) {
        heritability_result = readLines(file_path)
        ind = grep('^After merging with regression SNP LD', heritability_result)
        if(length(ind) == 0) {
            NA
        } else {
            as.numeric(str_split(heritability_result[ind]," ")[[1]][7])
        }
    } else {
        NA
    }
})

complete_results$z_score = complete_results$observed_h2/complete_results$observed_h2_se
complete_results$p_value = exp(-0.717*complete_results$z_score - 0.416*complete_results$z_score^2)
tmp = rep(FALSE,length(studies))
tmp[which(complete_results$z_score > heritability_threshold_zscore & complete_results$snps_used > snp_used_thr & complete_results$mean_chi_munged >= mean_chi_munged_thr )] = TRUE
complete_results$heritability_passed = tmp

a = apply_which_false(complete_results$z_score > heritability_threshold_zscore, heritability_threshold_zscore, "z_score")
b = apply_which_false(complete_results$snps_used > snp_used_thr, snp_used_thr, "snp_used")
c = apply_which_false(complete_results$mean_chi_munged, mean_chi_munged, "mean_chi_munged")

all_criteria = data.frame(x= a,y =b, z = c)

comment_heritability = unlist(apply(all_criteria, 1, function(x) {
    str = as.character(x[!is.na(x) & x!=""])
    xtr = paste0(str, collapse = " and ")
    if(!is.na(xtr) & xtr !="") {
        paste0(xtr," (heritability) ")
    } else {
        xtr
    }
}))

cm = sapply(c(1:length(comment_heritability)), function(x) {
    res=""
    if(!is.na(complete_results$qc_passed_comment[x]) &  complete_results$qc_passed_comment[x] != "") {
        res = complete_results$qc_passed_comment[x]
    }

    if(!is.na(comment_heritability[x]) &  comment_heritability[x] != "") {
        if(res != "") {
            res = paste0(res, "#", comment_heritability[x])
        } else {
            res = comment_heritability[x]
        }
    }
    res
})

complete_results$qc_passed_comment = cm
write.table(complete_results,paste0(studies_list_path,".heritability"), quote = F,col.names = T, row.names = F, sep = "\t")
