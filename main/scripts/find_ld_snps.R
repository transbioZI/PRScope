rm(list=ls())
gc()
library(hash)
library(stringr)
library(LDlinkR)
library(dplyr)
library(hash)
library(data.table)
library(tibble)
plink = "/zi/home/ersoy.kocak/Desktop/Tools/plink"
bfile = "/data/references/1000Genome/reference/EUR/eur_hg38.phase3"

filter_snps = function(snpslist, p_val) {
  snps = read.table(snpslist, sep = "\t",header = T)
  snps = snps[which(snps$P <=  p_val),]
  return(snps$SNP)
}

get_ld_proxies <- function(rsid, bfile, searchspace=NULL, tag_kb=1000, tag_nsnp=99999999, tag_r2=corr_th, threads=30, out=tempfile())
{
  searchspacename <- paste0(out, ".searchspace")
  targetsname <- paste0(out, ".targets")
  outname <- paste0(out, ".targets.ld")
  utils::write.table(rsid, file=targetsname, row.names = FALSE, col.names = FALSE, quote = FALSE)
  if(!is.null(searchspace))
  {
    stopifnot(is.character(searchspace))
    
    utils::write.table(unique(c(rsid, searchspace)), file=searchspacename, row.names = FALSE, col.names = FALSE, quote = FALSE)
    extract_param <- paste0(" --extract ", searchspacename)
  } else {
    extract_param <- " " 
  }
  cmd <- paste0(plink,
                " --bfile ", bfile, 
                extract_param,
                " --keep-allele-order ",
                " --r2 ",
                " --ld-snp-list ", targetsname,
                " --ld-window-kb ", tag_kb,
                " --ld-window-r2 ", tag_r2,
                " --ld-window ", tag_nsnp,
                " --out ", targetsname,
                " --threads ", threads,
                " 2>&1 > /dev/null"
  )
  message("Finding proxies...")
  system(cmd)
  
  if (!file.exists(outname)) {
    ld <- data.frame(CHR_A = integer(), BP_A = integer(), SNP_A = character(), MAF_A = double(), CHR_B = integer(), BP_B = integer(), 
                     SNP_B = character(), PHASE = character(), MAF_B = double(), R = double())
    message("Index SNP not found in the reference panel")
    return(ld)
  }
  ld <- data.table::fread(outname, header=TRUE) %>%
    dplyr::as_tibble(.name_repair="minimal") %>%
    dplyr::filter(.data[["R2"]] > tag_r2) #%>%
    #dplyr::filter(.data[["SNP_A"]] != .data[["SNP_B"]])
  
  unlink(searchspacename)
  unlink(targetsname)
  unlink(paste0(targetsname, c(".log", ".nosex")))
  #unlink(outname)
  if(nrow(ld) == 0)
  {
    message("No proxies found")
    return(ld)
  }
 
  
  return(ld)
}

find_overlap = function(proxy1,proxy2_hash) {
  
  overlap = apply(proxy1,1,function(x) {
    snps = x[[1]]
    found = FALSE
    for(snp in snps) {
      if(is.null(proxy2_hash[[snp]]) == FALSE) {
        found = TRUE
        break
      }
    }
    found
  })
  
  return(overlap)
}

find_overlap_by_list = function(snp_list,proxy2_hash) {
  
  overlap = lapply(snp_list,function(x) {
    found = FALSE
    if(is.null(proxy2_hash[[x]]) == FALSE) {
      found = TRUE
    }
    found
  })
  
  return(overlap)
}
#transbio031 scz
#GCST90029028 neuro
#GCST006952
#GCST006951

corr_th = 0.9
scz_ss <- as_tibble(fread("/data/sumstats/transformed/transbio031.transformed.h.tsv"))
neuro_ss <- as_tibble(fread("/data/sumstats/transformed/GCST90029028.transformed.h.tsv"))
pvalues = c("0.0000001","0.000001","0.00001","0.0001","0.001","0.01","0.05","0.1")
cmd_lines = c()
cmd_lines = c(cmd_lines,"rm /data/projects/on_going/multiprs_patient_stratification/cluster_analysis/overlap_bari/GCST90029028/*.all_score")
for(pp in pvalues ) {
  current_pval = as.numeric(pp)
  scz_ss_filtered <- scz_ss[ scz_ss$p_value <= current_pval,] %>% column_to_rownames('VARID')
  neuro_ss_filtered <- neuro_ss[neuro_ss$p_value <= current_pval,] %>% column_to_rownames('VARID')
  
  snp_list_scz = rownames(scz_ss_filtered)
  out_file_scz = "/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/scz_proxy_snps"
  ld_scz_proxy  = get_ld_proxies(snp_list_scz,bfile, out = out_file_scz)
  ld_scz_proxy_hash = hash(unique(ld_scz_proxy$SNP_B), 1)
  ld_scz_proxy = ld_scz_proxy %>% group_by(SNP_A) %>% summarise(proxies = list(SNP_B))
  ld_scz_proxy = ld_scz_proxy %>% column_to_rownames('SNP_A')
  
  snp_list_neuro = rownames(neuro_ss_filtered)
  out_file_neuro = "/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/neuro_proxy_snps"
  ld_neuro_proxy  = get_ld_proxies(snp_list_neuro,bfile, out = out_file_neuro)
  ld_neuro_proxy_hash = hash(unique(ld_neuro_proxy$SNP_B), 1)
  ld_neuro_proxy = ld_neuro_proxy %>% group_by(SNP_A) %>% summarise(proxies = list(SNP_B))
  ld_neuro_proxy = ld_neuro_proxy %>% column_to_rownames('SNP_A')
  
  scz_neuro_overlap = find_overlap(ld_scz_proxy, ld_neuro_proxy_hash)
  snps = rownames(ld_scz_proxy)
  overlap = snps[which(scz_neuro_overlap == TRUE)]
  non_overlap = snps[which(scz_neuro_overlap == FALSE)]
  
  overlap_scz = ld_scz_proxy[overlap,,drop = FALSE]
  non_overlap_scz = ld_scz_proxy[non_overlap,,drop = FALSE]
  overlap = unique(unlist(overlap_scz))
  non_overlap = unique(unlist(non_overlap_scz))
  
  inter = intersect(unique(unlist(overlap)),unique(unlist(non_overlap)))
  print(paste0("p_value: ", pp))
  print(paste0("length of intersect: ", length(inter)))
  print(paste0("length of overlap: ", length(overlap)))
  print(paste0("length of non-overlap: ", length(non_overlap)))
  print(head(overlap))
  
  overlap = overlap[which(overlap %in% inter == FALSE)]
  non_overlap = non_overlap[which(non_overlap %in% inter == FALSE)]
  total = c(overlap,non_overlap)
  #overlap_scz_hash = hash(overlap[which(overlap %in% inter == FALSE)], 1)
  #non_overlap_scz_hash = hash(non_overlap[which(non_overlap %in% inter == FALSE)], 1)
  
  #clumped_snps = read.table("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/scz.snp", sep ="\t", header = T)
  #clumped_snps = clumped_snps[which( clumped_snps$P <= current_pval),]
  #overlap = find_overlap_by_list(clumped_snps$SNP,overlap_scz_hash)
  #non_overlap = find_overlap_by_list(clumped_snps$SNP,non_overlap_scz_hash)
  
  #overlap = unlist(overlap)
  #non_overlap = unlist(non_overlap)
  
  #overlap = clumped_snps[overlap,]
  #non_overlap = clumped_snps[non_overlap,]
  
  #print(paste0("overlap: ",length(overlap$SNP)))
  #print(paste0("non_overlap: ",length(non_overlap$SNP)))
  #print(paste0("intersect ",length(intersect(overlap$SNP, non_overlap$SNP))))
  p_val_thres = pp
  
  writeLines(overlap,paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/overlap","_",p_val_thres,".txt"))
  writeLines(non_overlap,paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/non_overlap_all","_",p_val_thres,".txt"))
  writeLines(total,paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/overlap_non_overlap","_",p_val_thres,".txt"))
  
  for(i in 1:101) {
    
    if(length(overlap) > length(non_overlap)){
      print(paste0("length not equal, non_ov:",length(non_overlap)," ov:",length(overlap)," p_val:",p_val_thres))
    }
   
    sampling_size = min(length(overlap), length(non_overlap))
    exclude = sample(non_overlap,sampling_size)
    
    writeLines(setdiff(total,exclude),paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/non_overlap_",i,"_",p_val_thres,".txt"))
    writeLines(exclude,paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/excluded_",i,"_",p_val_thres,".txt"))
  }
  
  overlap_cmd = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_bari.yaml --jobs 1 -k --rerun-incomplete --config out_name=overlap_",p_val_thres," extract=/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/overlap","_",p_val_thres,".txt ", "p_val=",p_val_thres," &")
  
  overlap_non_overlap = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_bari.yaml --jobs 1 -k --rerun-incomplete --config out_name=overlap_non_overlap_",p_val_thres," extract=/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/overlap_non_overlap","_",p_val_thres,".txt ","p_val=",p_val_thres," &")
  non_overlap_all = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_bari.yaml --jobs 1 -k --rerun-incomplete --config out_name=non_overlap_all_",p_val_thres," extract=/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/non_overlap_all","_",p_val_thres,".txt ","p_val=",p_val_thres," &")
  
  
  cmd_lines = c(cmd_lines, overlap_cmd)
  cmd_lines = c(cmd_lines, non_overlap_all)
  cmd_lines = c(cmd_lines, overlap_non_overlap)
  
  for(i in 1:101) {
    non_ov = paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/non_overlap_",i,"_",p_val_thres,".txt")
    non_ov_cmd = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_bari.yaml --jobs 1 -k --rerun-incomplete --config extract=",non_ov, " p_val=",p_val_thres," out_name=non_",i,"_",p_val_thres, " &")
    cmd_lines = c(cmd_lines, non_ov_cmd)
    non_ov = paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/excluded_",i,"_",p_val_thres,".txt")
    non_ov_cmd = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_bari.yaml --jobs 1 -k --rerun-incomplete --config extract=",non_ov, " p_val=",p_val_thres," out_name=exclude_",i,"_",p_val_thres, " &")
    cmd_lines = c(cmd_lines, non_ov_cmd)
    if(i%%15 == 0) {
      cmd_lines = c(cmd_lines, "wait")
    }
  }
  
}

writeLines(cmd_lines,paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/overlap_bari/run_overlap_GCST90029028.sh"))

#cmd_lines = c()
#overlap_cmd = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_psycourse.yaml --jobs 1 -k --rerun-incomplete --config out_name=overlap extract=/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/overlap.txt p_val=",p_val_thres, " &")
#overlap_non_overlap = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_psycourse.yaml --jobs 1 -k --rerun-incomplete --config out_name=overlap_non_overlap extract=/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/overlap_non_overlap.txt p_val=",p_val_thres, " &")
#non_overlap_all = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_psycourse.yaml --jobs 1 -k --rerun-incomplete --config out_name=non_overlap_all extract=/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/non_overlap_all.txt p_val=",p_val_thres, " &")
#
#cmd_lines = c(cmd_lines,"rm /data/projects/on_going/multiprs_patient_stratification/cluster_analysis/overlap_psycourse/GCST90029028/*.all_score")
#cmd_lines = c(cmd_lines, overlap_cmd)
#cmd_lines = c(cmd_lines, non_overlap_all)
#cmd_lines = c(cmd_lines, overlap_non_overlap)
#
#for(i in 1:1001) {
#  non_ov = paste0("/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/down_sampling_snps/GCST90029028/non_overlap_",i,".txt")
#  non_ov_cmd = paste0("snakemake -s prs_calculation_prsice.snakefile --configfile config_psycourse.yaml --jobs 1 -k --rerun-incomplete --config extract=",non_ov, " p_val=",p_val_thres," out_name=non_",i, " &")
#  cmd_lines = c(cmd_lines, non_ov_cmd)
#  if(i%%30 == 0) {
#    cmd_lines = c(cmd_lines, "wait")
#  } 
#}
#
#writeLines(cmd_lines,"/data/projects/on_going/multiprs_patient_stratification/cluster_analysis/overlap_psycourse/run_overlap_GCST90029028.sh")