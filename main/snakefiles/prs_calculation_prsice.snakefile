import os
import pandas as pd

def studies_to_calculate():
    csvFile = pd.read_csv(config["studies_to_calculate"], sep='\t', engine='python').set_index("study_id")
    csvFile.dropna(subset=['genetic_correlation_passed'], inplace=True)
    csvFile = csvFile[csvFile['genetic_correlation_passed'] == True]
    return csvFile

def studies_to_calculate_list():
    return list(studies_to_calculate().index)

def get_path(wildcards):
    csvFile = studies_to_calculate()
    return str(csvFile.loc[wildcards.study].path)

rule all:
    input:
        config['results_path']+'/'+config['results_directory_name']+'/'+config['results_data_table_name'] + '_' + str(config['min_number_of_snps_included']) + '.tsv'

rule calculate_PRS: 
    input:
        bim_file = config['target_data_path'] + "/" + config['target_data_prefix'] + ".bim"
    conda: "../environment.yaml"
    output:
        config['results_path']+'/'+config['results_directory_name']+'/{study}.all_score'
    params:
        target = config['target_data_path'] + "/" + config['target_data_prefix'],
        out = config['results_path']+'/'+config['results_directory_name']+'/{study}',
        sumstat = get_path
    shell:
        """
        PRSice \
        --base {params.sumstat}/qced/{wildcards.study}.qced.h.tsv.gz \
        --snp VARID \
        --no-default \
        --chr CHR \
        --bp BP \
        --A1 A1 \
        --A2 A2 \
        --pvalue P \
        --bar-levels {config[p_values]} \
        --fastscore \
        --target {params.target} \
        --clump-kb {config[clump_kb]} \
        --clump-p {config[clump_p]} \
        --clump-r2 {config[clump_r2]} \
        --out {params.out} \
        --base-maf MAF:{config[base_maf]} \
        --thread 1 \
        --beta \
        --print-snp \
        --no-regress \
        --ld {config[clumping_reference_path]}/{config[clumping_reference_prefix]} \
        --stat BETA
        """

rule create_pgs_data_table:
    input:
        expand(config['results_path']+'/'+config['results_directory_name']+'/{study}.all_score' , study = studies_to_calculate_list())
    conda: "../environment.yaml"
    output:
        config['results_path']+'/'+config['results_directory_name']+'/'+config['results_data_table_name'] + '_' + str(config['min_number_of_snps_included']) + '.tsv'
    shell:
        """
        Rscript {config[repository]}/scripts/create_prs_datatable_prsice.R {config[results_path]}/{config[results_directory_name]} {config[min_number_of_snps_included]} {config[min_number_of_snps_included_of_the_highest_p_value]} {config[studies_to_calculate]} {config[results_data_table_name]} {config[output_prsice_gwas_list]}
        """
