import os

def read_studies(path):
    with open(path) as file:
        return(file.read().splitlines())

def studies_to_calculate():
    studies_map = read_studies(config['studies_to_calculate'])
    return set(studies_map)

rule all:
    input:
        config['results_path']+'/'+config['results_directory_name']+'/'+config['results_data_table_name'] + '_' + str(config['min_number_of_snps_included']) + '.tsv'

rule calculate_PRS: 
    input:
        sumstat = config['gwas_data_path'] + '/{study}.qced.h.tsv',
        bim_file = config['target_data_path'] + "/" + config['target_data_prefix'] + ".bim"
    conda: "../environment.yaml"
    output:
        config['results_path']+'/'+config['results_directory_name']+'/{study}.all_score'
    params:
        target = config['target_data_path'] + "/" + config['target_data_prefix'],
        out = config['results_path']+'/'+config['results_directory_name']+'/{study}'
    shell:
        """
        PRSice \
        --base {input.sumstat} \
        --snp VARID \
        --no-default \
        --chr CHR \
        --bp POS \
        --A1 Allele1 \
        --A2 Allele2 \
        --pvalue p_value \
        --bar-levels {config[p_values]} \
        --fastscore \
        --target {params.target} \
        --clump-kb {config[clump_kb]} \
        --clump-p {config[clump_p]} \
        --clump-r2 {config[clump_r2]} \
        --out {params.out} \
        --base-maf BASEMAF:{config[base_maf]} \
        --thread 1 \
        --beta \
        --print-snp \
        --no-regress \
        --ld {config[clumping_reference_path]}/{config[clumping_reference_prefix]} \
        --stat BETA || true
        touch {params.out}.all_score
        touch {params.out}.prsice
        """

rule create_pgs_data_table:
    input:
        expand(config['results_path']+'/'+config['results_directory_name']+'/{study}.all_score' , study = studies_to_calculate())
    conda: "../environment.yaml"
    output:
        config['results_path']+'/'+config['results_directory_name']+'/'+config['results_data_table_name'] + '_' + str(config['min_number_of_snps_included']) + '.tsv'
    shell:
        """
        Rscript {config[repository]}/scripts/create_prs_datatable_prsice.R {config[results_path]}/{config[results_directory_name]} {config[min_number_of_snps_included]} {config[min_number_of_snps_included_of_the_highest_p_value]} {config[studies_to_calculate]} {config[results_data_table_name]}
        """
