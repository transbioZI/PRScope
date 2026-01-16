import os
import pandas as pd

ldpred_path = config["ldpred_path"]

if ldpred_path is None:
    ldpred_path = config["repository"] + "/tools/ldpred2"

def studies_to_calculate():
    csvFile = pd.read_csv(config["studies_to_calculate_ldpred"], sep='\t', engine='python').set_index("study_id")
    csvFile.dropna(subset=['genetic_correlation_passed'], inplace=True)
    csvFile = csvFile[csvFile['genetic_correlation_passed'] == True]
    return csvFile

def studies_to_calculate_list():
    return list(studies_to_calculate().index)

def get_path(wildcards):
    csvFile = studies_to_calculate()
    return str(csvFile.loc[wildcards.study].path)

def get_effective_sample_size(wildcards):
    csvFile = studies_to_calculate()
    return float(csvFile.loc[wildcards.study].neff)

rule all:
    input:
        config['results_path_ldpred'] + '/results' + '.tsv'

rule convert_plink_file_rds:
    input:
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.bed'
    conda: "../environment.yaml"
    output:
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.rds',
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.bk'
    shell:
        """
        Rscript {ldpred_path}/createBackingFile.R --file-input {input} --file-output {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.rds
        """

rule impute:
    input: rules.convert_plink_file_rds.output
    output:
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.' + config['imputation_mode'] + '.nomiss.rds',
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.' + config['imputation_mode'] + '.nomiss.bk'
    conda: "../environment.yaml"
    shell:
        """
        cp {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.rds {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.{config[imputation_mode]}.nomiss.rds
        cp {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.bk {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.{config[imputation_mode]}.nomiss.bk
        Rscript {ldpred_path}/imputeGenotypes.R --impute-simple {config[imputation_mode]} --geno-file-rds {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.{config[imputation_mode]}.nomiss.rds
        """

rule calculate_LDSC:
    input:
        rds = config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.' + config['imputation_mode'] + '.nomiss.rds'
    output:
        config['results_path_ldpred'] + '/{study}' + '.' + config['mode']
    params:
        gwas = get_path,
        neff = get_effective_sample_size
    conda: "../environment.yaml"
    shell:
        """
        mkdir -p {config[results_path_ldpred]}
        mkdir -p {config[results_path_ldpred]}/tmp
        mkdir -p {config[results_path_ldpred]}/tmp/{wildcards.study}
        Rscript {ldpred_path}/ldpred2.R \
            --ldpred-mode {config[mode]} \
            --col-stat BETA \
            --col-stat-se SE \
            --stat-type BETA \
            --sumstats {params.gwas}/ldpred/{wildcards.study}.qced.h.tsv.gz \
            --out {config[results_path_ldpred]}/{wildcards.study}.{config[mode]} \
            --cores 5 \
            --genomic-build hg38 \
            --name-score {wildcards.study} \
            --col-pvalue P \
            --col-bp BP \
            --col-A1 A1 \
            --col-A2 A2 \
            --col-snp-id VARID \
            --col-chr CHR \
            --ld-meta-file {config[ldpred2_ref]}/map_hm3_plus.rds  \
            --ld-file {config[ldpred2_ref]}/ldref_hm3_plus/LD_with_blocks_chr@.rds \
            --geno-file-rds {input.rds} \
            --effective-sample-size {params.neff} \
            --tmp-dir {config[results_path_ldpred]}/tmp/{wildcards.study} 2>&1 | tee {config[results_path_ldpred]}/{wildcards.study}.{config[mode]}.log
        rm -f -r {config[results_path_ldpred]}/tmp/{wildcards.study}
        """

rule create_pgs_data_table:
    input:
        expand(config['results_path_ldpred'] + '/{study}' + '.' + config['mode'], study = studies_to_calculate_list())
    output:
        config['results_path_ldpred'] + '/results' + '.tsv'
    conda: "../environment.yaml"
    shell:
        """
        Rscript {config[repository]}/scripts/create_prs_datatable_ldpred.R {config[results_path_ldpred]} {config[mode]} {config[studies_to_calculate_ldpred]}
        """
