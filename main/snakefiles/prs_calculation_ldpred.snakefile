import os

ldpred_path = config["ldpred_path"]

if ldpred_path is None:
    ldpred_path = config["repository"] + "/tools/ldpred2"

def read_studies(path):
    with open(path) as file:
        return(file.read().splitlines())

def studies_to_calculate():
    studies_map = read_studies(config['studies_to_calculate_ldpred'])
    return set(studies_map)

rule all:
    input:
        config['results_path_ldpred']+'/'+config['results_directory_name_ldpred']+'/'+config['results_data_table_name_ldpred']+'.tsv'

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
        rds = config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.' + config['imputation_mode'] + '.nomiss.rds',
        gwas = ancient(config['gwas_data_path_ldpred'] + '/{study}.qced.h.tsv.gz')
    output:
        config['results_path_ldpred']+'/'+config['results_directory_name_ldpred'] + '/{study}' + '.' + config['mode']
    conda: "../environment.yaml"
    shell:
        """
        mkdir -p {config[results_path_ldpred]}/{config[results_directory_name_ldpred]}
        mkdir -p {config[results_path_ldpred]}/{config[results_directory_name_ldpred]}/tmp
        mkdir -p {config[results_path_ldpred]}/{config[results_directory_name_ldpred]}/tmp/{wildcards.study}
        Rscript {ldpred_path}/ldpred2.R \
            --ldpred-mode {config[mode]} \
            --col-stat BETA \
            --col-stat-se SE \
            --stat-type BETA \
            --sumstats {input.gwas} \
            --out {config[results_path_ldpred]}/{config[results_directory_name_ldpred]}/{wildcards.study}.{config[mode]} \
            --cores 1 \
            --genomic-build hg38 \
            --name-score {wildcards.study} \
            --col-n N \
            --col-pvalue P \
            --col-bp BP \
            --col-A1 A1 \
            --col-A2 A2 \
            --col-snp-id VARID \
            --col-chr CHR \
            --ld-meta-file {config[ldpred2_ref]}/map_hm3_plus.rds  \
            --ld-file {config[ldpred2_ref]}/ldref_hm3_plus/LD_with_blocks_chr@.rds \
            --geno-file-rds {input.rds} \
            --tmp-dir {config[results_path_ldpred]}/{config[results_directory_name_ldpred]}/tmp/{wildcards.study} \
            --hyper-p-max {config[hyper_p_max]} || true
        rm -f -r {config[results_path_ldpred]}/{config[results_directory_name_ldpred]}/tmp/{wildcards.study}
        touch {config[results_path_ldpred]}/{config[results_directory_name_ldpred]}/{wildcards.study}.{config[mode]}
        """

rule create_pgs_data_table:
    input:
        expand(config['results_path_ldpred']+'/'+config['results_directory_name_ldpred']+'/{study}'+ '.' + config['mode'] , study = studies_to_calculate())
    output:
        config['results_path_ldpred']+'/'+config['results_directory_name_ldpred']+'/'+config['results_data_table_name_ldpred']+'.tsv'
    conda: "../environment.yaml"
    shell:
        """
        Rscript {config[repository]}/scripts/create_prs_datatable_ldpred.R {config[results_path_ldpred]}/{config[results_directory_name_ldpred]} {config[mode]} {config[studies_to_calculate_ldpred]} {config[intercept_max]} {config[intercept_min]} {config[heritability]} {config[results_data_table_name_ldpred]}
        """
