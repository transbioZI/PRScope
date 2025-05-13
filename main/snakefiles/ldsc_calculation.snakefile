import os

def read_harmonised_list():
    harm = config['working_directory']+"/harmonised_list.txt"
    if not os.path.isfile(harm):
        urllib.request.urlretrieve(config['harmonised_list'], harm)
    with open(harm) as file:
        download_links = file.read().splitlines()
        study_download_link = dict()
        for link in download_links:
            study = link.split("/")[2]
            study_download_link[study] = link
    return study_download_link

def studies_harmonised():
    stds = studies_to_calculate()
    study_download_link = read_harmonised_list()
    intersection = list(set(study_download_link.keys()) & set(stds))
    return {key: study_download_link[key] for key in intersection}

def get_link(wildcards):
    studies = studies_harmonised()
    return "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"+ studies[wildcards.study][2:]

def get_keys():
    return list(studies_harmonised().keys())

def studies_to_calculate():
    with open(config['studies_to_calculate']) as file:
        studies = file.read().splitlines()
        return(studies)

rule all:
    input:
        expand(config['working_directory'] + '/' + config['results_name']+'/{study}' + '.her' , study = get_keys())


rule convert_plink_file_rds:
    input:
        config['preprocessed_target_data']+'.bed'
    output:
        config['preprocessed_target_data']+'.rds',
        config['preprocessed_target_data']+'.bk'
    threads:
        1
    shell:
        """
        Rscript {config[repository_path]}/tools/ldpred2/createBackingFile.R --file-input {input} --file-output {config[preprocessed_target_data]}.rds
        """

rule impute:
    input: rules.convert_plink_file_rds.output
    output:
        config['preprocessed_target_data'] + '.' + config['imputation_mode'] + '.nomiss.rds',
        config['preprocessed_target_data'] + '.' + config['imputation_mode'] + '.nomiss.bk'
    threads: config['max_threads']
    shell:
        """
        cp {config[preprocessed_target_data]}.rds {config[preprocessed_target_data]}.{config[imputation_mode]}.nomiss.rds
        cp {config[preprocessed_target_data]}.bk {config[preprocessed_target_data]}.{config[imputation_mode]}.nomiss.bk
        Rscript {config[repository_path]}/tools/ldpred2/imputeGenotypes.R --impute-simple {config[imputation_mode]} --geno-file-rds {config[preprocessed_target_data]}.{config[imputation_mode]}.nomiss.rds --cores {config[max_threads]}
        """

rule calculate_LDSC:
    input:
        rds = config['preprocessed_target_data'] + '.' + config['imputation_mode'] + '.nomiss.rds',
        gwas = ancient(config['transformed_path'] + '/{study}.transformed.h.tsv')
    output:
        config['working_directory'] + '/' + config['results_name'] + '/{study}' + '.her'
    threads:
        1
    shell:
        """
        mkdir -p {config[working_directory]}/{config[results_name]}
        mkdir -p {config[working_directory]}/{config[results_name]}/tmp
        mkdir -p {config[working_directory]}/{config[results_name]}/tmp/{wildcards.study}
        Rscript {config[repository_path]}/tools/ldpred2/calculate_heritability_intercept.R \
            --ldpred-mode {config[mode]} \
            --col-stat BETA \
            --col-stat-se standard_error \
            --stat-type BETA \
            --sumstats {input.gwas} \
            --out {config[working_directory]}/{config[results_name]}/{wildcards.study}.{config[mode]} \
            --cores 1 \
            --genomic-build hg38 \
            --set-seed 23545454 \
            --name-score {wildcards.study} \
            --col-n N \
            --col-pvalue p_value \
            --col-bp POS \
            --col-A1 Allele1 \
            --col-A2 Allele2 \
            --col-snp-id VARID \
            --col-chr CHR \
            --ld-meta-file {config[ldpred2_ref]}/map_hm3_plus.rds  \
            --ld-file {config[ldpred2_ref]}/ldref_hm3_plus/LD_with_blocks_chr@.rds \
            --geno-file-rds {input.rds} \
            --tmp-dir {config[working_directory]}/{config[results_name]}/tmp/{wildcards.study} \
            --hyper-p-max 0.9
        rm -f -r {config[working_directory]}/{config[results_name]}/tmp/{wildcards.study}
        touch {config[working_directory]}/{config[results_name]}/{wildcards.study}.her
        """
