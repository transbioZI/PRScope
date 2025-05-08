import os

def studies_to_calculate():
    with open(config['studies_to_calculate']) as file:
        studies = file.read().splitlines()
        return(studies)

rule all:
    input:
        expand(config['working_directory'] + '/' + config['results_name']+'/{study}' + '.weights' , study = studies_to_calculate())

rule calculate_LDpred_PRS:
    params:
        rds = config['preprocessed_target_data'],
        gwas = ancient(config['transformed_path'] + '/{study}.transformed.h.tsv')
    output:
        config['working_directory'] + '/' + config['results_name'] + '/{study}' + '.weights'
    threads:
        5
    shell:
        """
        mkdir -p {config[working_directory]}/{config[results_name]}
        mkdir -p {config[working_directory]}/{config[results_name]}/tmp
        mkdir -p {config[working_directory]}/{config[results_name]}/tmp/{wildcards.study}
        Rscript {config[repository_path]}/tools/ldpred2/calculateOnlyWeights.R \
            --ldpred-mode {config[mode]} \
            --col-stat BETA \
            --col-stat-se standard_error \
            --stat-type BETA \
            --sumstats {params.gwas} \
            --out {config[working_directory]}/{config[results_name]}/{wildcards.study}.weigths \
            --cores 5 \
            --genomic-build hg19 \
            --set-seed 43536435 \
            --name-score {wildcards.study} \
            --col-n N \
            --col-pvalue p_value \
            --col-bp POS \
            --col-A1 Allele1 \
            --col-A2 Allele2 \
            --col-snp-id SNPID \
            --col-chr CHR \
            --ld-meta-file {config[ldpred2_ref]}/map_hm3_plus.rds  \
            --ld-file {config[ldpred2_ref]}/ldref_hm3_plus/LD_with_blocks_chr@.rds \
            --geno-file-rds {params.rds} \
            --tmp-dir {config[working_directory]}/{config[results_name]}/tmp/{wildcards.study} \
            --hyper-p-max 0.9
        rm -f -r {config[working_directory]}/{config[results_name]}/tmp/{wildcards.study}
        touch {config[working_directory]}/{config[results_name]}/{wildcards.study}.{config[mode]}
        """
