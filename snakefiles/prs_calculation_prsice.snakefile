import os

def read_studies(path):
    with open(path) as file:
        return(file.read().splitlines())

def studies_to_calculate():
    studies_map = read_studies(config['working_directory']+"/"+config['studies_to_calculate'])
    return set(studies_map)

rule all:
    input:
        expand(config['working_directory']+'/'+config['results_name']+'/{study}.all_score' , study = studies_to_calculate())

rule calculate_PRS: 
    input:
        sumstat = ancient(config['transformed_path'] + '/{study}.transformed.h.tsv')
    output:
        config['working_directory']+'/'+config['results_name']+'/{study}.all_score'
    params:
        target = config['preprocessed_target_data'],
        out = config['working_directory']+'/'+config['results_name']+'/{study}'
    threads:
        1
    shell:
        """
        {config[repository_path]}/tools/PRSice/PRSice.R \
        --prsice {config[repository_path]}/tools/PRSice/PRSice_linux \
        --base {input.sumstat} \
        --snp VARID \
        --no-default \
        --chr CHR \
        --bp POS \
        --A1 Allele1 \
        --A2 Allele2 \
        --pvalue p_value \
        --bar-levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.1,0.5,1 \
        --fastscore \
        --target {params.target} \
        --clump-kb {config[clump_kb]} \
        --clump-p {config[clump_p]} \
        --clump-r2 {config[clump_r2]} \
        --out {params.out} \
        --memory 4Gb \
        --base-maf BASEMAF:{config[prsice_default_maf]} \
        --thread 1 \
        --beta \
        --print-snp \
        --no-regress \
        --ld {config[ld_reference]} \
        --seed 1640568366 \
        --stat BETA || true
        touch {params.out}.all_score
        touch {params.out}.prsice
        """
