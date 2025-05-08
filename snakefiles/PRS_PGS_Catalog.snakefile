import os

def studies_to_calculate():
    with open(config['studies_to_calculate']) as file:
        studies = file.read().splitlines()
    return studies

configfile: "config_PGS_Catalog.yaml"

rule all:
    input:
        expand(config['output_path']+'/'+config['test_name']+'/{study}.all_score' , study = studies_to_calculate())

rule transform_study:
    input:
        config['gwas_path'] + "/{study}_hmPOS_GRCh38.txt.gz"
    output:
        config['output_path'] + "/" + "transformed_pgs_catalog" + "/{study}.transformed.h.tsv"
    params:
        script = config['scripts_path'] + "/transform_pgs_catalogue.R",
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule calculate_PRS: 
    input:
        study = rules.transform_study.output
    output:
        config['output_path']+'/'+config['test_name']+'/{study}.all_score'
    params:
        target = config['preprocessed_target_data'],
        out = config['output_path']+'/'+config['test_name']+'/{study}'
    shell:
        """
        /zi/flstorage/group_transbio/ersoy.kocak/Tools/PRSice/PRSice.R \
        --prsice /zi/flstorage/group_transbio/ersoy.kocak/Tools/PRSice/PRSice_linux \
        --base {input.study} \
        --snp VARID \
        --chr CHR \
        --bp POS \
        --A1 Allele1 \
        --A2 Allele2 \
        --pvalue p.value \
        --bar-levels 1 \
        --fastscore \
        --target {params.target} \
        --out {params.out} \
        --thread 32 \
        --no-regress \
        --clump-kb 500 \
        --clump-p 1 \
        --clump-r2 0.2 \
        --beta \
        --ld /zi/flstorage/group_transbio/data/1000Genome/1kGP/EUR/all_hg38.EUR.ID.dedup \
        --seed 1640568366 \
        --stat BETA
        """

#--model dom \
#--missing SET_ZERO \
#--score std
#--base-maf MAF:0.01 \
#--clump-kb 500 \
#--clump-p 1 \
#--clump-r2 0.2 \