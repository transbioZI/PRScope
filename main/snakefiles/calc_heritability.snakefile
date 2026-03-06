import os
import pandas

ldsc_path = config["ldsc_path"]

if ldsc_path is None:
    ldsc_path = config["repository"] + "/tools/ldsc"

def studies_to_calculate():
    csvFile = pandas.read_csv(config["study_list_for_heritability"], sep='\t', engine='python').set_index("study_id")
    csvFile.dropna(subset=['sample_size'], inplace=True)
    csvFile["sample_size"] = csvFile['sample_size'].astype('int')
    csvFile = csvFile[csvFile['qc_passed_ldpred'] == True]
    csvFile = csvFile[csvFile['sample_size'] > 0]
    return csvFile

def studies_to_calculate_list():
    return list(studies_to_calculate().index)

def get_path():
    csvFile = studies_to_calculate()
    parent_path = list()
    for i in list(csvFile.index):
        p = csvFile.loc[i].path
        parent_path.append(str(p))
    return parent_path

def get_genome_build(wildcards):
    csvFile = studies_to_calculate()
    return str(csvFile.loc[str(wildcards.study)].genome_build)

rule all:
    input:
        config['output_heritability_gwas_list'] + ".heritability.tsv"

rule munge_study:
    conda: "../environment_for_ldsc.yaml"
    output:
        "{study_path}" + "/munged/{study}.sumstats.gz"
    params:
        genome = get_genome_build
    shell:
        """
        if [[ {params.genome}  == "hg38" ]]; then
            python2 {ldsc_path}/munge_sumstats.py --chunksize {config[chunksize]} --sumstats {wildcards.study_path}/ldpred/{wildcards.study}.qced.h.tsv.gz --N-col NEFF --out {wildcards.study_path}/munged/{wildcards.study} --merge-alleles {config[hm3_path_hg38]} --ignore VARID,OR,EAF,Z_SCORE,SE,N
        else
            python2 {ldsc_path}/munge_sumstats.py --chunksize {config[chunksize]} --sumstats {wildcards.study_path}/ldpred/{wildcards.study}.qced.h.tsv.gz --N-col NEFF --out {wildcards.study_path}/munged/{wildcards.study} --merge-alleles {config[hm3_path_hg37]} --ignore VARID,OR,EAF,Z_SCORE,SE,N
        fi
        """

rule calculate_heritability:
    input:
        rules.munge_study.output
    conda: "../environment_for_ldsc.yaml"
    output:
        "{study_path}"+"/heritability/{study}.log"
    params:
        genome = get_genome_build
    shell:
        """
        if [[ {params.genome}  == "hg38" ]]; then
            python2 {ldsc_path}/ldsc.py --h2 {wildcards.study_path}/munged/{wildcards.study}.sumstats.gz --ref-ld-chr {config[ld_ref_hg38]}/ --w-ld-chr {config[ld_ref_hg38]}/ --out {wildcards.study_path}/heritability/{wildcards.study}
        else
            python2 {ldsc_path}/ldsc.py --h2 {wildcards.study_path}/munged/{wildcards.study}.sumstats.gz --ref-ld-chr {config[ld_ref_hg37]}/ --w-ld-chr {config[ld_ref_hg37]}/ --out {wildcards.study_path}/heritability/{wildcards.study}
        fi
        """

rule filter_heritability:
    input:
        expand("{study_path}" + "/heritability/{study}.log", zip, study = studies_to_calculate_list(), study_path = get_path())
    conda: "../environment_for_ldsc.yaml"
    output:
        config['output_heritability_gwas_list'] + ".heritability.tsv"
    params:
        script = config['repository'] + "/scripts/filter_by_heritability.R"
    shell:
        """
        Rscript {params.script} {config[study_list_for_heritability]} {config[heritability_min_zscore]} {config[output_heritability_gwas_list]}
        """
