import os
import pandas

ldsc_path = config["ldsc_path"]

if ldsc_path is None:
    ldsc_path = config["repository"] + "/tools/ldsc"

def studies_to_calculate():
    csvFile = pandas.read_csv(config["study_list_for_heritability"], sep='\t', engine='python')
    csvFile.dropna(subset=['sample_size'], inplace=True)
    csvFile["sample_size"] = csvFile['sample_size'].astype('int')
    csvFile = csvFile[csvFile['qc_passed'] == True]
    csvFile = csvFile[csvFile['sample_size'] > 0]
    return csvFile["study_id"].tolist()

def get_sample_size(st):
    csvFile = pandas.read_csv(config["study_list_for_heritability"], sep='\t', engine='python')
    index_of_st = csvFile["study_id"].tolist().index(str(st))
    sample_sizes = csvFile["sample_size"].tolist()
    return int(sample_sizes[index_of_st])

def get_path(st):
    csvFile = pandas.read_csv(config["study_list_for_heritability"], sep='\t', engine='python')
    index_of_st = csvFile["study_id"].tolist().index(str(st))
    sample_sizes = csvFile["path"].tolist()
    return str(sample_sizes[index_of_st])

rule all:
    input:
        config['output_heritability_gwas_list'] + ".heritability.tsv"

rule munge_study:
    conda: "../environment_for_ldsc.yaml"
    output:
        config['gwas_data_path_heritability'] + "/{study}/munged/{study}.sumstats.gz"
    params:
        study_path = get_path
    shell:
        """
        python2 {ldsc_path}/munge_sumstats.py --chunksize {config[chunksize]} --sumstats {params.study_path}/{wildcards.study}.qced.h.tsv.gz --N-col N --out {config[gwas_data_path_heritability]}/{wildcards.study}/munged/{wildcards.study} --merge-alleles {config[hm3_path]} --ignore VARID,OR,EAF,Z_SCORE,SE,NEFF
        """

rule calculate_heritability:
    input:
        rules.munge_study.output
    conda: "../environment_for_ldsc.yaml"
    output:
        config['output_path_heritability'] + "/{study}.log"
    shell:
        """
        python2 {ldsc_path}/ldsc.py --h2 {config[gwas_data_path_heritability]}/{wildcards.study}/munged/{wildcards.study}.sumstats.gz --ref-ld-chr {config[ld_ref]}/ --w-ld-chr {config[ld_ref]}/ --out {config[output_path_heritability]}/{wildcards.study}
        """

rule filter_heritability:
    input:
        expand(config['output_path_heritability'] + "/{study}.log", study = studies_to_calculate())
    conda: "../environment_for_ldsc.yaml"
    output:
        config['output_heritability_gwas_list'] + ".heritability.tsv"
    params:
        script = config['repository'] + "/scripts/filter_by_heritability.R"
    shell:
        """
        Rscript {params.script} {config[output_path_heritability]} {config[study_list_for_heritability]} {config[heritability_min_zscore]} {config[output_path_heritability]} {config[gwas_data_path_heritability]} {config[output_heritability_gwas_list]}
        """
