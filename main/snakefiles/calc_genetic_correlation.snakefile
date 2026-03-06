import os
import pandas

ldsc_path = config["ldsc_path"]

if ldsc_path is None:
    ldsc_path = config["repository"] + "/tools/ldsc"

def studies_to_calculate():
    csvFile = pandas.read_csv(config["study_list_for_genetic_correlation"], sep='\t', engine='python').set_index("study_id")
    csvFile = csvFile[csvFile['heritability_passed'] == True]
    return csvFile

def studies_to_calculate_list():
    return list(studies_to_calculate().index)

def get_genetic_correlation_command(wildcards):
    all_studies = studies_to_calculate_list()
    all_paths = get_path()
    index_of_st = all_studies.index(wildcards.study)
    command_str = "--rg " + wildcards.study_path + "/munged/" + wildcards.study + ".sumstats.gz"
    if (index_of_st+1) == len(all_studies):
        return "--h2 " + wildcards.study_path + "/munged/" + wildcards.study + ".sumstats.gz"

    commands = list()
    commands.append(command_str)
    for x in range(index_of_st+1, len(all_studies)):
        commands.append(all_paths[x] + "/munged/" + all_studies[x] + ".sumstats.gz")

    return ",".join(commands)

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
        config["output_genetic_correlation_gwas_list"] + ".genetic_correlation.tsv"

rule munge_study:
    conda: "../environment_for_ldsc.yaml"
    output:
        "{study_path}" + "/munged/{study}.sumstats.gz"
    params:
        genome = get_genome_build
    shell:
        """
        if [[ {params.genome}   == "hg38" ]]; then
            python2 {ldsc_path}/munge_sumstats.py --chunksize {config[chunksize]} --sumstats {wildcards.study_path}/ldpred/{wildcards.study}.qced.h.tsv.gz --N-col NEFF --out {wildcards.study_path}/munged/{wildcards.study} --merge-alleles {config[hm3_path_hg38]} --ignore VARID,OR,EAF,Z_SCORE,SE,N
        else
            python2 {ldsc_path}/munge_sumstats.py --chunksize {config[chunksize]} --sumstats {wildcards.study_path}/ldpred/{wildcards.study}.qced.h.tsv.gz --N-col NEFF --out {wildcards.study_path}/munged/{wildcards.study} --merge-alleles {config[hm3_path_hg37]} --ignore VARID,OR,EAF,Z_SCORE,SE,N
        fi
        """

rule calculate_genetic_correlation:
    input:
        "{study_path}" + "/munged/{study}.sumstats.gz"
    conda: "../environment_for_ldsc.yaml"
    params:
        command_str = get_genetic_correlation_command,
        genome = get_genome_build
    output:
        "{study_path}" + "/genetic_correlation/{study}.log"
    shell:
        """
        if [[ {params.genome}   == "hg38" ]]; then
            python2 {ldsc_path}/ldsc.py {params.command_str} --ref-ld-chr {config[ld_ref_hg38]}/ --w-ld-chr {config[ld_ref_hg38]}/ --out {wildcards.study_path}/genetic_correlation/{wildcards.study}
        else
            python2 {ldsc_path}/ldsc.py {params.command_str} --ref-ld-chr {config[ld_ref_hg37]}/ --w-ld-chr {config[ld_ref_hg37]}/ --out {wildcards.study_path}/genetic_correlation/{wildcards.study}
        fi
        """

rule create_pairwise_correlation_matrix:
    input:
        expand("{study_path}" + "/genetic_correlation/{study}.log",zip, study = studies_to_calculate_list(), study_path=get_path())
    conda: "../environment_for_ldsc.yaml"
    output:
        config["output_genetic_correlation_gwas_list"] + ".genetic_correlation.tsv"
    shell:
        """
        Rscript {config[repository]}/scripts/create_genetic_correlation_table.R {config[study_list_for_genetic_correlation]} {config[rg_thr]} {config[output_genetic_correlation_gwas_list]}
        """
