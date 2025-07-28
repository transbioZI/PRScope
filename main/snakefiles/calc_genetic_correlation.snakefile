import os
import pandas

ldsc_path = config["ldsc_path"]

if ldsc_path is None:
    ldsc_path = config["repository"] + "/tools/ldsc"

def studies_to_calculate():
    csvFile = pandas.read_csv(config["study_list_for_genetic_correlation"], sep='\t', engine='python')
    csvFile = csvFile[csvFile['heritability_passed'] == True]
    return csvFile["study_id"].tolist()

def get_sample_size(st):
    csvFile = pandas.read_csv(config["study_list_for_genetic_correlation"], sep='\t', engine='python')
    index_of_st = csvFile["study_id"].tolist().index(str(st))
    sample_sizes = csvFile["sample_size"].tolist()
    return int(sample_sizes[index_of_st])

def get_genetic_correlation_command(st):
    st = str(st)
    all_studies = studies_to_calculate()
    index_of_st = all_studies.index(st)
    command_str = "--rg " + config['gwas_data_path_genetic_correlation'] + "/munged/" + st + ".sumstats.gz"
    if (index_of_st+1) == len(all_studies):
        return "--h2 " + config['gwas_data_path_genetic_correlation'] + "/munged/"+st+".sumstats.gz"

    commands = list()
    commands.append(command_str)
    for x in range(index_of_st+1, len(all_studies)):
        commands.append(config['gwas_data_path_genetic_correlation'] + "/munged/" + all_studies[x] + ".sumstats.gz")

    return ",".join(commands)

rule all:
    input:
        config["study_list_for_genetic_correlation"] + ".genetic_correlation"

rule munge_study:
    input:
        config['gwas_data_path_genetic_correlation'] + "/{study}.qced.h.tsv.gz"
    conda: "../environment_for_ldsc.yaml"
    output:
        config['gwas_data_path_genetic_correlation'] + "/munged/{study}.sumstats.gz"
    params:
        sample_size = get_sample_size
    shell:
        """
        python2 {ldsc_path}/munge_sumstats.py --chunksize {config[chunksize]} --sumstats {config[gwas_data_path_genetic_correlation]}/{wildcards.study}.qced.h.tsv.gz --N {params.sample_size} --out {config[gwas_data_path_genetic_correlation]}/munged/{wildcards.study} --merge-alleles {config[hm3_path]} --ignore VARID,OR,EAF,Z_SCORE
        """

rule calculate_genetic_correlation:
    input:
        expand(config['gwas_data_path_genetic_correlation'] + "/munged/{study}.sumstats.gz", study = studies_to_calculate())
    conda: "../environment_for_ldsc.yaml"
    params:
        command_str = get_genetic_correlation_command
    output:
        config['output_path_genetic_correlation'] + "/{study}.log"
    shell:
        """
        python2 {ldsc_path}/ldsc.py {params.command_str} --ref-ld-chr {config[ld_ref]}/ --w-ld-chr {config[ld_ref]}/ --out {config[output_path_genetic_correlation]}/{wildcards.study}
        """

rule create_pairwise_correlation_matrix:
    input:
        expand(config['output_path_genetic_correlation'] + "/{study}.log", study = studies_to_calculate())
    conda: "../environment_for_ldsc.yaml"
    output:
        config["study_list_for_genetic_correlation"] + ".genetic_correlation"
    shell:
        """
        Rscript {config[repository]}/scripts/create_genetic_correlation_table.R {config[output_path_genetic_correlation]} {config[study_list_for_genetic_correlation]} {config[rg_thr]}
        """
