import os
import urllib.request
import pandas

os.system("rm -f " +  config['output_path_qced_gwas'] + "/*_inprogress*")

def studies_to_calculate():
    csvFile = pandas.read_csv(config["study_list"], sep='\t', engine='python')
    csvFile["sample_size"] = csvFile['sample_size'].astype('int')
    return csvFile["study_id"].tolist()

def get_sample_size(wildcards):
    csvFile = pandas.read_csv(config["study_list"], sep='\t', engine='python')
    index_of_st = csvFile["study_id"].tolist().index(str(wildcards))
    sample_sizes = csvFile["sample_size"].tolist()
    return int(sample_sizes[index_of_st])

def read_harmonised_list():
    os.system("mkdir -p " + config['output_path_qced_gwas'])
    harm = config['output_path_qced_gwas']+"/harmonised_list.txt"
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

rule all:
    input:
        config['study_list'] + ".qced"

rule download_study:
    output:
        config['output_path_qced_gwas'] + "/{study}_inprogress.h.tsv.gz"
    params:
        link = get_link,
        output_inprogress = config['output_path_qced_gwas'] + "/{study}_inprogress.h.tsv.gz",
    shell:
        """
        wget {params.link} -O {params.output_inprogress}
        """

rule gzip_study:
    input:
        rules.download_study.output
    output:
        config['output_path_qced_gwas'] + "/{study}.to_qc.h.tsv"
    params:
        output_inprogress = config['output_path_qced_gwas'] + "/{study}_inprogress.h.tsv.gz",
        output_done = config['output_path_qced_gwas'] + "/{study}.to_qc.h.tsv.gz"
    shell:
        """
        mv {params.output_inprogress} {params.output_done}
        gzip -d {params.output_done}
        """

rule transform_study:
    input:
        rules.gzip_study.output
    conda: "../environment.yaml"
    output:
        config['output_path_qced_gwas'] + "/{study}.qced.h.tsv.gz"
    params:
        script = config['repository'] + "/scripts/qc_sumstats.R",
        maf_file = config['maf_file'],
        N = get_sample_size
    shell:
        """
        Rscript {params.script} {input} {output} {params.maf_file} {params.N}
        rm {input}
        """

rule create_studies_metadata:
    input:
        expand(config['output_path_qced_gwas'] + "/{study}.qced.h.tsv.gz", study = get_keys())
    conda: "../environment.yaml"
    output:
        config['study_list']+".qced"
    params:
        number_of_snps_script = config['repository'] + "/scripts/create_metadata_files_of_GWAS.R"
    shell:
        """
        Rscript {params.number_of_snps_script} {config[output_path_qced_gwas]} {config[study_list]} {config[number_of_snps_after_qc]}
        """

