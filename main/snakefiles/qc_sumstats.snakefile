import os
import urllib.request

def studies_to_calculate():
    with open(config['study_list']) as file:
        studies = file.read().splitlines()
    return studies

def read_harmonised_list():
    os.system("mkdir -p " + config['output_path_qced_gwas'])
    os.system("wget -O " + config['output_path_qced_gwas']+"/harmonised_list.txt "+ config["harmonised_list"])
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
        expand(config['output_path_qced_gwas'] + "/{study}.qced.h.tsv", study = get_keys())

rule download_study:
    output:
        config['output_path_qced_gwas'] + "/{study}.to_qc.h.tsv"
    params:
        link = get_link,
        output_inprogress = config['output_path_qced_gwas'] + "/{study}_inprogress.h.tsv.gz",
        output_done = config['output_path_qced_gwas'] + "/{study}.to_qc.h.tsv.gz"
    shell:
        """
        wget {params.link} -O {params.output_inprogress}
        mv {params.output_inprogress} {params.output_done}
        gzip -d {params.output_done}
        """

rule transform_study:
    input:
        rules.download_study.output
    output:
        config['output_path_qced_gwas'] + "/{study}.qced.h.tsv"
    params:
        script = config['repository'] + "/scripts/qc_sumstats.R",
        maf_file = config['maf_file']
    shell:
        """
        Rscript {params.script} {input} {output} {wildcards.study} NA {params.maf_file}
        rm {input}
        """
