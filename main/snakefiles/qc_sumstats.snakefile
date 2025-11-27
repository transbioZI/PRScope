import os
import urllib.request
import pandas

def studies_to_calculate():
    csvFile = pandas.read_csv(config["study_list"], sep='\t', engine='python').set_index("study_id")
    csvFile["sample_size"] = csvFile['sample_size'].astype('int')
    return csvFile

def studies_to_calculate_list():
    return list(studies_to_calculate().index)

def get_path(wildcards):
    csvFile = studies_to_calculate()
    return str(csvFile.loc[wildcards.study].path)

def get_sample_size(wildcards):
    csvFile = studies_to_calculate()
    return int(csvFile.loc[wildcards.study].sample_size)

def get_neff(wildcards):
    csvFile = studies_to_calculate()
    return int(csvFile.loc[wildcards.study].neff)

def get_genome_build(wildcards):
    csvFile = studies_to_calculate()
    return str(csvFile.loc[wildcards.study].genome_build)

def get_path_all():
    csvFile = studies_to_calculate()
    parent_path = list()
    for i in list(csvFile.index):
        p = csvFile.loc[i].path
        if pandas.isnull(p):
            parent_path.append(str(config["output_path_qced_gwas"] + "/" + i))
        else:
            parent_path.append(str(p))
    return parent_path

def get_gwas_to_qc(wildcards):
    csvFile = studies_to_calculate()
    p = csvFile.loc[wildcards.study].path
    if pandas.isnull(p):
        return str(config["output_path_qced_gwas"] + "/" + wildcards.study + "/" + wildcards.study + ".to_qc.h.tsv")
    else:
        return csvFile.loc[wildcards.study].path_to_qc

def read_harmonised_list():
    os.system("mkdir -p " + config['output_path_qced_gwas'])
    harm = config['output_path_qced_gwas'] + "/harmonised_list.txt"
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
    stds = studies_to_calculate_list()
    study_download_link = read_harmonised_list()
    intersection = list(set(study_download_link.keys()) & set(stds))
    return {key: study_download_link[key] for key in intersection}

studies = studies_harmonised()

def get_link(wildcards):    
    keys = list(studies.keys())
    if wildcards.study in keys:
        return "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"+ studies[wildcards.study][2:]
    else:
        return "downloaded"

def get_link_md5sum(wildcards):
    keys = list(studies.keys())
    if wildcards.study in keys:
        link1 = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"+ studies[wildcards.study][2:]
        md5sum = os.path.dirname(link1) + "/md5sum.txt"
        return md5sum
    else:
        return "downloaded"

rule all:
    input:
        config['output_qc_gwas_list'] + ".qced.tsv"

rule download_study:
    output:
        "{path_qc}" + "/raw/{study}.download.done"
    params:
        link = get_link,
        output_inprogress = "{path_qc}" +  "/raw/{study}.h.tsv.gz",
        md5sum_download = "{path_qc}" +  "/raw/{study}_md5sum_download.txt",
        download_done = "{path_qc}" + "/raw/{study}.download.done"
    shell:
        """
        if [[ {params.link} == "downloaded" ]]; then
            touch {params.md5sum_download} && touch {params.download_done} 
        else
            wget {params.link} -O {params.output_inprogress} --retry-connrefused --tries=10 && md5sum {params.output_inprogress} > {params.md5sum_download} && touch {params.download_done}
        fi
        """

rule download_md5sum:
    output:
        "{path_qc}" +  "/raw/{study}_md5sum_real.txt"
    params:
        md5sum = get_link_md5sum,
        md5sum_real = "{path_qc}" + "/raw/{study}_md5sum_real.txt"
    shell:
        """
        if [[ {params.md5sum}  == "downloaded" ]]; then
            touch {params.md5sum_real}
        else
            wget --spider --force-html -i {params.md5sum} && wget {params.md5sum} -O {params.md5sum_real} --retry-connrefused --tries=10 || touch {params.md5sum_real}
        fi
        """

rule gzip_study:
    input:
        rules.download_study.output
    output:
        "{path_qc}" + "/{study}.gzip.done"
    params:
        link = get_link,
        output_inprogress = "{path_qc}" +"/raw"+ "/{study}.h.tsv.gz",
        output_done = "{path_qc}" + "/{study}.to_qc.h.tsv",
        output_done_file = "{path_qc}" + "/{study}.gzip.done"
    shell:
        """
        if [[ {params.link} == "downloaded" ]]; then
            touch {params.output_done_file}
        else
            gzip -dkc {params.output_inprogress}  > {params.output_done} && touch {params.output_done_file}
        fi
        """

rule qc_gwas:
    input:
        rules.gzip_study.output
    conda: "../environment.yaml"
    output:
        "{path_qc}" + "/qced/{study}.qced.done"
    params:
        link = get_link,
        script = config['repository'] + "/scripts/qc_sumstats.R",
        maf_file = config['maf_file'],
        N = get_sample_size,
        Neff = get_neff,
        otput_other = "{path_qc}" + "/qced/{study}",
        inp = get_gwas_to_qc,
        qc_done = "{path_qc}"  + "/qced/{study}.qced.done",
        out_p = "{path_qc}" + "/qced/{study}.qced.h.tsv.gz"
    shell:
        """
        if [[ {params.link} == "downloaded" ]]; then
            Rscript {params.script} {params.inp} {params.out_p} {params.maf_file} {params.N} {params.otput_other} {params.Neff} && touch {params.qc_done}
        else
            Rscript {params.script} {params.inp} {params.out_p} {params.maf_file} {params.N} {params.otput_other} {params.Neff} && rm {params.inp} && touch {params.qc_done}
        fi
        """

rule qc_gwas_ldpred:
    input:
        rules.qc_gwas.output
    conda: "../environment.yaml"
    output:
        "{path_qc}" + "/ldpred/{study}.qced.done"
    params:
        script = config['repository'] + "/scripts/qc_sumstats_ldpred.R",
        inp = "{path_qc}"  + "/qced/{study}.qced.h.tsv.gz",
        metadata = "{path_qc}"  + "/qced/{study}.metadata.tsv",
        qc_done = "{path_qc}" + "/ldpred/{study}.qced.done",
        out_p = "{path_qc}"  + "/ldpred/{study}.qced.h.tsv.gz",
        genome_build = get_genome_build
    shell:
        """
        Rscript {params.script} {params.inp} {params.out_p} {config[ldpred2_ref]}/map_hm3_plus.rds {params.genome_build} {params.metadata} && touch {params.qc_done}
        """

rule create_studies_metadata:
    input:
        expand("{path_qc}"  + "/qced/{study}.qced.done",zip, study = studies_to_calculate_list(), path_qc=get_path_all()),
        expand("{path_qc}"  + "/ldpred/{study}.qced.done",zip, study = studies_to_calculate_list(), path_qc=get_path_all()),
        expand("{path_qc}"  + "/raw/{study}_md5sum_real.txt",zip, study = studies_to_calculate_list(), path_qc=get_path_all())
    conda: "../environment.yaml"
    output:
        config['output_qc_gwas_list']+".qced.tsv"
    params:
        number_of_snps_script = config['repository'] + "/scripts/create_metadata_files_of_GWAS.R"
    shell:
        """
        Rscript {params.number_of_snps_script} {config[output_path_qced_gwas]} {config[study_list]} {config[number_of_snps_after_qc]} {config[output_qc_gwas_list]}
        """
