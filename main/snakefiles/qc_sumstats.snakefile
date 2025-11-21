import os
import urllib.request
import pandas

def studies_to_calculate():
    csvFile = pandas.read_csv(config["study_list"], sep='\t', engine='python')
    csvFile["sample_size"] = csvFile['sample_size'].astype('int')
    return csvFile["study_id"].tolist()

def get_sample_size(wildcards):
    csvFile = pandas.read_csv(config["study_list"], sep='\t', engine='python')
    index_of_st = csvFile["study_id"].tolist().index(str(wildcards))
    sample_sizes = csvFile["sample_size"].tolist()
    return int(sample_sizes[index_of_st])

def get_neff(wildcards):
    csvFile = pandas.read_csv(config["study_list"], sep='\t', engine='python')
    index_of_st = csvFile["study_id"].tolist().index(str(wildcards))
    sample_sizes = csvFile["neff"].tolist()
    return int(sample_sizes[index_of_st])

def get_genome_build(wildcards):
    csvFile = pandas.read_csv(config["study_list"], sep='\t', engine='python')
    index_of_st = csvFile["study_id"].tolist().index(str(wildcards))
    sample_sizes = csvFile["genome_build"].tolist()
    return str(sample_sizes[index_of_st])

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

def get_link_md5sum(wildcards):
    studies = studies_harmonised()
    link1 = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"+ studies[wildcards.study][2:]
    md5sum = os.path.dirname(link1) + "/md5sum.txt"
    return md5sum

def get_keys():
    return list(studies_harmonised().keys())

rule all:
    input:
        config['output_qc_gwas_list'] + ".qced.tsv"

rule download_study:
    output:
        config['output_path_qced_gwas'] + "/{study}" + "/raw/{study}.download.done"
    params:
        link = get_link,
        output_inprogress = config['output_path_qced_gwas'] + "/{study}" + "/raw/{study}.h.tsv.gz",
        md5sum_download = config['output_path_qced_gwas'] + "/{study}" + "/raw/{study}_md5sum_download.txt",
        download_done = config['output_path_qced_gwas'] + "/{study}" + "/raw/{study}.download.done"
    shell:
        """
        wget {params.link} -O {params.output_inprogress} --retry-connrefused --tries=10 && md5sum {params.output_inprogress} > {params.md5sum_download} && touch {params.download_done}
        """

rule download_md5sum:
    output:
        config['output_path_qced_gwas'] + "/{study}"  + "/raw/{study}_md5sum_real.txt"
    params:
        md5sum = get_link_md5sum,
        md5sum_real = config['output_path_qced_gwas'] + "/{study}"  + "/raw/{study}_md5sum_real.txt"
    shell:
        """
        wget --spider --force-html -i {params.md5sum} && wget {params.md5sum} -O {params.md5sum_real} --retry-connrefused --tries=10 || touch {params.md5sum_real}
        """

rule gzip_study:
    input:
        rules.download_study.output
    output:
        config['output_path_qced_gwas'] + "/{study}" + "/{study}.gzip.done"
    params:
        output_inprogress = config['output_path_qced_gwas'] +"/{study}/raw"+ "/{study}.h.tsv.gz",
        output_done = config['output_path_qced_gwas'] + "/{study}/{study}.to_qc.h.tsv",
        output_done_file = config['output_path_qced_gwas'] + "/{study}/{study}.gzip.done"
    shell:
        """
        gzip -dkc {params.output_inprogress}  > {params.output_done} && touch {params.output_done_file}
        """

rule qc_gwas:
    input:
        rules.gzip_study.output
    conda: "../environment.yaml"
    output:
        config['output_path_qced_gwas'] + "/{study}" + "/qced/{study}.qced.done"
    params:
        script = config['repository'] + "/scripts/qc_sumstats.R",
        maf_file = config['maf_file'],
        N = get_sample_size,
        Neff = get_neff,
        genome_build = get_genome_build,
        otput_other = config['output_path_qced_gwas'] + "/{study}/qced/{study}",
        inp = config['output_path_qced_gwas'] + "/{study}" + "/{study}.to_qc.h.tsv",
        qc_done = config['output_path_qced_gwas'] + "/{study}" + "/qced/{study}.qced.done",
        out_p = config['output_path_qced_gwas'] + "/{study}" + "/qced/{study}.qced.h.tsv.gz"
    shell:
        """
        Rscript {params.script} {params.inp} {params.out_p} {params.maf_file} {params.N} {params.otput_other} {params.Neff} {params.genome_build} && rm {params.inp} && touch {params.qc_done}
        """
rule qc_gwas_ldpred:
    input:
        rules.qc_gwas.output
    conda: "../environment.yaml"
    output:
        config['output_path_qced_gwas'] + "/{study}" + "/ldpred/{study}.qced.done"
    params:
        script = config['repository'] + "/scripts/qc_sumstats_ldpred.R",
        inp = config['output_path_qced_gwas'] + "/{study}" + "/qced/{study}.qced.h.tsv.gz",
        qc_done = config['output_path_qced_gwas'] + "/{study}" + "/ldpred/{study}.qced.done",
        out_p = config['output_path_qced_gwas'] + "/{study}" + "/ldpred/{study}.qced.h.tsv.gz",
        genome_build = get_genome_build
    shell:
        """
        Rscript {params.script} {params.inp} {params.out_p} {config[ldpred2_ref]}/map_hm3_plus.rds {params.genome_build} && touch {params.qc_done}
        """

rule create_studies_metadata:
    input:
        expand(config['output_path_qced_gwas'] + "/{study}" + "/qced/{study}.qced.done", study = get_keys()),
        expand(config['output_path_qced_gwas'] + "/{study}" + "/ldpred/{study}.qced.done", study = get_keys()),
        expand(config['output_path_qced_gwas'] + "/{study}" + "/raw/{study}_md5sum_real.txt", study = get_keys())
    conda: "../environment.yaml"
    output:
        config['output_qc_gwas_list']+".qced.tsv"
    params:
        number_of_snps_script = config['repository'] + "/scripts/create_metadata_files_of_GWAS.R"
    shell:
        """
        Rscript {params.number_of_snps_script} {config[output_path_qced_gwas]} {config[study_list]} {config[number_of_snps_after_qc]} {config[output_qc_gwas_list]}
        """
