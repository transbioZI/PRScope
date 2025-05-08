import os

def get_raw_target_data():
    target_data = list()
    for file in os.listdir(config['target_data']+'/raw'):
            if file.endswith(".bim"):
                target_data.append(file.removesuffix('.bim'))
    return target_data

rule all:
    input:
        expand([
            config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.FINAL.snp_missingness.pdf',
            config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.FINAL.sample_missingness.pdf',
            config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.FINAL.frq.pdf'],
            target_data_name = get_raw_target_data())

rule remove_samples:
    input:
        multiext(config['target_data']+'/raw/{target_data_name}','.bim', ".bed", ".fam")
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.filtered', '.bim', ".bed", ".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/{wildcards.target_data_name}\
            --make-bed \
            --threads 20 \
            --remove {config[remove_samples]} \
            --out {config[target_data]}/raw/tmp/{wildcards.target_data_name}.filtered"

rule remove_snps_missingness:
    input:
        rules.remove_samples.output
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.snpsremoved', '.bim', ".bed", ".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/tmp/{wildcards.target_data_name}.filtered\
            --geno 0.2 \
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/raw/tmp/{wildcards.target_data_name}.snpsremoved"

rule remove_samples_missingness:
    input:
        rules.remove_snps_missingness.output
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.samplesremoved', '.bim', ".bed", ".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/tmp/{wildcards.target_data_name}.snpsremoved\
            --mind 0.2 \
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/raw/tmp/{wildcards.target_data_name}.samplesremoved"

rule remove_snps_missingness_deeply:
    input:
        rules.remove_samples_missingness.output
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.snpsremoved2', '.bim', ".bed", ".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/tmp/{wildcards.target_data_name}.samplesremoved\
            --geno 0.02 \
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/raw/tmp/{wildcards.target_data_name}.snpsremoved2"

rule remove_samples_missingness_deeply:
    input:
        rules.remove_snps_missingness_deeply.output
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.samplesremoved2', '.bim', ".bed", ".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/tmp/{wildcards.target_data_name}.snpsremoved2 \
            --mind 0.02 \
            --make-bed \
            --threads 20 \
            --allow-no-sex \
            --out {config[target_data]}/raw/tmp/{wildcards.target_data_name}.samplesremoved2"

rule remove_maf:
    input:
        rules.remove_samples_missingness_deeply.output
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.mafremoved', '.bim', ".bed", ".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/tmp/{wildcards.target_data_name}.samplesremoved2 \
            --maf {config[maf]} \
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/raw/tmp/{wildcards.target_data_name}.mafremoved"

rule remove_hwe:
    input:
        rules.remove_maf.output
    output:
        multiext(config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.QC',".fam",".bim",".bed")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/tmp/{wildcards.target_data_name}.mafremoved\
            --hwe 1e-6 \
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC"

rule exclude_high_ld_region:
    input:
        rules.remove_hwe.output
    output:
        multiext(config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.QC1', ".bim",".bed",".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
            --bfile {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC \
            --threads 20 \
            --make-bed \
            --allow-no-sex \
            --exclude {config[repository_path]}/resources/removeRegion.txt \
            --out {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC1"

rule remove_het:
    input:
        rules.exclude_high_ld_region.output
    output:
        multiext(config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.QC2', ".bim",".bed",".fam"),
        config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.FINAL.fail_het_imiss.sample'
    shell:
        """Rscript {config[repository_path]}/scripts/het.R {config[target_data]}/corrected/{config[processed_data_name]} {wildcards.target_data_name}.QC1 FINAL.fail_het_imiss.sample {config[repository_path]}/tools/plink
        {config[repository_path]}/tools/plink \
            --bfile {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC1 \
            --remove-fam {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL.fail_het_imiss.sample \
            --allow-no-sex \
            --threads 20 \
            --make-bed \
            --out {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC2
        """

rule exclude_related_samples:
    input:
        rules.remove_het.output
    threads: 20
    output:
        multiext(config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.QC3', ".bim",".bed",".fam"),
        config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}.FINAL.related.samples'
    shell:
        """Rscript {config[repository_path]}/scripts/relatedSamples.R {config[target_data]}/corrected/{config[processed_data_name]} {wildcards.target_data_name}.QC2 FINAL.related.samples {config[repository_path]}/tools/plink hg19
        {config[repository_path]}/tools/plink \
            --bfile {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC2 \
            --remove-fam {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL.related.samples\
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC3

        """

rule assessPopStratification:
    input:
        rules.exclude_related_samples.output
    output:
        multiext(config['target_data']+'/corrected/' + config['processed_data_name'] + '/{target_data_name}',".FINAL.bim",".FINAL.bed",".FINAL.fam")
    threads: 20
    shell:
        """
        Rscript {config[repository_path]}/scripts/convert_ids_and_calc_PCA.R {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC3 {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL {config[convert_rsid]} {config[repository_path]}/tools/plink {config[repository_path]}/tools/plink2 {config[repository_path]}/resources/removeRegion.txt
        """

rule plot:
    input:
        rules.assessPopStratification.output
    output:
        multiext(config['target_data']+'/corrected/'+config['processed_data_name']+'/{target_data_name}',".FINAL.sample_missingness.pdf",".FINAL.snp_missingness.pdf",".FINAL.frq.pdf")
    shell:
        """
        Rscript {config[repository_path]}/scripts/plotMisFrq.R {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL
        rm {config[target_data]}/corrected/{config[processed_data_name]}/{wildcards.target_data_name}.QC*
        find {config[target_data]}/corrected/{config[processed_data_name]}/ -type f ! -name "*.FINAL.*" -exec rm -rf {{}} \\;
        rm -rf {config[target_data]}/raw/tmp
        """
