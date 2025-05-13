import os

plink = config["plink"]
plink2 = config["plink2"]

if plink is None:
    plink = config["repository"] + "/tools/plink"

if plink2 is None:
    plink2 = config["repository"] + "/tools/plink2"

rule all:
    input:
        config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/' + config['genotype_data_name'] + '.FINAL.snp_missingness.jpeg',
        config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/' + config['genotype_data_name'] + '.FINAL.sample_missingness.jpeg',
        config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/' + config['genotype_data_name'] + '.FINAL.frq.jpeg'

rule remove_snps_missingness:
    input:
        multiext(config['target_data']+'/{target_data_name}','.bim', ".bed", ".fam")
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.snpsremoved', '.bim', ".bed", ".fam")
    shell:
        "{plink} \
        --bfile {config[target_data]}/{wildcards.target_data_name}\
            --geno {config[shallow_genotype_missingness]}  \
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved"

rule remove_samples_missingness:
    input:
        rules.remove_snps_missingness.output
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.samplesremoved', '.bim', ".bed", ".fam")
    shell:
        "{plink} \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved\
            --mind {config[shallow_sample_missingness]} \
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved"

rule remove_snps_missingness_deeply:
    input:
        rules.remove_samples_missingness.output
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.snpsremoved2', '.bim', ".bed", ".fam")
    shell:
        "{plink} \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved\
            --geno {config[deep_genotype_missingness]} \
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved2"

rule remove_samples_missingness_deeply:
    input:
        rules.remove_snps_missingness_deeply.output
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.samplesremoved2', '.bim', ".bed", ".fam")
    shell:
        "{plink} \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved2 \
            --mind {config[deep_sample_missingness]} \
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved2"

rule remove_maf:
    input:
        rules.remove_samples_missingness_deeply.output
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.mafremoved', '.bim', ".bed", ".fam")
    shell:
        "{plink} \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved2 \
            --maf {config[maf]} \
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}//tmp/{wildcards.target_data_name}.mafremoved"

rule remove_hwe:
    input:
        rules.remove_maf.output
    output:
        multiext(config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/{target_data_name}.QC',".fam",".bim",".bed")
    shell:
        "{plink} \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.mafremoved\
            --hwe {config[hwe]} \
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC"

rule exclude_high_ld_region:
    input:
        rules.remove_hwe.output
    output:
        multiext(config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/{target_data_name}.QC1', ".bim",".bed",".fam")
    shell:
        "{plink} \
            --bfile {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC \
            --make-bed \
            --allow-no-sex \
            --exclude {config[repository]}/resources/removeRegion.txt \
            --out {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC1"

rule remove_het:
    input:
        rules.exclude_high_ld_region.output
    output:
        multiext(config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/{target_data_name}.QC2', ".bim",".bed",".fam"),
        config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/{target_data_name}.FINAL.fail_het_imiss.sample'
    shell:
        """Rscript {config[repository]}/scripts/het.R {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]} {wildcards.target_data_name}.QC1 FINAL.fail_het_imiss.sample {plink} {config[heterozygosity]}
        {plink} \
            --bfile {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC1 \
            --remove-fam {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL.fail_het_imiss.sample \
            --allow-no-sex \
            --make-bed \
            --out {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC2
        """

rule exclude_related_samples:
    input:
        rules.remove_het.output
    output:
        multiext(config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/{target_data_name}.QC3', ".bim",".bed",".fam"),
        config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/{target_data_name}.FINAL.related.samples'
    shell:
        """Rscript {config[repository]}/scripts/relatedSamples.R {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]} {wildcards.target_data_name}.QC2 FINAL.related.samples {plink} hg19 {config[relatedness]}
        {plink} \
            --bfile {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC2 \
            --remove-fam {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL.related.samples\
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC3
        """

rule assessPopStratification:
    input:
        rules.exclude_related_samples.output
    output:
        multiext(config['output_path']+'/corrected_hg19/' + config['processed_data_directory_name'] + '/{target_data_name}',".FINAL.bim",".FINAL.bed",".FINAL.fam")
    shell:
        """
        Rscript {config[repository]}/scripts/convert_ids_and_calc_PCA.R {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC3 {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL {config[convert_rsid]} {plink} {plink2} {config[repository]}/resources/removeRegion.txt
        """

rule plot:
    input:
        rules.assessPopStratification.output
    output:
        multiext(config['output_path']+'/corrected_hg19/'+config['processed_data_directory_name']+'/{target_data_name}',".FINAL.sample_missingness.jpeg",".FINAL.snp_missingness.jpeg",".FINAL.frq.jpeg")
    shell:
        """
        Rscript {config[repository]}/scripts/plotMisFrq.R {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL
        rm {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC*
        find {config[output_path]}/corrected_hg19/{config[processed_data_directory_name]}/ -type f ! -name "*.FINAL.*" -exec rm -rf {{}} \\;
        rm -rf {config[output_path]}/tmp
        """
