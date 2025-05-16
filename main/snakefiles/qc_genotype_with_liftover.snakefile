import os

rule all:
    input:
        config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/' + config['genotype_data_name'] + '.FINAL.snp_missingness.jpeg',
        config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/' + config['genotype_data_name'] + '.FINAL.sample_missingness.jpeg',
        config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/' + config['genotype_data_name'] + '.FINAL.frq.jpeg'

rule sort_plink:
    input:
        multiext(config['target_data']+'/{target_data_name}','.bim', ".bed", ".fam")
    output:
        multiext(config['output_path'] + '/tmp/{target_data_name}.sort','.bim', ".bed", ".fam")
    conda: "../environment.yaml"
    params:
        prefix = config['target_data']+'/{target_data_name}',
        prefix_o = config['output_path'] + '/tmp/{target_data_name}.sort'
    shell:
        "plink --bfile {params.prefix} --make-bed --out {params.prefix_o}"

rule create_map_file:
    input:
        rules.sort_plink.output
    output:
        config['output_path'] + '/tmp/{target_data_name}.tab.map'
    conda: "../environment.yaml"
    params:
        prefix = config['output_path'] + '/tmp/{target_data_name}.sort',
        prefix_o = config['output_path'] + '/tmp/{target_data_name}.tab'
    shell:
        "plink --bfile {params.prefix} --recode tab --allow-extra-chr --out {params.prefix_o}"

rule liftover:
    input:
        rules.create_map_file.output
    output:
        config['output_path'] + '/lifted/{target_data_name}.lifted.map'
    conda: "../environment_for_liftover.yaml"
    params:
        chain = config['repository'] + '/resources/hg19ToHg38.over.chain.gz',
        lifted = config['output_path'] + '/lifted/{target_data_name}.lifted',
        pedfile = config['output_path'] + '/tmp/{target_data_name}.tab.ped'
    shell:
        "python2 {config[repository]}/tools/liftOverPlink/liftOverPlink.py -p {params.pedfile} -m {input} -c {params.chain} -o {params.lifted} -e {config[repository]}/tools/liftOver"

rule create_bim_bed_bam:
    input:
        rules.liftover.output
    output:
        config['output_path'] + '/lifted/{target_data_name}.lifted.bim',
        config['output_path'] + '/lifted/{target_data_name}.lifted.bed',
        config['output_path'] + '/lifted/{target_data_name}.lifted.fam'
    conda: "../environment.yaml"
    params:
        prefix = config['output_path'] + '/lifted/{target_data_name}.lifted'
    shell:
        "plink --file {params.prefix} --make-bed --allow-extra-chr --chr 1-22 --out {params.prefix}"

# for more information : https://shicheng-guo.github.io/bioinformatics/2017/08/01/hapmap3

rule remove_snps_missingness:
    input:
        rules.create_bim_bed_bam.output
    conda: "../environment.yaml"
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.snpsremoved', '.bim', ".bed", ".fam")
    shell:
        "plink \
        --bfile  {config[output_path]}/lifted/{wildcards.target_data_name}.lifted\
            --geno {config[shallow_genotype_missingness]} \
            --make-bed \
            --missing \
            --allow-no-sex \
            --freq \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved"

rule remove_samples_missingness:
    input:
        rules.remove_snps_missingness.output
    conda: "../environment.yaml"
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.samplesremoved', '.bim', ".bed", ".fam")
    shell:
        "plink \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved\
            --mind {config[shallow_sample_missingness]} \
            --make-bed \
            --missing \
            --allow-no-sex \
            --freq \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved"

rule remove_snps_missingness_deeply:
    input:
        rules.remove_samples_missingness.output
    conda: "../environment.yaml"
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.snpsremoved2', '.bim', ".bed", ".fam")
    shell:
        "plink \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved\
            --geno {config[deep_genotype_missingness]} \
            --make-bed \
            --missing \
            --allow-no-sex \
            --freq \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved2"

rule remove_samples_missingness_deeply:
    input:
        rules.remove_snps_missingness_deeply.output
    conda: "../environment.yaml"
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.samplesremoved2', '.bim', ".bed", ".fam")
    shell:
        "plink \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.snpsremoved2 \
            --mind {config[deep_sample_missingness]} \
            --make-bed \
            --missing \
            --allow-no-sex \
            --freq \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved2"

rule remove_maf:
    input:
        rules.remove_samples_missingness_deeply.output
    output:
        multiext(config['output_path']+'/tmp/{target_data_name}.mafremoved', '.bim', ".bed", ".fam")
    conda: "../environment.yaml"
    shell:
        "plink \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.samplesremoved2 \
            --maf {config[maf]} \
            --make-bed \
            --missing \
            --allow-no-sex \
            --freq \
            --out {config[output_path]}/tmp/{wildcards.target_data_name}.mafremoved"

rule remove_hwe:
    input:
        rules.remove_maf.output
    output:
        multiext(config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/{target_data_name}.QC', ".snplist",".fam")
    conda: "../environment.yaml"
    shell:
        "plink \
        --bfile {config[output_path]}/tmp/{wildcards.target_data_name}.mafremoved\
            --hwe {config[hwe]} \
            --missing \
            --freq \
            --write-snplist \
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC"

rule exclude_high_ld_region:
    input:
        rules.remove_hwe.output
    output:
        multiext(config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/{target_data_name}.QC1', ".bim",".bed",".fam")
    conda: "../environment.yaml"
    shell:
        "plink \
            --bfile {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC \
            --make-bed \
            --allow-no-sex \
            --exclude {config[repository]}/resources/removeRegion.txt \
            --out {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC1"

rule remove_het:
    input:
        rules.exclude_high_ld_region.output
    output:
        multiext(config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/{target_data_name}.QC2', ".bim",".bed",".fam"),
        config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/{target_data_name}.FINAL.fail_het_imiss.sample'
    conda: "../environment.yaml"
    shell:
        """Rscript {config[repository]}/scripts/het.R {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]} {wildcards.target_data_name}.QC1 FINAL.fail_het_imiss.sample $(which plink) {config[heterozygosity]}
        plink \
            --bfile {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC1 \
            --remove-fam {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL.fail_het_imiss.sample \
            --allow-no-sex \
            --make-bed \
            --out {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC2
        """

rule exclude_related_samples:
    input:
        rules.remove_het.output
    output:
        multiext(config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/{target_data_name}.QC3', ".bim",".bed",".fam"),
        config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/{target_data_name}.FINAL.related.samples'
    conda: "../environment.yaml"
    shell:
        """Rscript {config[repository]}/scripts/relatedSamples.R {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]} {wildcards.target_data_name}.QC2 FINAL.related.samples $(which plink) hg38 {config[relatedness]}
        plink \
            --bfile {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC2 \
            --remove-fam {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL.related.samples\
            --make-bed \
            --allow-no-sex \
            --out {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC3

        """

rule assessPopStratification:
    input:
        rules.exclude_related_samples.output
    output:
        multiext(config['output_path']+'/corrected_hg38/' + config['processed_data_directory_name'] + '/{target_data_name}',".FINAL.bim",".FINAL.bed",".FINAL.fam")
    conda: "../environment.yaml"
    shell:
        """
        Rscript {config[repository]}/scripts/convert_ids_and_calc_PCA.R {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC3 {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL {config[convert_rsid]} $(which plink) $(which plink2) {config[repository]}/resources/removeRegion.txt
        """

rule plot:
    input:
        rules.assessPopStratification.output
    output:
        multiext(config['output_path']+'/corrected_hg38/'+config['processed_data_directory_name']+'/{target_data_name}',".FINAL.sample_missingness.jpeg",".FINAL.snp_missingness.jpeg",".FINAL.frq.jpeg")
    conda: "../environment.yaml"
    shell:
        """
        Rscript {config[repository]}/scripts/plotMisFrq.R {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.FINAL
        rm {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/{wildcards.target_data_name}.QC*
        find {config[output_path]}/corrected_hg38/{config[processed_data_directory_name]}/ -type f ! -name "*.FINAL.*" -exec rm -rf {{}} \\;
        rm -rf {config[output_path]}/tmp
        """
