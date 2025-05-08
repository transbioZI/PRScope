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
            config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.FINAL.snp_missingness.pdf',
            config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.FINAL.sample_missingness.pdf',
            config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.FINAL.frq.pdf'],
            target_data_name = get_raw_target_data())

rule sort_plink:
    input:
        multiext(config['target_data']+'/raw/{target_data_name}','.bim', ".bed", ".fam")
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.sort','.bim', ".bed", ".fam")
    threads: 20
    params:
        prefix = config['target_data']+'/raw/{target_data_name}',
        prefix_o = config['target_data']+'/raw/tmp/{target_data_name}.sort'
    shell:
        "{config[repository_path]}/tools/plink --bfile {params.prefix} --make-bed --out {params.prefix_o} --threads 20"

rule create_map_file:
    input:
        rules.sort_plink.output
    threads: 20
    output:
        config['target_data'] + '/raw/tmp/{target_data_name}.tab.map'
    params:
        prefix = config['target_data']+'/raw/tmp/{target_data_name}.sort',
        prefix_o = config['target_data']+'/raw/tmp/{target_data_name}.tab'
    shell:
        "{config[repository_path]}/tools/plink --bfile {params.prefix} --recode tab --allow-extra-chr --out {params.prefix_o} --threads 20"

rule liftover:
    input:
        rules.create_map_file.output
    output:
        config['target_data'] + '/raw/lifted/{target_data_name}.lifted.map'
    conda: "../environment.yaml"
    params:
        chain = config['repository_path'] + '/resources/hg19ToHg38.over.chain.gz',
        lifted = config['target_data'] + '/raw/lifted/{target_data_name}.lifted',
        pedfile = config['target_data'] + '/raw/tmp/{target_data_name}.tab.ped'
    shell:
        "python2 {config[repository_path]}/tools/liftOverPlink/liftOverPlink.py -p {params.pedfile} -m {input} -c {params.chain} -o {params.lifted} -e {config[repository_path]}/tools/liftOver"

rule create_bim_bed_bam:
    input:
        rules.liftover.output
    threads: 20
    output:
        config['target_data'] + '/raw/lifted/{target_data_name}.lifted.bim',
        config['target_data'] + '/raw/lifted/{target_data_name}.lifted.bed',
        config['target_data'] + '/raw/lifted/{target_data_name}.lifted.fam'
    params:
        prefix = config['target_data'] + '/raw/lifted/{target_data_name}.lifted'
    shell:
        "{config[repository_path]}/tools/plink --file {params.prefix} --make-bed --allow-extra-chr --chr 1-22 --out {params.prefix} --threads 20"

# for more information : https://shicheng-guo.github.io/bioinformatics/2017/08/01/hapmap3

rule remove_samples:
    input:
        rules.create_bim_bed_bam.output
    output:
        multiext(config['target_data']+'/raw/tmp/{target_data_name}.filtered', '.bim', ".bed", ".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/lifted/{wildcards.target_data_name}.lifted\
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
            --missing \
            --allow-no-sex \
            --freq \
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
            --missing \
            --allow-no-sex \
            --freq \
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
            --missing \
            --allow-no-sex \
            --freq \
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
            --missing \
            --threads 20 \
            --allow-no-sex \
            --freq \
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
            --missing \
            --allow-no-sex \
            --freq \
            --threads 20 \
            --out {config[target_data]}/raw/tmp/{wildcards.target_data_name}.mafremoved"

rule remove_hwe:
    input:
        rules.remove_maf.output
    output:
        multiext(config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.QC', ".snplist",".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
        --bfile {config[target_data]}/raw/tmp/{wildcards.target_data_name}.mafremoved\
            --hwe 1e-6 \
            --missing \
            --freq \
            --write-snplist \
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC"

rule exclude_high_ld_region:
    input:
        rules.remove_hwe.output
    output:
        multiext(config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.QC1', ".bim",".bed",".fam")
    threads: 20
    shell:
        "{config[repository_path]}/tools/plink \
            --bfile {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC \
            --threads 20 \
            --make-bed \
            --allow-no-sex \
            --exclude {config[repository_path]}/resources/removeRegion.txt \
            --out {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC1"

rule remove_het:
    input:
        rules.exclude_high_ld_region.output
    output:
        multiext(config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.QC2', ".bim",".bed",".fam"),
        config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.FINAL.fail_het_imiss.sample'
    shell:
        """Rscript {config[repository_path]}/scripts/het.R {config[target_data]}/corrected_hg38/{config[processed_data_name]} {wildcards.target_data_name}.QC1 FINAL.fail_het_imiss.sample {config[repository_path]}/tools/plink
        {config[repository_path]}/tools/plink \
            --bfile {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC1 \
            --remove-fam {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL.fail_het_imiss.sample \
            --allow-no-sex \
            --threads 20 \
            --make-bed \
            --out {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC2
        """

rule exclude_related_samples:
    input:
        rules.remove_het.output
    threads: 20
    output:
        multiext(config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.QC3', ".bim",".bed",".fam"),
        config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}.FINAL.related.samples'
    shell:
        """Rscript {config[repository_path]}/scripts/relatedSamples.R {config[target_data]}/corrected_hg38/{config[processed_data_name]} {wildcards.target_data_name}.QC2 FINAL.related.samples {config[repository_path]}/tools/plink hg38
        {config[repository_path]}/tools/plink \
            --bfile {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC2 \
            --remove-fam {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL.related.samples\
            --make-bed \
            --allow-no-sex \
            --threads 20 \
            --out {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC3

        """

rule assessPopStratification:
    input:
        rules.exclude_related_samples.output
    output:
        multiext(config['target_data']+'/corrected_hg38/' + config['processed_data_name'] + '/{target_data_name}',".FINAL.bim",".FINAL.bed",".FINAL.fam")
    threads: 20
    shell:
        """
        Rscript {config[repository_path]}/scripts/convert_ids_and_calc_PCA.R {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC3 {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL {config[convert_rsid]} {config[repository_path]}/tools/plink {config[repository_path]}/tools/plink2 {config[repository_path]}/resources/removeRegion.txt
        """

rule plot:
    input:
        rules.assessPopStratification.output
    output:
        multiext(config['target_data']+'/corrected_hg38/'+config['processed_data_name']+'/{target_data_name}',".FINAL.sample_missingness.pdf",".FINAL.snp_missingness.pdf",".FINAL.frq.pdf")
    shell:
        """
        Rscript {config[repository_path]}/scripts/plotMisFrq.R {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.FINAL
        rm {config[target_data]}/corrected_hg38/{config[processed_data_name]}/{wildcards.target_data_name}.QC*
        find {config[target_data]}/corrected_hg38/{config[processed_data_name]}/ -type f ! -name "*.FINAL.*" -exec rm -rf {{}} \\;
        rm -rf {config[target_data]}/raw/tmp_hg19
        """
