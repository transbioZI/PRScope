rule all:
    input:
        config['output_path_gwas_search']+ "/" + config['output_name_gwas_search'] + ".txt"

rule download_harmonised_list:
    output:
        config['output_path_gwas_search']+ "/harmonised_list.txt"
    shell:
        """
        wget -O {output} {config[harmonised_list]}
        """

rule find_efo_studies:
    input:
        harm = rules.download_harmonised_list.output
    output:
        config['output_path_gwas_search']+ "/"+ config['output_name_gwas_search'] + ".txt"
    shell:
        """
        Rscript {config[repository]}/scripts/filter_studies_by_efo.R {input.harm} {config[efo_ids]} {config[output_path_gwas_search]} {config[output_name_gwas_search]} {config[sample_size]} {config[number_of_snps]} {config[publication_date]} {config[population]}
        """
