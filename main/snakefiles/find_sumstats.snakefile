rule all:
    input:
        config['output_path_gwas_search']+ "/" + config['output_name_gwas_search'] + ".txt.tsv",
        config['output_path_gwas_search']+ "/" + config['output_name_gwas_search'] + "_less_cases.txt.tsv",
        config['output_path_gwas_search']+ "/" + config['output_name_gwas_search'] + "_less_control.txt.tsv",
        config['output_path_gwas_search']+ "/" + config['output_name_gwas_search'] + "_not_harmonized.txt.tsv",
        config['output_path_gwas_search']+ "/" + config['output_name_gwas_search'] + "_other_ancestry_gwas.txt.tsv"

rule download_harmonised_list:
    output:
        config['output_path_gwas_search']+ "/harmonised_list.txt"
    shell:
        """
        wget -O {output} {config[harmonised_list]}
        """

rule find_efo_studies:
    input:
        rules.download_harmonised_list.output
    conda: "../environment.yaml"
    output:
        config['output_path_gwas_search']+ "/"+ config['output_name_gwas_search'] + ".txt"
    shell:
        """
        Rscript {config[repository]}/scripts/filter_studies_by_efo.R {input} {config[efo_ids]} {config[output_path_gwas_search]} {config[output_name_gwas_search]} {config[sample_size]} {config[number_of_snps]} {config[publication_date]} {config[population]} {config[pubmed_ids_to_exclude]} {config[number_of_cases]}
        """

rule create_studies_table:
    input:
        rules.find_efo_studies.output
    conda: "../environment.yaml"
    output:
        config['output_path_gwas_search']+ "/"+ config['output_name_gwas_search'] + ".txt.tsv"
    shell:
        """
        Rscript {config[repository]}/scripts/studyList_GWAS_Catalog.R {config[output_path_gwas_search]}/{config[output_name_gwas_search]}.txt
        """

rule create_studies_table_less_cases:
    input:
        rules.find_efo_studies.output
    conda: "../environment.yaml"
    output:
        config['output_path_gwas_search']+ "/"+ config['output_name_gwas_search'] + "_less_cases.txt.tsv"
    shell:
        """
        Rscript {config[repository]}/scripts/studyList_GWAS_Catalog.R {config[output_path_gwas_search]}/{config[output_name_gwas_search]}_less_cases.txt
        """

rule create_studies_table_less_control:
    input:
        rules.find_efo_studies.output
    conda: "../environment.yaml"
    output:
        config['output_path_gwas_search']+ "/"+ config['output_name_gwas_search'] + "_less_control.txt.tsv"
    shell:
        """
        Rscript {config[repository]}/scripts/studyList_GWAS_Catalog.R {config[output_path_gwas_search]}/{config[output_name_gwas_search]}_less_control.txt
        """

rule create_studies_table_other_ancestry_gwas:
    input:
        rules.find_efo_studies.output
    conda: "../environment.yaml"
    output:
        config['output_path_gwas_search']+ "/"+ config['output_name_gwas_search'] + "_other_ancestry_gwas.txt.tsv"
    shell:
        """
        Rscript {config[repository]}/scripts/studyList_GWAS_Catalog.R {config[output_path_gwas_search]}/{config[output_name_gwas_search]}_other_ancestry_gwas.txt
        """

rule create_studies_table_not_harmonized:
    input:
        rules.find_efo_studies.output
    conda: "../environment.yaml"
    output:
        config['output_path_gwas_search']+ "/"+ config['output_name_gwas_search'] + "_not_harmonized.txt.tsv"
    shell:
        """
        Rscript {config[repository]}/scripts/studyList_GWAS_Catalog.R {config[output_path_gwas_search]}/{config[output_name_gwas_search]}_not_harmonized.txt
        """
