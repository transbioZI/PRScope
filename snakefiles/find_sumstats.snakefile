rule all:
    input:
        config['working_directory']+ "/" + config['filter_study_output']

rule download_harmonised_list:
    output:
        config['working_directory']+ "/harmonised_list.txt"
    shell:
        """
        wget -O {output} {config[harmonised_list]}
        """

rule find_efo_studies:
    input:
        harm = rules.download_harmonised_list.output
    output:
        config['working_directory']+ "/"+ config['filter_study_output']
    params:
        efoIds = config['working_directory']+"/" + config['efoIds']
    shell:
        """
        Rscript {config[repository_path]}/scripts/filter_studies_by_efo.R {params.efoIds} {output} {input.harm}
        """
