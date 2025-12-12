import os
import pandas as pd

ldpred_path = config["ldpred_path"]

if ldpred_path is None:
    ldpred_path = config["repository"] + "/tools/ldpred2"

rule all:
    input:
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.' + config['imputation_mode'] + '.nomiss.rds'

rule convert_plink_file_rds:
    input:
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.bed'
    conda: "../environment.yaml"
    output:
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.rds',
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.bk'
    shell:
        """
        Rscript {ldpred_path}/createBackingFile.R --file-input {input} --file-output {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.rds
        """

rule impute:
    input: rules.convert_plink_file_rds.output
    output:
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.' + config['imputation_mode'] + '.nomiss.rds',
        config['target_data_path_ldpred'] + "/" + config['target_data_prefix_ldpred'] + '.' + config['imputation_mode'] + '.nomiss.bk'
    conda: "../environment.yaml"
    shell:
        """
        cp {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.rds {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.{config[imputation_mode]}.nomiss.rds
        cp {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.bk {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.{config[imputation_mode]}.nomiss.bk
        Rscript {ldpred_path}/imputeGenotypes.R --impute-simple {config[imputation_mode]} --geno-file-rds {config[target_data_path_ldpred]}/{config[target_data_prefix_ldpred]}.{config[imputation_mode]}.nomiss.rds
        """
