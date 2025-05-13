#!/bin/bash
# install snakemake : conda create -c conda-forge -c bioconda -n snakemake snakemake python=3.12.1
#module load R
#conda activate snakemake

REPOSITORY_PATH=/data/projects/on_going/PRScope

snakemake -s $REPOSITORY_PATH/main/snakefiles/find_sumstats.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs 1 --config repository=$REPOSITORY_PATH/main
snakemake -s $REPOSITORY_PATH/main/snakefiles/qc_sumstats.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs 5 --config repository=$REPOSITORY_PATH/main
snakemake -s $REPOSITORY_PATH/main/snakefiles/qc_genotype_with_liftover.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --use-conda --jobs 5 --conda-frontend conda --config repository=$REPOSITORY_PATH/main
snakemake -s $REPOSITORY_PATH/main/snakefiles/qc_genotype.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs 5 --config repository=$REPOSITORY_PATH/main
snakemake -s $REPOSITORY_PATH/main/snakefiles/prs_calculation_prsice.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs 5 -k --config repository=$REPOSITORY_PATH/main
