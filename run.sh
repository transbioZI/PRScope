#!/bin/bash
# install snakemake : conda create -c conda-forge -c bioconda -n snakemake snakemake python=3.12.1
#conda activate snakemake

REPOSITORY_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NUMBER_OF_THREADS=5
cd $REPOSITORY_PATH

snakemake -s $REPOSITORY_PATH/main/snakefiles/find_sumstats.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs $NUMBER_OF_THREADS --config repository=$REPOSITORY_PATH/main --use-conda --conda-frontend conda
snakemake -s $REPOSITORY_PATH/main/snakefiles/qc_sumstats.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs $NUMBER_OF_THREADS --config repository=$REPOSITORY_PATH/main --use-conda --conda-frontend conda
snakemake -s $REPOSITORY_PATH/main/snakefiles/qc_genotype_with_liftover.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --use-conda --jobs $NUMBER_OF_THREADS --conda-frontend conda --config repository=$REPOSITORY_PATH/main
snakemake -s $REPOSITORY_PATH/main/snakefiles/qc_genotype.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs 5 --config repository=$REPOSITORY_PATH/main --use-conda --conda-frontend conda
snakemake -s $REPOSITORY_PATH/main/snakefiles/prs_calculation_prsice.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs $NUMBER_OF_THREADS -k --config repository=$REPOSITORY_PATH/main --use-conda --conda-frontend conda
#snakemake -s $REPOSITORY_PATH/main/snakefiles/prs_calculation_ldpred.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs $NUMBER_OF_THREADS -k --config repository=$REPOSITORY_PATH/main --use-conda --conda-frontend conda
#snakemake -s $REPOSITORY_PATH/main/snakefiles/ldsc_heritability_calculation.snakefile --configfile $REPOSITORY_PATH/config/config_vanilla.yaml --jobs $NUMBER_OF_THREADS -k --config repository=$REPOSITORY_PATH/main --use-conda --conda-frontend conda
