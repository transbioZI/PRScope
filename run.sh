#!/bin/bash
# install snakemake : conda create -c conda-forge -c bioconda -n snakemake snakemake python=3.12.1
#module load R
#conda activate snakemake

REPOSITORY_PATH=/zi/home/ersoy.kocak/Desktop/Projects/PRSCalculator
SCRIPT_DIR="$(pwd)"

touch empty_remove_samples.txt

snakemake -s $REPOSITORY_PATH/snakefiles/find_sumstats.snakefile --configfile config.yaml --jobs 1
snakemake -s $REPOSITORY_PATH/snakefiles/qc_sumstats.snakefile --configfile config.yaml --jobs 10 --rerun-incomplete
snakemake -s $REPOSITORY_PATH/snakefiles/qc_genotype_with_liftover.snakefile --configfile config.yaml --use-conda --jobs 20 --conda-frontend conda
snakemake -s $REPOSITORY_PATH/snakefiles/qc_genotype.snakefile --configfile config.yaml --use-conda --jobs 20 --conda-frontend conda
snakemake -s $REPOSITORY_PATH/snakefiles/prs_calculation_prsice.snakefile --configfile config.yaml --jobs 10 -k

Rscript $REPOSITORY_PATH/scripts/create_prs_datatable_prsice.R $SCRIPT_DIR/test_prs 0 1000 $SCRIPT_DIR/test_studies.txt name_0.txt
Rscript $REPOSITORY_PATH/scripts/create_prs_datatable_prsice.R $SCRIPT_DIR/test_prs 100 1000 $SCRIPT_DIR/test_studies.txt name_100.txt
Rscript $REPOSITORY_PATH/scripts/create_prs_datatable_prsice.R $SCRIPT_DIR/test_prs 250 1000 $SCRIPT_DIR/test_studies.txt name_250.txt
Rscript $REPOSITORY_PATH/scripts/create_prs_datatable_prsice.R $SCRIPT_DIR/test_prs 500 1000 $SCRIPT_DIR/test_studies.txt name_500.txt
Rscript $REPOSITORY_PATH/scripts/create_prs_datatable_prsice.R $SCRIPT_DIR/test_prs 1000 1000 $SCRIPT_DIR/test_studies.txt name_1000.txt
