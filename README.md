# PRScope

## Outline

**PRScope** automatically generates all the polygenic scores (PGS) associated with selected ontology IDs (e.g., Experimental Factor Ontology, EFO) for a given genotype dataset.

### Inputs

- Ontology IDs defining the PGS of interest
- Genotype data to calculate the PGS

### Output

  - A dataset (multi-PGS matrix) containing:
  - Subject IDs from the genotype data
  - Values of all selected PGS

The setup is optimized for a minimal-effort "vanilla use case" but supports advanced configurations.

---

## System Requirements

### Hardware Requirements

PRScope requires only a standard computer with enough RAM (4GB) to support the in-memory operations, but for best performance, we suggest a computer with higher specifications:

RAM: 16+ GB <br>
CPU: 4+ cores, 3.3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores each 3.3 GHz) and internet of speed 100 Mbps.

### Software Requirements

PRScope requires the following:

- conda
- Python
- R
- Snakemake
- PLINK
- PRSice

Only conda must be installed manually. All other dependencies are managed via the Conda environment.

PRScope has been tested on the Ubuntu 22.04.5 and requires a Linux system.

## Setup

### 1. Setting Up

#### a. Cloning the Repository
```bash
git clone https://github.com/transbioZI/PRScope
```
#### b. Download the required files from the provided link

Place them into this folder `input/reference/`:

- The folder should now include:

   - `ldpred2_ref/`
   - `eur_hg38.phase3.bed`
   - `eur_hg38.phase3.bim`
   - `eur_hg38.phase3.fam`
   - `eur_hg38.phase3.frq`

## Installation Instructions

### 1. Install Conda

[Conda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

### 2. Create the Conda Environment

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake python=3.12.1
```
- Environment name: `snakemake`
- Wait for installation to complete (~15 minutes)

### 3. Activate the Environment

```bash
conda activate snakemake
```
---

## Running PRScope (see Demo below)

```bash
cd PRScope
./run.sh
```

This command initiates the PRScope pipeline.

- Wait for installation to complete conda environment
- May take up to an hour

## The pipeline tested with :

conda version : 23.1.0     <br>
snakemake version : 8.4.8       <br>
python : 3.12.1               <br>

### R-sessionInfo()

R version 4.4.1 (2024-06-14)    <br>
Platform: x86_64-pc-linux-gnu     <br>
Running under: Ubuntu 22.04.5 LTS

Matrix products: default   <br>
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3   <br>
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so; LAPACK version 3.10.0

locale: <br>LC_CTYPE=C.UTF-8, LC_NUMERIC=C, LC_TIME=C.UTF-8, LC_COLLATE=C.UTF-8, LC_MONETARY=C.UTF-8, LC_MESSAGES=C.UTF-8, LC_PAPER=C.UTF-8, LC_NAME=C, LC_ADDRESS=C, LC_TELEPHONE=C, LC_MEASUREMENT=C.UTF-8, LC_IDENTIFICATION=C

time zone: Europe/Berlin   <br>
tzcode source: system (glibc)

attached base packages:      <br>
parallel, stats, graphics, grDevices, utils, datasets, methods, base

other attached packages:     <br>
reshape2_1.4.4, cluster_2.1.7, xgboost_1.7.8.1, bigutilsr_0.3.4, reshape_0.8.9, ggpubr_0.6.0,doParallel_1.0.17, iterators_1.0.14, foreach_1.5.2, glmnet_4.1-8, Matrix_1.7-1, lubridate_1.9.4, forcats_1.0.0, purrr_1.0.2, readr_2.1.5, tidyr_1.3.1, tibble_3.2.1, tidyverse_2.0.0, data.table_1.16.4, dplyr_1.1.4, gwasrapidd_0.99.17, caret_7.0-1, lattice_0.22-6, ranger_0.17.0, stringr_1.5.1, ggplot2_3.5.1, fmsb_0.7.6, optparse_1.7.5, tidyselect_1.2.1, timeDate_4041.110, bigassertr_0.1.6, pROC_1.18.5, digest_0.6.37, rpart_4.1.23,timechange_0.3.0, lifecycle_1.0.4, survival_3.7-0, magrittr_2.0.3, compiler_4.4.1, rlang_1.1.4, tools_4.4.1, utf8_1.2.4, ggsignif_0.6.4, plyr_1.8.9, abind_1.4-8, withr_3.0.2, nnet_7.3-19, grid_4.4.1, stats4_4.4.1, fansi_1.0.6, colorspace_2.1-1, future_1.34.0, globals_0.16.3, scales_1.3.0, MASS_7.3-61, cli_3.6.3, generics_0.1.3, RSpectra_0.16-2, rstudioapi_0.17.1, future.apply_1.11.3, tzdb_0.4.0, getopt_1.20.4, splines_4.4.1, vctrs_0.6.5, hardhat_1.4.0, jsonlite_1.8.9, carData_3.0-5, car_3.1-3, hms_1.1.3, rstatix_0.7.2, Formula_1.2-5, listenv_0.9.1, gower_1.0.1, recipes_1.1.0, glue_1.8.0,parallelly_1.40.1, codetools_0.2-20, stringi_1.8.4, gtable_0.3.6, shape_1.4.6.1, munsell_0.5.1, pillar_1.9.0, ipred_0.9-15, lava_1.8.0, R6_2.5.1, backports_1.5.0, broom_1.0.7, class_7.3-22, Rcpp_1.0.13-1, nlme_3.1-166, prodlim_2024.06.25, ModelMetrics_1.2.2.2, pkgconfig_2.0.3

# Demo

#### a. Navigate to the repository path

```bash
cd PRScope
```

#### b. Folder Structure

- `config/` – For advanced parameter customization
- `input/` – The only folder requiring user modifications
- `main/` – Contains the main pipeline
- `output/` – Will contain output after pipeline execution

#### c. Preparing the Input

1. Navigate to `input/`
2. Edit `efo_ids.txt`:
   - Default content: `EFO_0003898`
   - Replace with your own EFO IDs as needed

3. In `input/genotype/`, you’ll find:

   - `EUR.bed`
   - `EUR.bim`
   - `EUR.fam`

   (Replace with your genotype data if desired)

4. Running PRScope

```bash
cd PRScope
./run.sh
```
