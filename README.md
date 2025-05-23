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

## Setup

### 1. Setting Up

#### a. Cloning the Repository

```bash
git clone https://github.com/transbioZI/PRScope
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

4. In `input/reference/`:

   - Download the required files from the provided link
   - Place them into this folder

   The folder should now include:

   - `ldpred2_ref/`
   - `eur_hg38.phase3.bed`
   - `eur_hg38.phase3.bim`
   - `eur_hg38.phase3.fam`
   - `eur_hg38.phase3.frq`

---

## Software Requirements

PRScope requires the following:

- Conda
- Python
- R
- Snakemake
- PLINK
- PRSice

Only Conda must be installed manually. All other dependencies are managed via the Conda environment.

---

## Installation Instructions

### 1. Install Conda

[Conda Installation Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

### 2. Create the Conda Environment

```bash
cd PRScope
conda env create -f main/environment.yaml
```

- Environment name: `prscope`
- May take up to an hour
- Wait for installation to complete

### 3. Activate the Environment

```bash
conda activate prscope
```

---

## Running PRScope

```bash
cd PRScope
./run.sh
```

This command initiates the PRScope pipeline.
