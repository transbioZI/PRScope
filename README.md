<!-- TABLE OF CONTENTS -->
## Table of Contents
<ol>
  <li>
    <a href="#getting-started">Getting Started</a>
    <ul>
      <li><a href="#installation--setup">Installation & Setup</a></li>
    </ul>
  </li>
  <li>
    <a href="#summary-statistics">Summary Statistics</a>
    <ul>
      <li><a href="#efo-ids">EFO IDs</a></li>
      <li><a href="#configuration-file">Configuration File</a></li>
    </ul>
  </li>
  <li><a href="#genetic-data">Genetic Data</a></li>
  <li><a href="#polygenic-risk-scores">Polygenic Risk Scores</a></li>
</ol>



<!-- GETTING STARTED -->
## Getting Started

This is a step-by-step guide on how to work with the gwas_base_utils repository.

### Installation & Setup

1. Clone the repo:-
   ```sh
   git clone https://github.com/transbioZI/PGScope.git
   ```
2. Install snakemake:
   ```sh
   conda create -c conda-forge -c bioconda -n snakemake snakemake python=3.12.1
   conda activate snakemake
   ```
3. Load R:
   ```sh
   module load R/4.2.3
   ```
4. Create Working directory e.g. `/Workspace/WorkDir`.
5. Copy [config.yaml](./config.yaml) and [run.sh](./run.sh) files to the working directory to edit them only locally.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>



<!-- SUMMARY STATISTICS -->
## Summary Statistics

The following is a detailed explanation of how to fetch the desired Summary Statistics from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/), extract its data and run it through preprocessing and quality control.


### EFO IDs
1. In your working directory create a blank `.txt` file
2. Go to the [Experimental Factor Ontology](https://www.ebi.ac.uk/ols4/ontologies/efo?tab=classes) website and search in the tree for the desired experiment.
3. You can either click on the desired experiment to copy its EFO ID, or you can simply hover on the experiment and write out the ID yourself to avoid restarting the tree search.

   **_NOTE:_** If you want all the experiments that are subgroups of one experiment it is enough to get the EFO ID from the parent group (the algorithm will automatically fetch all the subgroups).

4. Copy all obtained EFO IDs to your .txt file.

**_NOTE:_** If you are only working with a few summary statistics you can download their `.txt` files yourself from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/) to your working directory. You will need their paths in the following steps.

### Configuration File
Open the `config.yaml` file that you have copied to your local working directory and edit the paths as follows:
1. **_DO NOT MODIFY:_** `harmonised_list` is the URL that contains the harmonised list of all summary statistics.
2. `repository_path` is the path to the git repository you cloned during <a href="#installation-and-setup">Setup</a>.
3. `working_directory` is the path to the working directory you created during <a href="#installation-and-setup">Setup</a>.
4. **_OPTIONAL:_** `efoIds` is the name of the `.txt` file you created in your working directory for your fetched <a href="#efo-ids">EFO IDs</a>.
5. **_OPTIONAL:_**  `filter_study_output` is the name of a blank `.txt` file that you should create in your working directory so filtered summary statistics are stored in case you want the algorithm to fetch the summary statistics based on the EFO IDs.
6. `transform_study_input` is the same file name as in `filter_study_output`, if it was indicated, or the name of a file you created yourself, that contains the filtered summary statistics.
7. `transformed_path` is the path to an empty directory that you should create where the transformed summary statistics will be stored (e.g. `/Workspace/WorkDir/Transformed`). 

<!-- Eventuell hier .gif hochladen, um zu zeigen, wie man die EFO IDs findet (auf mein Windows User Profile: Workspace/tutorials/EFO_IDs.gif) -->

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>




<!-- Genetic Data -->
## Genetic Data

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>




<!-- Polygenic Risk Scores -->
## Polygenic Risk Scores


<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

