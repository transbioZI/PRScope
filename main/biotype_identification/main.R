rm(list=ls())
gc()
library(gplots)
library(cluster) 
library(fmsb)
library(ggplot2)
library(ggpubr)
library(stringr)
library(networkD3)
library(grid)
library(gridExtra)
library(tidytext)
library(EMCluster)
library(reshape)
library(caret)
library(ranger)
library(ggradar)
library(wesanderson)
library(gwasrapidd)
library(bigutilsr)
library(doParallel)
library(data.table)
library(plyr)

healthy_color = "#8B5742"
high_neuroticism_color = "#4876FF" 
low_neuroticism_color = "#27408B"
inverse_color = "#FA88BF"
group.colors <- c('Biotype 1' = low_neuroticism_color, 'Biotype 2'= high_neuroticism_color, Healthy =healthy_color, Inverse = inverse_color)

repository_path = "/data/projects/on_going/PRScope"
setwd(repository_path)
experiment_path = './results/ML' 
num_it = 5 # repeat for ML analysis
repeat_ = 5 # repeat for random forest models
features_sizes = c(1,3,6) # feature sizes
features_sizes = sort(features_sizes)

DISCOVERY_DATA = DiscoveryRES     # DELETE
DISCOVERY_GROUP = DiscoveryGroup  # DELETE
DISCOVERY_IDs = DiscoveryMeta[,2] # DELETE

VALIDATION_DATA = ValidationRES     # DELETE
VALIDATION_GROUP = ValidationGroup  # DELETE
VALIDATION_IDs = ValidationMeta[,1]  # DELETE

source("./biotype_identification/helper_scripts/helper_functions.R")
source("./biotype_identification/helper_scripts/ML_functions.R")
source("./biotype_identification/main_ml_code.R")
source("./biotype_identification/helper_scripts/Plot_helper.R")
source("./biotype_identification/helper_scripts/readResults.R")
source("./biotype_identification/plotting_or_pvalue_cluster.R")
