#!/bin/bash

echo "=========================================\n\n
phyloSI-RakaiAgeGender: begin installation\n\n"

# Create new conda environment called "phyloSI-RakaiAgeGender"
if [ -d $HOME/anaconda3/envs/phyloSI-RakaiAgeGender ]; then
    echo "###############################################"
    echo -e "\nphyloSI-RakaiAgeGender conda environment is already present"
    echo -e "\nIf you wish you re-install please remove the conda environment first with:"
    echo -e "\tconda remove -n phyloSI-RakaiAgeGender --all -y"
    echo -e "\n\n###############################################"
    exit 1
else
    echo -e "\nCreating Conda environment: phyloSI-RakaiAgeGender"
    conda create -n phyloSI-RakaiAgeGender -y
    source activate phyloSI-RakaiAgeGender
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
fi

# Install base R
echo "\n\n=========================================\n\n
Installing base R\n\n"
conda install -c conda-forge r r-base r-devtools

# Install R dependencies
echo "\n\n=========================================\n\n
Installing R dependencies\n\n"
R -e 'options(unzip = "internal");install.packages(c("data.table","dplyr","ggplot2","rstan","stringi","extraDistr","plyr","coda","LaplacesDemon","mvtnorm","R.utils","prodlim","ggpubr","broom","matrixStats","ggnewscale","doParallel","foreach","jcolors","lognorm","ggExtra","Hmisc","knitr","Matrix","scales","gridgraphics","cowplot","abind","ggrepel","truncnorm","invgamma","sna","lubridate","tidyverse","MASS","gridExtra","GGally","ggnetwork","binom","igraph","bh","bayesplot","loo","hexbin","purrr","viridis"),repos = "http://cran.us.r-project.org")'

echo "=========================================\n\n
phyloSI-RakaiAgeGender: completed installation.\n
For next steps see\n
https://github.com/MLGlobalHealth/phyloSI-RakaiAgeGender
"
