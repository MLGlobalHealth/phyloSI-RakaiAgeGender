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
    # following https://www.biostars.org/p/450316/ here:
    conda create -n phyloSI-RakaiAgeGender -y
    source activate phyloSI-RakaiAgeGender
    # https://stackoverflow.com/questions/53014306/error-15-initializing-libiomp5-dylib-but-found-libiomp5-dylib-already-initial
    # The Intel MKL functions (e.g. FFT, LAPACK, BLAS) are threaded with the OpenMP technology.
    # But on macOS you do not need MKL, because the Accelerate Framework comes with its own
    # optimization algorithms and already uses OpenMP
    # You should install all packages without MKL support:
    conda install -c anaconda nomkl -y
    conda install -c r r r-essentials r-textshaping r-ragg -y
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

    # Installation for some R libraries fails as the libgfortran is not
    # found at the Linker path, and I don t know what -L variable I need to change
    # for now just copy the libraries as follows
    cp $CONDA_PREFIX/lib/libgfortran.* $CONDA_PREFIX/lib/gcc/arm64-apple-darwin20.0.0/11.0.1
    conda install -c anaconda fribidi -y

    # test gfortran
    # https://fortran-lang.org/en/learn/quickstart/hello_world/
    # gfortran hello.f90 -o hello

    # Install Rcpp-Eigen
    # address error on MacOS M1 error: no member named 'Rlog1p' in namespace 'std';
    # https://github.com/conda-forge/r-base-feedstock/issues/163
    echo "\n\n=========================================\n\n
    Installing R RcppEigen\n\n"
    export PKG_CPPFLAGS="-DHAVE_WORKING_LOG1P"
    R -e 'options(unzip = "internal");
r <- getOption("repos");
r["CRAN"] <- "https://cloud.r-project.org";
options(repos=r);
if(!require(RcppEigen)){ install.packages("RcppEigen") };
'
    # Install cmdstanr
    echo "\n\n=========================================\n\n
    Installing R cmdstanr\n\n"

    R -e 'options(unzip = "internal");
r <- getOption("repos");
r["CRAN"] <- "https://cloud.r-project.org";
options(repos=r);
if(!require(pillar)){ install.packages("pillar", version="1.8.1") };
if(!require(distributional)){ install.packages("distributional") };
if(!require(cmdstanr)){ install.packages("cmdstanr",repos = c("https://mc-stan.org/r-packages/",getOption("repos"))) };
library("cmdstanr");
cmdstanr::check_cmdstan_toolchain();
cmdstanr::install_cmdstan();
'
    # Install all else
    R -e 'options(unzip = "internal");
r <- getOption("repos");
r["CRAN"] <- "https://cloud.r-project.org";
options(repos=r);
pkgs <- c("dplyr","ggplot2","stringi","extraDistr","plyr","coda","LaplacesDemon","mvtnorm","R.utils","prodlim","ggpubr","broom","Hmisc","matrixStats","ggnewscale","doPara\
llel","foreach","lognorm","ggExtra","knitr","Matrix","scales","gridGraphics","cowplot","abind","ggrepel","truncnorm","invgamma","sna","lubridate","tidyverse","M\
ASS","gridExtra","GGally","ggnetwork","binom","igraph","BH","bayesplot","loo","hexbin","purrr","viridis", "patchwork");
for(pkg in pkgs){ if(!require(pkg, character.only = TRUE)){ install.packages(pkg) } };
'  

    echo "=========================================\n\n
phyloSI-RakaiAgeGender: completed installation.\n
For next steps see\n
https://github.com/MLGlobalHealth/phyloSI-RakaiAgeGender
"
fi





