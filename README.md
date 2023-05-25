# phyloSI-RakaiAgeGender
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![medRxiv link](https://img.shields.io/badge/medRxiv-link%20to%20paper-blue)](https://doi.org/10.1101/2023.03.16.23287351)

**Welcome!** This repository contains the code and data for the analyses presented in the paper *Growing gender inequity in HIV infection in Africa: sources and policy implications* by M Monod, A Brizzi, R Galiwango, R Ssekubugu, Y Chen, X Xi et al.

## Table of Contentes
- [License](#license)
- [Warrenty](#warranty)
- [Citation](#cite)
- [Acknowledgements](#acknowledgements)
- [Funding](#funding)
- [Quick Start](#quick-start)
  - [System Requirements](#system-requirements)
  - [Installation](#installation)
  - [Reproducing our Analyses](#reproducing-our-analyses)
    - [RCCS Surveillance Analyses](#rccs-surveillance-analyses)
    - [Phylogenetic Analyses](#phylogenetic-analyses)
    - [Age-specific HIV incidence rates](#age-specific-hiv-incidence-rates)
    - [Transmission flows analysis](#transmission-flows-analysis)

## License
The code and data in this repository are licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) by Imperial College London. Copyright Imperial College London 2022. 

## Warranty

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.

## Citation
Please cite this work as 

Mélodie Monod, Andrea Brizzi, Ronald M Galiwango, Robert Ssekubugu, Yu Chen, Xiaoyue Xi, Edward Nelson Kankaka, Victor Ssempijja, Lucie Abeler-Dörner, Adam Akullian, Alexandra Blenkinsop, David Bonsall, Larry William Chang, Shozen Dan, Christophe Fraser, Tanya Golubchik, Ronald H Gray, Matthew Hall, Jade C. Jackson, Godfrey Kigozi, Oliver Laeyendecker, Lisa A Mills, Thomas Quinn, Steven J Reynolds, John Santelli, Nelson Sewankambo, Simon EF Spencer, Joseph Ssekasanvu, Laura Thomson, Maria Wawer, David Serwadda, Peter Godfrey-Faussett, Joseph Kagaayi, M Kate Grabowski, Oliver Ratmann, Rakai Health Sciences Program, PANGEA-HIV consortium; Growing gender disparity in HIV infection in Africa: sources and policy implications; medRxiv 2023.03.16.23287351; doi: https://doi.org/10.1101/2023.03.16.23287351; 2023

## Acknowledgements
We thank all contributors, program staff and participants to the Rakai Community Cohort Study; all members of the PANGEA-HIV consortium, the [Rakai Health Sciences Program](https://www.rhsp.org/index.php), and CDC Uganda for comments on an earlier version of the manuscript.

We also extend our gratitude to the [Imperial College Research Computing Service](https://doi.org/10.14469/hpc/2232) and the [Biomedical Research Computing Cluster](https://www.bdi.ox.ac.uk/about/biomedical-research-computing) at the University of Oxford for providing the computational resources to perform this study. Additionally, we thank the Office of Cyberinfrastructure and Computational Biology at the [National Institute for Allergy and Infectious Diseases](https://www.niaid.nih.gov/) for data management support.

## Funding

This study was supported by the following organizations:
- [Bill & Melinda Gates Foundation](https://www.gatesfoundation.org/) (OPP1175094, OPP1084362)
- [National Institute of Allergy and Infectious Diseases](https://www.niaid.nih.gov/) (U01AI051171, U01AI075115, UM1AI069530-16, R01AI087409, U01AI100031, R01AI110324, R01AI114438, K25AI114461, R01AI123002, K01AI125086, R01AI128779, R01AI143333, R21AI145682, R01AI155080)
- [National Institute of Mental Health](https://www.nimh.nih.gov/) (F31MH095649, R01MH099733, R01MH107275, R01MH115799)
- [National Institute of Child Health and Development](https://www.nichd.nih.gov/) (R01HD038883, R01HD050180, R01HD070769, R01HD091003)
- Division of Intramural Research of the [National Institute for Allergy and Infectious Diseases](https://www.niaid.nih.gov/) (K01AA024068)
- [National Heart, Lung, and Blood Institute](https://www.nhlbi.nih.gov/) (R01HL152813)
- [Fogarty International Center](https://www.fic.nih.gov/) (D43TW009578, D43TW010557)
- World Bank
- [Doris Duke Charitable Foundation](https://www.ddcf.org/)
- [Johns Hopkins University Center for AIDS Research](https://hopkinscfar.org/) (P30AI094189)
- [President’s Emergency Plan for AIDS Relief](https://www.state.gov/pepfar/) through the Centers for Disease Control and Prevention (NU2GGH000817)

## Quick Start

### System Requirements
- macOS or UNIX, the code was developed on macOS Big Sur 11.7
- [R](https://www.r-project.org/) version >= 4.1.2

### Installation 
Please use the following ```bash``` script to build a conda virtual environment and install all R dependencies:
```shell
git clone https://github.com/MLGlobalHealth/phyloSI-RakaiAgeGender.git
cd phyloSI-RakaiAgeGender
```
for UNIX, please use
```shell
bash phyloSI-RakaiAgeGender-install.sh
```
for MacOS Apple M1 chip, please use
```Shell
bash phyloSI-RakaiAgeGender-install-macosM1.sh
```
If not activated, activate the environment for use:
```shell
source activate phyloSI-RakaiAgeGender
```

### Reproducing our Analyses
We provide all pathogen genomic and epidemiologic input data to reproduce our analyses in non-identifiable aggregate form, or have anonymised individual-level sample identifiers and have randomized individual-level data entries throughout.

Our main analyses depend on estimates of population sizes, HIV prevalence, and HIV suppression, as well as outputs from two phylogenetic analyses: one to estimate the 'time since infection', and the other to detect source-recipient pairs.

To perform the data pre-processing and analysis, navigate to the root directory of the repository and follow the instructions in the next sections.

> **Note** Some scripts will generate figures and other output that is saved in a separate directory outside the repository under `phyloSI-RakaiAgeGender-outputs`. Detailed flowcharts of how data read and written by each R script can be found in `docs/README.md`.

#### RCCS Surveillance Analyses
> **Warning** **(Please read)** Several of the pre-processing code for the surveillance requires the running of computationally demanding Stan models which may take more than 24 hours to finish on a standard laptop computer. We provide summarized outputs for all Stan models and recommend users to skip Stage 1 and proceed directly to Stage 2. To run the full data pre-processing step, execute both Stage 1 and Stage 2.

> **Note** Runtime arguments for Stan models may be configured by editing the contents of `./stan_models/binomial_gp_config.yml`. If your computer has suffcient RAM, we recommend running 4 chains with 4 cores with sampling iterations of 2000 for each chain to reduce Stan runtime.

##### Estimate HIV status and prevalence
```shell
Rscript "./surveillance_pipeline_src/get_estimates_prevalence.R"
```

##### Estimate ART use
```shell
Rscript "./surveillance_pipeline_src/get_estimates_art_coverage_participants.R"
Rscript "./surveillance_pipeline_src/get_estimates_art_coverage_non_participants.R"
```

##### Estimate viral suppresion
```shell
Rscript "./surveillance_pipeline_src/get_estimates_unsuppressed_proportion_participants.R"
Rscript "./surveillance_pipeline_src/get_estimates_unsuppressed_proportion_non_participants.R"
```

##### Estimate treatment cascade
```shell
Rscript "./surveillance_pipeline_src/get_treatment_cascade_participants.R"
Rscript "./surveillance_pipeline_src/get_treatment_cascade_non_participants.R"
Rscript "./surveillance_pipeline_src/get_treatment_cascade_population.R"
```

#### Phylogenetic Analyses
The deep-sequence phylogenetic time since infection estimates were refined using exact patient meta-data that we do not share, and instead we provide the outputs in the `data` directory.  

To reproduce our analyses in part, we share randomized versions of the data that are suffixed by `randomized`. Use these as shown below in combination with the deep-sequence phylogenetic data from the Zenodo repository:
```shell
ZENODO="~/OneDrive/shifting-dynamics-zenodo/" # path to your zenodo directory
DATA_DIR_TSI=$ZENODO/"deep_sequence_phylogenies_primary/data_for_likely_transmission_pairs/phyloscanner-results"
DATA_DIR_PHYLOSCANNER=$ZENODO/"deep_sequence_phylogenies_primary/data_for_time_since_infection/phyloscanner-results"

# get linkage and direction scores
Rscript ./phylo_pipeline_src/find_chains_from_phylogenetics.R --confidential FALSE --phylo-pairs-dir $DATA_DIR_PHYLOSCANNER
# refine time since infection estimates for source-recipient pairs
Rscript ./scripts_for_confidential_data/get_infection_dates_for_phylopairs --confidential FALSE
```

The statistical models present in this repository are built on top of outputs from phylogenetic analyses.
In particular, the same phylogenies as in [Xi et al.](https://doi.org/10.1111/rssc.12544) are used to obtain potential source-recipient pairs, and these can be found in our Zenodo data repository.
Separate phylogenetic analyses were also performed to obtain individual-level estimates of time since infection using the [HIV-phylo-TSI algorithm](https://github.com/BDI-pathogens/HIV-phyloTSI) described in [Golubchik et al.](https://doi.org/10.1101/2022.05.15.22275117). The scripts to perform this analysis can be found in the subdirectory `phylo_pipeline_src`
After the analyses were over, we were able to date transmission events as per the script `confidential_data_src/get_infection_dates_for_phylopairs.R`.

#### Age-specific HIV incidence rates
To run the age-specific HIV incidence rates analysis, run: 
```shell
cd phyloSI-RakaiAgeGender
$ Rscript scripts/run_incidence_rates_estimation.R
```

#### Transmission flows analysis
For the transmission flows analysis, we provide a bash shell script that can be run on a laptop. 

Set the **absolute path** to the output directory where the results should be stored in **line 7** of `phyloSI-RakaiAgeGender-run_phyloflows-laptop.sh`: 
```bash
OUTDIR="absolute/path/to/output/directory"
```

Open up command line tools, cd into the repository, and run:
```bash
$ source activate phyloSI-RakaiAgeGender
$ bash phyloSI-RakaiAgeGender-run_phyloflows-laptop.sh
```

This script will create two shell scripts, `bash_gp_220108-cutoff_2014.sh` and `bash_gp_220108-cutoff_2014-postprocessing.sh`, in the directory specified by `OUTDIR`. The first script, `bash_gp_220108-cutoff_2014.sh`, fits the age- and gender-specific phyloSI model to the Rakai data with Stan. The second script, `bash_gp_220108-cutoff_2014-postprocessing.sh`, performs diagnostic checks on the fitted model, creates posterior summaries, and generates figures. You will need to execute both scripts one after the other:
```bash
$ source activate phyloSI-RakaiAgeGender
$ cd $OUTDIR
$ bash bash_gp_220108-cutoff_2014.sh
$ bash bash_gp_220108-cutoff_2014-postprocessing.sh
