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
The code in this repository is licensed under [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) by Imperial College London.

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
If not activated, activate the environment for use:
```shell
source activate phyloSI-RakaiAgeGender
```

### Reproducing our Analyses
We provide all pathogen genomic and epidemiologic input data to reproduce our analyses in non-identifiable aggregate form, or have anonymised individual-level sample identifiers and have randomized individual-level data entries throughout. Please download the data from our [Zenodo data repository](https:/zenodo.org/record/8412741). Note that for ethical considerations a small subset of the data are not shared.

To perform all data pre-processing and analysis, **the user must update `config.R`** which is located within the root directory of the repository. Locate the following code around line 20 of `config.R` and change `"your-user"` to your system username and `"your-path-to-zenodo-dir"` to the absolute path to the Zenodo directory you have downloaded.
```r
# 1) Change "your-user" to your system user name
# 2) Change "your-path-to-zenodo-dir" to the ABSOLUTE PATH to the Zenodo directory
dir.zenodo <- data.table::fcase(
  # ...
  user == "your-user", "your-path-to-zenodo-dir"
  # ...
)
```

#### Pre-processing steps

Our main analyses depend on estimates of population sizes, HIV prevalence, and HIV suppression, as well as outputs from phylogenetic reconstructions of source-recipient pairs and estimates of the time since infection from when the phylogenetically likely recipient was infected until their sequences were sampled. 

Several of the pre-processing code for the surveillance requires the running of computationally demanding Stan models which may take more than 24 hours to finish on a standard laptop computer. We provide summarized outputs in the Zenodo directory for all pre-processing steps and users may skip directly to the main analysis.

> **Note** Some scripts will generate figures and other output that is saved in a separate directory outside the repository under `phyloSI-RakaiAgeGender-outputs`. Detailed flowcharts of how data read and written by each R script can be found in `docs/README.md`.

> **Note** Runtime arguments for Stan models may be configured by editing the contents of `./stan_models/binomial_gp_config.yml`. If your computer has suffcient RAM, we recommend running 4 chains with 4 cores with sampling iterations of 2000 for each chain to reduce Stan runtime.

> **Note** Runtime arguments for Stan models may be configured by editing the contents of `./stan_models/binomial_gp_config.yml`. If your computer has suffcient RAM, we recommend running 4 chains with 4 cores with sampling iterations of 2000 for each chain to reduce Stan runtime.

> **Note** Deep-sequence phylogenetic trees are time consuming to generate, and provided as part of the data directory. These were generated using an in-house pipeline.  

> **Note** Phylogenetic time since infection estimates were obtained using the [HIV-phylo-TSI algorithm](https://github.com/BDI-pathogens/HIV-phyloTSI). Phylo-TSI estimates were refined using exact patient meta-data that are only available upon reasonable request to RHSP. We provide the outputs in the `data` directory. We share randomized versions of the data that are suffixed by `randomized`. Note using these randomised versions of the data will not exactly reproduce our analyses.

To run the full data pre-processing steps:

```shell
# Estimate HIV status and prevalence
Rscript "./surveillance_pipeline_src/get_estimates_prevalence.R"

# Estimate ART use
Rscript "./surveillance_pipeline_src/get_estimates_art_coverage_participants.R"
Rscript "./surveillance_pipeline_src/get_estimates_art_coverage_non_participants.R"

# Estimate viral suppresion
Rscript "./surveillance_pipeline_src/get_estimates_unsuppressed_proportion_participants.R"
Rscript "./surveillance_pipeline_src/get_estimates_unsuppressed_proportion_non_participants.R"

# Estimate treatment cascade
Rscript "./surveillance_pipeline_src/get_treatment_cascade_participants.R"
Rscript "./surveillance_pipeline_src/get_treatment_cascade_non_participants.R"
Rscript "./surveillance_pipeline_src/get_treatment_cascade_population.R"

# Phylogenetic linkage and direction scores
# Note: using pre-generated deep-sequence phylogenies 
Rscript "./phylo_pipeline_src/find_chains_from_phylogenetics.R" --confidential FALSE

# Refine time since infection estimates for source-recipient pairs
# Note: using randomized data
Rscript "./scripts_for_confidential_data/get_infection_dates_for_phylopairs.R" --confidential FALSE
```

#### Main analysis 
To run the age-specific HIV incidence rates analysis, run: 
```shell
cd phyloSI-RakaiAgeGender
Rscript "./scripts/run_incidence_rates_estimation.R"
```

For the transmission flows analysis, we provide a bash shell script that can be run on a laptop. Set the **absolute path** to the output directory where the results should be stored in **line 7** of `phyloSI-RakaiAgeGender-run_phyloflows-laptop.sh`: 
```bash
OUTDIR="absolute/path/to/output/directory"
```

Open up command line tools, cd into the repository, and run:
```bash
$ bash phyloSI-RakaiAgeGender-run_phyloflows-laptop.sh
```

This script will create two shell scripts, which you will need to execute one after the other:
```bash
source activate phyloSI-RakaiAgeGender
cd $OUTDIR
bash bash_gp_220108-cutoff_2014.sh
bash bash_gp_220108-cutoff_2014-postprocessing.sh
```
