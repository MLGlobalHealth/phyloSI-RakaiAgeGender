# phyloSI-RakaiAgeGender
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

**Welcome!** This repository contains the code and data for the analyses presented in the paper *Growing gender inequity in HIV infection in Africa: sources and policy implications* by M Monod, A Brizzi, R Galiwango, R Ssekubugu, Y Chen, X Xi et al.

- [License](#license)
- [Warrenty](#warranty)
- [Citation](#cite)
- [Acknowledgements](#acknowledgements)
- [Funding](#funding)
- [Quick Start](#quick-start)
  - [System Requirements](#system-requirements)
  - [Installation](#installation)
  - [Data preprocessing](#data-preprocessing)
  - [Age-specific HIV incidence rates](#age-specific-hiv-incidence-rates)
  - [Transmission flows analysis](#transmission-flows-analysis)
- [Data and Script Reference](#data-and-script-reference)
  - [Sample Data](#sample-data)
  - [Generated Data](#generated-data)

## License
The code and data in this repository are licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) by Imperial College London. Copyright Imperial College London 2022. 

## Warranty

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.

## Citation

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
A ```yml``` file is provided and can be used to build a conda virtual environment containing all R dependencies. Create the environment using:
```bash
$ cd phyloSI-RakaiAgeGender
$ conda env create -f phyloSI-RakaiAgeGender.yml
```
Then activate the environment for use:
```bash
$ source activate phyloSI-RakaiAgeGender
```

### Data preprocessing

We provide all pathogen genomic and epidemiologic input data to reproduce our analyses in non-identifiable aggregate form, or have anonymised individual-level sample identifiers and have randomized individual-level data entries throughout.

Our main analysis functions on estimating incidence dynamics and transmission flows depend on estimates of population sizes, prevalence, virus suppression.

To perform the data preprocessing required to run downstream models, navigate to the root directory of the repository and execute the following commands

#### Stage 1
```shell
Rscript "./misc/get_estimates_art_coverage_participants.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_participants.R"
Rscript "./misc/get_estimates_art_coverage_non_participants.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_non_participants.R"
Rscript "./misc/get_treatment_cascade_non_participants.R"
Rscript "./misc/get_treatment_cascade_participants.R"
```

#### Stage 2
```shell
Rscript "./misc/get_estimates_art_coverage_participants_vl200.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_participants.R"
Rscript "./misc/get_estimates_art_coverage_non_participants.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_non_participants.R"
Rscript "./misc/get_treatment_cascade_non_participants.R"
Rscript "./misc/get_treatment_cascade_participants.R"
```

#### Stage 3
```shell
Rscript "./misc/get_treatment_cascade_population.R"
Rscript "./misc/get_estimates_prevalence.R"
Rscript "./misc/get_unsuppressed_median_age.R"
Rscript "./misc/get_unsuppressed_share_sex.R"
Rscript "./misc/get_unsuppressed_ratio_sex.R"
Rscript "./misc/get_unsuppressed_prevalence_share_sex.R"
```

Alternatively, you can run the following shell script for convenience.
```shell
$ bash preproc.sh
```
`preproc.sh` will generate three new shell scripts: `preproc-s1.sh`, `preproc-s2.sh`, and `preproc-s3.sh`. Each script will invoke the required R scripts for each stage of the data preprocessing process. They must be executed in order:
```shell
$ bash preproc-s1.sh # stage 1
$ bash preproc-s2.sh # stage 2
$ bash preproc-s3.sh # stage 3
```
> :exclamation: Note that some preprocessing scripts will generate figures which is saved in a separate directory outside the repository under `phyloSI-RakaiAgeGender-outputs`.

Detailed flowcharts of how data read and written by each R script within each preprecessing stage can be found in `docs/README.md`.

### Age-specific HIV incidence rates
To run the age-specific HIV incidence rates analysis, navigate to the root directory of the repository. 
```shell
cd phyloSI-RakaiAgeGender
```
Then run the following R script: 
```shell
$ Rscript scripts/run_incidence_rates_estimation.R
```
 > :exclamation: The output will be save in a seperate directory outside the repository under the name `phyloSI-RakaiAgeGender-outputs`.

### Transmission flows analysis
For the transmission flows analysis, we provide a bash shell script that can be run on a laptop. 

Set the **absolute path** to the output directory where the results should be stored in **line 7** of `run_stan_laptop.sh`: 
```bash
OUTDIR="absolute/path/to/output/directory"
```

Open up command line tools, cd into the repository, and run:
```bash
$ source activate phyloSI-RakaiAgeGender
$ bash run_stan_laptop.sh
```

This script will create two shell scripts, `bash_gp_220108-cutoff_2014.sh` and `bash_gp_220108-cutoff_2014-postprocessing.sh`, in the directory specified by `OUTDIR`. The first script, `bash_gp_220108-cutoff_2014.sh`, fits the age- and gender-specific phyloSI model to the Rakai data with Stan. The second script, `bash_gp_220108-cutoff_2014-postprocessing.sh`, performs diagnostic checks on the fitted model, creates posterior summaries, and generates figures. You will need to execute both scripts one after the other:
```bash
$ source activate phyloSI-RakaiAgeGender
$ cd $OUTDIR
$ bash bash_gp_220108-cutoff_2014.sh
$ bash bash_gp_220108-cutoff_2014-postprocessing.sh
```


## Data and script reference

### Sample Data
The table below lists the data files within `/data` and a brief description of its contents

<details>
<summary><b>Click to show table</b></summary>
  <table>
    <thead>
      <tr>
        <th>File name</th>
        <th>Description</th>
      </tr>
    </thead>
  <tbody>
    <tr>
      <td><code>pairsdata_toshare_d1_w11_netfrompairs_postponessrm.rds</code></td>
      <td>HIV source-recipient pairs</td>
    </tr>
    <tr>
      <td><code>RCCS_census_eligible_individuals_221116.csv</code></td>
      <td>Count of census eligible population by age, gender and round</td>
    </tr>
      <tr>
      <td><code>RCCS_participation_221208.csv</code></td>
      <td>Participation rates to the RCCS survey by age, gender and round</td>
      </tr>
    <tr>
      <td><code>aggregated_count_hiv_positive.csv</code></td>
      <td>Count of participants by hiv status, age, gender and round used in <code>misc/get_estimates_prevalence.R</code> to estimate smooth proportion of hiv prevalence among population.</td>
    </tr>
    <tr>
      <td><code>aggregated_participants_count_art_coverage.csv</code></td>
      <td>count of participants by self-reported art use, age, gender and round, used in <code>misc/get_estimates_art_coverage_participants.R</code> to estimate smooth proportion of art coverage among participants.</td>
    </tr>
    <tr>
      <td><code>aggregated_participants_count_unsuppressed.csv</code></td>
      <td>count of participants by viremic viral loads, age, gender and round, used in <code>misc/get_estimates_unsuppressed_proportion_participants.R</code> to estimate smooth proportion of viral suppression among participants.</td>
    </tr>
    <tr>
      <td><code>aggregated_newlyregistered_count_art_coverage.csv</code></td>
      <td>count of first-time participants by self-reported art use, age, gender and round, used in <code>misc/get_estimates_art_coverage_non_participants.R</code> to estimate smooth proportion of art coverage among first-time participants.</td>
    </tr>
    <tr>
      <td><code>aggregated_newlyregistered_count_unsuppressed.csv</code></td>
      <td>count of first-time participants by viremic viral loads, age, gender and round, used in <code>misc/get_estimates_unsuppressed_proportion_non_participants.R</code> to estimate smooth proportion of viral suppression among first-time participants.</td>
    </tr>
    <tr>
      <td><code>seroconverter_cohort_R6R19.rds</code></td>
      <td>Individual-level seroconvert cohort data with information on rounds of enrollment, age, sex, hiv status</td>
      </tr>
      <tr>
      <td><code>Rakai_incpredictions_inland_221107.csv</code></td>
      <td>HIV incidence rates estimates by age, gender and round</td>
    </tr>
    <tr>
      <td><code>inland_R015_cntcts_rate_1130b.rds</code></td>
      <td>Sexual contact rate estimates by age and gender for round 15</td>
    </tr>
    </tbody>
  </table>
</details>

### Generated Data
<details>
<summary><b>Click to show table</b></summary>
<table>
  <thead>
    <tr>
      <th>File name</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>RCCS_prevalence_estimates_DATE.csv</code></td>
      <td>Smooth proportion of HIV prevalence in population by age, gender and round</td>
    </tr>
    <tr>
      <td><code>RCCS_treatment_cascade_population_estimates_DATE.csv</code></td>
      <td>Smooth estimates of ART coverage and viral suppression in population by age, gender and round</td>
    </tr>
    <tr>
      <td><code>RCCS_treatment_cascade_participants_estimates_DATE.csv</code></td>
      <td>Smooth estimates of ART coverage and viral suppression in participants by age, gender and round</td>
    </tr>
    <tr>
      <td><code>RCCS_treatment_cascade_nonparticipants_estimates_DATE.csv</code></td>
      <td>Smooth estimates of ART coverage and viral suppression in first-time participants by age, gender and round</td>
    </tr>
  </tbody>
</table>
</details>

### Scripts

The following table list the actions performed by each of the available scripts:
<details>
  <summary><b>Click to show table</b></summary>
  <table>
    <thead>
      <tr>
        <th>Action</th>
        <th>Script</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>Obtain smooth HIV prevalence by age, gender, and round</td>
        <td><code>misc/get_estimates_prevalence.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth ART coverage in participants by age, gender, and round</td>
        <td><code>misc/get_estimates_art_coverage_participants.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth ART coverage in first-time participants by age, gender, and round</td>
        <td><code>misc/get_estimates_art_coverage_non_participants.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth viral suppression in participants by age, gender, and round</td>
        <td><code>misc/get_estimates_unsuppressed_proportion_participants.R</code></td>
      </tr>
      <tr>
        <td>Obtain smooth viral suppression in first-time participants by age, gender, and round</td>
        <td><code>misc/get_estimates_unsuppressed_proportion_non_participants.R</code></td>
      </tr>
      <tr>
        <td>Combine ART coverage and viral suppression estimates to obtain treatment cascade in participants</td>
        <td><code>misc/get_treatment_cascade_participants.R</code></td>
      </tr>
      <tr>
        <td>Combine ART coverage and viral suppression estimates to obtain treatment cascade in first-time participants</td>
        <td><code>misc/get_treatment_cascade_non_participants.R</code></td>
      </tr>
      <tr>
        <td>Combine ART coverage and viral suppression estimates to obtain treatment cascade in population</td>
        <td><code>misc/get_treatment_cascade_population.R</code></td>
      </tr>
    </tbody>
  </table>
</details>