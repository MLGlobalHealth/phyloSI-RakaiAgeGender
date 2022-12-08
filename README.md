# phyloflows: phylo-SIR model for estimating shifting patterns in HIV transmission dynamics in Rakai, Uganda.

## Overview
This repository contains the code for the phylo-SIR model to estimate shifting patterns in HIV transmission dynamics in Rakai, Uganda.

## Data
The data in ```data/``` contain
* source-recipient pairs, in file ```pairsdata_toshare_d1_w11_netfrompairs_seropairs.rds```.

* count of census eligible population by age, gender and round, in file ```RCCS_census_eligible_individuals_221116.csv```
* participation rates to the RCCS survey by age, gender and round, in file ```RCCS_participation_221208.csv```
* count of participants by hiv status, age, gender and round, in file ```aggregated_count_hiv_positive.csv``` used in ```misc/get_estimates_prevalence.R``` to estimate smooth proportion of hiv prevalence among population.
* count of participants by self-reported art use, age, gender and round, in file ```aggregated_participants_count_art_coverage.csv``` used in ```misc/get_estimates_art_coverage_participants.R``` to estimate smooth proportion of art coverage among participants.
* count of participants by viremic viral loads, age, gender and round, in file ```aggregated_participants_count_unsuppressed.csv``` used in ```misc/get_estimates_unsuppressed_proportion_participants.R``` to estimate smooth proportion of viral suppression among participants.
* count of first-time participants by self-reported art use, age, gender and round, in file ```aggregated_newlyregistered_count_art_coverage.csv``` used in ```misc/get_estimates_art_coverage_non_participants.R``` to estimate smooth proportion of art coverage among first-time participants.
* count of first-time participants by viremic viral loads, age, gender and round, in file ```aggregated_newlyregistered_count_unsuppressed.csv``` used in ```misc/get_estimates_unsuppressed_proportion_non_participants.R``` to estimate smooth proportion of viral suppression among first-time participants.

* Individual-level seroconvert cohort data with information on rounds of enrollment, age, sex, hiv status, in file ```seroconverter_cohort_R6R19.rds```
* HIV incidence rates estimates by age, gender and round, in file ```Rakai_incpredictions_inland_221107.csv```

* Sexual contact rate estimates by age and gender for round 15, in file ```inland_R015_cntcts_rate_1130.rds```

## Fit
* smooth proportion of hiv prevalence in population by age, gender and round, in file ```RCCS_prevalence_estimates_DATE.csv```
* smooth estimates of art coverage and viral suppression in population by age, gender and round in file ```RCCS_treatment_cascade_population_estimates_DATE.csv```
* smooth estimates of art coverage and viral suppression in participants by age, gender and round in file ```RCCS_treatment_cascade_participants_estimates_DATE.csv```
* smooth estimates of art coverage and viral suppression in first-time participants by age, gender and round in file ```RCCS_treatment_cascade_nonparticipants_estimates_DATE.csv```

## License
- The code in this repository is licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) by Imperial College London. Copyright Imperial College London 2022. 
- The data ??

## Warranty

Imperial makes no representation or warranty about the accuracy or completeness of the data nor that the results will not constitute in infringement of third-party rights. Imperial accepts no liability or responsibility for any use which may be made of any results, for the results, nor for any reliance which may be placed on any such work or results.

## Cite

Please cite 
XXX

## Acknowledgements

XXX

## Funding

XXX

## System Requirements
- macOS or UNIX, the code was developed on macOS Big Sur 11.7
- [R](https://www.r-project.org/) version >= 4.1.2


## Installation 
A ```yml``` file is provided and can be used to build a conda virtual environment containing all R dependencies. Create the environment using:
```bash
$ cd phyloflows
$ conda env create -f phyloflows.yml
```
Then activate the environment for use:
```bash
$ source activate phyloflows
```

## Usage
### The files in ```src/```
The scripts in ```src/``` read the identifiable individual-level data, de-identifiable and aggregate them, and save them in ```data/```. 

### The files in ```misc/```
The scripts in ```misc/``` read the aggregated data and fit smoothing models to obtain hiv prevalence, art coverage and viral suppression for participants and first-time participants by age, gender and round.
* Obtain smooth hiv prevalence by age, gender and round with file ```get_estimates_prevalence.R```
* Obtain smooth art coverage in participants by age, gender and round with file ```get_estimates_art_coverage_participants.R```
* Obtain smooth art coverage in first-time participants by age, gender and round with file ```get_estimates_art_coverage_non_participants.R```
* Obtain smooth viral suppression in participants by age, gender and round with file ```get_estimates_unsuppressed_proportion_participants.R```
* Obtain smooth viral suppression in first-time participants by age, gender and round with file ```get_estimates_unsuppressed_proportion_non_participants.R```
* Combine art coverage and viral suppression estimates to obtain treatment cascade in participants with file ```get_treatment_cascade_participants.R```
* Combine art coverage and viral suppression estimates to obtain treatment cascade in first-time participants with file ```get_treatment_cascade_non_participants.R```
* Combine art coverage and viral suppression estimates to obtain treatment cascade in population with file ```get_treatment_cascade_population.R```

### Run the incidence rates analysis
* To run the age-specific HIV incidence rates analysis, run ```scripts/run_incidence_rates_estimation.R```

### Run the transmission flows analysis
We provide a script that can be run on a laptop. The following modifications need to be done to the start of the bash script
```bash
run_stan_laptop.sh
```
First, set the the directory in which the repository is and the the output directory (where the results should be stored):
```bash
INDIR="/rds/general/user/mm3218/home/git/phyloflows"
OUTDIR="/rds/general/user/mm3218/home/projects/2021/phyloflows"
```
Second, set as appropriatly the virtual environment 
```bash
module load anaconda3/personal
```
Lastly, from the repository directory, on the terminal console execute,
```bash
./run_stan_laptop.sh
```
This generates in the output directory two bash scripts. One for running the model and one for processing the results. Execute these bash scripts one after the other.



 
