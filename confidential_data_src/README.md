# Confidential Data Processing

This directory contains scripts for processing confidential data.
The directory will be removed from the repository prior to making it a public repository.
Scripts to run experiments or analysis should not depend on the code in this directory.

## Processing RCCS files

Scripts that load in sensitive data can be separated into two types:
* scripts that produce outputs stored in the repository's `data` folder. (some of these also store outputs on the HPC).
* scripts that process the original data and output sensitive summaries, stored on the HPC.

### Scripts outputing sensitive data.

| name                                             | outputs in `data`                                            | outputs on HPC                          |
| ------------------------------------------------ | ------------------------------------------------------------ | --------------------------------------- |
| `get_art_coverage_non_participant_count`         | `aggregated_newlyregistered_count_art_coverage.csv`          | NA                                      |
| `get_art_coverage_non_participant_count_vl200.R` | `aggregated_newlyregistered_count_art_coverage_vl200.csv`    | NA                                      |
| `get_art_coverage_participants_count.R`          | `aggregated_participants_count_art_coverage.csv`, `sensitivity_specificity_art.csv` | `table_sensitivity_specificity_art.rds` |
|                                                  |                                                              |                                         |

### Scripts outputing in the `data` folder

| name                        | outputs                                            |
| --------------------------- | -------------------------------------------------- |
| `create_hiv_dataset_R914.R` | `HIV_R09_R14.csv`                                  |
| `clean_quest_hiv.R`         | `HIV_R6_R18_221129.csv`, `Quest_R6_R18_221208.csv` |


(Do I want to summarise the original data? Probably not)
| file                  | code  | description                         |
|-----------------------|-------|-------------------------------------|
| `community_names.csv` | comm  | Questionnaire for rounds R09 to R14 |
| `quest_R09_R14.csv`   | q0914 | Questionnaire for rounds R09 to R14 |
| `HIV_`                | q0914 | Questionnaire for rounds R09 to R14 |







## Documentation

**Feel free to change the output directory for the below files. If not confidential, I decided to push stuff on the data subfolder**

The script `randomize_metadata_dates.R` loads the file `Rakai_Pangea2_RCCS_Metadata_20221128.RData` in the deepsequence data repositor, produces a randomized version: `Rakai_Pangea2_RCCS_Metadata_randomized.RData`, and stores it in the repository's data folder.
In particular, sex, community type are kept constant for each individual, visit dates are shuffled among participants of the same round, while dates of birth, dates of first positive and last negative tests are modified by adding a random number of days from -2.5*365 to 2.5*365. 
Additional operations guarantee that the date of first positive is successive to the date of last negative.

The script `get_sample_collection_dates_and_randomize.R` extracts the dates of collection for each sequenced blood sample. 
As the these data may be confidential, two versions are generated. The summary of the original data is stored in `confidential_data/sequences_collection_dates.rds`. A randomized version is instead stored in `data/sequences_collection_dates_randomized.rds`. Here, the visit dates are modified adding
a random number of days between -365 and 365.

The script `get_tsi_predictions.R` removes `visit_dt` and the `PANGEA_ID` columns, as well as estimates on the dates of infection. Instead, these are re-computed as needed when necessary, and true dates are only used in confidential execution. Results are stored in `data/TSI_estimates.csv`

The script `find_chains_from_phylogenetics.R` summarises the results from the phylogenetic analyses, extracting pairwise relationships from all the inferred phylogenies. Additionally, it loads metadata on hosts' serohistories to change phylogenetically inferred directions of transmission whenever these are inconsistent with the serohistory.
The results are then stored in `data/Rakai_phscnetworks_ruleo_sero.rda`.

The script `get_infection_dates_for_phylopairs.R` loads the pairwise relationships obtained above, builds the transmission network and uses these as well as the serohistory and HIV-phylo-TSI's estimate to obtain an estimate of the date of infection through an importance sampling algorithm.


## Execution

Each of the above script can be run from the command line, after checking that the `indir.deepsequencedata` and `indir.deepanalyses_xiaoyue` paths point to the relevant HPC directories.
```{bash}
Rscript confidential_data_src/randomize_metadata_dates.R
Rscript confidential_data_src/get_sample_collection_dates_and_randomize.R
Rscript confidential_data_src/get_tsi_predictions.R
```
or with an adequate substitute on Windows.
