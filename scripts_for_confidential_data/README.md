# Confidential Data Processing

This directory contains scripts for processing confidential data.
The directory will be removed from the repository prior to making it a public repository.
Scripts to run experiments or analysis should not depend on the code in this directory.

## Processing RCCS files

Scripts that load in sensitive data can be separated into two types:
* scripts that produce outputs stored in the repository's `data` folder. (some of these also store outputs on the HPC).
* scripts that process the original data and output sensitive summaries, stored on the HPC.



## Scripts outputing in the `data` folder


| name                                                         | outputs in `data`                                            | outputs on HPC                          |
| ------------------------------------------------------------ | ------------------------------------------------------------ | --------------------------------------- |
| `get_art_coverage_non_participant_count`                     | `aggregated_newlyregistered_count_art_coverage.csv`          | NA                                      |
| `get_art_coverage_non_participant_count_vl200.R`             | `aggregated_newlyregistered_count_art_coverage_vl200.csv`    | NA                                      |
| `get_art_coverage_participants_count.R`                      | `aggregated_participants_count_art_coverage.csv`, `sensitivity_specificity_art.csv` | `table_sensitivity_specificity_art.rds` |
| `get_art_coverage_participants_count_vl200.R`                | `aggregated_participants_count_art_coverage_vl200.csv`, `sensitivity_specificity_art_vl200.csv` |                                         |
| `get_unsuppressed_proportion_non_participants_count_vl200.R  ` | `aggregated_newlyregistered_count_unsuppressed_vl200.csv`    |                                         |
| `get_hiv_status_count.R`                                     | `aggregated_count_hiv_positive.csv`                          |                                         |
| `randomize_metadata_dates.R`                                 | `Rakai_Pangea2_RCCS_Metadata_randomized.RData`               | NA                                      |
| `get_sample_collection_dates_and_randomize.R`                | `sequences_collection_dates_randomized.rds`                  | `sequences_collection_dates.rds`        |

### Scripts outputing sensitive data in `$DEEPDATA`

| name                          | outputs                                            |
| ----------------------------- | -------------------------------------------------- |
| `create_hiv_dataset_R914.R`   | `HIV_R09_R14.csv`                                  |
| `clean_quest_hiv.R`           | `HIV_R6_R18_221129.csv`, `Quest_R6_R18_221208.csv` |
| `get_df_round.R`              | `RCCS_round_timeline_220905.RData`                 |
| `get_census_eligible_count.R` | `RCCS_census_eligible_individuals_221116.csv`      |
| `get_participant_count.R`     | `RCCS_participation_221208.csv`                    |
| `write_updated_serohistory.R` | `221128_requested_updated_serohistory.csv`         |



## Randomization details

The script `randomize_metadata_dates.R` loads the file `Rakai_Pangea2_RCCS_Metadata_20221128.RData` in the deepsequence data repositor, produces a randomized version: `Rakai_Pangea2_RCCS_Metadata_randomized.RData`, and stores it in the repository's data folder.
In particular, sex, community type are kept constant for each individual, visit dates are shuffled among participants of the same round, while dates of birth, dates of first positive test and last negative tests are modified by adding a random number of days from $-.25*365$ to $.25*365$ (so the dates are randomised within 6 months). 
Additional operations guarantee that the date of first positive is successive to the date of last negative.

The script `get_sample_collection_dates_and_randomize.R` extracts the dates of collection for each sequenced blood sample. 
As the these data may be confidential, two versions are generated. The summary of the original data is stored on the HPC at `sequences_collection_dates.rds`. A randomized version is instead stored in `data/sequences_collection_dates_randomized.rds`. Here, the visit dates are modified adding a random number of days between -$.25*365$ and $.25*365$.

### Execution

Each of the above script can be run from the command line, after checking that the `indir.deepsequencedata` and `indir.deepanalyses_xiaoyue` paths point to the relevant HPC directories.
```{bash}
Rscript confidential_data_src/randomize_metadata_dates.R
Rscript confidential_data_src/get_sample_collection_dates_and_randomize.R
```
or with an adequate substitute on indows.
