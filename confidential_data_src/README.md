# Confidential Data Processing

This directory contains scripts for processing confidential data.
The directory will be removed from the repository prior to making it a public repository.
Scripts to run experiments or analysis should not depend on the code in this directory.

## Documentation

The script `randomize_metadata_dates.R` loads the file `Rakai_Pangea2_RCCS_Metadata_20221128.RData` in the deepsequence data repositor, produces a randomized version: `Rakai_Pangea2_RCCS_Metadata_randomized.RData`, and stores it in the repository's data folder.
In particular, sex, community type are kept constant for each individual, visit dates are shuffled among participants of the same round, while dates of birth, dates of first positive and last negative tests are modified by adding a random number of days from -2.5*365 to 2.5*365. 
Additional operations guarantee that the date of first positive is successive to the date of last negative.

The script `get_sample_collection_dates_and_randomize.R` extracts the dates of collection for each sequenced blood sample. 
As the these data may be confidential, two versions are generated. The summary of the original data is stored in `confidential_data/sequences_collection_dates.rds`. A randomized version is instead stored in `data/sequences_collection_dates_randomized.rds`. Here, the visit dates are modified adding
a random number of days between -365 and 365.

## Execution

Each of the above script can be run from the command line, after checking that the `indir.deepsequencedata` and `indir.deepanalyses_xiaoyue` paths point to the relevant HPC directories.
```{bash}
Rscript confidential_data_src/randomize_metadata_dates.R
Rscript confidential_data_src/get_sample_collection_dates_and_randomize.R
```
or with an adequate substitute on Windows.
