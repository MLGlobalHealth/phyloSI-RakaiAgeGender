# Phylogenetic pipelines

This directory contains scripts documenting the phylogenetic analyses performed prior to our statistical modelling.
We performed two separate analyses:
- The source-recipient pair analyses, using the same pipeline as [Xi et al](https://doi.org/10.1111/rssc.12544)
- The Time Since Infection (TSI) analyses, requiring a modification in the above pipeline.

The most notable difference in the two is the algorithm to partition hosts sequences together.
In the first case, we want to maximise the likelihood of source-recipient pairs to be in the same group, and are interested expensive pairwise relationship analyses.
As such, we need split sequences in small groups with low-genetic distance among members.
On the other hand, in TSI analyses, we are interested in obtaining robust phylogenies summarising viral evolution.
From these, we extract host-level subtrees and summarise them through features such as largest root to tip distance. 
These features will be passed as input to the [HIV-phylo-TSI algorithm](https://github.com/BDI-pathogens/HIV-phyloTSI/tree/main/ExampleInputs).

## Retrieving source-recipient pairs.

As mentioned above, the source-recipient pairs were obtained through analysis of the same phylogenetic pipeline as in [Xi et al.]( https://doi.org/10.1111/rssc.12544)
Here, we modify the procedure to estimate transition networks, as per the script `find_chains_from_phylogenetics.R`.
This summarises the results from the phylogenetic analyses, extracting pairwise relationships from all the inferred phylogenies. 
In our analyses, it loads confidential metadata on hosts' serohistories to change phylogenetically inferred directions of transmission whenever these are inconsistent with the serohistory. 
The results are then stored in `data/Rakai_phscnetworks_ruleo_sero.rda`.


## Time Since Infection pipeline

This folder contains the scripts that were used to perform the phylogenetic analyses.
These require extensive computational resources, and in our case were run on Imperial College London HPC.
Note that these scripts cannot be run on any computer, but the code in the repository documents the different steps of the analyses. In particular, the code is currently maintained [in this repo](`https://github.com/olli0601/Phyloscanner.R.utilities/tree/master/misc_data_analysis_RCCS1519/software`).

The bash script `runall_TSI_pairs2.sh` describes the order in which the scripts are run, and is therefore denoted as the 'controller'.
The primary inputs are two: `\*bam` files and base frequency files `\*BaseFreqs_WithHXB2.csv` for each NGS sequence.
The second is only necessary to compute Minor Allele Frequencies necessary for the HIV-phylo-TSI algorithm.

| code | script                            | role                    | output                                |
| ---- | --------------------------------- | ----------------------- | ------------------------------------- |
| net  | `TSI_initialise.R`                | group sequences         | `clusters.rds`; `phscinput_runs_*rds` |
| ali  | `make_deep_sequence_alignments.R` | queue alignment jobs    | `readali*.sh`                         |
| btr  | `make_trees.R`                    | queue iqtree jobs       | `srx*job.sh`                          |
| ctr  | `check_trees.R`                   | check and re-queue btr  | `srx*job.sh`                          |
| atr  | `analyse_trees.R`                 | queue trees analyses    | `phsc_tsi*.sh`                        |
| tsi  | `TSI_run_predictions.R`           | run HIV-phylo-TSI       | `tsi*sh`                              |
| dti  | `TSI_estimate_dates.R`            | aggregate TSI estimates | `aggregated_TSI*.csv`                 |

Details on the shell jobs is summarised below:

| job             | program used | input                               | output                             | output by         |
| --------------- | ------------ | ----------------------------------- | ---------------------------------- | ----------------- |
| `readali*.sh`   | MAFFT        | `*bam` files + consensus sequences  | `In_Window_XX_toYY*.fasta`         | group , window   |
| `srx*job.sh`    | IQTREE       | `In_Window_XX_toYY*fasta`           | `*iqtree` , `*.treefile`   | group , window   |
| `phsc_tsi*.sh`       | Phyloscanner |`*iqtree`   , `*.treefile`| `PatStats*csv` among others |group|
| `tsi*sh`        | HIV-phyloTSI | `PatStats*csv`, `maf*csv`           | `ptyr*_tsi.csv`                    | group             |

### De-identification: is this necessary?

The script `get_tsi_predictions.R` removes `visit_dt` and the `PANGEA_ID` columns, as well as estimates on the dates of infection. Instead, these are re-computed as needed when necessary, and true dates are only used in confidential execution. Results are stored in `data/TSI_estimates.csv`
