#!/bin/sh

# Create stage 1 script
cat > "preproc-s1.sh" <<EOF
#!/bin/bash

Rscript "./surveillance_pipeline_src/get_estimates_prevalence.R"
Rscript "./surveillance_pipeline_src/get_estimates_art_coverage_participants.R"
Rscript "./surveillance_pipeline_src/get_estimates_unsuppressed_proportion_participants.R"
Rscript "./surveillance_pipeline_src/get_estimates_art_coverage_non_participants.R"
Rscript "./surveillance_pipeline_src/get_estimates_unsuppressed_proportion_non_participants.R"
EOF

# Create stage 2 script
cat > "preproc-s2.sh" <<EOF
#!/bin/bash

Rscript "./surveillance_pipeline_src/get_treatment_cascade_participants.R"
Rscript "./surveillance_pipeline_src/get_treatment_cascade_non_participants.R"
Rscript "./surveillance_pipeline_src/get_treatment_cascade_population.R"

EOF

cat > "preproc-s4.sh" <<EOF
#!/bin/bash

TSI_OUTPUTS= "" #TOFILL: TSI_outputs in our Zenodo data
PAIRS_OUTPUTS= "" #"TOFILL": PAIRS_outputs in our Zenodo data

# phylogenetic outputs processing:

# get individual level time since infection
Rscript "./phylo_pipeline_src/TSI_estimate_dates.R --confidential FALSE --tsi_out_dir $TSI_OUTPUTS"
#
# get posterior probabilities of linkage and direction
Rscript "./phylo_pipeline_src/find_chains_from_phylogenetics.R --confidential FALSE --phylo-pairs-dir $PAIRS_OUTPUTS"

# combine both and estimate date of infection for pairs
Rscript "./scripts_for_confidential_data/get_infection_dates_for_phylopairs --confidential FALSE"

EOF
