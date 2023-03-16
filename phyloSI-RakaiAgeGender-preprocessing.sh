#!/bin/sh

# Create stage 1 script
cat > "preproc-s1.sh" <<EOF
#!/bin/bash

Rscript "./misc/get_estimates_art_coverage_participants.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_participants.R"
Rscript "./misc/get_estimates_art_coverage_non_participants.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_non_participants.R"

Rscript "./misc/get_treatment_cascade_non_participants.R"
Rscript "./misc/get_treatment_cascade_participants.R"
EOF

# Create stage 2 script
cat > "preproc-s2.sh" <<EOF
#!/bin/bash

Rscript "./misc/get_estimates_art_coverage_participants_vl200.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_participants.R"
Rscript "./misc/get_estimates_art_coverage_non_participants.R"
Rscript "./misc/get_estimates_unsuppressed_proportion_non_participants.R"

Rscript "./misc/get_treatment_cascade_non_participants.R"
Rscript "./misc/get_treatment_cascade_participants.R"
EOF

# Create stage 3 script
cat > "preproc-s3.sh" <<EOF
#!/bin/bash

Rscript "./misc/get_treatment_cascade_population.R"
Rscript "./misc/get_estimates_prevalence.R"

Rscript "./misc/get_unsuppressed_median_age.R"
Rscript "./misc/get_unsuppressed_share_sex.R"
Rscript "./misc/get_unsuppressed_ratio_sex.R"
Rscript "./misc/get_unsuppressed_prevalence_share_sex.R"
EOF