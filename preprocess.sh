#!/bin/sh

# Rscript ./misc/get_estimates_prevalence.R

Rscript ./misc/get_estimates_unsuppressed_proportion_non_participants.R
Rscript ./misc/get_treatment_cascade_participants.R

Rscript ./misc/get_treatment_cascade_non_participants.R