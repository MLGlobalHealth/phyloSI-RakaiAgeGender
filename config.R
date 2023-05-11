# TODO: fill in path to zenodo dir and specify all the following files accordingly

config <- list(
    overwrite.existing.files  = FALSE,
    vl  = 200
)

usr <- Sys.info()[["user"]]

# set your user and path here:
dir.zenodo <- data.table::fcase(
    usr == "andrea", "~/OneDrive/shifting-dynamics-zenodo/",
    usr == "shozendan", "~/Imperial/phyloSI-RakaiAgeGender-data",
    usr == "your-user", "your-path-to-zenodo-dir",
    default = "")

zenodo.exists <- dir.exists(dir.zenodo)
if (zenodo.exists) {
    list.files(dir.zenodo, pattern = "")
    dir.zenodo.phyloprim <- file.path(dir.zenodo, "data_deep_sequence_phylogenies")
    dir.zenodo.pairs <- file.path(dir.zenodo.phyloprim, "data_for_likely_transmission_pairs")
    dir.zenodo.pairs.phsc <- file.path(dir.zenodo.pairs, "phyloscanner-results")
    dir.zenodo.tsi <- file.path(dir.zenodo.phyloprim, "data_for_time_since_infection")
    dir.zenodo.tsi.phsc <- file.path(dir.zenodo.tsi, "phyloscanner-results")

    dir.zenodo.phyloproc = file.path(dir.zenodo, "preprocessed_outputs_deep_sequence_phylogenies")

    # surveillance data (primary or processed_outputs)
    dir.zenodo.survprim <- file.path(dir.zenodo, "data_RCCS_surveillance")
    dir.zenodo.survproc <- file.path(dir.zenodo, "preprocessed_outputs_RCCS_surveillance")

    # something
    dir.zenodo.central <- file.path(dir.zenodo, "final_central_analysis")

    c(dir.zenodo.phyloprim,
      dir.zenodo.pairs,
      dir.zenodo.pairs.phsc,
      dir.zenodo.tsi,
      dir.zenodo.tsi.phsc,
      dir.zenodo.phyloproc,
      dir.zenodo.survprim,
      dir.zenodo.survproc,
      dir.zenodo.central
    ) |> dir.exists() |> all() |> stopifnot()

}

gitdir.data <- file.path(gitdir, "data")
gitdir.functions <- file.path(gitdir, "functions")
gitdir.misc <- file.path(gitdir, "surveillance_pipeline_src")

######################
# confidential paths #
######################


# defaults implies dir.exists(indir.deepsequencedata) == FALSE
indir.deepsequencedata <- ""
indir.deepsequence_analyses <- ""

if (usr == "andrea") {
    indir.deepsequencedata <- "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live"
    indir.deepsequence_analyses <- "/home/andrea/HPC/project/ratmann_deepseq_analyses/live"
    indir.deepanalyses.xiaoyue <- "/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live"
}

if (usr %in% c("melodie", "ratmann")) {
    indir.deepsequencedata <- ""
    indir.deepsequence_analyses <- "~/Box\ Sync/2021/ratmann_deepseq_analyses/live/"
}

if (dir.exists(indir.deepsequencedata)) {
    indir.deepsequencedata.r151r18 <- file.path(indir.deepsequencedata, "RCCS_R15_R18")
    file.community.keys <- file.path(indir.deepsequencedata.r151r18, "community_names.csv")

    path.meta.confidential <- file.path(indir.deepsequencedata.r151r18, "Rakai_Pangea2_RCCS_Metadata_20221128.RData")
    path.collection.dates.confidential <- file.path(indir.deepsequencedata.r151r18, "sequences_collection_dates.rds")
}

################################
# paths from/for phylogenetics #
################################

# get_infection_dates_for_phylopairs.R
path.meta.randomized <- file.path(dir.zenodo.survprim,"Rakai_Pangea2_RCCS_Metadata_randomized.RData" )
path.collection.dates.randomized <- file.path(dir.zenodo.survprim, "sequences_collection_dates_randomized.rds")
file.path.round.timeline <- file.path(dir.zenodo.survprim, "RCCS_round_timeline_220905.RData" )

if (zenodo.exists) {
    path.chains.data <- file.path(dir.zenodo.phyloproc, "Rakai_phscnetworks_ruleo_sero.rda")
    path.tsiestimates <- file.path(dir.zenodo.phyloproc, "aggregated_TSI.rds")
}

###################################
#  Paths for incidence estimation #
###################################

# primary data:
file.path.seroconverter_cohort <- file.path(dir.zenodo.survprim , "seroconverter_cohort_R6R19.rds")

# Sensitivity analyses:
# sero conversion cohort status (continuously surveyed communities)
file.path.seroconverter_cohort.30 <- file.path(dir.zenodo.survprim, "seroconverter_cohort_30comm_R6R19.rds")


# outputs
output.dir.incidence.estimation <- file.path(dir.zenodo.survproc, "run_incidence_rates_estimation")

file.incidence.inland	<- file.path(output.dir.incidence.estimation, "Rakai_incpredictions_inland_221107.csv")
file.incidence.30com.inland	<- file.path(output.dir.incidence.estimation, "Rakai_incpredictions_inland_221119.csv")
file.incidence.loess.inland	<- file.path(output.dir.incidence.estimation, "Rakai_incpredictions_loess_inland_221116.csv")
file.incidence.samples.inland <- file.path(output.dir.incidence.estimation, "Rakai_incpredictions_samples_inland_221107.csv")
file.incidence.30com.samples.inland <- file.path(output.dir.incidence.estimation, "Rakai_incpredictions_samples_inland_221119.csv")
file.incidence.loess.samples.inland	<- file.path(output.dir.incidence.estimation, "Rakai_incpredictions_loess_samples_inland_221116.csv")

######################
# scripts/run_stan.R #
######################

# obtained in phylo_pipeline-src
file.pairs <- file.path(
    dir.zenodo.phyloproc,
    "pairsdata_toshare_d1_w11_netfrompairs_postponessrem.rds"
)
file.pairs.nonrefined <- file.path(
    dir.zenodo.phyloproc,
    "pairsdata_toshare_d1_w11_netfrompairs_seropairs_sensnoref.rds"
)


# obtained in src/ for analysis
file.eligible.count <- file.path(
    dir.zenodo.survprim,
    "RCCS_census_eligible_individuals_221116.csv"
)
file.participation <- file.path(
    dir.zenodo.survprim,
    "RCCS_participation_221208.csv"
)
file.prevalence.prop <- file.path(
    dir.zenodo.survproc,
    "RCCS_prevalence_estimates_221116.csv"
)

# obtained in surveillance_pipeline_src/ for analysis

file.treatment.cascade.prop.participants <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_participants_estimates_221208.csv"
)
file.treatment.cascade.prop.nonparticipants <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_nonparticipants_estimates_221208.csv"
)
file.treatment.cascade.prop.participants.samples <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_participants_posterior_samples_221208.rds"
)
file.treatment.cascade.prop.nonparticipants.samples <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_nonparticipants_posterior_samples_221208.rds"
)

# Sensitivity analysis
file.treatment.cascade.prop.participants.vl200 <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_participants_estimates_vl200_221208.csv"
)
file.treatment.cascade.prop.nonparticipants.vl200 <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_nonparticipants_estimates_vl200_221208.csv"
)
file.treatment.cascade.prop.participants.vl200.samples <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_participants_posterior_samples_vl200_221208.rds"
)
file.treatment.cascade.prop.nonparticipants.vl200.samples <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_nonparticipants_posterior_samples_vl200_221208.rds"
)

# obtained in misc/ for plots
file.unsuppressed.share <- file.path(dir.zenodo.survproc, "RCCS_unsuppressed_share_sex_221208.csv")
file.unsuppressed_rate_ratio <- file.path(dir.zenodo.survproc, "RCCS_unsuppressed_ratio_sex_221208.csv")
file.prevalence.share <- file.path(dir.zenodo.survproc, "RCCS_prevalence_share_sex_221116.csv")
file.unsuppressed_median_age <-file.path(dir.zenodo.survproc, "RCCS_unsuppressed_median_age_221208.csv")

file.unsuppressed.agegroup <- file.path(dir.zenodo.survproc, "RCCS_propunsuppressed_age_group_221208.csv")
file.prevalence.agegroup <- file.path(dir.zenodo.survproc, "RCCS_prevalence_age_group_221116.csv")

# sexual partnerships  rates
file.number.sexual.partnerships <- file.path(dir.zenodo.survprim, "age-age-group-est-cntcts-r15.rds")
file.sexual.partnerships.rates <- file.path(dir.zenodo.survprim, "inland_R015_cntcts_rate_1130b.rds")

# obtained in script/ for plots

#############
# Files used by the code within surveillance_pipeline_src/
#############

path.newly.registered.art <- file.path(
    dir.zenodo.survprim,
    "aggregated_newlyregistered_count_art_coverage.csv"
)
path.newly.registered.art.vl200 <- file.path(
    dir.zenodo.survprim,
    "aggregated_newlyregistered_count_art_coverage_vl200.csv"
)
path.participant.art <- file.path(
    dir.zenodo.survprim,
    "aggregated_participants_count_art_coverage.csv"
)
path.participant.art.vl200 <- file.path(
    dir.zenodo.survprim,
    "aggregated_participants_count_art_coverage_vl200.csv"
)
path.count.hivpositive <- file.path(
    dir.zenodo.survprim,
    "aggregated_count_hiv_positive.csv"
)
path.count.newly.unsupp <- file.path(
    dir.zenodo.survprim,
    "aggregated_newlyregistered_count_unsuppressed.csv"
)
path.count.newly.unsupp.vl200 <- file.path(
    dir.zenodo.survprim,
    "aggregated_newlyregistered_count_unsuppressed_vl200.csv"
)
path.count.unsupp <- file.path(
    dir.zenodo.survprim,
    "aggregated_participants_count_unsuppressed.csv"
)
path.count.unsupp.vl200 <- file.path(
    dir.zenodo.survprim,
    "aggregated_participants_count_unsuppressed_vl200.csv"
)
file.treatment.cascade <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_population_posterior_samples_221208.rds"
)
file.prevalence <- file.path(
    dir.zenodo.survproc,
    "RCCS_prevalence_posterior_sample_221116.rds"
)
file.unsuppressedviralload.newly <- file.path(
    dir.zenodo.survproc,
    "RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds"
)
file.selfreportedart.newly <- file.path(
    dir.zenodo.survproc,
    "RCCS_art_posterior_samples_newlyregistered_221208.rds"
)
file.spec.sens.art <- file.path(
    dir.zenodo.survprim,
    "sensitivity_specificity_art.csv"
)
file.unsuppressedviralload.newly.vl200 <- file.path(
    dir.zenodo.survproc,
    "RCCS_nonsuppressed_proportion_posterior_samples_vl_200_newlyregistered_221121.rds"
)
file.selfreportedart.newly.vl200 <- file.path(
    dir.zenodo.survproc,
    "RCCS_art_posterior_samples_newlyregistered_vl200_221208.rds"
)
file.spec.sens.art.vl200 <- file.path(
    dir.zenodo.survprim,
    "sensitivity_specificity_art_vl200.csv"
)

# Files used by get_treatment_cascade_*
file.unsuppressedviralload <- file.path(
    dir.zenodo.survproc,
    "RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds"
)
file.selfreportedart <- file.path(
    dir.zenodo.survproc,
    "RCCS_art_posterior_samples_221208.rds"
)
file.unsuppressedviralload.vl200 <- file.path(
    dir.zenodo.survproc,
    "RCCS_nonsuppressed_proportion_posterior_samples_vl_200_221121.rds"
)
file.selfreportedart.vl200 <- file.path(
    dir.zenodo.survproc,
    "RCCS_art_posterior_samples_vl200_221208.rds"
)
file.treatment.cascade.population <- file.path(
    dir.zenodo.survproc,
    "RCCS_treatment_cascade_population_estimates_221208.csv"
)
