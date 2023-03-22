# TODO: fill in path to zenodo dir and specify all the following files accordingly
zenodo.dir <- "" # TODO

gitdir.data <- file.path(gitdir, 'data')
gitdir.fit <- file.path(gitdir, 'fit')
gitdir.functions <- file.path(gitdir, 'functions')
gitdir.misc <- file.path(gitdir, 'misc')

######################
# confidential paths #
######################

usr <- Sys.info()[['user']]

# defaults implies dir.exists(indir.deepsequencedata) == FALSE
indir.deepsequencedata <- ''
indir.deepsequence_analyses <- ''

if(usr == 'andrea')
{
    indir.deepsequencedata <- "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live"
    indir.deepsequence_analyses <- "/home/andrea/HPC/project/ratmann_deepseq_analyses/live"
}
if(usr == "melodie")
{
    indir.deepsequencedata <- ''
    indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
}
if(usr== "ratmann")
{
    indir.deepsequencedata <- ''
    indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
}

if(dif.exists(indir.deepsequencedata))
{
    indir.deepsequencedata.r151r18 <- file.path(indir.deepsequencedata, "RCCS_R15_R18")
    file.community.keys <- file.path(indir.deepsequence_analyses.r151r18, 'community_names.csv')
}


######################
# scripts/run_stan.R #
######################

# incidence predictions paths
file.incidence.inland	<- file.path(gitdir.data, "Rakai_incpredictions_inland_221107.csv")
file.incidence.30com.inland	<- file.path(gitdir.data, "Rakai_incpredictions_inland_221119.csv")
file.incidence.loess.inland	<- file.path(gitdir.data, "Rakai_incpredictions_loess_inland_221116.csv")

# obtained in src/ for analysis
file.path.round.timeline <- file.path(gitdir.data, 'RCCS_round_timeline_220905.RData')
file.eligible.count <- file.path(gitdir.data, 'RCCS_census_eligible_individuals_221116.csv')
file.participation <- file.path(gitdir.data, 'RCCS_participation_221208.csv')
file.prevalence.prop <- file.path(gitdir.fit, 'RCCS_prevalence_estimates_221116.csv')

# obtained in misc/ for analysis
file.pairs <- file.path(gitdir.data, 'pairsdata_toshare_d1_w11_netfrompairs_postponessrem.rds')
file.pairs.nonrefined <- file.path(gitdir.data, 'pairsdata_toshare_d1_w11_netfrompairs_seropairs_sensnoref.rds')

file.treatment.cascade.prop.participants <- file.path(gitdir.fit, "RCCS_treatment_cascade_participants_estimates_221208.csv")
file.treatment.cascade.prop.nonparticipants <- file.path(gitdir.fit, "RCCS_treatment_cascade_nonparticipants_estimates_221208.csv")
file.treatment.cascade.prop.participants.samples <- file.path(gitdir.fit, 'RCCS_treatment_cascade_participants_posterior_samples_221208.rds')
file.treatment.cascade.prop.nonparticipants.samples <- file.path(gitdir.fit, 'RCCS_treatment_cascade_nonparticipants_posterior_samples_221208.rds')

file.treatment.cascade.prop.participants.vl200 <- file.path(gitdir.fit, "RCCS_treatment_cascade_participants_estimates_vl200_221208.csv")
file.treatment.cascade.prop.nonparticipants.vl200 <- file.path(gitdir.fit, "RCCS_treatment_cascade_nonparticipants_estimates_vl200_221208.csv")
file.treatment.cascade.prop.participants.vl200.samples <- file.path(gitdir.fit, 'RCCS_treatment_cascade_participants_posterior_samples_vl200_221208.rds')
file.treatment.cascade.prop.nonparticipants.vl200.samples <- file.path(gitdir.fit, 'RCCS_treatment_cascade_nonparticipants_posterior_samples_vl200_221208.rds')

# obtained in misc/ for plots
file.unsuppressed.share <- file.path(gitdir.fit, 'RCCS_unsuppressed_share_sex_221208.csv')
file.unsuppressed_rate_ratio <- file.path(gitdir.fit, 'RCCS_unsuppressed_ratio_sex_221208.csv')
file.prevalence.share <- file.path(gitdir.fit, 'RCCS_prevalence_share_sex_221116.csv')
file.unsuppressed_median_age <-file.path(gitdir.fit, 'RCCS_unsuppressed_median_age_221208.csv')

# sexual partnerships  rates
file.number.sexual.partnerships <- file.path(gitdir.data, 'age-age-group-est-cntcts-r15.rds')
file.sexual.partnerships.rates <- file.path(gitdir.data, 'inland_R015_cntcts_rate_1130b.rds')

# obtained in script/ for plots
file.incidence.samples.inland <- file.path(gitdir.data, "Rakai_incpredictions_samples_inland_221107.csv")
file.incidence.30com.samples.inland <- file.path(gitdir.data, "Rakai_incpredictions_samples_inland_221119.csv")
file.incidence.loess.samples.inland	<- file.path(gitdir.data, "Rakai_incpredictions_loess_samples_inland_221116.csv")


############################################
# scripts/run_incidence_rates_estimation.R #
############################################

# Main analysis: sero conversion cohort status (all communities)
file.path.seroconverter_cohort <- file.path(gitdir.data, 'seroconverter_cohort_R6R19.rds')

# Sensitivity analyses:
# sero conversion cohort status (continuously surveyed communities)
# file.path.seroconverter_cohort <- file.path(gitdir.data, 'seroconverter_cohort_30comm_R6R19.rds')

# sero conversion cohort status (continuously surveyed communities < round 9 and all communities afterwards)
# file.path.seroconverter_cohort <- file.path(indir.deepsequencedata, 'RCCS_R9_R14/RCCS_data_estimate_incidence_inland_R6_R18_230218', 'seroconverter_cohort.rds')


#############
# misc/.... #
#############

path.newly.registered.art <- file.path(gitdir.data, 'aggregated_newlyregistered_count_art_coverage.csv')
path.newly.registered.art.vl200 <- file.path(gitdir.data, 'aggregated_newlyregistered_count_art_coverage_vl200.csv')

path.participant.art <- file.path(gitdir.data, 'aggregated_participants_count_art_coverage.csv')
path.participant.art.vl200 <- file.path(gitdir.data, 'aggregated_participants_count_art_coverage_vl200.csv')

path.count.hivpositive <- file.path(gitdir.data, 'aggregated_count_hiv_positive.csv')

path.count.newly.unsupp <- file.path(gitdir.data, 'aggregated_newlyregistered_count_unsuppressed.csv')
path.count.newly.unsupp.vl200 <- file.path(gitdir.data, 'aggregated_newlyregistered_count_unsuppressed_vl200.csv')

path.count.unsupp <- file.path(gitdir.data, 'aggregated_participants_count_unsuppressed.csv')
path.count.unsupp.vl200 <- file.path(gitdir.data, 'aggregated_participants_count_unsuppressed_vl200.csv')

file.treatment.cascade <- file.path(gitdir.fit, 'RCCS_treatment_cascade_population_posterior_samples_221208.rds')
file.prevalence <- file.path(gitdir.fit, 'RCCS_prevalence_posterior_sample_221116.rds')

file.unsuppressedviralload.newly <- file.path(gitdir.fit, 'RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds')
file.selfreportedart.newly <- file.path(gitdir.fit, 'RCCS_art_posterior_samples_newlyregistered_221208.rds')
file.spec.sens.art = file.path(gitdir.data, 'sensitivity_specificity_art.csv')

file.unsuppressedviralload.newly.vl200 <- file.path(gitdir, 'fit', 'RCCS_nonsuppressed_proportion_posterior_samples_vl_200_newlyregistered_221121.rds')
file.selfreportedart.newly.vl200 <- file.path(gitdir.fit, 'RCCS_art_posterior_samples_newlyregistered_vl200_221208.rds')
file.spec.sens.art.vl200 = file.path(gitdir.data, 'sensitivity_specificity_art_vl200.csv')

# get_treatment_cascade*
file.unsuppressedviralload <- file.path(gitdir.fit,'RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds')
file.selfreportedart <- file.path(gitdir.fit,'RCCS_art_posterior_samples_221208.rds')
file.unsuppressedviralload.vl200 <- file.path(gitdir.fit,'RCCS_nonsuppressed_proportion_posterior_samples_vl_200_221121.rds')
file.selfreportedart.vl200 <- file.path(gitdir.fit,'RCCS_art_posterior_samples_vl200_221208.rds')

# posterior samples non-participants
file.unsuppressedviralload.np <- file.path(gitdir.fit,'RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds')
file.selfreportedart.np <- file.path(gitdir.fit,'RCCS_art_posterior_samples_newlyregistered_221208.rds')
