# TODO: fill in path to zenodo dir and specify all the following files accordingly

config <- list(
    overwrite.existing.files  = FALSE
)

usr <- Sys.info()[["user"]]

# set your user and path here:
dir.zenodo <- data.table::fcase(
    usr == "andrea", "~/OneDrive/Imperial/shifting-dynamics-zenodo/",
    usr == "shozendan", "~/Imperial/phyloSI-RakaiAgeGender-data",
    usr == 'melodiemonod', '/Users/melodiemonod/Box Sync/2023/shifting-dynamics-zenodo',
    usr == 'mm3218', '/rds/general/user/mm3218/home/data/shifting-dynamics-zenodo/',
    usr == "Yu", "C:/Users/Yu/OneDrive - Imperial College London/shifting-dynamics-zenodo/",
    usr == 'alexb', '/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/shifting-dynamics-zenodo',
    usr == "your-user", "your-path-to-zenodo-dir",
    default = "")

zenodo.exists <- dir.exists(dir.zenodo)
if (zenodo.exists) {
    list.files(dir.zenodo, pattern = "")
  
  # data for deep sequence phylogenies
    dir.zenodo.phyloprim <- file.path(dir.zenodo, "data_deep_sequence_phylogenies")
    dir.zenodo.pairs <- file.path(dir.zenodo.phyloprim, "data_for_likely_transmission_pairs")
    dir.zenodo.pairs.phsc <- file.path(dir.zenodo.pairs, "phyloscanner-results")
    dir.zenodo.tsi <- file.path(dir.zenodo.phyloprim, "data_for_time_since_infection")
    dir.zenodo.tsi.phsc <- file.path(dir.zenodo.tsi, "phyloscanner-results")

    # processed output for deep sequence phylogenies
    dir.zenodo.phyloproc = file.path(dir.zenodo, "preprocessed_outputs_deep_sequence_phylogenies")

    # data, results, and processed output surveillance data (primary or processed_outputs)
    dir.zenodo.survprim <- file.path(dir.zenodo, "data_RCCS_surveillance")
    dir.zenodo.survfin <- file.path(dir.zenodo, "final_RCCS_surveillance")
    dir.zenodo.survproc <- file.path(dir.zenodo, "preprocessed_outputs_RCCS_surveillance")
    
    # data and results for incidence rate analysis
    dir.zenodo.dataincrate <- file.path(dir.zenodo, "data_incidence_rate")
    dir.zenodo.resincrate <- file.path(dir.zenodo, "final_incidence_rate")
    
    # results from sexual partnership analysis
    dir.zenodo.ressexpart <- file.path(dir.zenodo, "final_sexual_partnership")
    
    # results from transmission flows analysis
    dir.zenodo.transflow <- file.path(dir.zenodo, "final_transmission_flow")
    
    c(dir.zenodo.phyloprim,
      dir.zenodo.pairs,
      dir.zenodo.pairs.phsc,
      dir.zenodo.tsi,
      dir.zenodo.tsi.phsc,
      dir.zenodo.phyloproc,
      dir.zenodo.survprim,
      dir.zenodo.survproc
    ) %>% dir.exists() %>% all() %>% stopifnot()

}else{
  stop("Please specify path to zenodo data in config.R")
}

gitdir.outputs <- file.path(gitdir, "outputs")
gitdir.stan <- file.path(gitdir, "stan_models")
# gitdir.functions <- file.path(gitdir, "R") # DEPRECATED
gitdir.R <- file.path(gitdir, "R")
gitdir.R.flow <- file.path(gitdir.R, 'functions_transmission_flow')
gitdir.R.surv <- file.path(gitdir.R, 'functions_surveillance_pipeline')
gitdir.R.phylo <- file.path(gitdir.R, 'functions_phylo_pipeline')
gitdir.R.incid <- file.path(gitdir.R, 'functions_incidence_rate')
gitdir.R.conf <- file.path(gitdir.R, 'functions_confidential_data_pipeline')


######################
# confidential paths #
######################

# defaults implies dir.exists(indir.deepsequencedata) == FALSE
indir.deepsequencedata <- ""
indir.deepsequence_analyses <- ""

# set base paths
if (usr == "andrea" ) {
  indir.deepsequencedata <- "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live"
  indir.deepsequence_analyses <- "/home/andrea/HPC/project/ratmann_deepseq_analyses/live"
  indir.deepanalyses.xiaoyue <- "/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI"
}

if(usr == 'alexb') {
  indir.deepsequencedata <- '~/OneDrive - Imperial College London/PANGEA/ratmann_pangea_deepsequencedata/live'
  indir.deepsequence_analyses <- '~/OneDrive - Imperial College London/PANGEA/ratmann_deepseq_analyses/live'
  indir.deepanalyses.xiaoyue <- '~/OneDrive - Imperial College London/PANGEA/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
}

if(usr == "melodiemonod"){
  # main indir
  indir.deepsequencedata <- "~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/"
  indir.deepsequence_analyses <- "~/Box\ Sync/2021/ratmann_deepseq_analyses/live/"
  indir.deepanalyses.xiaoyue <- "~/Box\ Sync/2021/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI"
}

# define all other confidential paths
if( ! indir.deepsequencedata == "" ){

  indir.deepsequencedata.r151r18 <- file.path(indir.deepsequencedata, "RCCS_R15_R18")

  # some outdir
  output.dir.incidence.estimation <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'incidence_rate_inland')
  output.dir.incidence.estimation.30comms <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'incidence_rate_inland_30comms')
  
  # community keys
  file.community.keys <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS1519_UVRI", "community_names.csv")
  file.community.keys.aggregated <- file.path(indir.deepsequence_analyses,"PANGEA2_RCCS1519_UVRI", "community_id_index.csv")
  
  # geography
  file.Rakai_community_geography <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "Rakai_community_geography_R15.rda")
  
  # anonymisation keys
  file.anonymisation.keys <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS1519_UVRI", "important_anonymisation_keys_210119.csv")
  
  # hiv dataset round 15 to 18 
  file.path.hiv.1518 <- file.path(indir.deepsequencedata.r151r18, "HIV_R15_R18_VOIs_220129.csv")
  # hiv dataset round 6-14: 30 continuously surveyed communities
  file.path.hiv.614 <- file.path(indir.deepsequencedata, "RCCS_data_estimate_incidence_inland_R6_R18/220903", "hivincidence_1.dta")
  # hiv dataset round 9-14: all communities (created in create_hiv_dataset_from_quest_all_communities.R)
  file.hiv_R09_R14 = file.path(indir.deepsequencedata, "RCCS_R9_R14/HIV_R09_R14.csv")
  # hiv dataset combination of all rounds obtained in clean_quest_hiv
  file.path.hiv <- file.path(indir.deepsequencedata, "RCCS_data_estimate_incidence_inland_R6_R18/220903", "HIV_R6_R18_221129.csv")
  # hiv dataset round 19
  file.path.hiv_19 <- file.path(indir.deepsequencedata,"R019_VoIs","HIV_R019_VOIs_220607.csv")
  
  # quest dataset round 15 to 18 without 16
  file.path.quest.1518 <- file.path(indir.deepsequencedata.r151r18, "quest_R15_R18_VoIs_220129.csv")
  # quest dataset round 16
  file.path.quest.16 <- file.path(indir.deepsequencedata.r151r18, "quest_R15_R19_VoIs_Dec072022.csv")
  # quest dataset round 6-14: 30 continuously surveyed communities
  file.path.quest.614 <- file.path(indir.deepsequencedata, "RCCS_data_estimate_incidence_inland_R6_R18/220903", "quest_1.dta")
  # quest dataset round 9-14: all communities
  file.quest_R09_R14 = file.path(indir.deepsequencedata, "RCCS_R9_R14/quest_R09_R14.csv")
  # quest dataset combination of all rounds obtained in clean_quest_hiv 
  file.path.quest <- file.path(indir.deepsequencedata, "RCCS_data_estimate_incidence_inland_R6_R18/220903", "Quest_R6_R18_221208.csv")
  
  # flow dataset round <14: 30 continuously surveyed communities
  file.path.flow.614 <- file.path(indir.deepsequencedata, "RCCS_data_estimate_incidence_inland_R6_R18/220903", "verif_1.dta")
  # flow dataset round 9-14: all communities
  file.path.flow_914 <- file.path(indir.deepsequencedata, "RCCS_R9_R14", "FlowR09_R14.csv")
  # flow dataset round 15 to 18
  file.path.flow <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "FlowR15_R18_VoIs_221118.csv")
  # file.path.flow <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "FlowR15_R18_VoIs_220129.csv") # deprecated
  # flow dataset round 19
  file.path.flow_19 <- file.path(indir.deepsequencedata,"R019_VoIs","Flow_R019_VoIs_220607.csv")
  
  # Latest data from Rakai"s CCS (Joseph"s data from 2022-01-29)
  file.path.allhiv <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "All_HIVpcr_for_questR15_R18_220129.csv")
  
  # Latest update from Joseph concerning dates of infection
  file.path.update.first.positive <- file.path(indir.deepsequencedata.r151r18, "221128_requested_updated_serohistory.csv")
  
  # vl tests
  path.tests <- file.path(indir.deepsequencedata, "RCCS_R15_R20","all_participants_hivstatus_vl_220729.csv")
  
  # sequence files
  infile.sequence <- file.path(indir.deepsequencedata,"200422_PANGEA2_RCCSMRC_alignment.fasta")
  infile.seq.criteria <- file.path(indir.deepsequencedata,"PANGEA2_RCCS/221117_dct.rda")
  
  # characteristics sequenced participants 
  file.characteristics_sequenced_ind_R14_18 <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "participants_count_by_gender_loc_age", "characteristics_sequenced_ind_R14_18_221206.rds")
  file.characteristics_ever_sequenced <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "participants_count_by_gender_loc_age", "characteristics_ever_sequenced.rds")
  file.characteristics_sequenced_R14_18 <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "participants_count_by_gender_loc_age", "characteristics_sequenced_R14_18.rds")
  
  # metadata: Latest data from Rakai"s CCS (Kate"s data from 2022-03-08)
  file.path.metadata <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "Rakai_Pangea2_RCCS_Metadata__12Nov2019.csv")
  file.path.neuro.metadata <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "Pangea_Rakai_NeuroStudy_Metadata_11Dec2015.csv")
  
  # combined meta data (from preprocess meta_data)
  path.meta.confidential <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "Rakai_Pangea2_RCCS_Metadata_20221128.RData")
  
  # seroconvert individuals restrained to 30 continuously surveyed communities in round < 14 but not afterwards
  file.anonymised.id <- file.path(indir.deepsequencedata,"RCCS_data_estimate_incidence_inland_R6_R18","220903","anonymized_id_for_incidence_estimate_221129.csv")
  
  # seroconvert individuals not restrained to 30 continuously surveyed communities in round < 14
  file.anonymised.id.all_comm <- file.path(indir.deepsequencedata, "RCCS_R9_R14","RCCS_data_estimate_incidence_inland_R6_R18_230218","anonimized_id_for_incidence_estimate.csv")
  file.seroconvert.cohort.all.comm<- file.path(indir.deepsequencedata,"RCCS_R9_R14", "RCCS_data_estimate_incidence_inland_R6_R18_230218","seroconverter_cohort.rds")

  # sequence collection date
  path.collection.dates.confidential <- file.path(indir.deepsequencedata.r151r18, "sequences_collection_dates.rds")

  # phyloscanner samples
  path.selected.samples <- file.path(indir.deepanalyses.xiaoyue,"210120_RCCSUVRI_phscinput_samples.rds" )
  
  # pangea_db
  path.sdates.rccs <- file.path(indir.deepsequencedata, "PANGEA2_RCCS", "200316_pangea_db_sharing_extract_rakai.csv")
  path.sdates.mrc <- file.path(indir.deepsequencedata, "PANGEA2_MRC","200319_pangea_db_sharing_extract_mrc.csv")
}



################################
# PRIMARY DATA 
################################

# used in src/phylo_pipeline
path.meta.randomized <- file.path(dir.zenodo.survprim,"Rakai_Pangea2_RCCS_Metadata_randomized.RData" )
path.collection.dates.randomized <- file.path(dir.zenodo.survprim, "sequences_collection_dates_randomized.rds")
file.path.round.timeline <- file.path(dir.zenodo.survprim, "RCCS_round_timeline_220905.RData" )
path.chains.data <- file.path(dir.zenodo.phyloproc, "Rakai_phscnetworks_ruleo_sero.rda")
path.tsiestimates <- file.path(dir.zenodo.phyloproc, "aggregated_TSI.rds")

# used in src/incidence_rate
file.path.seroconverter_cohort <- file.path(dir.zenodo.dataincrate , "seroconverter_cohort_R6R19.rds")
file.path.seroconverter_cohort.30 <- file.path(dir.zenodo.dataincrate, "seroconverter_cohort_30comm_R6R19.rds")

# obtained in another repository sexual partnerships  rates
file.number.sexual.partnerships <- file.path(dir.zenodo.ressexpart, "age-age-group-est-cntcts-r15.rds")
file.sexual.partnerships.rates <- file.path(dir.zenodo.ressexpart, "inland_R015_cntcts_rate_1130b.rds")
file.age_dist_cntct_area <- file.path(dir.zenodo.ressexpart, "inland-R015_age-dist_cntct_area_1549.rds")
file.age_dist_ma_cntct_area <- file.path(dir.zenodo.ressexpart, "inland-R015_age-dist_ma_cntct_area_1549.rds")

# final flow data for Extended data fig 7
file.unsuppressed.share.fig <- file.path(dir.zenodo.survproc, "RCCS_unsuppressed_share_sex_221208_fig.csv")

file.expected_contribution_age_source <- file.path(dir.zenodo.transflow, "gp_221201d-central3-output-log_lambda_latentby_direction_round_age_transmission.sourcestandardisedby_direction_round.rds")
file.expected_contribution_sliced_age_source <- file.path(dir.zenodo.transflow, "gp_221201d-central3-output-log_lambda_latentby_direction_round_age_transmission.source_age_infection.recipientstandardisedby_direction_round.rds")

###################################################
#  OUTPUTS GENERATED IN scripts_for_confidential_data #
###################################################

file.eligible.count <- file.path(dir.zenodo.survprim, "RCCS_census_eligible_individuals_221116.csv")
file.participation <- file.path(dir.zenodo.survprim, "RCCS_participation_221208.csv")

path.count.hivpositive <- file.path(dir.zenodo.survprim,"aggregated_count_hiv_positive.csv")

path.newly.registered.art <- file.path(dir.zenodo.survprim,"aggregated_newlyregistered_count_art_coverage.csv")
path.newly.registered.art.vl200 <- file.path(dir.zenodo.survprim,"aggregated_newlyregistered_count_art_coverage_vl200.csv")
path.participant.art <- file.path(dir.zenodo.survprim,"aggregated_participants_count_art_coverage.csv")
path.participant.art.vl200 <- file.path(dir.zenodo.survprim,"aggregated_participants_count_art_coverage_vl200.csv")

path.count.newly.unsupp <- file.path(dir.zenodo.survprim,"aggregated_newlyregistered_count_unsuppressed.csv")
path.count.newly.unsupp.vl200 <- file.path(dir.zenodo.survprim,"aggregated_newlyregistered_count_unsuppressed_vl200.csv")
path.count.unsupp <- file.path(dir.zenodo.survprim,"aggregated_participants_count_unsuppressed.csv")
path.count.unsupp.vl200 <- file.path(dir.zenodo.survprim,"aggregated_participants_count_unsuppressed_vl200.csv")

file.spec.sens.art <- file.path(dir.zenodo.survprim,"sensitivity_specificity_art.csv")
file.spec.sens.art.vl200 <- file.path(dir.zenodo.survprim,"sensitivity_specificity_art_vl200.csv")

############################################
# OUTPUTS GENERATED IN confidential_data_src
############################################

file.pairs <- file.path( dir.zenodo.phyloproc, "pairsdata_toshare_d1_w11_netfrompairs_postponessrem.rds")
file.pairs.nonrefined <- file.path( dir.zenodo.phyloproc, "pairsdata_toshare_d1_w11_netfrompairs_seropairs_sensnoref.rds")

###########################################
#  OUTPUTS GENERATED IN src/incidence_rate #
###########################################

file.incidence.fits	<- file.path(dir.zenodo.resincrate, "fit_incidence_rates_221109.RData")
file.incidence.30com.fits	<- file.path(dir.zenodo.resincrate, "fit_incidence_rates_221109.RData")

file.incidence.inland	<- file.path(dir.zenodo.resincrate, "Rakai_incpredictions_inland_221107.csv")
file.incidence.samples.inland <- file.path(dir.zenodo.resincrate, "Rakai_incpredictions_samples_inland_221107.csv")

file.incidence.loess.inland	<- file.path(dir.zenodo.resincrate, "Rakai_incpredictions_loess_inland_221116.csv")
file.incidence.loess.samples.inland	<- file.path(dir.zenodo.resincrate, "Rakai_incpredictions_loess_samples_inland_221116.csv")

file.incidence.30com.inland	<- file.path(dir.zenodo.resincrate, "Rakai_incpredictions_inland_221119.csv")
file.incidence.30com.samples.inland <- file.path(dir.zenodo.resincrate, "Rakai_incpredictions_samples_inland_221119.csv")

##########################################
# OUTPUTS GENERATE IN src/surveillance_pipeline #
##########################################

file.prevalence.prop <- file.path(dir.zenodo.survfin,"RCCS_prevalence_estimates_221116.csv")
file.prevalence <- file.path(dir.zenodo.survfin,"RCCS_prevalence_posterior_sample_221116.rds")

file.unsuppressedviralload <- file.path(dir.zenodo.survfin,"RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds")
file.unsuppressedviralload.newly <- file.path(dir.zenodo.survfin,"RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds")

file.unsuppressedviralload.vl200 <- file.path(dir.zenodo.survfin,"RCCS_nonsuppressed_proportion_posterior_samples_vl_200_221121.rds")
file.unsuppressedviralload.newly.vl200 <- file.path(dir.zenodo.survfin,"RCCS_nonsuppressed_proportion_posterior_samples_vl_200_newlyregistered_221121.rds")

file.selfreportedart <- file.path(dir.zenodo.survfin,"RCCS_art_posterior_samples_221208.rds")
file.selfreportedart.newly <- file.path(dir.zenodo.survfin,"RCCS_art_posterior_samples_newlyregistered_221208.rds")

file.selfreportedart.vl200 <- file.path(dir.zenodo.survfin,"RCCS_art_posterior_samples_vl200_221208.rds")
file.selfreportedart.newly.vl200 <- file.path(dir.zenodo.survfin,"RCCS_art_posterior_samples_newlyregistered_vl200_221208.rds")

file.treatment.cascade.prop.participants <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_participants_estimates_221208.csv")
file.treatment.cascade.prop.nonparticipants <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_nonparticipants_estimates_221208.csv")
file.treatment.cascade.prop.participants.samples <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_participants_posterior_samples_221208.rds")
file.treatment.cascade.prop.nonparticipants.samples <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_nonparticipants_posterior_samples_221208.rds")

file.treatment.cascade.prop.participants.vl200 <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_participants_estimates_vl200_221208.csv")
file.treatment.cascade.prop.nonparticipants.vl200 <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_nonparticipants_estimates_vl200_221208.csv")
file.treatment.cascade.prop.participants.vl200.samples <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_participants_posterior_samples_vl200_221208.rds")
file.treatment.cascade.prop.nonparticipants.vl200.samples <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_nonparticipants_posterior_samples_vl200_221208.rds")

file.treatment.cascade.population <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_population_estimates_221208.csv")
file.treatment.cascade <- file.path(dir.zenodo.survfin,"RCCS_treatment_cascade_population_posterior_samples_221208.rds")

file.unsuppressed.share <- file.path(dir.zenodo.survproc, "RCCS_unsuppressed_share_sex_221208.csv")
file.unsuppressed_rate_ratio <- file.path(dir.zenodo.survproc, "RCCS_unsuppressed_ratio_sex_221208.csv")
file.prevalence.share <- file.path(dir.zenodo.survproc, "RCCS_prevalence_share_sex_221116.csv")
file.unsuppressed_median_age <-file.path(dir.zenodo.survproc, "RCCS_unsuppressed_median_age_221208.csv")
file.unsuppressed.agegroup <- file.path(dir.zenodo.survproc, "RCCS_propunsuppressed_age_group_221208.csv")
file.prevalence.agegroup <- file.path(dir.zenodo.survproc, "RCCS_prevalence_age_group_221116.csv")

##########################################
# OUTPUTS GENERATE IN transmission/flows #
##########################################

file.lambda.est <- file.path(dir.zenodo.transflow, 'gp_221201d-central3-output-log_lambda_latentby_direction_round_age_transmission.source_age_infection.recipientstandardisedby_direction_round.rds')

##########################################
# PATH TO STAN MODELS #
##########################################

path_stan_binomialgp <- file.path(gitdir.stan, "binomial_gp.stan")
path_binomialgp_model_config <- file.path(gitdir.stan, "binomial_gp_config.yml")

