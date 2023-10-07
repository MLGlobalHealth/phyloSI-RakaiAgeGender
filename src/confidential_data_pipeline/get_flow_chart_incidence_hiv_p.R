library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library(haven)
library(here)
library(tidyverse)
library(tidyr)
library(mgcv)
library(zoo)
library(remotes)
library(metR)
library(ggpubr)

# directory of the repository
gitdir <- here()

# load file paths
source(file.path(gitdir, 'config.R'))

file.exists(c(
  file.community.keys ,
  file.eligible.count,
  path.tests,
  file.characteristics_sequenced_ind_R14_18,
  file.path.hiv,
  file.path.quest,
  file.path.metadata))  |> all() |> stopifnot()

c(  file.path.flow.614,
    file.path.hiv.614,
    file.hiv_R09_R14,
    file.path.flow,
    file.path.flow_914,
    file.path.hiv.1518,
    file.path.flow_19,
    file.path.hiv_19) |> file.exists() |> all() |> stopifnot()

#
# DO NOT CHANGE. 
# COPY-PASTE CODE FROM OTHER SCRIPTS TO OBTAIN INCIDENCE COHORT AND PARTICIPANTS
#
#

if(1){
  
  # load files
  community.keys <- as.data.table(read.csv(file.community.keys))
  
  # rounds of interest
  df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 10:18)), 
                    data.table(COMM = 'fishing', ROUND = paste0('R0', c(15, '15S', 16:18))))
  
  # community
  community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
  
  
  ###############################################
  
  # GET META DATA 
  
  ###############################################
  
  #
  # Meta data
  #
  
  meta_data <- as.data.table(read.csv(file.path.metadata)) #additional meta_data
  
  # find age
  meta_data[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
  meta_data[, AGEYRS := round(lubridate::time_length(difftime(sample_date, date_birth),"years"))]
  meta_data[is.na(AGEYRS), AGEYRS := round(lubridate::time_length(difftime(firstposvd, date_birth),"years"))]
  meta_data[is.na(AGEYRS)]
  
  # restrict age
  meta_data <- meta_data[AGEYRS > 14 & AGEYRS < 50]
  
  # find community
  meta_data[, COMM := 'inland']
  meta_data[LakeVictoria_FishingCommunity == 'yes', COMM := 'fishing']
  
  # find hiv status
  meta_data[, HIV := ifelse(is.na(firstposvd), 'N', 'P')]
  
  # find art use
  meta_data[, ART := artslfuse == 'yes']
  
  # subset round 14 to 30 continuously surveyed communities
  meta_data <- meta_data[!(round=='14' & !COMM %in% c(1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 89, 94, 106, 107, 108, 120, 391, 602, 754))]
  
  # keep variable of interest
  meta_data[, round := paste0('R0', round)]
  colnames(meta_data) <- toupper(colnames(meta_data))
  meta_data <- meta_data[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV, ART)]
  
  # set 15.1 to be 15S
  meta_data[ROUND == 'R015.1', ROUND := 'R015S']
  
  
  #
  # Quest
  #
  
  quest <- as.data.table(read.csv(file.path.quest))
  
  # keep variable of interest
  rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, arvmed, cuarvmed)]
  
  # find  community
  rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')
  
  # to upper
  colnames(rinc) <- toupper(colnames(rinc))
  
  # restrict age
  rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]
  
  # get ART status
  rinc[!ROUND %in% c('R016', 'R017', 'R018'), ART := ARVMED ==1]
  rinc[!ROUND %in% c('R016', 'R017', 'R018') & is.na(ARVMED), ART := F]
  rinc[ROUND == 'R016', ART := ARVMED ==1 | CUARVMED ==1]
  rinc[ROUND == 'R016' & (is.na(ARVMED) | is.na(CUARVMED)), ART := F]
  rinc[ROUND %in% c('R017', 'R018'), ART := CUARVMED ==1]
  rinc[ROUND %in% c('R017', 'R018') & is.na(CUARVMED), ART := F]
  
  # art was not reported in round 10
  rinc[ROUND == 'R010', ART := NA]
  
  # add meta data from Kate
  tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], rinc[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
  tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
  rinc <- rbind(rinc[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)])
  
  # SET ROUND 15S IN INLAND AS 15
  rincp <- copy(rinc)
  rincp[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
  rincp[ROUND =='R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
  rincp <- rincp[!(ROUND =='R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]
  
  
  #
  # HIV
  #
  
  hiv <- as.data.table(read.csv(file.path.hiv))
  hiv[, round := gsub(' ', '', round)] # remove space in string
  
  # get hiv status
  rhiv <- hiv[, .(study_id, round, hiv)]
  rhiv[, round := gsub(" ", '', round, fixed = T)]
  colnames(rhiv) <- toupper(colnames(rhiv))
  
  # add meta data from Joseph 
  hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))
  
  # add meta data from Kate
  tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], hivs[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
  tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
  hivs <- rbind(hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)])
  
  # SET ROUND 15S IN INLAND AS 15
  hivs[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
  hivs[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
  hivs <- hivs[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]
  
  # restric age
  hivs <- hivs[AGEYRS > 14 & AGEYRS < 50]
  
  
  #####################################
  
  #   FIND INCIDENCE COHORT
  
  ####################################

  # function
  source(file.path(gitdir.R.incid, "incidence_rate_estimation_functions.R"))
  
  # utils
  rounds_group_1 <- c("R006","R007", "R008", "R009", "R010",
                      "R011", "R012", "R013", "R014")
  rounds_group_2 <- c("R015", "R015S", "R016", "R017", "R018")
  rounds_group_3 <- c("R019")
  
  # make df round
  df_round <- make_df_round(rounds_group_1, rounds_group_2, rounds_group_3)
  
  rounds_numeric_group_1 <- df_round[visit %in% rounds_group_1, round_numeric]
  rounds_numeric_group_2 <- df_round[visit %in% rounds_group_2, round_numeric]
  rounds_numeric_group_3 <- df_round[visit == rounds_group_3, round_numeric]
  
  
  ##############
  
  # LOAD DATA #
  
  ##############
  
  # load HIV status data
  hivstatus_vlcopies_raw <- read_hiv_data(file.path.hiv.614,
                                          file.path.hiv.1518,
                                          file.path.hiv_19)
  
  # load HIV verification data set
  verif_raw <- read_flow_data(file.path.flow.614,
                              file.path.flow,
                              file.path.flow_19)
  
  
  ################
  
  # FORMAT DATA  #
  
  ################
  
  # format and clean hiv data
  hivstatus_vlcopies_1 <- format_hiv_data(hivstatus_vlcopies_raw)
  rm(hivstatus_vlcopies_raw)
  
  # format and clean flow data
  verif_1 <- format_flow_data(verif_raw)

  
  ####################################
  
  # RESTRICT ANALAYSIS BY CRITERIA  #
  
  ####################################
  
  # Use Verification data to restrict analysis based on following criteria
  ## 1. Age 15-49
  ## 2. residency status: resident==1
  ## 3. Resides in inland communities
  
  # central analysis
  verif_el <- verif_1[resident == 1 &
                        ageyrs >= 15 &
                        ageyrs <= 49 &
                        !comm_num %in% c(38, 770, 771, 774), ]
  verif_el[, table(round)]
  
  # only include HIV data as determined
  verified_visit_ids <- unique(verif_el[, visit_id])
  hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[visit_id %in% verified_visit_ids, ]
  table(hivstatus_vlcopies_1$round)
  
  
  ###################################################
  
  # PREPARE HIV STATUS FOR EVERY AGE AND YEAR #
  
  ###################################################
  
  hivstatus_vlcopies_1_inc <- prepare_hiv_status(hivstatus_vlcopies_1, verif_el)
  
  ###############################
  
  # FIND SEROCONVERSION STATUS #
  
  ###############################
  
  status_df <- find_seroconvert_status(hivstatus_vlcopies_1_inc)
  
  #
  # subset serially negative or seroconvert 
  #
  
  # number of followed-up participants by age group, sex and round
  part <- as.data.table(hivstatus_vlcopies_1_inc)
  part[, age := floor(as.numeric(date-birthdate)/365.25)]
  part[age < 15, age := 15]; part[age > 49, age := 49] # discrepancies in age
  part <- unique(part[, list(age = min(age), sex = unique(sex), 
                             hivstatus_imputed = min(hivstatus_imputed),
                             missed_visit = is.na(visit_id)), by = c('research_id', 'round')])
  df_ne <- unique(part[, .(research_id, hivstatus_imputed, round)])
  df_ne[, index_round_after_sero := 0] # round after seroconversion
  df_ne[hivstatus_imputed == 1, index_round_after_sero := 1:length(hivstatus_imputed), by = 'research_id']
  df_ne[, index_round := 1:length(round), by = 'research_id'] # index round observed
  df_ne[, seroconvert_at_first_round := any(index_round == 1 & hivstatus_imputed == 1), by = 'research_id']
  part <- merge(part, df_ne, by = c('research_id', 'hivstatus_imputed', 'round'))

  # flag individuals who have mistake in their date of first infection, i.e., first pos < last neg
  ids_rm <- c('A007828', 'B009659', 'B012398', 'C010218', 'D011727', 'D022911', 'E026988', 'E057461', 'F009328', 'F009865', 
              'F023120', 'F038827', 'F055446', 'G039734', 'G046253', 'H038942', 'J023937')
  part[, erroneous_first_pos_date := F]
  part[research_id %in% ids_rm, erroneous_first_pos_date := T]
  
  # keep included indiv
  part1 <- part[seroconvert_at_first_round == F] # remove hivp at first round 
  part2 <- part1[index_round_after_sero <= 1]   #keep rounds before and on seroconversion
  part2 <- part2[erroneous_first_pos_date == F] 
  df_age_aggregated <- data.table(age = 15:49, age_group = c(rep('15-24', 10),  rep("25-34", 10), rep("35-49", 15)))
  part2 <- merge(part2, df_age_aggregated, by = 'age')
  part_all <- part2[, list(N = length(unique(research_id))), by = c('round', 'age_group', 'sex')]
  part_s <- part2[, list(N = length(unique(research_id))), by = c('round', 'sex')]
  part_r <- part2[, list(N = length(unique(research_id))), by = c('round')]
}

# get table
tmp <- part[hivstatus_imputed == 1 &  missed_visit == F]
tmp <- tmp[, list(research_id_before = length(unique(research_id)),
                  seroconvert_at_first_round = sum(seroconvert_at_first_round == T), 
                  after_seroversion = sum((index_round_after_sero > 1)[seroconvert_at_first_round == F]),
                  erroneous_first_pos_date = length(unique(research_id[erroneous_first_pos_date == T & seroconvert_at_first_round == F & index_round_after_sero <= 1]))), 
                  by = 'round']

tmp2 <- part2[hivstatus_imputed == 1 &  missed_visit == F]
tmp2 <- tmp2[, list(research_id_after = length(unique(research_id))), by = 'round']

tmp <- merge(tmp, tmp2, by = 'round')
tmp[research_id_after != research_id_before - seroconvert_at_first_round - after_seroversion - erroneous_first_pos_date]
tmp <- merge(tmp, df_round, by.x = 'round', by.y = 'round_numeric' )

hivsp <- hivs[HIV == 'P' & COMM == 'inland']
hivsp <- merge(hivsp, unique(verif_1[, .(visit, research_id, resident)]), 
               by.x = c('STUDY_ID', 'ROUND'), by.y = c('research_id', 'visit'), all.x = T)
hivsp[is.na(resident), resident := 10]
# hivsp[FROM_KATE_DATA == T, table(resident)]

tmp3 <- hivsp[, list(ALL_HIV = length(unique(STUDY_ID)),
                     ALL_HIV_RESIDENT1 = length(unique(STUDY_ID[resident == 1])),
                     ALL_HIV_RESIDENTNOT1 = length(unique(STUDY_ID[resident != 1]))
                     # ALL_HIV_FROM_KATE = length(unique(STUDY_ID[FROM_KATE_DATA == T ]))
), by = 'ROUND']
tmp3[, resident_does_not_equal_1 := ALL_HIV_RESIDENTNOT1]
# tmp3[, from_Kate_analysis := ALL_HIV_FROM_KATE]

tmp <- merge(tmp, tmp3, by.x = 'visit', by.y = 'ROUND')
setnames(tmp, 'research_id_after', 'ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS')
tmp[ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS != ALL_HIV - seroconvert_at_first_round - after_seroversion - resident_does_not_equal_1 - erroneous_first_pos_date]

fin <- tmp[, .(visit, ALL_HIV, seroconvert_at_first_round, after_seroversion, resident_does_not_equal_1, erroneous_first_pos_date, ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS)]
setnames(fin, 'visit', 'ROUND')
fin <- fin[ROUND %in% paste0('R0', 10:18)]

# save
file.name <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "participants_count_by_gender_loc_age","flow_chart_hivp_included_in_incidence_analysis_230926.csv")
write.csv(fin, file = file.name, row.names = F)

