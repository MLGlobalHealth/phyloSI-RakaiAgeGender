# Rakai data analysis for EMOD inputs
# August 12, 2022
# Adam Akullian/Kate Grabowski/Melodie Monod;
# Shozen Dan (Minor edits, March 3, 2023)

library(tidyverse)
library(data.table)
library(haven)
library(tidyr)
library(mgcv)
library(zoo)
library(remotes)
library(metR)
library(lubridate)
library(ggpubr)

# set to directory
indir <- "~/git/phyloflows"

# set up path to data
indir_deepsequencedata <- file.path("~/Box\ Sync",
                                    "2019",
                                    "ratmann_pangea_deepsequencedata",
                                    "live")
indir_deepsequence_analyses <- file.path("~/Box\ Sync",
                                         "2021",
                                         "ratmann_deepseq_analyses",
                                         "live")

# set path to save results
outdir <- file.path(indir_deepsequence_analyses,
                    "PANGEA2_RCCS",
                    "incidence_rate_inland")

# round 6 to 18: only use round 6 to 14
path_flow_614 <- file.path(indir_deepsequencedata,
                           "RCCS_data_estimate_incidence_inland_R6_R18",
                           "220903",
                           "verif_1.dta")
path_hiv_614 <- file.path(indir_deepsequencedata,
                               "RCCS_data_estimate_incidence_inland_R6_R18",
                               "220903",
                               "hivincidence_1.dta")

# round 15 to 18
path_flow_1518 <- file.path(indir_deepsequencedata,
                            "RCCS_R15_R18",
                            "FlowR15_R18_VoIs_220129.csv")
path_hiv_1518 <- file.path(indir_deepsequencedata,
                           "RCCS_R15_R18",
                           "HIV_R15_R18_VOIs_220129.csv")

# round 19
path_flow_19 <- file.path(indir_deepsequencedata,
                               "R019_VoIs",
                               "Flow_R019_VoIs_220607.csv")
path_hiv_19 <- file.path(indir_deepsequencedata,
                         "R019_VoIs",
                         "HIV_R019_VOIs_220607.csv")

# function
source(file.path(indir, "functions", "incidence_rate_estimation_functions.R"))

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
hivstatus_vlcopies_raw <- read_hiv_data(path_hiv_614,
                                        path_hiv_1518,
                                        path_hiv_19)

# load HIV verification data set
verif_raw <- read_flow_data(path_flow_614,
                            path_flow_1518,
                            path_flow_19)


################

# FORMAT DATA  #

################

# format and clean hiv data
hivstatus_vlcopies_1 <- format_hiv_data(hivstatus_vlcopies_raw)
rm(hivstatus_vlcopies_raw)

# format and clean flow data
verif_1 <- format_flow_data(verif_raw)
rm(verif_raw)


################

# FIND ANONIMIZED ID  #

################

# anonimized id
anonimized_id <- data.table(research_id = verif_1[, sort(unique(research_id))])
anonimized_id <- anonimized_id[research_id != ""]

# randomly assign new id
set.seed(12)
anonimized_id[, anonimized_id := paste0("RkINC-",
                                        sample(str_pad(1:nrow(anonimized_id),
                                               width = 5,
                                               pad = "0"),
                                        replace = FALSE))]


####################################

# RESTRICT ANALAYSIS BY CRITERIA  #

####################################

# Use Verification data to restrict analysis based on following criteria
## 1. Age 15-49
## 2. residency status: resident==1
## 3. Resides in inland communities

# Subset data based inclusion criteria
inland_comms <- c("1", "2", "4", "5", "6", "7", "8",
                  "16", "19", "22", "24", "29", "33",
                  "34", "40", "56", "57", "58", "62",
                  "74", "77", "89", "94", "106", "107",
                  "108", "120", "391", "602", "754")

# sensitivity analysis
# verif_el <- verif_1[resident == 1 &
#                     ageyrs >= 15 &
#                     ageyrs <= 49 &
#                     comm_num %in% inland_comms, ]

# central analysis
verif_el <- verif_1[resident == 1 &
                    ageyrs >= 15 &
                    ageyrs <= 49 &
                    !comm_num %in% c(38, 770, 771, 774), ]
verif_el[, table(round)]
rm(verif_1)

# only include HIV data as determined
verified_visit_ids <- unique(verif_el[, visit_id])
hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[visit_id %in% verified_visit_ids, ]
table(hivstatus_vlcopies_1$round)


###################################################

# PREPARE HIV STATUS FOR EVERY AGE AND YEAR #

###################################################

hivstatus_vlcopies_1_inc <- prepare_hiv_status(hivstatus_vlcopies_1, verif_el)
rm(hivstatus_vlcopies_1)
rm(verif_el)


###############################

# FIND SEROCONVERSION STATUS #

###############################

status_df <- find_seroconvert_status(hivstatus_vlcopies_1_inc)

################

# FIND STATISTICS COHORT

################

stats <- save_statistics_incidence_cohort(hivstatus_vlcopies_1_inc, status_df)


###############################

# FIND SEROCONVERSION DATE #

###############################

N <- 50
seed <- 12

seroconverter_cohort_list <- vector(mode = "list", length = N)

set.seed(seed)
for (i in 1:N) {
  cat("\niteration", i)
  # generate random date of infection and find perosn year
  seroconverter_cohort <- find_seroconvert_cohort(status_df,
                                                  hivstatus_vlcopies_1_inc)

  # add anonimized id
  seroconverter_cohort <- merge(seroconverter_cohort,
                                anonimized_id,
                                by = "research_id")

  # keep only necessary covariate
  seroconverter_cohort <- seroconverter_cohort[, .(anonimized_id,
                                                   age,
                                                   sex,
                                                   number_missing_visits,
                                                   round, py, hivinc)]
  # prepare
  seroconverter_cohort$iterations <- i
  seroconverter_cohort_list[[i]] <- seroconverter_cohort
}


################

# SAVE

################

# anonimized id
file <- file.path(indir_deepsequencedata,
                  "RCCS_data_estimate_incidence_inland_R6_R18",
                  "220903",
                  "anonymized_id_for_incidence_estimate_221129.csv")
write.csv(anonimized_id, file = file, row.names = FALSE)

# statistics
file_name <- file.path(outdir,
                       "incidence_inland_statistics_cohort_for_paper_221129.rds")
saveRDS(stats, file_name)

# seroconvert cohort central analysis
file <- file.path(indir, "data", "seroconverter_cohort_R6R19.rds")
saveRDS(seroconverter_cohort_list, file = file)

# seroconvert cohort sensitivity analysis (continuously surveyed commm)
# file <- file.path(indir,
#                   "data",
#                   "seroconverter_cohort_30comm_R6R19.rds")
# saveRDS(seroconverter_cohort_list, file = file)
