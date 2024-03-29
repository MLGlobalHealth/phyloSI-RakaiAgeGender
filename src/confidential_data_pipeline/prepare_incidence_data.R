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
library(here)

# do you want to restrict the analysis to the 30 communities continuously surveyed?
restrict_to_30_comms <- F

# paths
gitdir <- here()
source(file.path(gitdir, "config.R"))

c(  file.path.flow.614,
    file.path.hiv.614,
    file.path.flow,
    file.path.hiv.1518,
    file.path.flow_19,
    file.path.hiv_19) |> file.exists() |> all() |> stopifnot()

# save intermediary results not on Zolodo
outdir <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "incidence_rate_inland")

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

if(restrict_to_30_comms){
  # sensitivity analysis
  inland_comms <- c("1", "2", "4", "5", "6", "7", "8",
                    "16", "19", "22", "24", "29", "33",
                    "34", "40", "56", "57", "58", "62",
                    "74", "77", "89", "94", "106", "107",
                    "108", "120", "391", "602", "754")
  
  verif_el <- verif_1[resident == 1 &
                      ageyrs >= 15 &
                      ageyrs <= 49 &
                      comm_num %in% inland_comms, ]

}

# central analysis
verif_el <- verif_1[resident == 1 &
                    ageyrs >= 15 &
                    ageyrs <= 49 &
                    !comm_num %in% c(38, 770, 771, 774), ]
verif_el[, table(round)]
rm(verif_1)

# only include HIV data as determined
verified_visit_ids <- unique(verif_el[, visit_id])
verified_visit_ids_r1018 <- unique(verif_el[visit %in% paste0('R0', 10:18), visit_id])
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
file.name <- file.anonymised.id
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  write.csv(anonimized_id, file = file.name, row.names = FALSE)
}else{
  cat("\n Output file", file.name, "already exists\n")
}

# statistics
file.name <- file.path(outdir,"incidence_inland_statistics_cohort_for_paper_230715.rds")
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(stats, file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}

# seroconvert cohort central analysis
if(!restrict_to_30_comms){
  file.name <- file.path.seroconverter_cohort
}else{
  file.name <- file.path.seroconverter_cohort.30
}

if( !file.exists(file.name))
{
  cat('\n Careful: This data should already exist exist in ', file.name  )
  cat('\n check that your Zenodo path is correctly specified in config.R ' )
  cat('\nIf you wish to proceed, and save this file anyway run the commented line below')
  #   saveRDS(seroconverter_cohort_list, file = file)
}else{
  cat('\n Output file', file.name,'already exists.\n')
}
