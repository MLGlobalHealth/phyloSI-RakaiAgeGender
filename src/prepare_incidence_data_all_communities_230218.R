#################################################################################################
# Rakai data analysis for EMOD inputs
# August 12, 2022
# Adam Akullian/Kate Grabowski/Melodie Monod;
#################################################################################################

library("ggplot2")
library("data.table")
library("dplyr")
library("haven")
library("tidyr")
library("tidyverse")
library("mgcv")
library("zoo")
library("remotes")
library("metR")
library("lubridate")
library(ggpubr)

# set to directory
indir <- '~/git/phyloSI-RakaiAgeGender'

# set up path to data
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'

# set path to save results
outdir <- file.path(indir.deepsequence.analyses, 'PANGEA2_RCCS', 'incidence_rate_inland')

# round 6 to 18: only use round 6 to 8
file.path.flow.68 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'verif_1.dta')
file.path.hiv.68 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'hivincidence_1.dta')

# round 9 to 14:
file.path.flow.914 <- file.path(indir.deepsequencedata, 'RCCS_R9_R14', 'FlowR09_R14.csv')
file.path.hiv.914 <- file.path(indir.deepsequencedata, 'RCCS_R9_R14', 'HIV_R09_R14.csv')

# round 15 to 18
file.path.flow.1518 <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'FlowR15_R18_VoIs_220129.csv')
file.path.hiv.1518 <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'HIV_R15_R18_VOIs_220129.csv')

# round 19
file.path.flow.19 <- file.path(indir.deepsequencedata, 'R019_VoIs', 'Flow_R019_VoIs_220607.csv')
file.path.hiv.19 <- file.path(indir.deepsequencedata, 'R019_VoIs', 'HIV_R019_VOIs_220607.csv')

# function
source(file.path(indir, 'functions', 'incidence_rate_estimation_functions.R'))

# utils
rounds_group_0 <- c("R006","R007", "R008")
rounds_group_1 <- c("R006","R007", "R008", "R009", "R010", "R011", "R012", "R013", "R014")
rounds_group_2 <- c("R015", "R015S", "R016", "R017", "R018")
rounds_group_3 <- c('R019')

# make df round
df_round <- make_df_round(rounds_group_1, rounds_group_2, rounds_group_3)

rounds_numeric_group_1 <- df_round[visit %in% rounds_group_1, round_numeric]
rounds_numeric_group_2 <- df_round[visit %in% rounds_group_2, round_numeric]
rounds_numeric_group_3 <- df_round[visit == rounds_group_3, round_numeric]



##############

# LOAD DATA #

##############

# load HIV status data 
hivstatus_vlcopies_raw <- read_hiv_data_230218(file.path.hiv.68, file.path.hiv.914, file.path.hiv.1518, file.path.hiv.19)

# load HIV verification data set 
verif_raw <- read_flow_data_230218(file.path.flow.68, file.path.flow.914, file.path.flow.1518, file.path.flow.19)


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
anonimized_id <- anonimized_id[research_id != '']

# randomly assign new id
set.seed(12)
anonimized_id[, anonimized_id := paste0('RkINC-', sample(str_pad(1:nrow(anonimized_id), 5, pad = "0"), replace = F))]


####################################

# RESTRICT ANALAYSIS BY CRITERIA  #

####################################

# Use Verification data to restrict analysis based on following criteria
## 1. Age 15-49
## 2. residency status: resident==1 
## 3. Resides in inland communities 

# Subset data based inclusion criteria 
verif_el<-verif_1[resident==1 &ageyrs>=15 & ageyrs<=49 & !comm_num%in%c(38,770,771,774),]
verif_el[, table(round)]
rm(verif_1)

# only include HIV data as determined 
hivstatus_vlcopies_1<-hivstatus_vlcopies_1[visit_id%in%verif_el[, visit_id],]
table(hivstatus_vlcopies_1$ROUND)


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


###############################

# FIND SEROCONVERSION DATE #

###############################

N <- 50
seed = 12

seroconverter_cohort.list = vector(mode = 'list', length = N)

set.seed(seed)
for(i in 1:N){
  
  cat('\niteration', i)
  
  # generate random date of infection and find perosn year
  seroconverter_cohort <- find_seroconvert_cohort(status_df, hivstatus_vlcopies_1_inc)
  
  # add anonimized id
  seroconverter_cohort <- merge(seroconverter_cohort, anonimized_id, by = 'research_id')
  
  # keep only necessary covariate
  seroconverter_cohort <- seroconverter_cohort[, .(anonimized_id, age, sex, number_missing_visits, 
                                                   round, py, hivinc)]
  
  # prepare
  seroconverter_cohort$iterations <- i
  seroconverter_cohort.list[[i]] <- seroconverter_cohort
  
}


################

# SAVE

################

# anonimized id
file <- file.path(indir.deepsequencedata, 'RCCS_R9_R14/RCCS_data_estimate_incidence_inland_R6_R18_230218','anonimized_id_for_incidence_estimate.csv')
write.csv(anonimized_id, file = file, row.names = F)

# seroconvert cohort central analysis
file <- file.path(indir.deepsequencedata, 'RCCS_R9_R14/RCCS_data_estimate_incidence_inland_R6_R18_230218', 'seroconverter_cohort.rds')
saveRDS(seroconverter_cohort.list, file = file)


