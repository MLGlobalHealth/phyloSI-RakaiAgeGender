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
indir <- '~/git/phyloflows'

# set up path to data
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'

# set path to save results
outdir <- file.path(indir.deepsequence.analyses, 'PANGEA2_RCCS', 'incidence_rate_inland')

# round 6 to 18: only use round 6 to 14
file.path.flow.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'verif_1.dta')
file.path.hiv.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'hivincidence_1.dta')

# round 15 to 18
file.path.flow.1518 <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'FlowR15_R18_VoIs_220129.csv')
file.path.hiv.1518 <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'HIV_R15_R18_VOIs_220129.csv')

# round 19
file.path.flow.19 <- file.path(indir.deepsequencedata, 'R019_VoIs', 'Flow_R019_VoIs_220607.csv')
file.path.hiv.19 <- file.path(indir.deepsequencedata, 'R019_VoIs', 'HIV_R019_VOIs_220607.csv')

# function
source(file.path(indir, 'functions', 'incidence_rate_estimation_functions.R'))

# utils
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
hivstatus_vlcopies_raw <- read_hiv_data(file.path.hiv.614, file.path.hiv.1518, file.path.hiv.19)

# load HIV verification data set 
verif_raw <- read_flow_data(file.path.flow.614, file.path.flow.1518, file.path.flow.19)


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

# ADD ANONIMIZED ID AND KEEP NECESSARY COVARIATE  #

################

# anonimized id
anonimized_id <- data.table(research_id = verif_1[, sort(unique(research_id))])
anonimized_id <- anonimized_id[research_id != '']
anonimized_id[, anonimized_id := paste0('RkINC-', str_pad(1:nrow(anonimized_id), 5, pad = "0"))]

# verif 1
verif_1 <- merge(verif_1, anonimized_id, by = 'research_id')
verif_1 <- verif_1[, .(anonimized_id, round, ROUND, visit, sex, comm_num, ageyrs, resident, birthdat)]
verif_1[, visit_id:=paste(anonimized_id, round, sep=":")]
setnames(verif_1, 'anonimized_id', 'research_id')

# hivstatus_vlcopies_1
hivstatus_vlcopies_1 <- merge(hivstatus_vlcopies_1, anonimized_id, by = 'research_id')
hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[, .(anonimized_id, round, ROUND, hiv, hivstatus, hivdate, year)]
hivstatus_vlcopies_1[, visit_id:=paste(anonimized_id, round, sep=":")]
setnames(hivstatus_vlcopies_1, 'anonimized_id', 'research_id')

################

# SAVE

################

# anonimized id
file <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/','anonimized_id_for_incidence_estimate.csv')
write.csv(anonimized_id, file = file, row.names = F)

# verif 1
file <- file.path(indir, 'data', 'metadata_for_incidence_estimation_R6R19.rds')
saveRDS(verif_1, file = file)

# hivstatus_vlcopies_1
file <- file.path(indir, 'data', 'hivstatus_for_incidence_estimation_R6R19.rds')
saveRDS(hivstatus_vlcopies_1, file = file)


