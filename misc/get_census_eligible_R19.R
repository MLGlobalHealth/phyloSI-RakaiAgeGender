library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.repository <- '~/git/phyloflows'

file.path.flow <- file.path(indir.deepsequencedata, 'R019_VoIs', 'Flow_R019_VoIs_220607.csv')

flow <- as.data.table(read.csv(file.path.flow))

# source(file.path(indir.repository, 'functions', 'utils.R'))
# source(file.path(indir.repository, 'misc', 'functions', 'preprocess_meta_data-functions.R'))


###############################

# CENSUS ELIGIBLE INDIVIDUALS #

###############################

# Code for ineligibility

flow[, reason_ineligible := NA_character_]

flow[locate1==10 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==2 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==13 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==5 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==10, reason_ineligible := "Out_migrated"]

flow[locate1==7 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==2 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==6 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==3 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==5 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==17 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==88, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==2, reason_ineligible := "Already_seen"]

flow[locate1==11 & locate2==8, reason_ineligible := "Dead"]
flow[locate2==11, reason_ineligible := "Dead"]

flow[resident==0, reason_ineligible := "not_resident"]

flow[ageyrs<15 | ageyrs > 49, reason_ineligible := "Not_within_eligible_age_range"]

flow[is.na(reason_ineligible), reason_ineligible := 'none']

# find count eligible
re <- flow[, list(count = .N), by = c('reason_ineligible', 'round', 'comm_num', 'ageyrs', 'sex')]
re <- dcast.data.table(re, round + comm_num + ageyrs + sex ~ reason_ineligible, value.var = 'count')
re[is.na(re)] = 0
re[, ELIGIBLE := round(none + Out_migrated / 2)]
re <- re[ELIGIBLE != 0]

# fill missing count with 0
tmp <- data.table(expand.grid(comm_num = re[, (unique(comm_num))], ageyrs = 15:49, sex = c('F', 'M')))
re <- merge(re, tmp, by = c('comm_num', 'ageyrs', 'sex'), all.y = T)
re[is.na(re)] = 0

# additional variable
colnames(re) <- toupper(colnames(re))
re[, ROUND := substring(ROUND, 3)]

# save
write.csv(re, file.path(indir.deepsequencedata, 'R019_VoIs', 'RCCS_census_eligible_individuals_220621.csv'), row.names = F)


