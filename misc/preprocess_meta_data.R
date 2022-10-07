library(data.table)
library(lubridate)
library(dplyr)

# change as appropriate
make.hiv.history.plots <- 0

if(dir.exists('/Users/melodiemonod'))
{
  indir.repository <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
}

if(dir.exists('/home/andrea'))
{
  indir.repository <-'~/git/phyloflows'
  make.hiv.history.plots <- 1
  indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata/'
  outdir <- '/home/andrea/Documents/Box/2021/phyloflows/testing_history/'
}

if(dir.exists('/rds/general/user/'))
{
  indir.repository <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
}

# file paths keys
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')

# Latest data from Rakai's CCS (Joseph's data from 2022-01-29)
file.path.allhiv <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'All_HIVpcr_for_questR15_R18_220129.csv')
file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'HIV_R15_R18_VOIs_220129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'quest_R15_R18_VoIs_220129.csv')

# Latest data from Rakai's CCS (Kate's data from 2022-03-08)
file.path.metadata <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata__12Nov2019.csv')

# Latest data from Rakai's CCS (Kate's data from 2022-03-08)
file.path.neuro.metadata <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Pangea_Rakai_NeuroStudy_Metadata_11Dec2015.csv')


#
# LOAD FUNCTIONS
#

source(file.path(indir.repository, 'functions', 'utils.R'))
source(file.path(indir.repository, 'misc', 'functions', 'preprocess_meta_data-functions.R'))

#
# LOAD DATA
#

# Load keys
aik <- .read(file.anonymisation.keys); aik$X <- NULL
colnames(aik) <- tolower(colnames(aik))
community.keys <-.read(file.community.keys)

# Load Joseph's data from 2022-01-29
allhiv <- .read(file.path.allhiv) 
hiv <- .read(file.path.hiv)
quest <- .read(file.path.quest)

# Load Kate's data from 2022-03-08
raw_metadata <- .read(file.path.metadata)

# Load Neuro cohort's data from 2022-03-17
raw_neuro_metadata <- .read(file.path.neuro.metadata)
setnames(raw_neuro_metadata,  'studyid', 'study_id')


#
# PROCESS RAW DATA
#

# Add 'RK-' to all data.tables study_ids
invisible(lapply(list(hiv, allhiv, quest, raw_metadata, raw_neuro_metadata),
                 function(dt) {
                   dt[!grepl('RK-', study_id), study_id := paste0('RK-', study_id)]
                 }))

# process quest and make date.birth
quest <- process.quest(quest)

# make date of first positive and last negative test with allhiv
date.first.positive <- make.date.first.positive(allhiv)

# process hiv and find date first and last visit 
hiv <- process.hiv(hiv)

date.first.last.visit <- make.date.first.last.visit(hiv)

# what s the difference between hiv and allhiv?
# hiv probably contains all tests, while allhiv only those from people ever tested positive?
# There are 6 entries with different first positive diagnoses between hiv and allhivl
# THink it makes more sense to take values from hiv, BUT be careful to those that
# turned negative after a first positive.
date.first.positive <- compare.hiv.allhiv.firstpositivedates(hiv, date.first.positive)
stopifnot(.vars.with.multiple.values(date.first.positive, 'study_id')[, .N == 0])

# process Kate's meta data
meta_data_2 <- process.meta.data(raw_metadata, aik, community.keys)

# process Neuro's meta data
meta_data_neuro <- process.neuro.meta.data(copy(raw_neuro_metadata), aik)

#
# MAKE META DATA
#

# get Joseph's meta data
meta_data <- get.meta.data(quest, date.first.positive, date.first.last.visit, aik, community.keys)

# add Kate's data for missing individuals
meta_data <- rbind(meta_data, meta_data_2[!study_id %in% meta_data[, study_id]])

# get round by study_id
# tmp <- merge(quest[!is.na(round), .(study_id, round)], raw_metadata[!is.na(round), .(study_id, round)], by = 'study_id', all.x = T, all.y = T)
# tmp[, round.x := gsub('R0(.*)', '\\1', round.x)]
# tmp <- tmp[, list(round = paste0('R0', sort(unique(na.omit(c(round.x, round.y)))), collapse = '_')), by = 'study_id']
# meta_data <- merge(select(meta_data, -round), tmp[, .(study_id, round)], by = 'study_id')

# add Neuro's data for missing individuals
meta_data <- rbind(meta_data, meta_data_neuro[!study_id %in% meta_data[, study_id]], fill=TRUE)

# compare to missing
# missing <- as.data.table(read.csv(file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'missing_study_id_RCCS_R15_R18_220113.csv')))
# stillmissing <- missing[!RCCS_studyid %in% meta_data[, gsub('RK-(.+)', '\\1', study_id)]]
# stillmissing[study_id_in_neuro ==F] # only Discarded study_id
# stillmissing[study_id_in_neuro ==T]


#
# SAVE META DATA
#
save(meta_data, file = file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_20220329.RData'), row.names = F)
