library(data.table)
library(lubridate)
library(dplyr)

# change as appropriate
make.hiv.history.plots <- 0

# paths
gitdir <- here()
source(file.path(gitdir, "config.R"))

# make sure all files exist
file.exists(c(
  file.anonymisation.keys ,
  file.community.keys ,
  file.path.allhiv,
  file.path.hiv.1518,
  file.path.quest.1518,
  file.path.metadata,
  file.path.neuro.metadata,
  file.path.update.first.positive))  |> all() |> stopifnot()


#
# LOAD FUNCTIONS
#

source(file.path(gitdir.functions, 'utils.R'))
source(file.path(gitdir.functions, "functions_confidential_data_pipeline", "preprocess_meta_data-functions.R"))


#
# LOAD DATA
#

# Load keys
aik <- .read(file.anonymisation.keys); aik$X <- NULL
colnames(aik) <- tolower(colnames(aik))
community.keys <-.read(file.community.keys)

# Load Joseph's data from 2022-01-29
allhiv <- .read(file.path.allhiv) 
hiv <- .read(file.path.hiv.1518)
quest <- .read(file.path.quest.1518)

# Load Kate's data from 2022-03-08
raw_metadata <- .read(file.path.metadata)

# Load Neuro cohort's data from 2022-03-17
raw_neuro_metadata <- .read(file.path.neuro.metadata)
setnames(raw_neuro_metadata,  'studyid', 'study_id')

# update with joseph's first positive frp, 2022-11-28
firstpos_update <- fread(file.path.update.first.positive )
setnames(firstpos_update, c('firstposdat', 'lastnegdat'), c('date_first_positive', 'date_last_negative'))


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
meta_data[, date_last_negative := as.character(date_last_negative)]
meta_data_2[, date_last_negative := as.character(date_last_negative)]
meta_data[, date_first_positive := as.character(date_first_positive)]
meta_data_2[, date_first_positive := as.character(date_first_positive)]
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
file.name <- path.meta.confidential

if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  save(meta_data, file = , row.names = F)
}else{
  cat("\n Output file", file.name, "already exists\n")
}
cat("\n Done \n")



