library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

file.census.count <- file.path(indir.repository, 'data', 'RCCS_census_eligible_individuals_221116.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_220909.csv')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))
ncen <- as.data.table(read.csv(file.census.count))
quest <- as.data.table(read.csv(file.path.quest))


################################

# FIND COUNT OF PARTICIPANT

################################

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate)]

# Set to date format
rin[, intdate := as.Date(intdate, format = '%d-%b-%y')]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# SET ROUND 15S IN INLAND AS 15
rinc[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
rinc[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
rinc <- rinc[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]

# GET PARTICIPANT
rinc <- rinc[, list(PARTICIPANT  = .N), by = c('AGEYRS', 'ROUND', 'SEX', 'COMM')]

# GET PARTICIPANT SMOOTH with loess smooth
rinc <- rinc[order(COMM, ROUND, SEX, AGEYRS)]
rinc <- rinc[, {
  loessMod50 <- loess(PARTICIPANT ~ AGEYRS, span=0.5)
  smoothed50 <- predict(loessMod50, new_data = AGEYRSPREDICT) 

  list(AGEYRS = AGEYRS, PARTICIPANT_SMOOTH = smoothed50, PARTICIPANT = PARTICIPANT)
}, by = c('COMM', 'SEX', 'ROUND')]


################################

# GET PROPORTION OF PARITCIPATION

################################

tmp <- select(ncen, c('AGEYRS', 'ROUND', 'SEX', 'COMM', 'ELIGIBLE_NOT_SMOOTH', 'ELIGIBLE'))
rinc[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
rpr <- merge(rinc,tmp , by =  c('AGEYRS', 'ROUND', 'SEX', 'COMM'))
rpr[, PARTICIPATION := PARTICIPANT / ELIGIBLE_NOT_SMOOTH]
rpr[, PARTICIPATION_SMOOTH := PARTICIPANT_SMOOTH / ELIGIBLE]
rpr[PARTICIPATION > 1, PARTICIPATION := 1]
rpr[PARTICIPATION_SMOOTH > 1, PARTICIPATION_SMOOTH := 1]


################################

# SAVE

################################

tmp <- rpr[!ROUND %in% c("06", "07", "08", "09"), list(MEAN = paste0(round(mean(PARTICIPATION)*100))), by = c('SEX', 'COMM')]
saveRDS(tmp, file.path(outdir, 'Participation.rds'))

file <- file.path(indir.repository, 'data', 'RCCS_participation_221116.csv')
tmp <- rpr[, .(AGEYRS, ROUND, SEX, COMM, PARTICIPANT, PARTICIPATION_SMOOTH)]
setnames(tmp, 'PARTICIPATION_SMOOTH', 'PARTICIPATION')
write.csv(tmp, file = file, row.names = F)

