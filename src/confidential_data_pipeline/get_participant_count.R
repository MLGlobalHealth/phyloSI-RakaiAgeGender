library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(here)

gitdir <- here()
source(file.path(gitdir, "config.R"))

c(  file.eligible.count,
    file.community.keys,
    file.path.quest) |> file.exists() |> all() |> stopifnot()

# load files
community.keys <- fread(file.community.keys)
ncen <- fread(file.eligible.count)
quest <- fread(file.path.quest)

# path for intermediary results
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')


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

# participation rate aggregated by age and round but not by sex, comm
tmp <- rpr[!ROUND %in% c("06", "07", "08", "09"), list(MEAN = paste0(round(mean(PARTICIPATION)*100))), by = c('SEX', 'COMM')]
tmp1 <- rpr[!ROUND %in% c("06", "07", "08", "09") , list(MEAN = paste0(round(mean(PARTICIPATION)*100))), by = 'COMM']
tmp1[, `:=` (SEX = 'Total')]
tmp <- rbind(tmp, tmp1)

# min and max participation rate aggregated by age, sex but not by comm, round
tmp1 <- rpr[!ROUND %in% c("06", "07", "08", "09"), list(MEAN = paste0(round(mean(PARTICIPATION)*100, 1))), by = c('ROUND', 'COMM')]
tmp1 <- tmp1[, list(MIN = min(MEAN), MAX= max(MEAN)), by = 'COMM']

# participation rate non-aggregated
tmp <- rpr[, .(AGEYRS, ROUND, SEX, COMM, PARTICIPANT, PARTICIPATION_SMOOTH)]
setnames(tmp, 'PARTICIPATION_SMOOTH', 'PARTICIPATION')

# save
file.name <- file.participation
if( !file.exists(file.name))
{
  cat('\n Careful: This data should already exist exist in ', file.name  )
  cat('\n check that your Zenodo path is correctly specified in config.R ' )
  cat('\nIf you wish to proceed, and save this file anyway run the commented line below')
  #  write.csv(tmp, file = file, row.names = F)
}else{
  cat('\n Output file', file.name,'already exists.\n')
}

