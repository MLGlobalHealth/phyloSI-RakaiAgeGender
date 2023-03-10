library(data.table)
library(dplyr)
library(scales)
library(lubridate)
library("haven")
library(here)

usr <- Sys.info()[['user']]
indir.repository <- here()
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'

if(usr=="andrea")
{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live'
}

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')
file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'HIV_R6_R18_221129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_221208.csv')

c(  file.community.keys,
    file.path.hiv,
    file.path.quest) |> file.exists() |> all() |> stopifnot()

# load files
community.keys <- fread(file.community.keys)
quest <- fread(file.path.quest)
hiv <- fread(file.path.hiv)


#################################

# HIV TESTS USING HIV DATA SET #

#################################

if(0)
{ # check percentage with hiv tests
  
  for(Round in hiv[, sort(unique(round))]){
    hiv_n <- hiv[round == Round, length(unique(study_id))]
    participant_n <- quest[round == gsub(' ', '', Round), length(unique(study_id))]
    cat('There is ', participant_n, 'participants in round', Round, ', ')
    cat(hiv_n, 'of them have an hiv test result (', round(hiv_n  / participant_n, 4), '%)\n')
  }
  
}

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

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# SET ROUND 15S IN INLAND AS 15
hivs[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
hivs[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
hivs <- hivs[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]

# find HIV prevalence rate for participant
rprev <- hivs[, list(COUNT = sum(HIV == 'P'),
                     TOTAL_COUNT = length(HIV)), by = c('ROUND', 'SEX', 'COMM', 'AGEYRS')]


#########

# SAVE #

#########

filename = file.path(indir.repository, "data", "aggregated_count_hiv_positive.csv")
if(! file.exists(filename))
{
    cat("\n Saving", filename, "...\n")
    write.csv(rprev, filename, row.names = F)
}else{
    cat("\n Output file", filename, "already exists\n")
}
