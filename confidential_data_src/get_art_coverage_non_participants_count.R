library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library("haven")
library(here)

usr <- Sys.info()[['user']]
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
    indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
    indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
}else if(usr == 'andrea')
{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live'
}

indir.repository <- here()

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'HIV_R6_R18_221129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_221208.csv')
path.tests <- file.path(indir.deepsequencedata, 'RCCS_R15_R20',"all_participants_hivstatus_vl_220729.csv")

c(  file.community.keys,
    file.path.hiv,
    file.path.quest,
    path.tests) |> file.exists() |> all() |> stopifnot()

# load files
community.keys <- fread(file.community.keys)
quest <- fread(file.path.quest)
hiv <- fread(file.path.hiv)


#################################

# FIND SELF-REPORTED ART USE  #

#################################

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate, arvmed, cuarvmed)]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# find index of round
rinc <- rinc[order(STUDY_ID, ROUND)]
rinc[, INDEX_ROUND := 1:length(ROUND), by = 'STUDY_ID']

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# keep HIV positive
rprev <- hivs[HIV == 'P']

# get ART status
rprev[!ROUND %in% c('R016', 'R017', 'R018'), ART := ARVMED ==1]
rprev[!ROUND %in% c('R016', 'R017', 'R018') & is.na(ARVMED), ART := F]
rprev[ROUND == 'R016', ART := ARVMED ==1 | CUARVMED ==1]
rprev[ROUND == 'R016' & (is.na(ARVMED) | is.na(CUARVMED)), ART := F]
rprev[ROUND %in% c('R017', 'R018'), ART := CUARVMED ==1]
rprev[ROUND %in% c('R017', 'R018') & is.na(CUARVMED), ART := F]
rprev[, table(ROUND)]


#################################

# ADD VIRAL LOAD DATA  #

#################################

# for round with suppressed set art to true if indiv is suppressed

# tuning
VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 1000 # WHO standards

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND %in% c(15:18, 15.5)]

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# remove HIV+ individuals with missing VLs  
DT <- subset(dall, HIV_STATUS==0 | HIV_AND_VL==1)

# set ARVMED to 0 for HIV-
set(DT, DT[, which(ARVMED==1 & HIV_STATUS==0)], 'ARVMED', 0) 

# define VLC as VL_COPIES for HIV+ and as 0 for HIV-
set(DT, NULL, 'VLC', DT$VL_COPIES)
set(DT, DT[,which(HIV_STATUS==0)], 'VLC', 0)

# define detectable VL as VLD and undetectable VL as VLU (machine-undetectable)
set(DT, NULL, 'VLU', DT[, as.integer(VLC<VL_DETECTABLE)])
set(DT, NULL, 'VLD', DT[, as.integer(VLC>=VL_DETECTABLE)])

# define suppressed VL as VLS and unsuppressed as VLNS (according to WHO criteria)	
set(DT, NULL, 'VLS', DT[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
set(DT, NULL, 'VLNS', DT[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])

# find HIV+ and detectable
set(DT, NULL, 'HIV_AND_VLD', DT[, as.integer(VLD==1 & HIV_AND_VL==1)])

# reset undetectable to VLC 0
set(DT, DT[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
setkey(DT, ROUND, FC, SEX, AGEYRS)

# keep within census eligible age
DT <- subset(DT, AGEYRS <= 50)

# keep infected
DT <- DT[HIV_STATUS ==1]

# merge to self-reported data
tmp <- DT[, .(STUDY_ID, ROUND, SEX, FC, VLNS, AGEYRS)]
tmp[, ROUND := paste0('R0', ROUND)]
setnames(tmp, 'AGEYRS', 'AGEYRS2')
rprev <- merge(rprev, tmp, by.x = c('STUDY_ID', 'ROUND', 'SEX', 'COMM'), by.y = c('STUDY_ID', 'ROUND', 'SEX', 'FC'), all.x = T, all.y = T)

# set ageyrs to the viral load data if available
rprev[!is.na(AGEYRS2), AGEYRS := AGEYRS2]
set(rprev, NULL, 'AGEYRS2', NULL)

# set art to true if viremic viral load
rprev[VLNS == 0, ART := T]
rprev <- rprev[!is.na(ART)]

# remove na vlns for round 16 onwards otherwise it leads to % art < % suppressed
nrow(rprev[ROUND == 'R016' & is.na(VLNS)]) / nrow(rprev[ROUND == 'R016' & !is.na(VLNS)])
nrow(rprev[ROUND == 'R017' & is.na(VLNS)]) / nrow(rprev[ROUND == 'R017' & !is.na(VLNS)])
nrow(rprev[ROUND == 'R018' & is.na(VLNS)]) / nrow(rprev[ROUND == 'R018' & !is.na(VLNS)])
rprev <- rprev[!(ROUND %in% c('R016', 'R017', 'R018') & is.na(VLNS))]


#################################

# KEEP INDIVIDUALS SEEN FOR THE FIRST TIME  
# THAT ARE THE CLOSEST TO NON-PARTICIPANTS

#################################

# get proportion of participants seen for the first time 
newpart <- rprev[HIV == 'P' & COMM == 'inland', round(sum(INDEX_ROUND == 1) / length(INDEX_ROUND) * 100, 2), by = 'ROUND']
newpart <- newpart[!ROUND %in% c('R006', 'R007', 'R008', 'R009')]
print(newpart[V1 == min(V1)])
print(newpart[V1 == max(V1)])

#KEEP INDIVIDUALS SEEN FOR THE FIRST TIME  
rprev <- rprev[INDEX_ROUND == 1]


#################################

# SET ROUND 15S IN INLAND AS 15

#################################

rprev[ROUND == 'R015S' & COMM == 'inland', ROUND := 'R015']


#################################

# AGGREGATE BY ROUND, SEX, COMM AND AGE  #

#################################

# find self reported under art for participant
rart <- rprev[, list(COUNT = sum(ART == T), TOTAL_COUNT = length(ART)), by = c('ROUND', 'SEX', 'COMM', 'AGEYRS')]


#################################

# SAVE DE-IDENTIFIED DATA  #

#################################

file.name=file.path(indir.repository, 'data', 'aggregated_newlyregistered_count_art_coverage.csv')
if( !file.exists(file.name))
{
    cat('\n Saving', file.name,'...\n')
    write.csv(rart, file =file.name , row.names = F)
}else{
    cat('\n Output file', file.name,'already exists.\n')
}
