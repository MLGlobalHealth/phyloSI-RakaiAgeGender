library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library(haven)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

file.eligible.count <- file.path(indir.repository, 'data', 'RCCS_census_eligible_individuals_221116.csv')
path.tests <- file.path(indir.deepsequencedata, 'RCCS_R15_R20',"all_participants_hivstatus_vl_220729.csv")
file.seq.count <- file.path(outdir, 'characteristics_sequenced_ind_R14_18.rds')

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'HIV_R6_R18_221129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_220909.csv')

# Latest data from Rakai's CCS (Kate's data from 2022-03-08)
file.path.metadata <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata__12Nov2019.csv')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 10:18)), 
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(15, '15S', 16:18))))

# community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']


###############################################

# GET META DATA 

###############################################

#
# Meta data
#

meta_data <- as.data.table(read.csv(file.path.metadata)) #additional meta_data

# find age
meta_data[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
meta_data[, AGEYRS := round(lubridate::time_length(difftime(sample_date, date_birth),"years"))]
meta_data[is.na(AGEYRS), AGEYRS := round(lubridate::time_length(difftime(firstposvd, date_birth),"years"))]
meta_data[is.na(AGEYRS)]

# restrict age
meta_data <- meta_data[AGEYRS > 14 & AGEYRS < 50]

# find community
meta_data[, COMM := 'inland']
meta_data[LakeVictoria_FishingCommunity == 'yes', COMM := 'fishing']

# find hiv status
meta_data[, HIV := ifelse(is.na(firstposvd), 'N', 'P')]

# find art use
meta_data[, ART := artslfuse == 'yes']

# keep variable of interest
meta_data[, round := paste0('R0', round)]
colnames(meta_data) <- toupper(colnames(meta_data))
meta_data <- meta_data[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV, ART)]

# set 15.1 to be 15S
meta_data[ROUND == 'R015.1', ROUND := 'R015S']


#
# Quest
#

quest <- as.data.table(read.csv(file.path.quest))

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, arvmed)]

# find  community
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restrict age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# get ART status
rinc[, ART := ARVMED ==1]
rinc[is.na(ARVMED), ART := F]

# add meta data from Kate
tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], rinc[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
rinc <- rbind(rinc[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)])

# SET ROUND 15S IN INLAND AS 15
rincp <- copy(rinc)
rincp[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
rincp[ROUND =='R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
rincp <- rincp[!(ROUND =='R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]


#
# HIV
#

hiv <- as.data.table(read.csv(file.path.hiv))
hiv[, round := gsub(' ', '', round)] # remove space in string

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))

# add meta data from Joseph 
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# add meta data from Kate
tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], hivs[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
hivs <- rbind(hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)])

# SET ROUND 15S IN INLAND AS 15
hivs[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
hivs[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
hivs <- hivs[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]

# restric age
hivs <- hivs[AGEYRS > 14 & AGEYRS < 50]


#################################

# GET CENSUS ELIGIBLE #

#################################

# load census eligible count
eligible_count <- as.data.table(read.csv(file.eligible.count))

# find census eligible
census <- eligible_count[, list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Total'), by = c('COMM', 'ROUND')]
census <- rbind(census, eligible_count[SEX == 'F', list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Female'), by = c('COMM', 'ROUND')])
census <- rbind(census, eligible_count[SEX == 'M', list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Male'), by = c('COMM', 'ROUND')])
census <- rbind(census, eligible_count[SEX == 'F' & AGEYRS < 25, list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
census <- rbind(census, eligible_count[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
census <- rbind(census, eligible_count[SEX == 'F' & AGEYRS > 34, list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
census <- rbind(census, eligible_count[SEX == 'M' & AGEYRS < 25, list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
census <- rbind(census, eligible_count[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
census <- rbind(census, eligible_count[SEX == 'M' & AGEYRS > 34, list(ELIGIBLE = sum(ELIGIBLE_NOT_SMOOTH),  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])
census[, ROUND := paste0('R0', ROUND)]

# keep round of interest
census <- merge(census, df_round, by = c('COMM', 'ROUND'))


#################################

# GET PARTICIPANTS #

#################################

# find participant
part <- rincp[, list(PARTICIPANT = .N,  TYPE = 'Total'), by = c('COMM', 'ROUND')]
part <- rbind(part, rincp[SEX == 'F', list(PARTICIPANT = .N,  TYPE = 'Female'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M', list(PARTICIPANT = .N,  TYPE = 'Male'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'F' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'F' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
part <- merge(part, df_round, by = c('COMM', 'ROUND'))


########################################

# GET HIV-POSITIVE AMONG PARTICIPANTS #

########################################

# keep only positive
rprev <- hivs[, list(COUNT = sum(HIV == 'P')), by = c('ROUND', 'SEX', 'COMM', 'AGEYRS')]

# get hiv table
hivp <- rprev[, list(HIV = sum(COUNT),  TYPE = 'Total'), by = c('COMM', 'ROUND')]
hivp <- rbind(hivp, rprev[SEX == 'F', list(HIV = sum(COUNT),  TYPE = 'Female'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M', list(HIV = sum(COUNT),  TYPE = 'Male'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'F' & AGEYRS < 25, list(HIV = sum(COUNT),  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(HIV = sum(COUNT),  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'F' & AGEYRS > 34, list(HIV = sum(COUNT),  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M' & AGEYRS < 25, list(HIV = sum(COUNT),  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(HIV = sum(COUNT),  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M' & AGEYRS > 34, list(HIV = sum(COUNT),  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
hivp <- merge(hivp, df_round, by = c('COMM', 'ROUND'))


######################################################

# GET HIV-POSITIVE AND ART NAIVE AMONG PARTICIPANTS #

######################################################

# keep HIV positive  
sart <- hivs[HIV == 'P']

# find participant
sartp <- sart[, list(SELF_REPORTED_ART = sum(ART == F),  TYPE = 'Total'), by = c('COMM', 'ROUND')]
sartp <- rbind(sartp, sart[SEX == 'F', list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M', list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'F' & AGEYRS < 25, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'F' & AGEYRS > 34, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M' & AGEYRS < 25, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M' & AGEYRS > 34, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
sartp <- merge(sartp, df_round, by = c('COMM', 'ROUND'))


######################################################################

# GET HIV-POSITIVE WITH UNSUPPRESSED VIRAL LOADS AMONG PARTICIPANTS #

######################################################################

# tuning
VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 1000 # WHO standards

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND %in% c(15:18)]
# dall <- dall[ROUND == round]

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# keep within census eligible age
DT <- subset(dall, AGEYRS > 14 & AGEYRS < 49)

# remove HIV+ individuals with missing VLs  
DT <- subset(DT, HIV_STATUS==0 | HIV_AND_VL==1)

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

# use age from quest
DT[, round := paste0('R0', ROUND)]
set(DT, NULL, 'AGEYRS', NULL)
# DT <- merge(DT, hiv, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
DT <- merge(DT, quest, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
setnames(DT, 'ageyrs', 'AGEYRS')

# get count for every categories
tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
tmp1 <- DT[, sort(unique(ROUND))]
vla <- as.data.table(expand.grid(ROUND=tmp1,
                                 FC=c('fishing','inland'),
                                 SEX=c('M','F'),
                                 AGEYRS=tmp))
vla <- vla[, {		
  z <- which(DT$ROUND==ROUND & DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
  list(N          = length(z), # number of participants
       HIV_N      = sum(DT$HIV_STATUS[z]==1), # number of HIV+
       VLNS_N     = sum(DT$VLNS[z]==1), # number of unsuppressed from viral load
       ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z])) # number of unsuppressed from self-reporting
  )				
}, by=names(vla)]
setnames(vla, 'FC', "COMM")

# get unsuppressed table
uns <- vla[, list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Total'), by = c('COMM', 'ROUND')]
uns <- rbind(uns, vla[SEX == 'F', list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Female'), by = c('COMM', 'ROUND')])
uns <- rbind(uns, vla[SEX == 'M', list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Male'), by = c('COMM', 'ROUND')])
uns <- rbind(uns, vla[SEX == 'F' & AGEYRS < 25, list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
uns <- rbind(uns, vla[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
uns <- rbind(uns, vla[SEX == 'F' & AGEYRS > 34, list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
uns <- rbind(uns, vla[SEX == 'M' & AGEYRS < 25, list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
uns <- rbind(uns, vla[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
uns <- rbind(uns, vla[SEX == 'M' & AGEYRS > 34, list(UNSUPPRESSED = sum(VLNS_N), INFECTED_TESTED = sum(HIV_N),  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])
uns[, ROUND := paste0('R0', ROUND)]



########################################################################

# GET PARTICIPANTS CONSIDERED FOR TRANSMISSION NETWORK RECONSTRUCTION

########################################################################

# load seq count
sequ <- as.data.table(readRDS(file.seq.count))
sequ <- unique(sequ[, .(PT_ID, ROUND)])

# get first round when sequenced
sem <- sequ[, list(MIN_ROUND = min(ROUND)), by = c('PT_ID')]
sem[, STUDY_ID := gsub('RK-(.+)', '\\1', PT_ID)]
set(sem, NULL, 'PT_ID', NULL)

# merge to meta_data
semt <- merge(rincp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS)], sem, by = 'STUDY_ID')
stopifnot(semt[, length(unique(STUDY_ID))] == sem[, length(unique(STUDY_ID))])

# keep only rounds after first sequencing
semt <- semt[ROUND >= MIN_ROUND]
stopifnot(semt[, length(unique(STUDY_ID))] == sem[, length(unique(STUDY_ID))])

# get sequenced table
seqs <- semt[, list(SEQUENCE = .N,  TYPE = 'Total'), by = c('COMM', 'ROUND')]
seqs <- rbind(seqs, semt[SEX == 'F', list(SEQUENCE = .N,  TYPE = 'Female'), by = c('COMM', 'ROUND')])
seqs <- rbind(seqs, semt[SEX == 'M', list(SEQUENCE = .N,  TYPE = 'Male'), by = c('COMM', 'ROUND')])
seqs <- rbind(seqs, semt[SEX == 'F' & AGEYRS < 25, list(SEQUENCE = .N,  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
seqs <- rbind(seqs, semt[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(SEQUENCE = .N, TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
seqs <- rbind(seqs, semt[SEX == 'F' & AGEYRS > 34, list(SEQUENCE = .N, TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
seqs <- rbind(seqs, semt[SEX == 'M' & AGEYRS < 25, list(SEQUENCE = .N, TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
seqs <- rbind(seqs, semt[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(SEQUENCE = .N, TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
seqs <- rbind(seqs, semt[SEX == 'M' & AGEYRS > 34, list(SEQUENCE = .N, TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])



########################

# MAKE TABLE

########################

# merge all information
tab <- merge(census, part, by = c('TYPE', 'COMM', 'ROUND'))
tab <- merge(tab, hivp, by = c('TYPE', 'COMM', 'ROUND'))
tab <- merge(tab, sartp, by = c('TYPE', 'COMM', 'ROUND'))
tab <- merge(tab, uns, by = c('TYPE', 'COMM', 'ROUND'), all.x = T)
tab <- merge(tab, seqs, by = c('TYPE', 'COMM', 'ROUND'), all.x = T)
tab[is.na(tab)] = '--'

# make factor for population categories
tab[, unique(TYPE)]
tab[, TYPE := factor(TYPE, levels = c('Total', 'Female', 'Female, 15-24', "Female, 25-34", "Female, 35-49", 
                                         "Male",  "Male, 15-24", "Male, 25-34", "Male, 35-49"))]
tab <- tab[order(COMM, ROUND, TYPE)]

# check that the counts make sense (e.g., there cannot be more participant than census eligible)
stopifnot(nrow(tab[ELIGIBLE  < PARTICIPANT ]) == 0)
stopifnot(nrow(tab[PARTICIPANT  < HIV ]) == 0)
stopifnot(nrow(tab[HIV  < SELF_REPORTED_ART ]) == 0)
stopifnot(nrow(tab[HIV  < as.numeric(SEQUENCE) & COMM == 'inland']) == 0)
stopifnot(nrow(tab[HIV  < as.numeric(INFECTED_TESTED) ]) == 0)

# add comma thousands separator
comma_thousands <- function(x) format(x, big.mark=",")
tab[, ELIGIBLE := comma_thousands(ELIGIBLE)]
tab[, PARTICIPANT := comma_thousands(PARTICIPANT)]
tab[, HIV := comma_thousands(HIV)]
tab[, INFECTED_TESTED := comma_thousands(INFECTED_TESTED)]
tab[, SELF_REPORTED_ART := comma_thousands(SELF_REPORTED_ART)]
tab[, UNSUPPRESSED := comma_thousands(UNSUPPRESSED)]
tab[, SEQUENCE := comma_thousands(SEQUENCE)]

# save
tab <- tab[, .(COMM, TYPE, ROUND, ELIGIBLE, PARTICIPANT, HIV, INFECTED_TESTED, SELF_REPORTED_ART, UNSUPPRESSED, SEQUENCE)]
saveRDS(tab, file.path(outdir, 'characteristics_study_population.rds'))



