library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library("haven")
library(here)

gitdir <- here()
source(file.path(gitdir, "config.R"))

c(  file.eligible.count,
    file.community.keys ,
    file.eligible.count,
    path.tests,
    file.characteristics_sequenced_ind_R14_18,
    file.path.hiv,
    file.path.quest,
    file.path.metadata) |> file.exists() |> all() |> stopifnot()

# load files
community.keys <- as.data.table(read.csv(file.community.keys))

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 10:18)),
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(15, '15S', 16:18))))

# community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']


## count census eligible and ineligible ----

# load census eligible count
dt <- fread(file.eligible.count)

dt <- subset(dt,COMM=='inland')
dt <- dt[, list(N_census = sum(ELIGIBLE_NOT_SMOOTH + ALREADY_SEEN + DEAD + NOT_WITHIN_ELIGIBLE_AGE_RANGE + OUT_MIGRATED + NOT_RESIDENT),
                N_census_eligible = sum(round(ELIGIBLE_NOT_SMOOTH + OUT_MIGRATED)),
                N_eligible = sum(round(ELIGIBLE_NOT_SMOOTH + (OUT_MIGRATED/2) )),
                N_ineligible = sum(ALREADY_SEEN + DEAD + NOT_WITHIN_ELIGIBLE_AGE_RANGE + OUT_MIGRATED + NOT_RESIDENT),
                N_ineligible_age = sum(NOT_WITHIN_ELIGIBLE_AGE_RANGE),
                N_ineligible_not_res = sum(NOT_RESIDENT),
                N_ineligible_outmig = sum(round(OUT_MIGRATED/2))), by=c('ROUND')]

# NOTE: ELIGIBLE_NOT_SMOOTH already has half the out-migrated added!
dt <- dt[, list(N_census_eligible = round(sum(ELIGIBLE_NOT_SMOOTH + OUT_MIGRATED/2)),
                N_ineligible_outmig = round(sum(OUT_MIGRATED/2)),
                N_eligible = round(sum(ELIGIBLE_NOT_SMOOTH))),
           by=c('ROUND')]
dt[, Tot_census:= N_eligible + N_ineligible_outmig]

####
## get testing data ----

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
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, arvmed, cuarvmed)]

# find  community
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restrict age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# get ART status
rinc[!ROUND %in% c('R016', 'R017', 'R018'), ART := ARVMED ==1]
rinc[!ROUND %in% c('R016', 'R017', 'R018') & is.na(ARVMED), ART := F]
rinc[ROUND == 'R016', ART := ARVMED ==1 | CUARVMED ==1]
rinc[ROUND == 'R016' & (is.na(ARVMED) | is.na(CUARVMED)), ART := F]
rinc[ROUND %in% c('R017', 'R018'), ART := CUARVMED ==1]
rinc[ROUND %in% c('R017', 'R018') & is.na(CUARVMED), ART := F]

# art was not reported in round 10
rinc[ROUND == 'R010', ART := NA]

# add meta data from Kate
tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], rinc[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
rinc <- rbind(rinc[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)])

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

tmp <- hivs[COMM=='inland', list(N_test=sum(.N)),by=c('ROUND','HIV')]
tmp[HIV=='', HIV:='M']
tmp <- dcast(tmp, ROUND~HIV, value.var='N_test')


###### ###### ###### ###### ###### ###### ######
# Get viral load data ----
###### ###### ###### ###### ###### ###### ######

# tuning
VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 1000 # WHO standards

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND %in% c(15:18)] ## AB removed

# just keep inland
dall <- subset(dall,COMM=='inland')

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# count missing VL and not ART naive at any visit from R15+ (think ARVMED means ever on ART)
dall[, HIV_AND_NOVL := ifelse( HIV_STATUS == 1 & is.na(VL_COPIES), 1, 0)]
dall[, HIV_AND_NOVL_AND_NOART := ifelse( HIV_STATUS == 1 & is.na(VL_COPIES) & is.na(ARVMED), 1, 0)]

# count those positive with a missing a VL who are ART naive
vl <- dall[, list(miss_VL=sum(HIV_AND_NOVL_AND_NOART==1)),by='ROUND']

# remove HIV- individuals with missing VLs
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

DT <- subset(DT,FC=='inland')

# count number of individuals who were ever unsuppressed
tmp <- DT[, list(rounds_unsupp=sum(VLNS)),by=c('STUDY_ID')]
tmp[, always_supp:= ifelse(rounds_unsupp==0,1,0)]
DT <- merge(DT,subset(tmp,select=c('STUDY_ID','always_supp')),by=c('STUDY_ID'),all.x=T)

# find HIV+ and detectable
set(DT, NULL, 'HIV_AND_VLD', DT[, as.integer(VLD==1 & HIV_AND_VL==1)])

# reset undetectable to VLC 0
set(DT, DT[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
setkey(DT, ROUND, FC, SEX, AGEYRS)

# keep within census eligible age
DT <- subset(DT, AGEYRS > 14 & AGEYRS < 50)

# keep infected
DT <- DT[HIV_STATUS ==1]

# get number always suppressed R15-18 by round
supp <- DT[, list(N_always_supp=length(unique(STUDY_ID[always_supp==1]))),by='ROUND']


# get transmission cohort (sequenced) ----
# from get_participant_sequenced_count.R
###############################################

# GET SEQUENCES DATA

###############################################

alignment <- read.fasta(file = infile.sequence)
unique(do.call(c, lapply(alignment, unique)))
nsequence <- length(alignment)
npos <- unique(lengths(alignment))

# map alignments to studyid
dinfo <- data.table(pangea_id=names(alignment))
id.dt <- fread(path.sdates.rccs)
id.dt <- subset(id.dt,select = c("pt_id","pangea_id", 'cd4_count', 'visit_dt'))
id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
tmp <- fread(path.sdates.mrc)
tmp <- subset(tmp,select = c("pt_id","pangea_id", 'cd4_count', 'visit_dt'))
tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
id.dt <- rbind(id.dt,tmp)
id.dt <- unique(id.dt)
dinfo <- merge(dinfo, id.dt, by="pangea_id", all.x=T)
tmp <- dinfo[is.na(pt_id)]
dinfo <- dinfo[!is.na(pt_id)]
cat('No personal information found for ',nrow(dinfo[is.na(pt_id)]), ' sequences \n')

# remove mrc cohort
dinfo <- dinfo[!grepl('MRCUVRI', pangea_id)]

# remove neuro data
neuro.metadata <- as.data.table(read.csv(file.path.neuro.metadata))
dinfo <- dinfo[!pt_id %in% neuro.metadata[, paste0('RK-', studyid)]]

# find meta data
hivs[, pt_id := paste0('RK-', STUDY_ID)]
dm <- merge(dinfo, hivs, by = c('pt_id'))
colnames(dm) <- toupper(colnames(dm))
stopifnot(dm[, length(unique(PT_ID))] == dinfo[, length(unique(pt_id))])

# keep sequences which meet minimum criteria
load(infile.seq.criteria)
dct[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
dm <- merge(dm,dct,by='PANGEA_ID',all.x=T)

# AB - just keep inland
dm <- subset(dm,COMM=='inland')

saveRDS(dm,file='sequenced_met_criteria.RDS')

# AB: tabulate seqs not meeting criteria
tmp <- dm[V1=='FALSE', list(N=length(unique(STUDY_ID))),by=c('ROUND','V1')]

# just keep those meeting min criteria
dm <- subset(dm,V1=='TRUE')


# keep meta info closer to sample date
dm[, VISIT_DT := as.Date(VISIT_DT)]
dm[, SAMPLE_DATE := as.Date(SAMPLE_DATE)]
dm[, DIFF_DATE := abs(VISIT_DT - SAMPLE_DATE), by = 'PANGEA_ID']
dm[, IS_MIN := DIFF_DATE == min(na.omit(DIFF_DATE)), by = 'PANGEA_ID']
dcount <- dm[IS_MIN == 1]
stopifnot(nrow(dcount) == dm[, length(unique(PANGEA_ID))])
dcount[, table(ROUND, COMM)]

dcount[HIV == 'N'] ## negative but sequenceD?

# set round to 15 if inland 15S
dcount[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'pt_id']
dcount[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
dcount <- dcount[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]
set(dcount, NULL, 'PARTICIPATED_TO_ROUND_RO15', NULL)

# keep round of interest
dcount <- merge(dcount, df_round, by = c('COMM', 'ROUND'))

# create age groups
dcount[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,49),include.lowest=T,right=F,
                     labels=c('15-24','25-34','35-49'))]

# keep necessary variable
dcount <- dcount[, .(COMM, ROUND, PANGEA_ID, PT_ID, STUDY_ID, SEX, AGEYRS, AGEGP, DIFF_DATE)]

# save sequenced id
file.name <- file.characteristics_sequenced_ind_R14_18
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(dcount, file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}

########################################################################

# GET PARTICIPANTS CONSIDERED FOR TRANSMISSION NETWORK RECONSTRUCTION

########################################################################

# load seq count
sequ <- as.data.table(readRDS(file.characteristics_sequenced_ind_R14_18))

# keep age within 15-49
sequ <- sequ[AGEYRS > 14 & AGEYRS < 50]

# merge to meta_data
semt <- merge(hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV)], unique(sequ[, .(STUDY_ID)]), by = 'STUDY_ID')
stopifnot(semt[, length(unique(STUDY_ID))] == sequ[, length(unique(STUDY_ID))])

# keep only positive participants
semt <- semt[HIV == 'P']

# get sequenced table
seqs <- semt[, list(SEQUENCE = .N,  TYPE = 'Total'), by = c('COMM', 'ROUND')]
