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

# NOTE: ELIGIBLE_NOT_SMOOTH already has half the out-migrated added so only add half to get census eligible
dt <- dt[, list(N_census_eligible = round(sum(ELIGIBLE_NOT_SMOOTH + OUT_MIGRATED/2)),
                N_ineligible_outmig = round(sum(OUT_MIGRATED/2)),
                N_eligible = round(sum(ELIGIBLE_NOT_SMOOTH))),
           by=c('ROUND')]
dt[, Tot_census:= N_eligible + N_ineligible_outmig]

# census eligible population
dt

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

# save community IDs
comm_ids_r <- subset(meta_data,select=c('study_id','round','community_number'))

# find community
meta_data[, COMM := 'inland']
meta_data[LakeVictoria_FishingCommunity == 'yes', COMM := 'fishing']

# find hiv status
meta_data[, HIV := ifelse(is.na(firstposvd), 'N', 'P')]

# find art use
meta_data[, ART := artslfuse == 'yes']

# subset round 14 to 30 continuously surveyed communities
idr14not30comm <- meta_data[(round=='14' & !COMM %in% c(1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 89, 94, 106, 107, 108, 120, 391, 602, 754)), unique(study_id)]
meta_data <- meta_data[!(round=='14' & !community_number %in% c(0, 1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 89, 94, 106, 107, 108, 120, 391, 602, 754))]

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

## count sequenced #################################

# patient-level data on positive participants
hpos <- subset(hivs,COMM=='inland')
hpos[, list(PARTICIPANT = length(unique(STUDY_ID[HIV=='P']))), by = c('COMM', 'ROUND')]
# NB ART naive comes from metadata

# merge in VL data
DT[, ROUND:=paste0('R0',ROUND)]
dat <- merge(subset(hpos,HIV=='P'),subset(DT,select=c('STUDY_ID','ROUND','ARVMED','VLC','VLU','VLD','VLNS','HIV_AND_VLD')),
             by=c('STUDY_ID','ROUND'),all=T)
# find those who were suppressed at every round
tmp <- dat[!is.na(VLC), list(rounds_unsupp=sum(VLNS,na.rm=T)),by=c('STUDY_ID')]
tmp[, always_supp:= ifelse(rounds_unsupp==0,1,0)] # 1==individual was never unsuppressed

dat <- merge(dat,subset(tmp,select=c('STUDY_ID','always_supp')),by=c('STUDY_ID'),all=T)

# count number with (no VL and not ART naive) at R15 or later
#tmp <- dat[, list(N_VL=length(ROUND[!is.na(VLC)]),
#                  ART_R15_onw=(sum(ARVMED,na.rm=T))),by=c('STUDY_ID')]
#tmp[ART_R15_onw>0, ART_R15_onw:= 1]
#dat <- merge(dat,tmp,by=c('STUDY_ID'),all=T)
#dat[, HIV_AND_NOVL_AND_ART := ifelse( HIV == 'P' & is.na(VLC) & ART_R15_onw==1 & always_supp!=1, 1, 0)]
# generate variable for no VL measurement and not ART naive for EVERY round
# generate a flag if this is 1 at any round 15-18 PER study ID
#dat[, HIVP_NOVL_AND_ART:= ifelse( HIV == 'P' & N_VL==0 & ARVMED==1, 1, 0)]
tmp <- dat[ROUND %in% c('R015','R016','R017','R018'),
           list(N_VL_AND_NOART=sum(is.na(VLC) & ART==T),
                PART_R15_ONW=1),
           by=c('STUDY_ID')]
dat <- merge(dat,tmp,by=c('STUDY_ID'),all.x=T)
dat[N_VL_AND_NOART>0, N_VL_AND_NOART:=1]
dat[N_VL_AND_NOART==1 & always_supp==1, always_supp:=0]

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
semt <- subset(semt,COMM=='inland')
semt[, trsm_cohort:= 1]

dt <- merge(dat,subset(semt,select=c('STUDY_ID','ROUND','trsm_cohort')),by=c('STUDY_ID','ROUND'),all.x=T)
dt[is.na(trsm_cohort), trsm_cohort:=0]
dt <- merge(dt, df_round, by = c('COMM', 'ROUND'))

# summarise transmission cohort and exclusions
dt[, list(POSITIVE = length(na.omit(unique(STUDY_ID))),
          N_not_seen_since_R15 = length(na.omit(unique(STUDY_ID[is.na(PART_R15_ONW) & trsm_cohort==0 & (always_supp!=1 | is.na(always_supp)) & (N_VL_AND_NOART!=1 | is.na(N_VL_AND_NOART))]))),
          N_always_supp=length(na.omit(unique(STUDY_ID[always_supp==1 & trsm_cohort!=1]))),
          N_noVL_and_not_ARTnaive=length(na.omit(unique(STUDY_ID[N_VL_AND_NOART==1 & trsm_cohort!=1]))),
          N_not_met_criteria=length(na.omit(unique(STUDY_ID[trsm_cohort==0 & (always_supp!=1 | is.na(always_supp)) & N_VL_AND_NOART!=1]))),
          N_in_trsm_cohort=length(na.omit(unique(STUDY_ID[trsm_cohort==1])))), by = c('COMM', 'ROUND')]

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


file.name <- file.path(dirname(file.characteristics_sequenced_ind_R14_18),'transm_cohort_230926.RDS')
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(semt, file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}

########################################################################

# GET PARTICIPANTS IN TRANSMISSION NETWORKS AND INCIDENCE COHORT

########################################################################

file.name <- file.path(dirname(file.characteristics_sequenced_ind_R14_18),'transm_cohort_230926.RDS')

semt <- readRDS(file.name)
semt[, transm:=1]
semt <- subset(semt,COMM=='inland')

comm_ids_r[, ROUND:= paste0('R0',round)]
setnames(comm_ids_r,'study_id','STUDY_ID')
semt <- merge(semt,subset(comm_ids_r,ROUND=='R014'),by=c('STUDY_ID','ROUND'),all.x=T)

# subset round 14 to 30 continuously surveyed communities
semt <- semt[!(ROUND=='R014' & !is.na(community_number) & !community_number %in% c(0, 1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 89, 94, 106, 107, 108, 120, 391, 602, 754))]


file.name <- file.path(dirname(file.characteristics_sequenced_ind_R14_18),'flow_chart_IDs_in_incidence_cohort_230926.csv')

incid <- fread(file.name)
setnames(incid,c('visit','ROUND','research_id'),c('ROUND','round_long','STUDY_ID'))

dm <- merge(semt,
              subset(incid,select=c('STUDY_ID','ROUND','ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS'),ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS==1),
              by=c('STUDY_ID','ROUND'),all=T)
dm[,incid_not_transm:= ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS==1 & is.na(transm)]
dm[,transm_not_incid:= transm==1 & is.na(ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS)]
dm[,transm_and_incid:= transm==1 & ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS==1]

dm <- subset(dm, ROUND %in% c('R010','R011','R012','R013','R014','R015',
                              'R016','R017','R018'))
tmp <- dm[, list(incidence_cohort = length(na.omit(unique(STUDY_ID[ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS==1]))),
                  transmission_cohort = length(na.omit(unique(STUDY_ID[transm==1]))),
                  incid_not_transm = length(na.omit(unique(STUDY_ID[incid_not_transm==T]))),
                  transm_not_incid = length(na.omit(unique(STUDY_ID[transm_not_incid==T]))),
                  transm_and_incid = length(na.omit(unique(STUDY_ID[transm_and_incid==T])))), by = 'ROUND']

tmp[, ROUND:= factor(ROUND,levels=c('R010','R011','R012','R013','R014','R015',
                                    'R016','R017','R018'))]
                       
tmp <- tmp[order(ROUND),]

file.name <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "participants_count_by_gender_loc_age","flow_charts_incidence_transm_cohort_round_250926.csv")
write.csv(tmp,file=file.name)

tmp2 <- dm[, list(incidence_cohort = length(na.omit(unique(STUDY_ID[ALL_HIV_INCLUDED_IN_INCIDENCE_ANALYSIS==1]))),
                 transmission_cohort = length(na.omit(unique(STUDY_ID[transm==1]))),
                 incid_not_transm = length(na.omit(unique(STUDY_ID[incid_not_transm==T]))),
                 transm_not_incid = length(na.omit(unique(STUDY_ID[transm_not_incid==T]))),
                  transm_and_incid = length(na.omit(unique(STUDY_ID[transm_and_incid==T]))))]

file.name <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "participants_count_by_gender_loc_age","RCCS_incidence_transm_cohorts_unique_230926.csv")
write.csv(tmp2,file=file.name)
saveRDS(tmp2,file=gsub('.csv','.rds',file.name))

