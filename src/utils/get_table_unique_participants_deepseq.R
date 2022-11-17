library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library(haven)

indir.deepsequencedata <- '~/OneDrive - Imperial College London/PANGEA/ratmann_pangea_deepsequencedata/live'
indir.deepsequence_analyses <- '~/OneDrive - Imperial College London/PANGEA/ratmann_deepseq_analyses/live'
indir.repository <- '~/Documents/GitHub/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

file.seq.count <- file.path(outdir, 'characteristics_sequenced_R14_18.rds')
file.seq.count.ind <- file.path(outdir, 'characteristics_sequenced_ind_R14_18.rds')

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903', 'HIV_R6_R18_220909.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903', 'Quest_R6_R18_220909.csv')

# Latest data from Rakai's CCS (Kate's data from 2022-03-08)
file.path.metadata <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata__12Nov2019.csv')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 14:18)),
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(14, '15S', 16:18))))

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
rincp[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
rincp <- rincp[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]


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

# keep rounds of interest and only positives
hivsp <- merge(hivs, df_round, by = c('COMM', 'ROUND'))
hivsp <- subset(hivsp,HIV=='P')

# group into analysis age groups
hivsp[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                   labels=c('15-24','25-34','35-49'))]

# save individual-level data
hivi <- copy(hivsp)

# keep first round only so age is unique per participant
hivsp[, r:= as.numeric(gsub('R','',gsub('S','.1',ROUND)))]
setkey(hivsp,STUDY_ID,r)
hivsp <- hivsp[,.SD[1],by = STUDY_ID]

# count unique participants with positive test during rounds
hivp <- hivsp[,list(HIV = length(unique(STUDY_ID))), by = c('COMM','SEX','AGEGP')]

tot1 <- hivsp[,list(SEX = 'Total', AGEGP = 'Total', HIV = length(unique(STUDY_ID))), by = c('COMM')]
tot2 <- hivsp[,list(AGEGP = 'Total', HIV = length(unique(STUDY_ID))), by = c('COMM','SEX')]
hivp <- rbind(hivp,tot1,tot2)

hivp[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
hivp[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
hivp[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
hivp <- hivp[order(COMM,SEX,AGEGP),]


######################################################

# GET HIV-POSITIVE AND USING ART AMONG PARTICIPANTS #

######################################################

# keep HIV positive and find art status
sart <- hivs[HIV == 'P']

# age groups
sart[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                   labels=c('15-24','25-34','35-49'))]

# keep round of interest
sart <- merge(sart, df_round, by = c('COMM', 'ROUND'))

# save individual-level record
arti <- copy(sart)

# keep first visit only
sart[, r:= as.numeric(gsub('R','',gsub('S','.1',ROUND)))]
setkey(sart,STUDY_ID,r)
sart <- sart[,.SD[1],by = STUDY_ID]

art <- sart[,list(NO_ART = length(unique(STUDY_ID[ART==F]))), by = c('COMM','SEX','AGEGP')]

tot1 <- sart[,list(SEX = 'Total', AGEGP = 'Total', NO_ART = length(unique(STUDY_ID[ART==F]))), by = c('COMM')]
tot2 <- sart[,list(AGEGP = 'Total', NO_ART = length(unique(STUDY_ID[ART==F]))), by = c('COMM','SEX')]
art <- rbind(art,tot1,tot2)

art[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
art[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
art[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]


########################################################################

# GET PARTICIPANTS CONSIDERED FOR TRANSMISSION NETWORK RECONSTRUCTION

########################################################################

# load seq count
sequ <- as.data.table(readRDS(file.seq.count.ind))
sequ[, STUDY_ID:= gsub('RK-','',PT_ID)]

# merge three datasets
do <- merge(arti,sequ,by=c('STUDY_ID','ROUND','COMM','SEX','AGEYRS','AGEGP'),all=T)

# keep round of interest
do <- merge(do, df_round, by = c('COMM', 'ROUND'))

# count unique participants with positive test during rounds
tab <- do[,list(HIV = length(unique(STUDY_ID[HIV=='P'])),
               NO_ART = length(unique(STUDY_ID[ART==F])),
               SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX','AGEGP')]

tot1 <- do[,list(SEX = 'Total', AGEGP = 'Total',
                 HIV = length(unique(STUDY_ID[HIV=='P'])),
                 NO_ART = length(unique(STUDY_ID[ART==F])),
                 SEQUENCE = length(unique(PT_ID))), by = c('COMM')]
tot2 <- do[,list(AGEGP = 'Total', HIV = length(unique(STUDY_ID[HIV=='P'])),
                 NO_ART = length(unique(STUDY_ID[ART==F])),
                 SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX')]
tab <- rbind(tab,tot1,tot2)

tab[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
tab[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
tab[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
tab <- tab[order(COMM,SEX,AGEGP),]
tab[, c("COMM", "SEX", "AGEGP") := lapply(list(COMM, SEX, AGEGP), as.character)]

tab[, pct:= round(SEQUENCE/HIV*100,0)]

# save
saveRDS(tab, file.path(outdir, 'RCCS_transmission_cohort_characteristics_R14_18.rds'))

# check
`%notin%` <- Negate(`%in%`)
tmp <- sequ$STUDY_ID[sequ$STUDY_ID %notin% hivs$STUDY_ID & sequ$COMM=='inland']
tmp # no individual missing
tmp <- sequ$STUDY_ID[sequ$STUDY_ID %notin% hivi$STUDY_ID & sequ$COMM=='inland']
hivs[STUDY_ID %in% tmp] # negative during round of interest
