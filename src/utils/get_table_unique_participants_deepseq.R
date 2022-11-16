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

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903', 'HIV_R6_R18_220909.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903', 'Quest_R6_R18_220909.csv')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 14:18)),
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(14, '15S', 16:18))))

#################################

# GET PARTICIPANTS #

#################################

# load datasets
quest <- as.data.table(read.csv(file.path.quest))

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate)]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# find participant
part <- rinc[, list(PARTICIPANT = .N,  TYPE = 'Total'), by = c('COMM', 'ROUND')]
part <- rbind(part, rinc[SEX == 'F', list(PARTICIPANT = .N,  TYPE = 'Female'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M', list(PARTICIPANT = .N,  TYPE = 'Male'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'F' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'F' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
part <- merge(part, df_round, by = c('COMM', 'ROUND'))

########################################

# GET HIV-POSITIVE AMONG PARTICIPANTS #

########################################

# load datasets
hiv <- as.data.table(read.csv(file.path.hiv))
hiv[, round := gsub(' ', '', round)] # remove space in string

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# keep rounds of interest and only positives
hivs <- merge(hivs, df_round, by = c('COMM', 'ROUND'))
hivs <- subset(hivs,HIV=='P')

# keep first round only so age is unique per participant
hivs[, r:= as.numeric(gsub('R','',gsub('S','.1',ROUND)))]
setkey(hivs,STUDY_ID,r)
hivs <- hivs[,.SD[1],by = STUDY_ID]

# count unique participants with positive test during rounds
hivs[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                     labels=c('15-24','25-34','35-49'))]
hivp <- hivs[,list(HIV = length(unique(STUDY_ID))), by = c('COMM','SEX','AGEGP')]

tot1 <- hivs[,list(SEX = 'Total', AGEGP = 'Total', HIV = length(unique(STUDY_ID))), by = c('COMM')]
tot2 <- hivs[,list(AGEGP = 'Total', HIV = length(unique(STUDY_ID))), by = c('COMM','SEX')]
hivp <- rbind(hivp,tot1,tot2)

hivp[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
hivp[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
hivp[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
hivp <- hivp[order(COMM,SEX,AGEGP),]



########################################################################

# GET PARTICIPANTS CONSIDERED FOR TRANSMISSION NETWORK RECONSTRUCTION

########################################################################

# load seq count
sequ <- as.data.table(readRDS(file.seq.count))

# merge datasets
do <- merge(hivp,sequ)

# order rows
do <- do[order(COMM,SEX,AGEGP),]
do[, c("COMM", "SEX", "AGEGP") := lapply(list(COMM, SEX, AGEGP), as.character)]

# save
saveRDS(do, file.path(outdir, 'RCCS_transmission_cohort_characteristics_R14_18.rds'))
