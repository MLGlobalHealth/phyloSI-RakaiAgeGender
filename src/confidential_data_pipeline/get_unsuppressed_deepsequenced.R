{
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(lubridate)
  library(rstan)
  library(haven)
} |> suppressPackageStartupMessages()


usr <- Sys.info()[['user']]

if(usr == 'alexb'){
  
  indir.deepsequencedata <- '~/OneDrive - Imperial College London/PANGEA/ratmann_pangea_deepsequencedata/live'
  indir.deepsequence_analyses <- '~/OneDrive - Imperial College London/PANGEA/ratmann_deepseq_analyses/live'
  indir.repository <- '~/Documents/GitHub/phyloflows'
  
}else if(usr == 'melodie'){
  
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
  indir.repository <- '~/git/phyloflows'
  
}else if(usr == 'andrea'){
  
  indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
  indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live'
  indir.repository <- '~/git/phyloflows'
  
}

gitdir <- here::here()

source(file.path(gitdir, "config.R"))

# load functions
source(file.path(gitdir.R.flow, 'summary_functions.R'))
source(file.path(gitdir.R.flow, 'plotting_functions.R') )
naturemed_reqs()

# indicators for sensitivity analyses: default values
nonparticipants.treated.like.participants <- FALSE
nonparticipants.not.treated <- FALSE
nonparticipants.male.relative.infection <- 1
nonparticipants.female.relative.infection <- 1

# check files exist
file.exists(c(
  file.treatment.cascade.prop.participants,
  file.treatment.cascade.prop.nonparticipants,
  file.eligible.count,
  file.participation,
  file.prevalence.prop,
  file.path.round.timeline,
  file.community.keys,
  file.characteristics_sequenced_ind_R14_18,
  file.path.hiv,
  file.path.quest,
  file.path.metadata,
  file.characteristics_sequenced_ind_R14_18))  |> all() |> stopifnot()

# load files
community.keys <- as.data.table(read.csv(file.community.keys))

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 10:18)), 
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(15, '15S', 16:18))))

# community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']

# Get metadata and HIV+ ----

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


## calculate denominators (number infected and unsuppressed) ----

# load census eligible ount
eligible_count_smooth <- fread(file.eligible.count)

# load participation (% of census eligible population)
participation <- fread(file.participation)

# load proportion prevalence
proportion_prevalence <- fread(file.prevalence.prop)

# load non-suppressed proportion 
treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants, file.treatment.cascade.prop.nonparticipants)

# load round timeline
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]

# make df round
df_round <- make.df.round(df_round_inland)

# get number of unsuppressed in the poulation
eligible_count_round <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence, participation, 
                                                 nonparticipants.male.relative.infection, 
                                                 nonparticipants.female.relative.infection)
eligible_count_round <- add_infected_unsuppressed(eligible_count_round, treatment_cascade, participation, 
                                                  nonparticipants.treated.like.participants, 
                                                  nonparticipants.not.treated)

# aggregate counts by age groups
eligible_count_round[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                   labels=c('15-24','25-34','35-49'))]
eligible_count_round[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]

eligible_count_round <- eligible_count_round[, list(
  INFECTED_NON_SUPPRESSED = sum(INFECTED_NON_SUPPRESSED)),
  by=c('COMM','AGEGP','SEX','ROUND')]

## get numerator (ever deep sequenced) ----

# load seq count
sequ <- as.data.table(readRDS(file.characteristics_sequenced_ind_R14_18))

# keep age within 15-49
sequ <- sequ[AGEYRS > 14 & AGEYRS < 50]

sequ[, round:= gsub('S','',ROUND)]
sequ[, round:= as.integer(gsub('R0','',round))]

# last sequenced
ls <- sequ[, list(N=length(ROUND),MAXROUND=max(round)),by=c('STUDY_ID')]
sequ <- merge(sequ,ls,by=c('STUDY_ID'))

# merge to meta_data
semt <- merge(hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV)], unique(sequ[, .(STUDY_ID,MAXROUND)]), by = 'STUDY_ID')
stopifnot(semt[, length(unique(STUDY_ID))] == sequ[, length(unique(STUDY_ID))])

# keep only positive participants
semt <- semt[HIV == 'P']

# only count participants who were sequenced at current or future round
semt[, round:= gsub('S','',ROUND)]
semt[, round:= as.integer(gsub('R0','',round))]
semt <- subset(semt,round<=MAXROUND)

# keep only rounds 10-18
semt <- semt[ROUND %in% c('R010','R011','R012','R013','R014','R015','R016','R017','R018')]

# make age groups
semt[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                   labels=c('15-24','25-34','35-49'))]
semt[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]

seqs <- semt[COMM=='inland', list(N_EVERSEQ= .N), by=c('COMM','ROUND','SEX','AGEGP')]

# merge together

tmp <- merge(eligible_count_round,seqs,by=c('COMM','ROUND','SEX','AGEGP'))

# save
saveRDS(tmp,file=file.path(outdir,'prop_unsupp_eventuallydeepseq_byroundagesex.RDS'))

# plot

tmp[, (c('P_EVERSEQ', 'CL', 'CU')) := Hmisc::binconf(x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED, return.df=TRUE) ]
tmp[, ROUND_LAB := gsub("R0", "Round ", ROUND)]
tmp[, SEX_LAB := fifelse(SEX=="F", "Women", "Men")]
tab_seq_unsup <- copy(tmp)

# Make Figure

p_everseq_givenunspp <- ggplot(tmp, aes(x=ROUND_LAB, color=SEX_LAB, pch=AGEGP, linetype=AGEGP)) + 
  geom_point(aes(y=P_EVERSEQ), position=position_dodge(width=.6) ) +
  geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width =.6)) +
  # facet_grid(ROUND ~ .) +
  scale_y_continuous(limits = c(0,1), expand=c(0,0),labels=scales::percent) +
  scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
  labs( 
    x="Interview round",
    y="Proportion of unsuppressed participants eventually deep-sequenced", 
    linetype = "Age group", 
    pch = "Age group", 
    color = "Gender"
  ) +
  theme_bw() + 
  theme(legend.position = "bottom", strip.background = element_blank()) + 
  reqs

filename <- file.path(outdir, "prop_unsupp_eventuallydeepseq_byroundagesex.pdf") 
filename <- file.path(outdir, "prop_unsupp_eventuallydeepseq_byroundagesex.png") 
ggsave_nature(filename=filename, p=p_everseq_givenunspp, w=12, h=10)

