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

# specify outdir
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

# load files
community.keys <- fread(file.community.keys)

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 10:18)), 
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(15, '15S', 16:18))))

# community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']


# 
# Helpers 
# 


.make.plot <- function(DT, .y, .ylab )
{
    .y <- enquo(.y)

    add_linerange <- FALSE
    if('CL' %in% names(DT))
        add_linerange <- TRUE

    ggplot(DT, aes(x=ROUND_LAB, color=SEX_LAB, pch=AGEGP, linetype=AGEGP, y=!!.y)) + 
        geom_point(position=position_dodge(width=.6) ) + {
            if(add_linerange)
                geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width =.6))
        } +
        scale_y_continuous(limits = c(0,1), expand=c(0,0),labels=scales::percent) +
        scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
        labs( 
            x="Interview round",
            y=.ylab,
            linetype = "Age group", 
            pch = "Age group", 
            color = "Gender"
        ) +
        theme_bw() + 
        theme(legend.position = "bottom", strip.background = element_blank()) + 
        reqs
}


# Get metadata and HIV+ ----

#
# Meta data
#

meta_data <- fread(file.path.metadata) #additional meta_data

# find age
meta_data[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
meta_data[, AGEYRS := round(lubridate::time_length(difftime(sample_date, date_birth),"years"))]
meta_data[is.na(AGEYRS), AGEYRS := round(lubridate::time_length(difftime(firstposvd, date_birth),"years"))]
meta_data[is.na(AGEYRS)]

# restrict age
meta_data <- meta_data[AGEYRS %between% c(15,49)]

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

quest <- fread(file.path.quest)

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

hiv <- fread(file.path.hiv)
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
hivs <- hivs[AGEYRS %between% c(15, 49)]


## calculate denominators (number infected and unsuppressed) ----


# load participation (% of census eligible population)
participation <- fread(file.participation)

# load non-suppressed proportion 
treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants, file.treatment.cascade.prop.nonparticipants)

# load round timeline
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]

# make df round
df_round <- make.df.round(df_round_inland)


# get counts among census eligible
eligible_count_round <- add_susceptible_infected(
    eligible_count_smooth = fread(file.eligible.count),
    proportion_prevalence = fread(file.prevalence.prop),
    participation = participation,
    nonparticipants.male.relative.infection = nonparticipants.male.relative.infection,
    nonparticipants.female.relative.infection = nonparticipants.female.relative.infection 
) |> add_infected_unsuppressed(
        treatment_cascade,
        participation,
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
sequ <- sequ[AGEYRS %between% c(15,49) ]

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

tmp <- merge( eligible_count_round, seqs,  by=c('COMM','ROUND','SEX','AGEGP'))

# save saveRDS(tmp,file=file.path(outdir,'prop_unsupp_eventuallydeepseq_byroundagesex.RDS'))

# plot

tmp[, (c('P_EVERSEQ', 'CL', 'CU')) := Hmisc::binconf(x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED, return.df=TRUE) ]
tmp[, ROUND_LAB := gsub("R0", "Round ", ROUND)]
tmp[, SEX_LAB := fifelse(SEX=="F", "Women", "Men")] 
tab_seq_unsup <- copy(tmp)


# Make Figure

p_everseq_givenunspp <- .make.plot(tmp, .y=P_EVERSEQ, .ylab = "Proportion of unsuppressed participants eventually deep-sequenced")

filename <- file.path(outdir, "prop_unsupp_eventuallydeepseq_byroundagesex.pdf") 
filename <- file.path(outdir, "prop_unsupp_eventuallydeepseq_byroundagesex.png") 
ggsave_nature(filename=filename, p=p_everseq_givenunspp, w=12, h=10)


# what is the diff with art prop?
# -------------------------------

hivps <- hivs |> subset(HIV == 'P')
hivps[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F, labels=c('15-24','25-34','35-49'))]
hivps[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
hivps <- hivps[,.(
    N_HIV = .N,
    N_noART = sum(! ART, na.rm=TRUE)
), by=c("ROUND", 'COMM',"SEX","AGEGP")]
tmp_2 <- merge(
    tmp,
    hivps,
    by=c("ROUND", 'COMM',"SEX","AGEGP")
)
tmp_2[,  P_EVERSEQ_OLD := P_EVERSEQ, ]
tmp_2[, (c('P_EVERSEQ', 'CL', 'CU')) := Hmisc::binconf(x=N_EVERSEQ, n=N_noART, return.df=TRUE) ]
p_everseq_givennoart <- .make.plot(tmp_2, .y=P_EVERSEQ, .ylab="Proportion of ART naive participants eventually deep-sequenced")

###################
# get known unsup #
###################

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
DT <- subset(DT, AGEYRS > 14 & AGEYRS < 50)

# keep infected
DT <- DT[HIV_STATUS ==1]

# use age from quest
DT[, round := paste0('R0', ROUND)]
set(DT, NULL, 'AGEYRS', NULL)
# DT <- merge(DT, hiv, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
DT <- merge(DT, quest, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
setnames(DT, 'ageyrs', 'AGEYRS')

# get count for every categories
tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
tmp1 <- DT[, sort(unique(ROUND))]

# all those for which we know they were supp
unsup_id <- subset(DT, VLNS==1, select=c(round, STUDY_ID, SEX, FC,  HIV_AND_VL, VLNS, AGEYRS))  |>
    setnames(c('round','FC'), c('ROUND', 'COMM'))

unsup_id[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                   labels=c('15-24','25-34','35-49'))]
unsup_id[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
unsup_id[, AGEYRS := NULL ]

unsup_and_seq <- merge(unsup_id, sequ, by=c('STUDY_ID', 'ROUND', 'SEX', 'COMM', 'AGEGP'), all.x=TRUE, all.y=TRUE)

unsup_and_seq[, table(is.na(VLNS))]
unsup_and_seq[, table(ROUND, is.na(VLNS))]

n_unsup_and_seq <- unsup_and_seq[!is.na(VLNS), list(
    UNSUP=.N,
    EVENTUALLY_SUPPRESSED = sum(!is.na(PANGEA_ID))
), by=c('COMM','ROUND', 'SEX', 'AGEGP')] |> 
   subset(!is.na(COMM))

n_unsup_and_seq[, (c('P_EVERSEQ', 'CL', 'CU')) := Hmisc::binconf(x=EVENTUALLY_SUPPRESSED, n=UNSUP, return.df=TRUE) ]
n_unsup_and_seq[, SEX_LAB := fifelse(SEX=="F", "Women", "Men")] 
n_unsup_and_seq[, ROUND_LAB := gsub('R0','Round ',ROUND)] 
# setnames(n_unsup_and_seq, c('ROUND'), c('ROUND_LAB'), skip_absent =TRUE)

p_everseq_givenunspp_participants <- .make.plot(n_unsup_and_seq[COMM=='inland'], .y=P_EVERSEQ, .ylab='Proportion of known unsuppressed who are eventually deep-sequenced')
filename <- file.path(outdir, "prop_participant_unsupp_eventuallydeepseq_byroundagesex.png") 
filename <- file.path(outdir, "prop_participant_unsupp_eventuallydeepseq_byroundagesex.pdf") 
ggsave_nature(filename=filename, p=p_everseq_givenunspp_participants, w=12, h=11)


by_cols <- c('COMM', 'ROUND_LAB', 'SEX_LAB', 'AGEGP')

props_parts <- merge(
    tmp |> subset(select=c(by_cols, 'INFECTED_NON_SUPPRESSED')),
    n_unsup_and_seq |> subset(select=c(by_cols, 'UNSUP')),
    by=by_cols
)
# props_parts[, P := UNSUP/INFECTED_NON_SUPPRESSED,]
props_parts[, (c('P', 'CL', 'CU')) := Hmisc::binconf(x=UNSUP, n=INFECTED_NON_SUPPRESSED, return.df=TRUE) ]
p_prop_unsup_particpant <- .make.plot(props_parts, .y=P, .ylab='Proportion of known unsuppressed out of estimated unsuppressed size')
filename <- file.path(outdir, "prop_knwonunsup_among_estunsup_byroundagesex.pdf") 
filename <- file.path(outdir, "prop_knwonunsup_among_estunsup_byroundagesex.png") 
ggsave_nature(filename=filename, p=p_prop_unsup_particpant, w=12, h=11)


merge(
    merge( tmp_2, props_parts, by=by_cols)[, P_EVERSEQ * P , by=by_cols],
    tmp,
    by=by_cols
) |> ggplot(aes(x=ROUND, color=SEX_LAB, pch=AGEGP)) + 
    geom_point(aes(y=V1), position=position_dodge(width=.6) ) +
    geom_point(aes(y=P_EVERSEQ), position=position_dodge(width=.6), alpha=.5 ) +
    scale_y_continuous(limits = c(0,1), expand=c(0,0),labels=scales::percent) +
    theme_bw() 

.make.plot(tmp[ROUND %like% '15|16|17|18|19'], .y=P_EVERSEQ, .ylab = "Proportion of unsuppressed participants eventually deep-sequenced")


tmp-2
p_pro
