{
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggpattern)
    library(ggpubr)
    library(scales)
    library(lubridate)
    library(rstan)
    library(haven)
    # single imports
    binconf <- Hmisc::binconf
} |> suppressPackageStartupMessages()


usr <- Sys.info()[['user']]

gitdir <- here::here()
source(file.path(gitdir, "config.R"))
source(file.path(gitdir.R, 'utils.R'))

# load functions
source(file.path(gitdir.R.flow, 'summary_functions.R'))
source(file.path(gitdir.R.flow, 'plotting_functions.R') )
naturemed_reqs()

# indicators for sensitivity analyses: default values
nonparticipants.treated.like.participants <- FALSE
nonparticipants.not.treated <- FALSE
nonparticipants.male.relative.infection <- 1
nonparticipants.female.relative.infection <- 1
sequenced.at.round.or.later <- FALSE

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

# load round timeline
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]
df_round <- make.df.round(df_round_inland)

# community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']


# 
# Helpers 
# 


ageyrs2agegp = function(x){
    cut(x,breaks=c(15,25,35,50),include.lowest=T,right=F, labels=c('15-24','25-34','35-49')) |> 
        factor(levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))
}

load.viralsuppression.test.results <- function(QUEST=quest)
{
    # VL_DETECTABLE = 400 # Not needed
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
    DT[,  VLC := fifelse(HIV_STATUS == 0, yes=0, no=VL_COPIES)]

    # define suppressed VL as VLS and unsuppressed as VLNS (according to WHO criteria)	
    DT[, `:=` (
        VLS = as.integer(VLC<VIREMIC_VIRAL_LOAD),
        VLNS = as.integer(VLC>=VIREMIC_VIRAL_LOAD)
    )]

    # reset undetectable to VLC 0
    setkey(DT, ROUND, FC, SEX, AGEYRS)

    # keep within census eligible age
    DT <- subset(DT, AGEYRS > 14 & AGEYRS < 50)

    # keep infected
    DT <- DT[HIV_STATUS ==1]

    # use age from quest
    DT[, round := paste0('R0', ROUND)]
    set(DT, NULL, 'AGEYRS', NULL)
    # DT <- merge(DT, hiv, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
    DT <- merge(DT, QUEST, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
    setnames(DT, 'ageyrs', 'AGEYRS')
}

prettify_labs <- function(DT){
    nms <- names(DT)

    if (! 'ROUND_LAB' %in% nms & 'ROUND' %in% nms )
        DT[, ROUND_LAB := gsub("R0", "Round ", ROUND)]

    if (! 'SEX_LAB' %in% nms)
        DT[, SEX_LAB := fifelse(SEX=="F", "Women", "Men")] 

    return(DT)
}

.make.plot.with.binconf <- function(DT, x, n, .ylab, ylims = c(0,1) )
{
    x <- enexpr(x)
    n <- enexpr(n)

    dplot <- copy(DT)
    dplot[, (c('P', 'CL', 'CU')) := binconf(x=eval(x), n=eval(n), return.df=TRUE) ]
    dplot[, ROUND_LAB := gsub('^R0','Round', ROUND_LAB)]

    ggplot(dplot, aes(x=ROUND_LAB, color=SEX_LAB, pch=AGEGP, linetype=AGEGP, y=P )) + 
        geom_hline(yintercept = 1, linetype='dashed', color='grey50') +
        geom_point(position=position_dodge(width=.8) ) + 
        geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width =.8)) +
        scale_y_continuous(limits = ylims, expand=c(0,0),labels=scales::percent) +
        scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
        labs( 
            x = NULL,
            y = .ylab,
            linetype = "Age", 
            pch = "Age", 
            color = NULL,
        ) +
        theme_bw() + 
        theme(legend.position = "bottom", strip.background = element_blank()) + 
        reqs
}

.binconf.ratio.plot <- function(DT, x, n, .ylab, xaes=ROUND_LAB)
{
    x <- enexpr(x)
    n <- enexpr(n)

    dplot <- copy(DT)
    dplot[, (c('P', 'CL', 'CU')) := binconf(x=eval(x), n=eval(n), return.df=TRUE) ]

    if('ROUND_LAB' %in% names(DT))
        dplot[, ROUND_LAB := gsub('^R0','Round', ROUND_LAB)]

    by_col <- fifelse('PERIOD' %in% names(DT), yes='PERIOD', no='ROUND' )
    dplot[, (c('LOG_RATIO_P', 'LOG_RATIO_CL', 'LOG_RATIO_CU')):= {
        P.both <- P[which(SEX == 'Both')]
        list( 
            LOG_RATIO_P = log(P/P.both),
            LOG_RATIO_CL = log(CL/P.both),
            LOG_RATIO_CU = log(CU/P.both)
        )
    }, by=by_col]

    dplot |> subset(SEX_LAB != 'Both') |> 
        ggplot(aes(x=!!xaes, color=SEX_LAB, pch=AGEGP, linetype=AGEGP, y=LOG_RATIO_P )) + 
        geom_hline(yintercept = 0, linetype='dashed', color='grey50') +
        geom_point(position=position_dodge(width=.8) ) + 
        geom_linerange(aes(ymin=LOG_RATIO_CL, ymax=LOG_RATIO_CU), position=position_dodge(width =.8)) +
        # scale_y_continuous(expand=c(0,0),labels=scales::percent) +
        scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
        # coord_cartesian(ylim = ylims) +
        labs( 
            x = NULL,
            y = .ylab,
            linetype = "Age", 
            pch = "Age", 
            color = NULL,
        ) +
        theme_bw() + 
        theme(legend.position = "bottom", strip.background = element_blank()) + 
        rotate_x_axis(60) +
        reqs
}

plot.hist.numerators.denominators <- function(DT, DRANGE)
{
    # can show agegroups as different bar fills! How to do this? 
    dplot <- copy(DT)
    dplot[, variable_lab := fifelse(variable == 'N_EVERSEQ', 'Ever-sequenced', 'Never deep-sequenced') ]
    lims <- dplot[ , levels(interaction(variable_lab, SEX_LAB))[c(1,3)]] 

    ggplot(dplot, aes(x=ROUND_LAB, pch=AGEGP, y=value, fill=interaction(variable_lab, SEX_LAB) )) +
        geom_col(data=dplot[variable == 'INFECTED_NON_SUPPRESSED_M'], position=position_dodge(width=.9), width=.9, color='black') + 
        geom_linerange(data=DRANGE, position = position_dodge(width =.9),
            aes(y=NULL, ymin=INFECTED_NON_SUPPRESSED_CL, ymax=INFECTED_NON_SUPPRESSED_CU, fill=NULL)) + 
        geom_col_pattern(
            data=dplot[variable == 'N_EVERSEQ'],
            aes(pattern=AGEGP),
            position=position_dodge(width=.9), color='black', width=.9,
            pattern_fill = 'white', pattern_color='white', pattern_density = .2, pattern_spacing=.01) +
        facet_grid(SEX_LAB~.) + 
        theme_bw() + 
        scale_y_continuous(expand = expansion(c(0,0.1)) ) +
        scale_pattern_manual(values = c('15-24' = 'circle' , '25-34'= 'none', '35-49' = 'stripe' )) + 
        guides(pattern = guide_legend(override.aes = list(fill = "#BDC5D0"), pattern_density = .1)) +
        guides(fill = guide_legend(override.aes = list(pattern = "none"))) +
        scale_fill_manual(
            values=c("#85D4E3", "#F4B5BD",  'white', 'white'), 
            labels=c('Men', 'Women', NA_character_, NA_character_),
            limits = lims,
            na.value = 'white'
        ) +  
        theme(legend.position = "bottom", strip.background = element_blank()) + 
        labs(
            x=NULL,
            y="Estimated number of individuals with\nunsuppressed HIV in the population",
            fill=NULL,
            pattern='Age',
            color=NULL
        )
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

# get unsuppressed
supptests <- load.viralsuppression.test.results()

#
# Quest
#

quest <- fread(file.path.quest)

# keep variable of interest, merge communities
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, arvmed, cuarvmed)]
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')
colnames(rinc) <- toupper(colnames(rinc))
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
hivs <- hivs[AGEYRS %between% c(15, 49)]


## calculate denominators (number infected and unsuppressed) ----


# load participation (% of census eligible population)
participation <- fread(file.participation)

# extract posterior samples for prevalence of suppression
ps <- c(M=0.5, CL=0.025, CU=0.975)
by_cols <- c('SEX', 'COMM', 'ROUND')
samples <- readRDS(file.treatment.cascade) |> 
    subset(COMM == 'inland', select=c(by_cols, 'AGEYRS','iterations', 'PROP_SUPPRESSED_POSTERIOR_SAMPLE'))
samples[, ROUND := gsub('R0', '', ROUND)]
samples <- merge(
    samples,
    n_partsuscinf[, .SD, .SDcols=c(by_cols, 'AGEYRS','INFECTED')],
    by=c(by_cols, 'AGEYRS'))
samples[, N_UNSUPPRESSED := INFECTED * (1-PROP_SUPPRESSED_POSTERIOR_SAMPLE) ]
samples[, AGEGP := ageyrs2agegp(AGEYRS)]
samples <- samples[, .(
    N_UNSUPPRESSED = sum(N_UNSUPPRESSED)),
    by=c(by_cols, 'AGEGP','iterations')]
samples <- samples[, {
    z <- as.list(quantile(N_UNSUPPRESSED, probs = ps))
    names(z) <- paste0('INFECTED_NON_SUPPRESSED_',names(ps))
    z
}, by=c(by_cols, 'AGEGP') ]
samples[, ROUND := paste0('R0', ROUND)]


## get numerator (ever deep sequenced) ----

# load seq count
sequ <- as.data.table(readRDS(file.characteristics_sequenced_ind_R14_18)) |> 
    subset(AGEYRS %between% c(15,49))
sequ[, round:= as.integer(gsub('R0|S','',ROUND))]

last_sequenced <- sequ[, list(N=length(ROUND),MAXROUND=max(round)),by=c('STUDY_ID')]
sequ <- merge(sequ,last_sequenced,by=c('STUDY_ID'))

# merge to meta_data
semt <- merge(
    hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV)],
    unique(sequ[, .(STUDY_ID,MAXROUND)]),
    by = 'STUDY_ID')
stopifnot(semt[, uniqueN(STUDY_ID)] == sequ[, uniqueN(STUDY_ID)])

# keep only positive participants and rounds 10 -> 18
semt <- semt[HIV == 'P']
semt[, round:= as.integer(gsub('R0|S','',ROUND))]
semt <- semt[round %between% c(10, 18)]

if(sequenced.at.round.or.later){
    # only count participants who were sequenced at current or future round
    semt <- subset(semt,round<=MAXROUND)
}

# exclude those that are unsuppressed at a given round 
keys <- c('STUDY_ID', 'ROUND')
suppressed <- supptests |> 
    subset(VLNS == 0) |>
    setnames('round', 'ROUND') |> 
    subset(select=keys) |> 
    setkeyv(keys)
suppressed[, SUP := TRUE ]
setkeyv(semt, keys)
semt <- merge(semt, suppressed, by=keys, all.x=TRUE)
semt[is.na(SUP), SUP := FALSE ]

# make age groups
semt[, AGEGP:= ageyrs2agegp(AGEYRS)]

# get subset who were sequenced
seqs <- semt[COMM=='inland', list(
    N_EVERSEQ_OLD = .N,
    N_EVERSEQ= sum(!SUP)
), by=c('COMM','ROUND','SEX','AGEGP')]

# merge together
by_cols <- c('COMM','ROUND','SEX','AGEGP')
tab_seq_unsup <- merge( samples, seqs, by=by_cols) 


######################
cat('\n Saving... \n')
######################

new_cols <- paste0('PROP_UNSUP_EVERDEEPSEQ', c('', '_CL', '_CU'))
tosave <- tab_seq_unsup |> 
    subset(select=c(by_cols, "N_EVERSEQ", 'INFECTED_NON_SUPPRESSED_M'))
tosave[, (new_cols) := binconf(x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED_M, return.df =TRUE) ]

file.name <- path.prop.unsuppressed.deepsequenced
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
}else{
  cat("File:", file.name, "already exists...\n")
}

################################
cat('\nDot-linearange plots \n')
################################


ylab1 <- fifelse(sequenced.at.round.or.later == TRUE,
    yes = "Proportion of individuals with unsuppressed HIV in the population\nwho were eventually deep-sequenced",
    no = "Proportion of individuals with unsuppressed HIV in the population\nwho were ever deep-sequenced",
) 

p_everseq_givenunspp <- .make.plot.with.binconf(
    tab_seq_unsup, 
    x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED_M, 
    .ylab = ylab1)


############################
cat('\nFilled histogram \n')
############################

# hard to show uncertainty pop sizes as we previous uncertainty is specified in one-year bands 
dplot <- melt(tab_seq_unsup,
    id.vars=c('ROUND_LAB', 'SEX_LAB', 'AGEGP') ,
    measure.vars = c('N_EVERSEQ', 'INFECTED_NON_SUPPRESSED_M' )) |> suppressWarnings()
p_hist <- plot.hist.numerators.denominators(dplot, tab_seq_unsup)

#####################
cat('\nRatio plot\n')
#####################
# take aggregate stratifying over round only
by_cols <- c('ROUND')
daggregates <- tab_seq_unsup[, list(
    COMM=unique(COMM),
    ROUND_LAB = unique(ROUND_LAB),
    SEX="Both",
    SEX_LAB="Both",
    AGEGP='All',
    INFECTED_NON_SUPPRESSED_M = sum(INFECTED_NON_SUPPRESSED_M), 
    INFECTED_NON_SUPPRESSED_CL = sum(INFECTED_NON_SUPPRESSED_CU), 
    INFECTED_NON_SUPPRESSED_CU = sum(INFECTED_NON_SUPPRESSED_CL), 
    N_EVERSEQ = sum(N_EVERSEQ)
) , by='ROUND']


p_ratio <- rbind(daggregates, tab_seq_unsup[, -c("N_EVERSEQ_OLD")]) |> 
    .binconf.ratio.plot(
        x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED_M, 
        .ylab='Log ratio of sampling probabilities in each population strata\nrelative to entire population'
    )

######################
cat('\nFINAL PLOTS\n')
######################

filenames <- paste0('unsupeverdeepseq_', c('hist', 'range', 'ratio' ), '.png')
filenames <- file.path(outdir, filenames)
ggsave_nature(filename = filenames[1], p=p_hist, w = 13, h = 11 )
ggsave_nature(filename = filenames[2], p=p_everseq_givenunspp, w = 13, h = 11 )
ggsave_nature(filename = filenames[3], p=p_ratio, w = 13, h = 11 )

#####################
cat('\nxi_j plots\n')
#####################

cat('\n  * get pairs * \n')

pairs <- read_pairs(file.pairs) |> 
    select.pairs.for.analysis(
    only.one.community='inland', 
    use_30com_pairs=FALSE, 
    only.transmission.after.start.observational.period=TRUE, 
    only.transmission.before.stop.observational.period=TRUE, 
    remove.pairs.from.rounds=NULL
)

# keep only pairs with source-recipient with a time of infection
pairs <- pairs[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
pairs[, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT := DATE_INFECTION.RECIPIENT < start_second_period_inland]



cat('\n  * get incidence * \n')

# prepare round data.tables
with(df_round_inland, {
    start_first_period_inland   <<- min_sample_date[round == 'R010']    # [1] "2003-09-26"
    stop_first_period_inland    <<- max_sample_date[round == 'R015']    # [1] "2013-07-05"
    start_second_period_inland  <<- min_sample_date[round == 'R016']    # [1] "2013-07-08"
    stop_second_period_inland   <<- max_sample_date[round == 'R018']    # [1] "2018-05-22"
})
stopifnot(start_first_period_inland < stop_first_period_inland)
stopifnot(stop_first_period_inland < start_second_period_inland)
stopifnot(start_second_period_inland < stop_second_period_inland)
df_round <- make.df.round(df_round_inland)
df_period <- make.df.period(
    start_first_period_inland, 
    stop_first_period_inland, 
    start_second_period_inland, 
    stop_second_period_inland, 
    df_round)


# get incidence
treatment_cascade <- read_treatment_cascade(
    file.treatment.cascade.prop.participants, 
    file.treatment.cascade.prop.nonparticipants)
eligible_count_round <- add_susceptible_infected(
    eligible_count_smooth = fread(file.eligible.count),
    proportion_prevalence = fread(file.prevalence.prop),
    participation = fread(file.participation),
    nonparticipants.male.relative.infection  =  1, 
    nonparticipants.female.relative.infection  =  1
) |> add_infected_unsuppressed(
    treatment_cascade = treatment_cascade, 
    participation = fread(file.participation), 
    nonparticipants.treated.like.participants = FALSE, 
    nonparticipants.not.treated = FALSE
)

incidence_cases <- get_incidence_cases_round(
    incidence.inland = fread(file.incidence.inland),
    eligible_count_round ) |> 
    summarise_incidence_cases_period(df_period)

# finally get proportions
proportion_sampling <- get_proportion_sampling(pairs, incidence_cases, 'a_nonexisting_directory_name')

# make plots

detectionprob <- subset(proportion_sampling, 
    select=c('COMM','SEX.RECIPIENT', 'AGEYRS.RECIPIENT', 'PERIOD', 'count', 'INCIDENT_CASES')
) |> unique()
names(detectionprob) <- gsub('.RECIPIENT','', names(detectionprob)) 
detectionprob[, AGEGP := ageyrs2agegp(AGEYRS)]
by_cols <- c('COMM', 'SEX', 'PERIOD', 'AGEGP')
detectionprob <- detectionprob[, list(count=sum(count), INCIDENT_CASES=sum(INCIDENT_CASES)), by=by_cols]

# a
dplot <- detectionprob[,  {
    z <- binconf(x=count, n=INCIDENT_CASES, return.df = TRUE)
    names(z) <- c('M', 'CL', 'CU')
    z
}, by=by_cols ]
dplot <- prettify_labs(dplot)

p2b <- ggplot(dplot, aes( 
    x=PERIOD, y=M, ymin=CL, ymax=CU, color=SEX_LAB, linetype=AGEGP, pch=AGEGP) ) + 
    geom_linerange(position=position_dodge(width=.7) ) +
    geom_point(position=position_dodge(width=.7) ) +
    scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
    labs( 
        x = NULL,
        y = "Detection probability",
        linetype = "Age", 
        pch = "Age", 
        color = NULL,
    ) +
    scale_y_continuous(limits=c(0, .3), expand=c(0,0), labels=scales::percent) +
    theme_bw() + 
    theme(legend.position = "bottom", strip.background = element_blank()) + 
    NULL 
p2b


dplot <- melt( detectionprob, id.vars = by_cols, measure.vars = c('count', 'INCIDENT_CASES')) |> suppressWarnings()
dplot <- prettify_labs(dplot)
fill_lims <- with(dplot,levels(interaction(variable, SEX_LAB)))[c(1,3)]

p_hist2 <- ggplot(dplot, aes(x=PERIOD, pch=AGEGP, y=value, fill=interaction(variable, SEX_LAB) )) +
    geom_col(data=dplot[variable == 'INCIDENT_CASES'], position=position_dodge(width=.9), width=.9, color='black') + 
    # geom_linerange(data=DRANGE, position = position_dodge(width =.9),
    #     aes(y=NULL, ymin=INFECTED_NON_SUPPRESSED_CL, ymax=INFECTED_NON_SUPPRESSED_CU, fill=NULL)) + 
    geom_col_pattern(
        data=dplot[variable == 'count'],
        aes(pattern=AGEGP),
        position=position_dodge(width=.9), color='black', width=.9,
        pattern_fill = 'white', pattern_color='white', pattern_density = .3, pattern_spacing=.01) +
    facet_grid(SEX_LAB~.) + 
    theme_bw() + 
    scale_y_continuous(expand = expansion(c(0,0.1)) ) +
    scale_pattern_manual(values = c('15-24' = 'circle' , '25-34'= 'none', '35-49' = 'stripe' )) + 
    guides(pattern = guide_legend(override.aes = list(fill = "#BDC5D0"), pattern_density = .1)) +
    guides(fill = guide_legend(override.aes = list(pattern = "none"))) +
    scale_fill_manual(
        values=c("#85D4E3", "#F4B5BD",  'white', 'white'), 
        labels=c('Men', 'Women', NA_character_, NA_character_),
        limits = fill_lims,
        na.value = 'white'
    ) +  
    theme(legend.position = "bottom", strip.background = element_blank()) + 
    labs(
        x=NULL,
        y="Detection probability among incident cases",
        fill=NULL,
        pattern='Age',
        NULL
    )
p_hist2

daggregates2 <- tmp[, list(
    COMM=unique(COMM),
    SEX="Both",
    AGEGP='All',
    count = sum(count),
    INCIDENT_CASES = sum(INCIDENT_CASES)
) , by='PERIOD']

p_ratio2 <- rbind(daggregates2, tmp) |> 
    prettify_labs() |>
    .binconf.ratio.plot(
        x=count, n=INCIDENT_CASES, 
        .ylab='Log ratio of detection probabilities in each population strata\nrelative to estimated incident cases',
        xaes = expr(PERIOD)
    )


filename <- '~/Downloads/tmp_xii_plot.png'
p <- ggarrange(p_hist2 + reqs, p2b + reqs ,p_ratio2 + reqs , ncol=1)
ggsave_nature(filename = filename, p=p, w = 13, h = 30 )



# ALL TOGETHER: 

# arrange b and c because they have same legend:
p2_bc <- ggarrange(p2b + reqs, p_ratio2 + reqs, labels = c('b', 'c'), ncol=1, common.legend = TRUE, legend='bottom' )
p1_bc <- ggarrange(p_everseq_givenunspp + reqs, p_ratio + reqs, labels = c('e','f'), ncol=1, common.legend = TRUE, legend='bottom' )

p_all <- ggarrange(

    p_hist2 + reqs,
    p_hist + reqs, 

    p2_bc, 
    p1_bc,

    nrow=  2, ncol = 2 , heights = c(1,2), labels = c('a', '', 'c', '')
)

filename <- '~/Downloads/tmp_all.png'
ggsave_nature(filename = filename, p=p_all, w = 20, h = 28 )
