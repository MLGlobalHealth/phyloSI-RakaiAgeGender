{
    library(data.table)
    library(dplyr)
    library(ggplot2)
    # library(ggpattern)
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
source(file.path(gitdir.R.conf, 'helpers_sampling_proportions.R'))

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

 


## calculate denominators (number infected and unsuppressed) ----

# extract posterior dunsuppressed for prevalence of suppression
participation <- fread(file.participation)

# set DT=... with datatable containing n_infected to (slightly) speed up
dunsuppressed <- extract_posterior_suppression_samples()

quest <- fread(file.path.quest)

# get known unsuppressed in r > 15
supptests <- load.viralsuppression.test.results(QUEST=quest)



## calculate numerators: # of HIV positive at round eventually deep sequenced 
#   if no viral rebound, unsuppressed at sequencing time were unsuppressed test

hivs <- get_hivids()


## get numerator (ever deep sequenced) ----

# load seq count
sequ <- as.data.table(readRDS(file.characteristics_sequenced_ind_R14_18)) |> 
    subset(AGEYRS %between% c(15,49)) |> 
    within(round <- as.integer(gsub('R0|S', '', ROUND)) )
last_sequenced <- sequ[, list(N=length(ROUND),MAXROUND=max(round)),by=c('STUDY_ID')]
sequ <- merge(sequ,last_sequenced,by=c('STUDY_ID'))

# merge to hiv positives 
semt <- merge(
    hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV)],
    unique(sequ[, .(STUDY_ID,MAXROUND)]),
    by = 'STUDY_ID')
stopifnot(semt[, uniqueN(STUDY_ID)] == sequ[, uniqueN(STUDY_ID)])

# keep only positive participants and rounds 10 -> 18
semt <- semt[HIV == 'P']
semt[, round:= as.integer(gsub('R0|S','',ROUND))]
semt <- semt[round %between% c(10, 18)]

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
tab_seq_unsup <- merge( dunsuppressed, seqs, by=by_cols) 


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
    saveRDS(object=tosave,file=file.name)
}else{
    cat("File:", file.name, "already exists...\n")
}


#############################
cat('\n### MAKE PLOTS ###\n')
#############################

prettify_labs(tab_seq_unsup)

cat('\nDot-linearange plots \n')
#________________________________


ylab1 <- fifelse(sequenced.at.round.or.later == TRUE,
    yes = "Proportion of individuals with unsuppressed HIV\nin the population who were eventually deep-sequenced",
    no = "Proportion of individuals with unsuppressed HIV\nin the population who were ever deep-sequenced",
) 

p_everseq_givenunspp <- .make.plot.with.binconf(
    tab_seq_unsup, 
    x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED_M, 
    .ylab = ylab1)

cat('\nFilled histogram \n')
# __________________________

# hard to show uncertainty pop sizes as we previous uncertainty is specified in one-year bands 
dplot <- melt( tab_seq_unsup,
    id.vars=c('ROUND_LAB', 'SEX_LAB', 'AGEGP') ,
    measure.vars = c('N_EVERSEQ', 'INFECTED_NON_SUPPRESSED_M' )) |> suppressWarnings()
p_hist <- plot.hist.numerators.denominators(dplot, tab_seq_unsup)
p_hist

cat('\nRatio plot\n')
# ___________________

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
        .ylab='Log ratio of sampling probabilities in each strata\nrelative to entire population'
    )

cat('\nxi_j plots\n')
# ___________________

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

cat('\nDot-linearange plots \n')
#________________________________

p2b <- .make.plot.with.binconf(detectionprob, xvar=PERIOD,
    x=count, n=INCIDENT_CASES, ylims = c(0, .35),
    .ylab="Detection Probability")



cat('\nFilled histogram \n')
# __________________________

dplot <- melt( detectionprob, id.vars = by_cols, measure.vars = c('count', 'INCIDENT_CASES')) |> suppressWarnings()
dplot <- prettify_labs(dplot)
fill_lims <- with(dplot,levels(interaction(variable, SEX_LAB)))[c(1,3)]

p_hist2 <- ggplot(dplot, aes(x=PERIOD, pch=AGEGP, y=value, fill=interaction(variable, SEX_LAB) )) +
    geom_col(data=dplot[variable == 'INCIDENT_CASES'], position=position_dodge(width=.9), width=.9, color='black') + 
    geom_col(
        data=dplot[variable == 'count'],
        aes(alpha=AGEGP),
        position=position_dodge(width=.9), color='black', width=.9
        )+
    facet_grid(SEX_LAB~.) + 
    theme_bw() + 
    scale_y_continuous(expand = expansion(c(0,0.1)) ) +
    scale_alpha_manual(values = c('15-24' = .33, '25-34'= .66, '35-49' = 1 )) + 
    guides(alpha = guide_legend(override.aes = list(fill = "#BDC5D0"))) +
    guides(fill = guide_legend(override.aes = list(alpha= 1))) +
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
        alpha='Age',
        NULL
    ) +
    rotate_x_axis(30)

cat('\nRatio plot\n')
# ___________________

daggregates2 <- detectionprob[, list(
    COMM=unique(COMM),
    SEX="Both",
    AGEGP='All',
    count = sum(count),
    INCIDENT_CASES = sum(INCIDENT_CASES)
) , by='PERIOD']

p_ratio2 <- rbind(daggregates2, detectionprob) |> 
    prettify_labs() |>
    .binconf.ratio.plot(
        x=count, n=INCIDENT_CASES, 
        .ylab='Log ratio of detection probabilities in each strata\nrelative to estimated incident cases',
        xaes = expr(PERIOD)
    )


##########################
cat('\n ALL TOGETHER: \n')
##########################


library(patchwork)
t <- function(lab) labs(tag=lab)
r2 <- reqs + theme(legend.position='none', plot.tag = element_text(face='bold'))

# get top row
leg <- cowplot::get_legend(p_hist + reqs)
combined <- (p_hist2 + labs(tag= 'a')  + p_hist + labs(tag = 'b')) & r2
combined <- combined/leg + plot_layout(heights = c(20,1))

# get last two rows
leg2 <- cowplot::get_legend(p2b + reqs)
combined2 <- (p2b +t('c') + p_everseq_givenunspp + t('d') ) / (p_ratio2 + t('e') + p_ratio + t('f')) & r2
combined2

p_all2 <- ggarrange(combined, combined2, ncol=1, heights = c(1,2))

filename <- '~/Downloads/tmp_all.png'
ggsave_nature(filename = filename, p=p_all2, w = 20, h = 28 )
