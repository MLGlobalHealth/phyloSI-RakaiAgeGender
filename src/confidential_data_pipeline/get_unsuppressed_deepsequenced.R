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


# TODO: finalise plots:
# - 1 : histograms consistent with Alex's plot (unsuppressed among participants) 
# - 2 : redo point-range plot 
# - 3 : relative ratios plots: are there differences among age groups

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

    if (! 'ROUND_LAB' %in% nms)
        DT[, ROUND_LAB := gsub("R0", "Round ", ROUND)]

    if (! 'SEX_LAB' %in% nms)
        DT[, SEX_LAB := fifelse(SEX=="F", "Women", "Men")] 

    return(DT)
}

# plot.eligible.count <- function(DT=eligible_count_round, var){
#     var <- enexpr(var)
#     dplot <- copy(DT)
#     dplot[, AGEGP := ageyrs2agegp(AGEYRS)]
#     dplot[, AGEYRS := NULL ]
#     dplot <- dplot[, lapply(.SD, sum), by=c('COMM', 'ROUND','SEX', 'AGEGP')]
#     dplot |> prettify_labs()
#     ggplot(dplot, aes( x=ROUND_LAB, y=!!var, color=SEX_LAB, fill=SEX_LAB, pch=AGEGP))  +
#         # geom_point(data=dplot[SEX %like% "M"], position=position_dodge(.8)) +
#         # geom_point(data=dplot[SEX %like% "F"], position=position_dodge(.8)) +
#         geom_line(aes(group = interaction(SEX, AGEGP), linetype=AGEGP )) +
#         theme_bw() + 
#         theme(legend.position = "bottom", strip.background = element_blank()) + 
#         labs(x=NULL) +
#         reqs
# }


.make.plot.with.binconf <- function(DT, x, n, .ylab, ylims = c(0,1) )
{
    x <- enexpr(x)
    n <- enexpr(n)

    dplot <- copy(DT)
    dplot[, (c('P', 'CL', 'CU')) := binconf(x=eval(x), n=eval(n), return.df=TRUE) ]
    dplot[, ROUND_LAB := gsub('^R0','Round', ROUND_LAB)]

    add_hline = fifelse(max(ylims) < 1, yes=TRUE, no=FALSE)

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

.binconf.ratio.plot <- function(DT, x, n, .ylab)
{
    x <- enexpr(x)
    n <- enexpr(n)

    dplot <- copy(DT)
    dplot[, (c('P', 'CL', 'CU')) := binconf(x=eval(x), n=eval(n), return.df=TRUE) ]
    dplot[, ROUND_LAB := gsub('^R0','Round', ROUND_LAB)]

    add_hline = TRUE

    dplot[, (c('RATIO_P', 'RATIO_CL', 'RATIO_CU')):= {
        P.both <- P[which(SEX == 'Both')]
        list( 
            RATIO_P = P/P.both ,
            RATIO_CL = CL/P.both,
            RATIO_CU = CU/P.both
        )
    }, by=ROUND]

    dplot |> subset(SEX_LAB != 'Both') |> 
        ggplot(aes(x=ROUND_LAB, color=SEX_LAB, pch=AGEGP, linetype=AGEGP, y=RATIO_P )) + 
        geom_hline(yintercept = 1, linetype='dashed', color='grey50') +
        geom_point(position=position_dodge(width=.8) ) + 
        geom_linerange(aes(ymin=RATIO_CL, ymax=RATIO_CU), position=position_dodge(width =.8)) +
        scale_y_continuous(expand=c(0,0),labels=scales::percent) +
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
        reqs
}

plot.hist.numerators.denominators <- function(DT)
{
    # can show agegroups as different bar fills! How to do this? 
    dplot <- copy(DT)
    dplot[, variable_lab := fifelse(variable == 'N_EVERSEQ', 'Ever-sequenced', 'Never deep-sequenced') ]
    lims <- dplot[ , levels(interaction(variable_lab, SEX_LAB))[c(1,3)]] 
    ggplot(dplot, aes(x=ROUND_LAB, pch=AGEGP, y=value, fill=interaction(variable_lab, SEX_LAB) )) +
        geom_col(data=dplot[variable == 'INFECTED_NON_SUPPRESSED'], position=position_dodge(width=.9), width=.9, color='black') + 
        geom_col_pattern(
            data=dplot[variable == 'N_EVERSEQ'],
            aes(pattern=AGEGP),
            position=position_dodge(width=.9), color='black', width=.9,
            pattern_fill = 'white', pattern_color='white', pattern_density = .2, pattern_spacing=.01) +
        facet_grid(SEX_LAB~.) + 
        theme_bw() + 
        scale_y_continuous(expand = expansion(c(0,0.1)) ) +
        scale_pattern_manual(values = c('15-24' = 'circle' , '25-34'= 'none', '35-49' = 'stripe' )) + 
        guides(pattern = guide_legend(override.aes = list(fill = "purple"), pattern_density = .1)) +
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
            y="Number of unsuppressed individuals among census eligible",
            fill="Ever deep-sequenced",
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

#
# Quest
#

quest <- fread(file.path.quest)

# get unsuppressed
supptests <- load.viralsuppression.test.results()

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
hivs <- hivs[AGEYRS %between% c(15, 49)]


## calculate denominators (number infected and unsuppressed) ----


# load participation (% of census eligible population)
participation <- fread(file.participation)

# load non-suppressed proportion 
treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants, file.treatment.cascade.prop.nonparticipants)

# make a second version with no unsuppressed out-of-study for ease in later calculations
nms <- grep('PROP_UNSUPPRESSED_NONPARTICIPANTS',names(treatment_cascade), value=TRUE)
treatment_cascade_2 <- copy(treatment_cascade)
treatment_cascade_2[, (nms) := lapply(.SD, function(x) rep(0, length(x))), .SDcols = nms]


# get counts among census eligible
tmp <- add_susceptible_infected(
    eligible_count_smooth = fread(file.eligible.count),
    proportion_prevalence = fread(file.prevalence.prop),
    participation = participation,
    nonparticipants.male.relative.infection = nonparticipants.male.relative.infection,
    nonparticipants.female.relative.infection = nonparticipants.female.relative.infection 
) 
eligible_count_round <- add_infected_unsuppressed(
    tmp,
    treatment_cascade,
    participation,
    nonparticipants.treated.like.participants, 
    nonparticipants.not.treated)

# repeat the same process to get inf not suppressed amoang participants only
cols <- setdiff(names(eligible_count_round), c('COMM', 'ROUND', 'SEX', 'AGEYRS'))
new_cols <- paste0(cols, '_2')
eligible_count_round_onlyparts <- add_infected_unsuppressed(
    tmp,
    treatment_cascade_2,
    participation,
    nonparticipants.treated.like.participants, 
    nonparticipants.not.treated) |> 
    setnames(cols, new_cols)


if(0)   # check my understanding of the function above is correct 
{
    by_cols <- c('ROUND','COMM','SEX', 'AGEYRS')
    part_tmp <- copy(participation) |> subset(COMM=='inland')
    part_tmp[, ROUND := paste0('R0', ROUND)]
    # part_tmp[, AGEGP := ageyrs2agegp(AGEYRS) ]

    treat_tmp <- treatment_cascade[COMM=='inland']
    # treat_tmp[, AGEGP := ageyrs2agegp(AGEYRS) ]

    contrib_participants_unsuppression <- merge(part_tmp, treat_tmp)[, {
        weights.parts <- PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS_M 
        weights.nonparts <- (1-PARTICIPATION) * PROP_UNSUPPRESSED_NONPARTICIPANTS_M 
        list(MULTIPLIER = weights.parts/(weights.parts + weights.nonparts))
    }, by=by_cols]

    tmp <- merge(eligible_count_round, contrib_participants_unsuppression, by=by_cols) |>
        merge(
            subset(eligible_count_round_onlyparts, select=c(by_cols, 'INFECTED_NON_SUPPRESSED_2')),
            by = by_cols
        )
    # so yes, these are almost the same...
    tmp[, INFECTED_NON_SUPPRESSED * (1-MULTIPLIER) - INFECTED_NON_SUPPRESSED_2]
}

# aggregate counts by age groups
eligible_count_round[ , AGEGP := ageyrs2agegp(AGEYRS) ] 
eligible_count_round <- eligible_count_round[, list(
  INFECTED_NON_SUPPRESSED = sum(INFECTED_NON_SUPPRESSED)),
  by=c('COMM','AGEGP','SEX','ROUND')]

eligible_count_round_onlyparts[ , AGEGP := ageyrs2agegp(AGEYRS) ] 
eligible_count_round_onlyparts <- eligible_count_round_onlyparts[, list(
  INFECTED_NON_SUPPRESSED_2 = sum(INFECTED_NON_SUPPRESSED_2)),
  by=c('COMM','AGEGP','SEX','ROUND')]


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
if(TRUE){
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
}

# make age groups
semt[, AGEGP:= ageyrs2agegp(AGEYRS)]

seqs <- semt[COMM=='inland', list(
    N_EVERSEQ_OLD = .N,
    N_EVERSEQ= sum(!SUP)
), by=c('COMM','ROUND','SEX','AGEGP')]
seqs

# merge together

tmp  <- merge( eligible_count_round, seqs,  by=c('COMM','ROUND','SEX','AGEGP')) |> prettify_labs()
tmp1 <- merge( eligible_count_round_onlyparts, seqs,  by=c('COMM','ROUND','SEX','AGEGP')) |> prettify_labs()

# save saveRDS(tmp,file=file.path(outdir,'prop_unsupp_eventuallydeepseq_byroundagesex.RDS'))

# plot

tab_seq_unsup <- copy(tmp)

################################
cat('\nDot-linearange plots \n')
################################


ylab1 <- fifelse(sequenced.at.round.or.later == TRUE,
    yes = "Proportion of unsuppressed census eligible\neventually deep-sequenced",
    no = "Proportion of unsuppressed census eligible\nwho were ever deep-sequenced"
) 
p_everseq_givenunspp <- .make.plot.with.binconf(
    tab_seq_unsup, 
    x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED, 
    .ylab = ylab1)

# filename <- file.path(outdir, "prop_unsupp_eventuallydeepseq_byroundagesex.pdf") 
# if(sequenced.at.round.or.later) 
#     filename <- gsub('_eventuallydeepseq_', '_eventuallydeepseq_nopast_', filename)
# ggsave_nature(filename=filename, p=p_everseq_givenunspp, w=12, h=10)

############################
cat('\nFilled histogram \n')
############################

# hard to show uncertainty pop sizes as we previous uncertainty is specified in one-year bands 
dplot <- melt(tab_seq_unsup,
    id.vars=c('ROUND_LAB', 'SEX_LAB', 'AGEGP') ,
    measure.vars = c('N_EVERSEQ', 'INFECTED_NON_SUPPRESSED')) |> suppressWarnings()
p_hist <- plot.hist.numerators.denominators(dplot)

# filename <- file.path(outdir, "hist_unsupp_eventuallydeepseq_byroundagesex.pdf") 
# filename <- file.path(outdir, "hist_unsupp_eventuallydeepseq_byroundagesex.png") 
# if(sequenced.at.round.or.later) 
#     filename <- gsub('_eventuallydeepseq_', '_eventuallydeepseq_nopast_', filename)
# ggsave_nature(filename=filename, p=p_hist, w=13, h=11)

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
    INFECTED_NON_SUPPRESSED = sum(INFECTED_NON_SUPPRESSED), 
    N_EVERSEQ = sum(N_EVERSEQ)
) , by='ROUND']

dall <- rbind(daggregates, tab_seq_unsup[, -"N_EVERSEQ_OLD"])

p_ratio <- dall |> .binconf.ratio.plot(
    x=N_EVERSEQ, 
    n=INFECTED_NON_SUPPRESSED, 
    .ylab='Ratio of sequencing success relative to round-aggregates'
    )

# filename <- file.path(outdir, "ratio_unsupp_eventuallydeepseq_byroundagesex.pdf") 
# filename <- file.path(outdir, "ratio_unsupp_eventuallydeepseq_byroundagesex.png") 
# if(sequenced.at.round.or.later) 
#     filename <- gsub('_eventuallydeepseq_', '_eventuallydeepseq_nopast_', filename)
# ggsave_nature(filename=filename, p=p_hist, w=13, h=11)

######################
cat('\nFINAL PLOTS\n')
######################

filenames <- paste0('unsupeverdeepseq_', c('hist', 'range', 'ratio' ), '.png')
filenames <- file.path(outdir, filenames)
ggsave_nature(filename = filenames[1], p=p_hist, w = 13, h = 11 )
ggsave_nature(filename = filenames[2], p=p_everseq_givenunspp, w = 13, h = 11 )
ggsave_nature(filename = filenames[3], p=p_ratio, w = 13, h = 11 )
