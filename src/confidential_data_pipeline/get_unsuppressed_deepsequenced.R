{
    library(data.table)
    library(dplyr)
    library(ggplot2)
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

file.name <- file.prop.unsuppressed.deepsequenced
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

detectionprob <- readRDS(file.detection.probability.round)
names(detectionprob) <- gsub('.RECIPIENT', '' ,names(detectionprob)) 
detectionprob <- subset(detectionprob, select=c('SEX', 'COMM', 'AGEYRS', 'ROUND','INCIDENT_CASES', 'count', 'prop_sampling'))

detectionprob[, `:=` (AGEGP = ageyrs2agegp(AGEYRS))]
detectionprob <- detectionprob[, 
    lapply(.SD, sum),
    .SDcols = c('INCIDENT_CASES', 'count'),
    by=c('ROUND','SEX', 'COMM', 'AGEGP')]
prettify_labs(detectionprob)


cat('\nDot-linearange plots \n')
#________________________________

p2b <- .make.plot.with.binconf(detectionprob, xvar=ROUND_LAB,
    x=count, n=INCIDENT_CASES,
    .ylab="Proportion of estimated incident cases\nappearing as recipient")


cat('\nFilled histogram \n')
# __________________________

by_cols <- c('ROUND', 'SEX', 'COMM', 'AGEGP')
dplot <- melt( detectionprob, id.vars = by_cols, measure.vars = c('count', 'INCIDENT_CASES')) |> suppressWarnings()
dplot <- prettify_labs(dplot)
fill_lims <- with(dplot,levels(interaction(variable, SEX_LAB)))[c(1,3)]

p_hist2 <- ggplot(dplot, aes(x=ROUND_LAB, pch=AGEGP, y=value, fill=interaction(variable, SEX_LAB) )) +
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
        y="Estimated number of incident cases\nappearing as recipient",
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
) , by='ROUND']

p_ratio2 <- rbind(daggregates2, detectionprob[, -c("ROUND_LAB", "SEX_LAB")]) |> 
    prettify_labs() |>
    .binconf.ratio.plot(
        x=count, n=INCIDENT_CASES, 
        .ylab='Log ratio of detection probabilities in each strata\nrelative to estimated incident cases',
        xaes = expr(ROUND_LAB),
        .yrange = c(-5, 5)
    )


##########################
cat('\n ALL TOGETHER: \n')
##########################

# make the legend
cowplot::get_legend(
    ggplot(dplot, aes(x=value, y=value,  fill=SEX_LAB,  pch=AGEGP, linetype=AGEGP)) + 
        geom_point(color='black') +
        geom_ribbon(aes(ymin=value, ymax=value, alpha=AGEGP)) +
        geom_linerange(aes(ymin=value, ymax=value)) +
        scale_color_manual(values=c("black", "black" )) + 
        scale_fill_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
        scale_alpha_manual(values=c(`15-24`=.3, `25-34`=.6, `35-49`=.8)) +
        guides(
            # linetype = guide_legend(override.aes = list(color='black')),
            # pch = guide_legend(override.aes = list(color='white')),
        ) +
        labs(fill = 'Gender', linetype=NULL, pch=NULL, group=NULL, alpha=NULL, color=NULL) +
        theme(legend.position='bottom')  +
    reqs
) -> legend



library(patchwork)
r2 <- reqs + 
    theme(legend.position='none', plot.tag = element_text(face='bold'))

t <- function(lab){ labs(tag=lab) + r2 }


p_all <- ggarrange( ncol=2, nrow=3,
    p_hist2 + r2, p_hist + r2,
    p2b +r2 , p_everseq_givenunspp + r2  ,
    p_ratio2 + r2 , p_ratio + r2 ,
    labels = 'auto', font.label = list(size=8) 
    )


p_all2 <- ggarrange( p_all, legend, ncol=1, heights = c(40,1))

filename <- '~/Downloads/tmp_all.png'
ggsave_nature(filename = filename, p=p_all2, w = 18, h = 24 )
