{
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(scales)
    library(lubridate)
    library(ggpattern)
    library(patchwork)
    library(here)
    # single imports
    binconf <- Hmisc::binconf
} |> suppressPackageStartupMessages()


usr <- Sys.info()[['user']]

if(usr == 'ab1820'){
    gitdir <- '/rds/general/user/ab1820/home/git/phyloflows'
}else{
    gitdir <- here::here()
}



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
  file.characteristics_sequenced_ind_R14_18))  |> 
    all() |> 
    stopifnot()

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
p_hist <- plot.hist.numerators.denominators(dplot, tab_seq_unsup, filltext="Ever deep sequenced")

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
        .yrange = c(-1,1),
        x=N_EVERSEQ, n=INFECTED_NON_SUPPRESSED_M, 
        .ylab='Log ratio of sampling probabilities in each strata\nrelative to entire population'
    )

cat('\nxi_j plots\n')
# ___________________

# load detection and subset out susceptibles

file.name <- file.detection.probability.round.custom.agegroups
if(! file.exists(file.name) | config$overwrite.existing.files )
{

    detectionprob <- readRDS(file.detection.probability.round)
    names(detectionprob) <- gsub('.RECIPIENT', '' ,names(detectionprob)) 
    detectionprob[, `:=` (AGEGP = ageyrs2agegp(AGEYRS))]
    dsusc <- subset(detectionprob, select=c('ROUND','SEX','COMM','AGEYRS','AGEGP','SUSCEPTIBLE'))
    dcount <- subset(detectionprob, select=c('ROUND','SEX','COMM','AGEYRS','AGEGP','count'))

    # find uncertainty around incident cases
    cat('\n\n (loading samples...) \n\n')

    dincsamples <- load_incidence_rates_samples(file.incidence.samples.inland)
    dincsamples <-  merge(dincsamples, dsusc, by=c("ROUND", 'SEX', 'AGEYRS', 'COMM'))
    dincsamples[, `:=` (YEAR_LENGTH = .year.diff(MAX_SAMPLE_DATE, MIN_SAMPLE_DATE))]
    tmp <- dincsamples[, list(  DRAW = sum(INCIDENCE.DRAW * SUSCEPTIBLE *  YEAR_LENGTH)), by=c("ROUND", "SEX", "AGEGP", 'iterations' ) ]
    tmp <- tmp[, list( 
      Q = quantile(DRAW, probs = c(.025, .5, .975)),
      Q_LEVEL = c("CL", "M", "CU")
      ), by=c("ROUND", "SEX", "AGEGP" ) ] |> 
      dcast(ROUND + SEX + AGEGP ~ Q_LEVEL, value.var='Q')
    setnames(tmp, c("M", "CL", "CU"), paste0('INCIDENT_CASES', c("", "_LB", "_UB")) )

    # rename detectionprob with correct aggregation
    detectionprob <- merge(dcount, tmp, by =c("ROUND", "SEX", "AGEGP" ) )
    prettify_labs(detectionprob)

    cat("Saving file:", file.name, '\n')
    saveRDS(object=detectionprob,file=file.name)
}else{

    cat("Loading previously obtained file:", file.name, '\n')
    detectionprob <- readRDS(file.name)
}

detectionprob <- detectionprob[, .(
    INCIDENT_CASES = unique(INCIDENT_CASES),
    INCIDENT_CASES_LB = unique(INCIDENT_CASES_LB),
    INCIDENT_CASES_UB = unique(INCIDENT_CASES_UB),
    count = sum(count)
), by=c('ROUND', 'SEX', "AGEGP", "COMM", "ROUND_LAB", "SEX_LAB" )]


cat('\nDot-linearange plots \n')
#________________________________

p2b <- .make.plot.with.binconf(detectionprob, xvar=ROUND_LAB,
    x=count, n=INCIDENT_CASES,
    .ylab="Proportion of estimated infection events\nappearing as recipient")

cat('\nFilled histogram \n')
# __________________________

p_hist2 <- plot.hist.numerators.denominators.2(detectionprob, range=TRUE, filltext = 'hello')

cat('\nRatio plot\n')
# ___________________

daggregates2 <- detectionprob[, list(
    COMM=unique(COMM),
    SEX="Both",
    AGEGP='All',
    count = sum(count),
    INCIDENT_CASES = sum(INCIDENT_CASES)
) , by='ROUND']

tmp <- subset(detectionprob, select=names(daggregates2))
p_ratio2 <- rbind(daggregates2, tmp) |> 
    prettify_labs() |>
    .binconf.ratio.plot(
        x=count, n=INCIDENT_CASES, 
        .ylab='Log ratio of detection probabilities in each strata\nrelative to estimated infection events',
        xaes = expr(ROUND_LAB),
        .yrange = c(-5, 5)
    )


##########################
cat('\n ALL TOGETHER: \n')
##########################

r2 <- reqs + theme(
    legend.position='none',
    plot.tag = element_text(face='bold'),
    strip.background = element_blank(),
    panel.grid.minor= element_blank(),
    panel.grid.major= element_blank()
)
t <- function(lab){ labs(tag=lab, caption=lab) + r2 }
naturemed_reqs()
patchwork.way <- TRUE
if(! patchwork.way ) # the ggarrange way
{
    top_row_2 <- ggarrange(p_hist2 + labs( tag='a'), p_hist + labs( tag='d'), common.legend = TRUE, legend = "bottom")
    bottow_rows_2 <- ggarrange( ncol =2 , nrow=2,
        p2b + labs( tag='b') , p_everseq_givenunspp + labs( tag='e')  ,
        p_ratio2 + labs( tag='c') , p_ratio + labs( tag='f') ,
        common.legend = TRUE, legend = "bottom")
    all_rows_2 <- ggarrange( ncol =1 , nrow=2, heights = c(1,2),
        top_row_2, bottow_rows_2,
        common.legend = TRUE, legend = "bottom")
    all_rows <- copy(all_rows_2)

}else{

    # top_row_3 <-(p_hist2 + r2 + labs( tag='a')) + (p_hist + r2 + labs( tag='d')) + plot_layout(guides = 'collect') & theme(legend.position ='bottom', legend.box.margin = margin(r=-15, unit='pt'))
    # bottow_rows_3 <- (p2b + r2 + labs( tag='b')) + (p_everseq_givenunspp + r2 + labs( tag='e')) + (p_ratio2 + r2 + labs( tag='c')) + (p_ratio + r2 + labs( tag='f')) + plot_layout(guides = 'collect') & theme(legend.position ='bottom')
    # all_rows_3 <- top_row_3 / bottow_rows_3 + plot_layout(heights = c(.94,2))
    # all_rows_3 <- all_rows_3 & reqs 
    # all_rows <- copy(all_rows_3)

    all_rows_3b <- (p_hist2 + labs(tag = 'a') + r2 | p_hist + labs(tag = 'd') + r2) / (p2b + labs(tag = 'b') + r2 | p_everseq_givenunspp + labs(tag = 'e') + r2) / (p_ratio2 + labs(tag = 'c') + r2 | p_ratio + labs(tag = 'f') + r2) + plot_layout(guides='collect') & theme(legend.position = 'bottom') + reqs
    all_rows <- copy(all_rows_3b)

}

filename <- file.path(outdir, 'extendeddatafigure_samplingproportions.pdf')
ggsave_nature(filename = filename, p=all_rows, w = 18, h = 24 )
filename <- file.path(outdir, 'extendeddatafigure_samplingproportions.png')
ggsave_nature(filename = filename, p=all_rows, w = 18, h = 24 )
