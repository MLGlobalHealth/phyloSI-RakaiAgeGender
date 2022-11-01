library(data.table)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
library(grid)
library(ggtree)
library(ggnet)
require(lubridate)
library(rstan)
library(gridExtra)

# laptop
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
  outdir <- '~/Box\ Sync/2021/phyloflows/'

  jobname <- 'new_period'
  stan_model <- 'gp_221101a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/home/andrea'))
{
  indir <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata/'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'test'
  stan_model <- 'gp_220911a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/rds/general'))
{
  indir.deepsequence_analyses   <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-stan_model')
  stopifnot(args_line[[7]]=='-jobname')
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  stan_model <- args_line[[6]]
  jobname <- args_line[[8]]
}

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))
outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
if(!dir.exists(dirname(outfile.figures))) dir.create(dirname(outfile.figures))
if(!dir.exists(dirname(outdir.table))) dir.create(dirname(outdir.table))

# indicators
include.only.heterosexual.pairs <- T
include.pairs.uncleardirection.disparateviralloads.complextopology <- F
threshold.likely.connected.pairs <- 0.5
use.tsi.estimates <- F
use.network.derived.infection.dates <- T
use.tsi.oneyear.before.first.positive <- F
remove.inconsistent.infection.dates <- F
remove.young.individuals <- T
remove.missing.community.recipient <- T
remove.neuro.individuals <- T
only.transmission.after.start.observational.period <- T
only.transmission.before.stop.observational.period <- T
use.diagonal.prior <- F
use.informative.prior <- F
only.transmission.same.community <- F
only.participant.treated <- T
remove.pairs.from.rounds <- NULL
only.one.community <- 'inland'

# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
# file.path.tsiestimates <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', 'TSI_estimates_220119.csv')
file.path.tsiestimates <- file.path(dirname(indir.deepsequence_analyses), 'PANGEA2_RCCS_UVRI_TSI2','2022_07_26_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T', 'aggregated_adjusted_TSI_with_estimated_dates.csv')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')

# from EMODO_RAKAI repo
file.incidence.inland	<- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', "Rakai_incpredictions_inland_220930.csv")
file.incidence.fishing	<- file.path(indir.deepsequencedata, 'RCCS_R15_R18', "Rakai_incpredictions_fishing_220930.csv")

# from misc/
file.path.meta <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_20220329.RData')
file.path.round.timeline <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'RCCS_round_timeline_220905.RData')

file.eligible.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220830.csv')
file.participation <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'RCCS_participation_220915.csv')
file.unsuppressed.prop <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', "RCCS_artcoverage_estimates_220906.csv")
file.unsuppressed.share <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_artcoverage_share_sex_220906.csv'))
file.prevalence.prop <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_prevalence_estimates_220811.csv')
file.prevalence.share <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('RCCS_prevalence_share_sex_220830.csv'))
file.unsuppressed_rate_ratio <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_artcoverage_ratio_sex_220926.csv'))
file.reported.sexual.partnerships <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('cont_age-R015.rds'))

path.to.stan.model <- file.path(indir, 'stan_models', paste0(stan_model, '.stan'))

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'statistics_functions.R'))
source(file.path(indir, 'functions', 'stan_utils.R'))

# check args
stopifnot( use.tsi.estimates + use.network.derived.infection.dates + use.tsi.oneyear.before.first.positive == 1)

# load anonymous aid
aik <- .read(file.anonymisation.keys); aik$X <- NULL

# load chains
load(file.path.chains.data)
dchain <- as.data.table(dchain)

# load meta data and round timeline
load(file.path.meta)
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]

# load Tanya's estimate time since infection using phylogenetic data
if(use.tsi.estimates)
        time.since.infection <- make.time.since.infection(fread(file.path.tsiestimates))

# load census eligible ount
eligible_count_smooth <- fread(file.eligible.count)

# load participation (% of census eligible population)
participation <- fread(file.participation)

# load proportion prevalence
proportion_prevalence <- fread(file.prevalence.prop)

# load non-suppressed proportion 
proportion_unsuppressed <- fread(file.unsuppressed.prop)

# load incidence estimates from Adam
incidence.inland <- fread(file.incidence.inland)

#for plots
unsuppressed_rate_ratio <- fread(file.unsuppressed_rate_ratio) # sex ratio of unsuppression rate
df_reported_contact <- as.data.table(readRDS(file.reported.sexual.partnerships)) # reported sexual contacts
unsuppressed_share <- fread(file.unsuppressed.share) # share of unsuppressed count by sex
infected_share <- fread(file.prevalence.share) # share of infected count by sex


#
# Define start time, end time and cutoff
#

start_observational_period_inland <- df_round_inland[round == 'R010', min_sample_date] # "2003-09-26"
stop_observational_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

cutoff_date <- df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"

stopifnot(start_observational_period_inland <= cutoff_date & stop_observational_period_inland >= cutoff_date)

df_period <- make.df.period(start_observational_period_inland, stop_observational_period_inland, 
                            cutoff_date)

df_round <- make.df.round(df_round_inland, df_period)


#
# Find count eligible susceptible / infected / infected unsuppressed 
# 


# by round
eligible_count_round <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence)
eligible_count_round <- add_infected_unsuppressed(eligible_count_round, proportion_unsuppressed, participation, only.participant.treated)
eligible_count_round[, table(ROUND, COMM)]


#
# Find incidence cases
#

# by round
incidence_cases_round <- get_incidence_cases_round(incidence.inland, eligible_count_round)
incidence_cases_round[, table(ROUND, COMM)]

# summarise by time period
incidence_cases <- summarise_incidence_cases_period(incidence_cases_round, cutoff_date, df_period)
incidence_cases[, table(PERIOD, COMM)]


#
# Find phylo pairs 
#

# get time of infection (using Tanya's estimate if use.tsi.estimates == T)
meta_data <- find.time.of.infection(meta_data, time.since.infection, use.tsi.estimates)

# get likely transmission pairs
chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)

# merge meta data to source and recipient
pairs.all <- pairs.get.meta.data(chain, meta_data, aik)

if(use.tsi.estimates){
  plot.coherent.tsi.estimates.with.seroconversion(outfile.figures)
  plot.tsi.relationships.among.source.recipient.pairs(outfile.figures)
}

if(include.only.heterosexual.pairs)
{
  cat('Keep only heterosexual pairs\n')
  cat('Removing ', nrow(pairs.all[! ((SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M'))]), ' pairs\n')
  pairs.all <- pairs.all[(SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M')]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')

  # remove from chain too for later step.
  setkey(chain, SOURCE, RECIPIENT)
  chain <- chain[pairs.all[, .(SOURCE, RECIPIENT)]]
}

# this below changes the values of metadata, chain & pairs.all
if(use.network.derived.infection.dates)
{
        cat('Assign infection dates according to network relationships...\n')

        # helper functions here are in helpers.
        filename <- paste0('network_attributed_doi',
                           '_chains', .assign.code.meta(file.path.meta),
                           '_meta',  .assign.code.chain(file.path.chains.data),
                           '_onlyhetero', as.integer(include.only.heterosexual.pairs),
                           '.rds')
        if(filename %like% 'NA')
                stop('Updated path for chains/meta? Need to update .assign.code.meta/meta\n')

        filename <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', filename)

        out <- update.meta.pairs.after.doi.attribution(path=filename, outfile.figures, overwrite = F)
        stopifnot(nrow(meta_data) == nrow(out$meta_data))
        

        cat('\t Because of plausible dates of infection inconsistencies:\n')
        .f <- function(DT){DT[, .(SOURCE, RECIPIENT)]}
        tmp <- setdiff(.f(out$pairs.all), .f(pairs.all))[, .N]
        cat('\t - switched direction of transmission for', tmp,'couples.\n')
        tmp <- nrow(pairs.all) - nrow(out$pairs.all)
        cat('\t - Removed', tmp,' couples.\n\n')

        pairs.all <- out$pairs.all
        chain <- out$chain
        meta_data <- out$meta_data
}

if(remove.inconsistent.infection.dates){
  cat('Remove infections for which estimated date at infection of source is after the estimated date at infection of the recipient.\n ')
  cat('Removing ', nrow(pairs.all[ DATE_INFECTION.SOURCE >= DATE_INFECTION.RECIPIENT ]), ' pairs\n')
  pairs.all <- pairs.all[! DATE_INFECTION.SOURCE >= DATE_INFECTION.RECIPIENT ]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(remove.young.individuals){
  # exclude young indivis
  cat('\nExcluding sources and recipients younger than 15\n')
  cat('Removing ', nrow(pairs.all[AGE_TRANSMISSION.SOURCE < 15 | AGE_INFECTION.RECIPIENT < 15]), ' pairs\n')
  pairs.all <- pairs.all[AGE_TRANSMISSION.SOURCE >= 15 & AGE_INFECTION.RECIPIENT >= 15]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}

# plot time of infection
plot_hist_time_infection(copy(pairs.all), cutoff_date, outfile.figures)

if(!is.null(only.one.community)){
  cat('\nExcluding sources and recipients in ',   pairs.all[COMM.RECIPIENT != only.one.community, unique(COMM.RECIPIENT)] ,'\n')
  cat('Removing ', nrow(pairs.all[COMM.RECIPIENT != only.one.community]), ' pairs\n')
  pairs.all <- pairs.all[COMM.RECIPIENT == only.one.community]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(only.transmission.after.start.observational.period){
  cat('\nFor inland excluding recipients infected before ', as.character(start_observational_period_inland), '\n')
  cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland')]

  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}

if(only.transmission.before.stop.observational.period){
  cat('\nFor inland excluding recipients infected after ', as.character(stop_observational_period_inland), '\n')
  cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland')]
  
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(remove.missing.community.recipient){
  cat('\nExcluding recipients without community \n')
  cat('Removing ', nrow(pairs.all[is.na(COMM.RECIPIENT)]), ' pairs\n')
  pairs.all <- pairs.all[!is.na(COMM.RECIPIENT)]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
} 
if(only.transmission.same.community ){
  cat('\nExcluding transmission events between communities (I->F or F->I) \n')
  cat('Removing ', nrow(pairs.all[COMM.SOURCE != COMM.RECIPIENT]), ' pairs\n')
  pairs.all <- pairs.all[COMM.SOURCE == COMM.RECIPIENT]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(!is.null(remove.pairs.from.rounds)){
  cat('\nExcluding pairs in inland community from round', remove.pairs.from.rounds, '\n')
  tmp <- df_round_inland[round %in% remove.pairs.from.rounds, list(min_exclusion = min(min_sample_date), 
                                                            max_exclusion = max(max_sample_date))]
  cat('Removing ', nrow(pairs.all[COMM.RECIPIENT == 'inland' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ]]), ' pairs\n')
  pairs.all <- pairs.all[!(COMM.RECIPIENT == 'inland' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ])]
  
  cat('\nExcluding pairs in fishing community from round', remove.pairs.from.rounds, '\n')
  tmp <- df_round_fishing[round %in% remove.pairs.from.rounds, list(min_exclusion = min(min_sample_date), 
                                                                   max_exclusion = max(max_sample_date))]
  cat('Removing ', nrow(pairs.all[COMM.RECIPIENT == 'fishing' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ]]), ' pairs\n')
  pairs.all <- pairs.all[!(COMM.RECIPIENT == 'fishing' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ])]
}

print.which.NA(pairs.all)
print.statements.about.pairs(copy(pairs.all))

# keep only pairs with source-recipient with a time of infection
pairs <- pairs.all[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
pairs[, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT := DATE_INFECTION.RECIPIENT < cutoff_date]
tab <- pairs[, list(count = .N), by = c('DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT', 'COMM.RECIPIENT', 'SEX.RECIPIENT')]
print_table(tab[order(DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT, COMM.RECIPIENT, SEX.RECIPIENT)])


#
# Find probability of observing a transmissing event
#

proportion_sampling <- get_proportion_sampling(pairs, incidence_cases, outfile.figures)


#
# PREPARE MAPS
#

# prepare age map
df_age <- get.age.map(pairs, age_bands_reduced = 4)
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))

# prepare direciton and commuity
df_direction <- get.df.direction()
df_community <- get.df.community()


#
# PREPARE STAN DATA
#

stan_data <- add_stan_data_base()
stan_data <- add_phylo_data(stan_data, pairs)
stan_data <- add_incidence_cases(stan_data, incidence_cases_round)
stan_data <- add_incidence_rates(stan_data, incidence_cases_round)
stan_data <- add_incidence_rates_lognormal_parameters(stan_data, incidence_cases_round)
stan_data <- add_offset(stan_data, eligible_count_round)
stan_data <- add_offset_time(stan_data, eligible_count_round)
stan_data <- add_offset_susceptible(stan_data, eligible_count_round)
stan_data <- add_probability_sampling(stan_data, proportion_sampling)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3,
                                      n_knots_rows = 6, n_knots_columns = 6,
                                      X = unique(df_age$AGE_TRANSMISSION.SOURCE),
                                      Y = unique(df_age$AGE_INFECTION.RECIPIENT))
stan_init <- add_init(stan_data)



#
# MAKE EXPLANATORY PLOTS
#

if(1){
  
  # find color palette of rounds
  find_palette_round()
  
  # plot count eligible susceptible / infected / infected unsuppressed and incident cases
  plot_data_by_round(eligible_count_round, proportion_unsuppressed, proportion_prevalence, incidence_cases_round, outfile.figures)
  plot_data_by_period(incidence_cases, outfile.figures)
  
  # plot tansmission events over time
  plot_transmission_events_over_time(pairs, outfile.figures)
  save_statistics_transmission_events(pairs, outdir.table)
  
  # plot incident rates and cases over time
  plot_incident_cases_over_time(incidence_cases_round, participation, outfile.figures)
  plot_incident_rates_over_time(incidence_cases_round, eligible_count_round, outfile.figures, outdir.table)
  plot_incident_cases_to_unsuppressed_rate_ratio(incidence_cases_round, unsuppressed_rate_ratio, outfile.figures, outdir.table)
    
  # plot offset
  plot_offset(stan_data, outfile.figures)
  
  # plot pair from chains
  # phsc.plot.transmission.network(copy(as.data.table(dchain)), copy(as.data.table(dc)), pairs,outdir=outfile, arrow=arrow(length=unit(0.02, "npc"), type="open"), edge.size = 0.1)
  plot_hist_age_infection(copy(pairs), outfile.figures)
  plot_age_infection_source_recipient(pairs[SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F'], 'Male -> Female', 'MF', outfile.figures)
  plot_age_infection_source_recipient(pairs[SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M'], 'Female -> Male', 'FM', outfile.figures)
  plot_pairs(pairs, outfile.figures)
  plot_CI_age_infection(pairs, outfile.figures)
  plot_CI_age_transmission(pairs, outfile.figures)

}


## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=paste0(outfile, "-stanin_",jobname,".RData"))

#
# RUN STAN DATA
#

# make stan model
model = rstan::stan_model(path.to.stan.model)

# sample
if(0){
  fit <- sampling(model, data = stan_data, iter = 10, warmup = 5, chains=1, thin=1, init = rep(list(stan_init), 1))
}else{
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- sampling(model, data = stan_data,
                  iter = 3000, warmup = 500, chains=4, thin=1, seed = 5,
                  verbose = FALSE, control = list(adapt_delta = 0.99,max_treedepth=15), 
                  init = rep(list(stan_init), 4))
}

file = paste0(outfile, "-stanout_", jobname, ".rds")
cat("Save file ", file)
saveRDS(fit,file = file)
