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

  jobname <- 'test_new'
  stan_model <- 'gp_220901a'
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
  stan_model <- 'gp_220108'
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
if(!dir.exists(dirname(outfile.figures))) dir.create(dirname(outfile.figures))

# indicators
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.5
use.tsi.estimates <- F
remove.inconsistent.infection.dates <- F
remove.young.individuals <- T
remove.missing.community.recipient <- T
remove.neuro.individuals <- T
only.transmission.after.start.observational.period <- T
only.transmission.before.stop.observational.period <- T
use.diagonal.prior <- F
use.informative.prior <- F
only.transmission.same.community <- F


# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.path.meta <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_20220329.RData')
file.path.tsiestimates <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', 'TSI_estimates_220119.csv')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')

file.incidence.inland	<- file.path(indir.deepsequencedata, 'RCCS_R15_R18', "Rakai_incpredictions_220524.csv")
file.incidence.fishing	<- file.path(indir.deepsequencedata, 'RCCS_R15_R18', "Rakai_incpredictions_fishing_220825.csv")

file.eligible.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220830.csv')
# file.unsuppressed.prop <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', "RCCS_nonsuppressed_proportion_arvmed_220801.csv")
file.unsuppressed.prop <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', "RCCS_nonsuppressed_proportion_vl_1000_220803.csv")
file.unsuppressed.share <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_share_sex_vl_1000_220830.csv'))
#file.partnership.rate <- file.path(indir.deepsequence_analyses,'RCCS_partnership_rate_220422.csv')
file.prevalence.prop <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_prevalence_estimates_220811.csv')
file.prevalence.share <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('RCCS_prevalence_share_sex_220830.csv'))

path.to.stan.model <- file.path(indir, 'stan_models', paste0(stan_model, '.stan'))

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'stan_utils.R'))
source(file.path(indir, 'functions', 'check_potential_TNet.R'))

# load chains
load(file.path.chains.data)
dchain <- as.data.table(dchain)

# load meta data
load(file.path.meta)

# load Tanya's estimate time since infection using phylogenetic data
time.since.infection <- make.time.since.infection(as.data.table(read.csv(file.path.tsiestimates)))

# load census eligible ount
eligible_count_smooth <- as.data.table(read.csv(file.eligible.count))

# load proportion prevalence
proportion_prevalence <- as.data.table(read.csv(file.prevalence.prop))

# load non-suppressed proportion 
proportion_unsuppressed <- as.data.table(read.csv(file.unsuppressed.prop))

# load anonymous aid
aik <- .read(file.anonymisation.keys); aik$X <- NULL

# load incidence estimates from Adam
incidence.inland <- as.data.table(read.csv(file.incidence.inland))
incidence.fishing <- as.data.table(read.csv(file.incidence.fishing))

#
# Define start time, end time and cutoff
#

start_observational_period_inland <- df_round_inland[round == 'R014', min_sample_date] #"2010-01-20" 
stop_observational_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

start_observational_period_fishing <- start_observational_period_inland #"2010-01-20" 
stop_observational_period_fishing <- df_round_inland[round == 'R018', max_sample_date] #  "2017-08-14"

cutoff_date <- df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"

stopifnot(start_observational_period_inland <= cutoff_date & stop_observational_period_inland >= cutoff_date)
stopifnot(start_observational_period_fishing <= cutoff_date & stop_observational_period_fishing >= cutoff_date)

df_period <- make.df.period(start_observational_period_inland, stop_observational_period_inland, 
                            start_observational_period_fishing, stop_observational_period_fishing, 
                            cutoff_date)

df_round <- make.df.round(df_round_inland, df_round_fishing, df_period)


#
# Find count eligible susceptible / infected / infected unsuppressed 
# 

# by round
eligible_count_round <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence)
eligible_count_round <- add_infected_unsuppressed(eligible_count_round, proportion_unsuppressed)
eligible_count_round[, table(ROUND, COMM)]

# summarise by time period
eligible_count <- summarise_eligible_count_period(eligible_count_round, cutoff_date, df_period)
eligible_count[, table(PERIOD, COMM)]


#
# Find incidence cases
#

# by round
incidence_cases_round <- get_incidence_cases_round(incidence.inland, incidence.fishing, eligible_count_round)
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

if(include.only.heterosexual.pairs){
  cat('Keep only heterosexual pairs\n')
  cat('Removing ', nrow(pairs.all[! ((SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M'))]), ' pairs\n')
  pairs.all <- pairs.all[(SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M')]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
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

if(only.transmission.after.start.observational.period){
  cat('\nFor inland excluding recipients infected before ', as.character(start_observational_period_inland), '\n')
  cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland')]
  
  cat('\nFor fishing excluding recipients infected before ', as.character(start_observational_period_fishing), '\n')
  cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT < start_observational_period_fishing & COMM.RECIPIENT == 'fishing']), ' pairs\n')
  pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT < start_observational_period_fishing & COMM.RECIPIENT == 'fishing')]
  
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(only.transmission.before.stop.observational.period){
  cat('\nFor inland excluding recipients infected after ', as.character(stop_observational_period_inland), '\n')
  cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland')]
  
  cat('\nFor fishing excluding recipients infected after ', as.character(stop_observational_period_fishing), '\n')
  cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT > stop_observational_period_fishing & COMM.RECIPIENT == 'fishing']), ' pairs\n')
  pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT > stop_observational_period_fishing & COMM.RECIPIENT == 'fishing')]
  
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

print.which.NA(pairs.all)
print.statements.about.pairs(copy(pairs.all))

# which base frequency files we have on the HPC
# atm gives error: maybe TODO when I understand more about PHSC pipeline
# missing_bff <- print.statements.about.basefreq.files(pairs.all)

# keep only pairs with source-recipient with proxy for the time of infection
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
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3,
                                      n_knots_rows = 6, n_knots_columns = 6,
                                      X = unique(df_age$AGE_TRANSMISSION.SOURCE),
                                      Y = unique(df_age$AGE_INFECTION.RECIPIENT))
stan_data <- add_incidence_cases(stan_data, incidence_cases_round)
stan_data <- add_offset(stan_data, eligible_count_round)
stan_data <- add_probability_sampling(stan_data, proportion_sampling)
stan_init <- add_init(stan_data)


#
# MAKE EXPLANATORY PLOTS
#

if(1){
  # plot count eligible susceptible / infected / infected unsuppressed and incident cases
  plot_data_by_round(eligible_count_round, proportion_unsuppressed, proportion_prevalence, incidence_cases_round, outfile.figures)
  plot_data_by_period(incidence_cases, outfile.figures)
  
  # plot tansmission events over time
  plot_transmission_events_over_time(eligible_count_round, incidence_cases_round, pairs, outfile.figures)
    
  # plot offset
  plot_offset(stan_data, outfile.figures)
  
  # plot pair from chains
  # plot_pairs_infection_dates(pairs.all, outfile.figures)
  # phsc.plot.transmission.network(copy(as.data.table(dchain)), copy(as.data.table(dc)), pairs,outdir=outfile, arrow=arrow(length=unit(0.02, "npc"), type="open"), edge.size = 0.1)
  plot_hist_age_infection(copy(pairs), outfile.figures)
  plot_age_infection_source_recipient(pairs[SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F'], 'Male -> Female', 'MF', outfile.figures)
  plot_age_infection_source_recipient(pairs[SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M'], 'Female -> Male', 'FM', outfile.figures)
  plot_pairs(pairs, outfile.figures)
  plot_CI_age_infection(pairs, outfile.figures)
  plot_CI_age_transmission(pairs, outfile.figures)

}


# for now ignore fishing
if(0){
  stan_data[['y']][,,1,] =  stan_data[['y']][,,2,]
  stan_data[['z']][,,1,] =  stan_data[['z']][,,2,]
  stan_data[['sampling_index_y']][,,1,] =  stan_data[['sampling_index_y']][,,2,]
  stan_data[['log_offset']][,1,,] =  stan_data[['log_offset']][,2,,]
  stan_data[['log_prop_sampling']][,1,,] =  stan_data[['log_prop_sampling']][,2,,]
  stan_data[['n_sampling_index_y']][,1,] =stan_data[['n_sampling_index_y']][,2,] 
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


