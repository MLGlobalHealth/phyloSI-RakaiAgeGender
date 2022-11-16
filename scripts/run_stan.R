library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
library(grid)
require(lubridate)
library(rstan)
library(gridExtra)

# laptop
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir <- '~/git/phyloflows'
  outdir <- '~/Box\ Sync/2021/phyloflows/'

  jobname <- 'new_offset'
  stan_model <- 'gp_221115a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/home/andrea'))
{
  indir <-'~/git/phyloflows'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'test'
  stan_model <- 'gp_220911a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
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
only.transmission.after.start.observational.period <- T
only.transmission.before.stop.observational.period <- T
nonparticipants.treated.like.participants <- F
nonparticipants.not.treated <- F
remove.pairs.from.rounds <- NULL
only.one.community <- 'inland'
use_number_susceptible_offset <- F
use_loess_inc_estimates <- F

# obtained in EMODO_RAKAI repo
file.incidence.inland	<- file.path(indir, 'data', "Rakai_incpredictions_inland_221107.csv")
file.incidence.loess.inland	<- file.path(indir, 'data', "Rakai_incpredictions_loess_inland_221116.csv")

# obtained in src/ for analysis
file.path.round.timeline <- file.path(indir, 'data', 'RCCS_round_timeline_220905.RData')
file.eligible.count <- file.path(indir, 'data', 'RCCS_census_eligible_individuals_220830.csv')
file.participation <- file.path(indir, 'data', 'RCCS_participation_220915.csv')
file.prevalence.prop <- file.path(indir, 'fit', 'RCCS_prevalence_estimates_220811.csv')

# obtained in misc/ for analysis
file.pairs <- file.path(indir, 'data', 'pairsdata_toshare_d1_w11.rds')
file.treatment.cascade.prop.participants <- file.path(indir, 'fit', "RCCS_treatment_cascade_participants_estimates_221101.csv")
file.treatment.cascade.prop.nonparticipants <- file.path(indir, 'fit', "RCCS_treatment_cascade_nonparticipants_estimates_221101.csv")

# obtained in misc/ for plots
file.unsuppressed.share <- file.path(indir, 'fit', paste0('RCCS_unsuppressed_share_sex_221101.csv'))
file.unsuppressed_rate_ratio <- file.path(indir, 'fit', paste0('RCCS_unsuppressed_ratio_sex_221101.csv'))
file.prevalence.share <- file.path(indir, 'fit', paste0('RCCS_prevalence_share_sex_220830.csv'))
file.reported.sexual.partnerships <- file.path(indir, 'data', paste0('age-age-group-est-cntcts-r15.rds'))

# obtained in EMODO_RAKAI repo for plots
file.incidence.samples.inland	<- file.path(indir, 'data', "Rakai_incpredictions_samples_inland_221107.csv")
file.incidence.loess.samples.inland	<- file.path(indir, 'data', "Rakai_incpredictions_loess_samples_inland_221116.csv")

path.to.stan.model <- file.path(indir, 'stan_models', paste0(stan_model, '.stan'))

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'statistics_functions.R'))
source(file.path(indir, 'functions', 'stan_utils.R'))

# load pairs
pairs.all <- read_pairs(file.pairs)

# load round timeline
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]

# load census eligible ount
eligible_count_smooth <- fread(file.eligible.count)

# load participation (% of census eligible population)
participation <- fread(file.participation)

# load proportion prevalence
proportion_prevalence <- fread(file.prevalence.prop)

# load non-suppressed proportion 
treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants, 
                                                             file.treatment.cascade.prop.nonparticipants)

# load incidence estimates 
if(use_loess_inc_estimates){
  incidence.inland <- fread(file.incidence.loess.inland)
  file.incidence.samples.inland	<- file.incidence.loess.samples.inland	
}else{
  incidence.inland <- fread(file.incidence.inland)
}

#for plots
unsuppressed_rate_ratio <- fread(file.unsuppressed_rate_ratio) # sex ratio of unsuppression rate
df_reported_contact <- as.data.table(readRDS(file.reported.sexual.partnerships)) # reported sexual contacts
unsuppressed_share <- fread(file.unsuppressed.share) # share of unsuppressed count by sex
infected_share <- fread(file.prevalence.share) # share of infected count by sex


#
# Define start time, end time and cutoff
#

start_first_period_inland <- df_round_inland[round == 'R010', min_sample_date] # "2003-09-26"
stop_first_period_inland <- df_round_inland[round == 'R015', max_sample_date] # "2013-07-05"
start_second_period_inland <-df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"
stop_second_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

stopifnot(start_first_period_inland < stop_first_period_inland)
stopifnot(stop_first_period_inland < start_second_period_inland)
stopifnot(start_second_period_inland < stop_second_period_inland)

df_round <- make.df.round(df_round_inland)

df_period <- make.df.period(start_first_period_inland, stop_first_period_inland, 
                            start_second_period_inland, stop_second_period_inland, 
                            df_round)


#
# Find count eligible susceptible / infected / infected unsuppressed 
# 


# by round
eligible_count_round <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence)
eligible_count_round <- add_infected_unsuppressed(eligible_count_round, treatment_cascade, participation, 
                                                  nonparticipants.treated.like.participants, 
                                                  nonparticipants.not.treated)
eligible_count_round[, table(ROUND, COMM)]


#
# Find incidence cases
#

# by round
incidence_cases_round <- get_incidence_cases_round(incidence.inland, eligible_count_round)
incidence_cases_round[, table(ROUND, COMM)]

# summarise by time period
incidence_cases <- summarise_incidence_cases_period(incidence_cases_round, df_period)
incidence_cases[, table(PERIOD, COMM)]


#
# Find phylo pairs 
#
pairs <- copy(pairs.all)
if(1)
{
  cat('Keep only heterosexual pairs\n')
  cat('Removing ', nrow(pairs[! ((SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M'))]), ' pairs\n')
  pairs <- pairs[(SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M')]
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
  
  cat('Keep only RCCS participants\n')
  cat('Removing ', nrow(pairs[(COMM.SOURCE == 'neuro' | COMM.RECIPIENT == "neuro")]), ' pairs\n')
  pairs <- pairs[COMM.SOURCE != 'neuro' & COMM.RECIPIENT != "neuro"]
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(!is.null(only.one.community)){
  cat('\nExcluding sources and recipients in ',   pairs[COMM.RECIPIENT != only.one.community, unique(COMM.RECIPIENT)] ,'\n')
  cat('Removing ', nrow(pairs[COMM.RECIPIENT != only.one.community]), ' pairs\n')
  pairs <- pairs[COMM.RECIPIENT == only.one.community]
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(only.transmission.after.start.observational.period){
  cat('\nFor inland excluding recipients infected before ', as.character(start_first_period_inland), '\n')
  cat('Removing ', nrow(pairs[DATE_INFECTION.RECIPIENT < start_first_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs <- pairs[!(DATE_INFECTION.RECIPIENT < start_first_period_inland & COMM.RECIPIENT == 'inland')]

  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(only.transmission.before.stop.observational.period){
  cat('\nFor inland excluding recipients infected after ', as.character(stop_second_period_inland), '\n')
  cat('Removing ', nrow(pairs[DATE_INFECTION.RECIPIENT > stop_second_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs <- pairs[!(DATE_INFECTION.RECIPIENT > stop_second_period_inland & COMM.RECIPIENT == 'inland')]
  
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(!is.null(remove.pairs.from.rounds)){
  cat('\nExcluding pairs in inland community from round', remove.pairs.from.rounds, '\n')
  tmp <- df_round_inland[round %in% remove.pairs.from.rounds, list(min_exclusion = min(min_sample_date), 
                                                            max_exclusion = max(max_sample_date))]
  cat('Removing ', nrow(pairs[COMM.RECIPIENT == 'inland' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ]]), ' pairs\n')
  pairs <- pairs[!(COMM.RECIPIENT == 'inland' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ])]
  
  cat('\nExcluding pairs in fishing community from round', remove.pairs.from.rounds, '\n')
  tmp <- df_round_fishing[round %in% remove.pairs.from.rounds, list(min_exclusion = min(min_sample_date), 
                                                                   max_exclusion = max(max_sample_date))]
  cat('Removing ', nrow(pairs[COMM.RECIPIENT == 'fishing' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ]]), ' pairs\n')
  pairs <- pairs[!(COMM.RECIPIENT == 'fishing' & DATE_INFECTION.RECIPIENT <= tmp[, max_exclusion] & DATE_INFECTION.RECIPIENT >= tmp[,min_exclusion ])]
}

print.which.NA(pairs)
print.statements.about.pairs(copy(pairs))

# keep only pairs with source-recipient with a time of infection
pairs <- pairs[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
pairs[, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT := DATE_INFECTION.RECIPIENT < start_second_period_inland]
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
stan_data <- add_offset(stan_data, eligible_count_round, use_number_susceptible_offset)
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
  plot_data_by_round(eligible_count_round, treatment_cascade, proportion_prevalence, outfile.figures)
  plot_data_by_period(incidence_cases, outfile.figures)
  
  # plot incident rates and cases over time
  plot_incident_cases_over_time(incidence_cases_round, participation, outfile.figures)
  incidence_rates_round.samples <- load_incidence_rates_samples(file.incidence.samples.inland)# need to load incidence rates sample to compute statistics such as ratio
  plot_incident_rates_over_time(incidence_cases_round, incidence_rates_round.samples, eligible_count_round, outfile.figures, outdir.table)
  plot_incident_cases_to_unsuppressed_rate_ratio(incidence_cases_round, unsuppressed_rate_ratio, outfile.figures, outdir.table)
    
  # plot offset
  plot_offset(stan_data, outfile.figures)
  
  # plot pair from chains
  plot_hist_age_infection(copy(pairs), outfile.figures)
  plot_hist_time_infection(copy(pairs), start_second_period_inland, outfile.figures)
  plot_age_infection_source_recipient(pairs[SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F'], 'Male -> Female', 'MF', start_second_period_inland, outfile.figures)
  plot_age_infection_source_recipient(pairs[SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M'], 'Female -> Male', 'FM', start_second_period_inland, outfile.figures)
  plot_CI_age_infection(pairs, start_second_period_inland, outfile.figures)
  plot_CI_age_transmission(pairs, start_second_period_inland, outfile.figures)
  plot_pairs(pairs, outfile.figures)
  plot_pairs_all(pairs.all, outfile.figures)
  plot_transmission_events_over_time(pairs, outfile.figures)
  save_statistics_transmission_events(pairs, outdir.table)
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
