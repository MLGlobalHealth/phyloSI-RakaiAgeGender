library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
require(lubridate)
library(rstan)
library(gridExtra)
library(lognorm)
library(ggExtra)
library(Hmisc)
library(patchwork)

# laptop
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  gitdir <- '~/git/phyloSI-RakaiAgeGender'
  outdir <- '~/Box\ Sync/2021/phyloflows/'

  jobname <- 'firstrun'
  stan_model <- 'gp_230602'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

# hpc
if(dir.exists('/rds/general/user/mm3218/'))
{
  gitdir <- '~/git/phyloSI-RakaiAgeGender'
  outdir <- '/rds/general/user/mm3218/home/projects/2021/phyloSI-RakaiAgeGender'
  
  jobname <- 'test'
  stan_model <- 'gp_221201d'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

# andrea
if(dir.exists('/home/andrea'))
{
  gitdir <-'~/git/phyloSI-RakaiAgeGender'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'central'
  stan_model <- 'gp_221201d'
  # outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  # dir.create(outdir)
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-stan_model')
  stopifnot(args_line[[7]]=='-jobname')
  gitdir <- args_line[[2]]
  outdir <- args_line[[4]]
  stan_model <- args_line[[6]]
  jobname <- args_line[[8]]
}

source(file.path(gitdir, 'config.R'))

# load functions
source(file.path(gitdir.R, 'utils.R'))
source(file.path(gitdir.R.flow, 'summary_functions.R'))
source(file.path(gitdir.R.flow, 'plotting_functions.R'))
source(file.path(gitdir.R.flow, 'statistics_functions.R'))
source(file.path(gitdir.R.flow, 'stan_utils.R'))

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))
outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
if(!dir.exists(dirname(outfile.figures))) dir.create(dirname(outfile.figures))
if(!dir.exists(dirname(outdir.table))) dir.create(dirname(outdir.table))

# indicators -- fixed
only.transmission.after.start.observational.period <- TRUE
only.transmission.before.stop.observational.period <- TRUE
use_number_susceptible_offset <- TRUE
only.one.community <- 'inland'
use_contact_rates_prior <- TRUE

# indicators for sensitivity analyses: default values
nonparticipants.treated.like.participants <- FALSE
nonparticipants.not.treated <- FALSE
nonparticipants.male.relative.infection <- 1
nonparticipants.female.relative.infection <- 1
remove.pairs.from.rounds <- NULL
use_loess_inc_estimates <- FALSE
pairs_replicates.seed <- NULL
viremic_viral_load_200ml <- FALSE
use_30com_inc_estimates <- FALSE
use_30com_pairs <- FALSE
use_tsi_non_refined <- FALSE

# modify values based on jobname for sensivity analyses 
set.sensitivity.indicators.from.jobname(jobname)

# stan model
path.to.stan.model <- file.path(gitdir.stan, paste0(stan_model, '.stan'))

# load pairs
if(use_tsi_non_refined){
  pairs.all <- read_pairs(file.pairs.nonrefined)
} else{
  pairs.all <- read_pairs(file.pairs)
}

# load round timeline
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]

# load census eligible ount
eligible_count_smooth <- fread(file.eligible.count)

# load participation (% of census eligible population)
participation <- fread(file.participation)

# load proportion prevalence
proportion_prevalence <- fread(file.prevalence.prop)

# load proportion unsuppressed deep sequenced
proportion_unsuppressed_deepsequenced <- as.data.table(readRDS(file.prop.unsuppressed.deepsequenced))

# load non-suppressed proportion 
if(viremic_viral_load_200ml){
  treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants.vl200, 
                                              file.treatment.cascade.prop.nonparticipants.vl200)
  if(file.exists(file.treatment.cascade.prop.participants.vl200.samples) & file.exists(file.treatment.cascade.prop.nonparticipants.vl200.samples)){
    treatment_cascade_samples <- read_treatment_cascade_samples(file.treatment.cascade.prop.participants.vl200.samples, 
                                                                file.treatment.cascade.prop.nonparticipants.vl200.samples)
  }
}else{
  treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants, 
                                              file.treatment.cascade.prop.nonparticipants)
  if(file.exists(file.treatment.cascade.prop.participants.samples) & file.exists(file.treatment.cascade.prop.nonparticipants.samples)){
    treatment_cascade_samples <- read_treatment_cascade_samples(file.treatment.cascade.prop.participants.samples, 
                                                                file.treatment.cascade.prop.nonparticipants.samples)
  }
}

# load incidence estimates 
if(use_loess_inc_estimates){
  incidence.inland <- fread(file.incidence.loess.inland)
  file.incidence.samples.inland	<- file.incidence.loess.samples.inland	
}else if(use_30com_inc_estimates){
  incidence.inland <- fread(file.incidence.30com.inland)
  file.incidence.samples.inland	<- file.incidence.30com.samples.inland 	
} else{
  incidence.inland <- fread(file.incidence.inland)
}

# for offset
df_estimated_contact_rates <- as.data.table(readRDS(file.sexual.partnerships.rates))#estimated secxual contact rate
if('cntct.rates' %in% colnames(df_estimated_contact_rates))
  setnames(df_estimated_contact_rates, 'cntct.rates', 'cntct.rate')

#for plots
unsuppressed_rate_ratio <- fread(file.unsuppressed_rate_ratio) # sex ratio of unsuppression rate
df_reported_contact <- as.data.table(readRDS(file.number.sexual.partnerships)) # estimated number of sexual contacts
unsuppressed_share <- fread(file.unsuppressed.share) # share of unsuppressed count by sex
infected_share <- fread(file.prevalence.share) # share of infected count by sex
unsuppressed_median_age <-fread(file.unsuppressed_median_age) # median age of unsuppressed

#
# Define start time, end time and cutoff
#

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

df_period <- make.df.period(start_first_period_inland, stop_first_period_inland, 
                            start_second_period_inland, stop_second_period_inland, 
                            df_round)


#
# Find count eligible susceptible / infected / infected unsuppressed 
# 


# by round
eligible_count_round <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence, participation, 
                                                 nonparticipants.male.relative.infection, 
                                                 nonparticipants.female.relative.infection)
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

pairs <- copy(pairs.all)
pairs <- select.pairs.for.analysis(pairs.all, 
    only.one.community=only.one.community, 
    use_30com_pairs=use_30com_pairs, 
    only.transmission.after.start.observational.period=only.transmission.after.start.observational.period, 
    only.transmission.before.stop.observational.period=only.transmission.before.stop.observational.period, 
    remove.pairs.from.rounds=remove.pairs.from.rounds
)
#
print.which.NA(pairs)
print.statements.about.pairs(copy(pairs))

# keep only pairs with source-recipient with a time of infection
pairs <- pairs[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
pairs[, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT := DATE_INFECTION.RECIPIENT < start_second_period_inland]
tab <- pairs[, list(count = .N), by = c('DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT', 'COMM.RECIPIENT', 'SEX.RECIPIENT')]
print_table(tab[order(DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT, COMM.RECIPIENT, SEX.RECIPIENT)])

# replace pairs with a bootstrap sample 
if(!is.null(pairs_replicates.seed)){
  cat('\nSeed for replicate ', pairs_replicates.seed)
  set.seed(pairs_replicates.seed)
  stopifnot(all(duplicated(pairs) == F))
  pairs <- pairs[sample(nrow(pairs), replace = T)]
  stopifnot(any(duplicated(pairs) == T))
}

#
# Find detection probability of observing a transmissing event
#

proportion_sampling <- get_proportion_sampling(pairs, incidence_cases, incidence_cases_round, outfile.figures)


#
# PREPARE MAPS
#

# prepare age map
df_age <- get.age.map(age_bands_reduced = 4)
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
stan_data <- add_offset(stan_data, eligible_count_round, df_estimated_contact_rates,
                        use_number_susceptible_offset, use_contact_rates_prior)
stan_data <- add_offset_time(stan_data, df_round)
stan_data <- add_offset_susceptible(stan_data, eligible_count_round)
stan_data <- add_probability_sampling(stan_data, proportion_sampling)
stan_data <- add_probability_detection(stan_data, proportion_sampling)
stan_data <- add_probability_sampling_source(stan_data, proportion_unsuppressed_deepsequenced)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3,
                                      n_knots_rows = 10, n_knots_columns = 10,
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
  incidence_rates_round.samples <- load_incidence_rates_samples(file.incidence.samples.inland)# need to load incidence rates sample to compute statistics such as ratio
  plot_incident_cases_over_time(incidence_cases_round, participation, outfile.figures)
  plot_incident_rates_over_time(incidence_cases_round, incidence_rates_round.samples, eligible_count_round, outfile.figures, outdir.table)
  plot_incident_cases_to_unsuppressed_rate_ratio(incidence_cases_round, unsuppressed_rate_ratio, outfile.figures, outdir.table)
  make_time_trends_HIV_incidence(eligible_count_round, incidence_cases_round, incidence_rates_round.samples,outfile.figures)

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
  plot_date_collection_pairs(pairs, df_round_inland, outfile.figures)
  save_statistics_transmission_events(pairs, pairs.all, outdir.table)
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
