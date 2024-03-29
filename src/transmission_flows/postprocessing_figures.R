cat("Start of postprocessing_figures.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)
library(ggnewscale)

usr <- Sys.info()[['user']]
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

# paths
path.to.stan.output = paste0(outfile, "-stanout_", jobname, ".rds")
.outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
.outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))

# load data
path.to.stan.data <- paste0(outfile, "-stanin_",jobname,".RData")
load(path.to.stan.data)
outfile.figures <- .outfile.figures
outdir.table <- .outdir.table

# load functions
gitdir <- copy(indir)
source(file.path(gitdir, 'config.R'))
source(file.path(gitdir.R.flow, 'postprocessing_summary_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_plot_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_utils_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_statistics_functions.R'))

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)


#
# offset
#

log_offset_round <- find_log_offset_by_round(stan_data, eligible_count_round, df_estimated_contact_rates, 
                                             use_number_susceptible_offset, use_contact_rates_prior)

# log offset formula (per year)
log_offset_formula <- 'log_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED'
if(!use_number_susceptible_offset)
  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED'
if(use_contact_rates_prior)
  log_offset_formula = paste0(log_offset_formula, ' + log_CONTACT_RATES')

# log offset formula (per year per unsuppressed)
log_offset_formula_perunsuppressed <- 'log_SUSCEPTIBLE'
if(!use_number_susceptible_offset)
  log_offset_formula_perunsuppressed = 'log_PROP_SUSCEPTIBLE'
if(use_contact_rates_prior)
  log_offset_formula_perunsuppressed = paste0(log_offset_formula_perunsuppressed, ' + log_CONTACT_RATES')

# log offset formula (per year per susceptible)
log_offset_formula_persusceptible <- 'log_INFECTED_NON_SUPPRESSED'
if(use_contact_rates_prior)
  log_offset_formula_persusceptible = paste0(log_offset_formula_persusceptible, ' + log_CONTACT_RATES')


#
# Summarise data and merge to maps for figures
#

count_data <- prepare_count_data(stan_data)
incidence_cases_recipient_round <- prepare_incidence_cases(incidence_cases_round)
unsuppressed_share_sex <- prepare_unsuppressed_share(unsuppressed_share, c('SEX'), outdir.table)
unsuppressed_share_sex_age <- prepare_unsuppressed_share(unsuppressed_share, c('SEX', 'AGEYRS'), outdir.table)
unsuppressed_share_sex_age_stand <- prepare_unsuppressed_share(unsuppressed_share, c('SEX', 'AGEYRS'), outdir.table, 'SEX')
prevalence_prop_sex <- prepare_infected_share(infected_share, 'SEX', outdir.table)
reported_contact <- clean_reported_contact(df_reported_contact, outdir.table)
df_unsuppressed_median_age<-prepare_unsuppressed_median_age(unsuppressed_median_age, outdir.table)


#
## PPC
#

cat("\nPlot PPC\n")

# rate of observed transmission
intensity_PP_sampled <- find_summary_output(samples, 'log_lambda', c('INDEX_DIRECTION', 'INDEX_TIME', 'INDEX_AGE'), transform = 'exp')
plot_intensity_PP(intensity_PP_sampled, count_data, outfile.figures)

# rate of transmission
intensity_PP <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'INDEX_AGE'),
                                             transform = 'exp',
                                             log_offset_round = log_offset_round,
                                             log_offset_formula = log_offset_formula)
plot_intensity_PP_by_round(intensity_PP, outfile.figures)

# predicted observed transmission
predict_y_source <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'))
predict_y_recipient <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'))
predict_y_source_recipient <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'))
plot_PPC_observed_source(predict_y_source, count_data, outfile.figures)
plot_PPC_observed_recipient(predict_y_recipient, count_data, outfile.figures)

# predicted transmission
predict_z_source <- find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'))
predict_z_recipient_round <- find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'))
plot_PPC_augmented_recipient_round(predict_z_recipient_round, incidence_cases_recipient_round, outfile.figures)

# predicted incidence rate
predict_incidence_rate_round <- find_summary_output_by_round(samples, 'ir_predict', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                             names = c('INDEX_AGE_INFECTION.RECIPIENT', 'INDEX_DIRECTION', 'INDEX_ROUND'))
plot_PPC_incidence_rate_round(predict_incidence_rate_round, incidence_cases_recipient_round,outfile.figures)
save_statistics_PPC(predict_y_source_recipient, count_data, predict_incidence_rate_round, incidence_cases_recipient_round, outdir.table)

# predicted observed transmission vs all transmission
plot_observed_to_augmented(predict_y_source, predict_z_source, outfile.figures)


#
## Probability of detecting a transmission
#

cat("\nPlot probability of detecting a transmisison\n")

if('log_thinning_prob' %in% names(samples)){
  # estimated probability of sampling a recipient given that the source is sampled
  cond_prob_sampling_recipient <- find_summary_output_by_round(samples, 'log_conditional_sampling_recipient_prob', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), transform = 'exp')
  plot_cond_prob_sampling_recipient(cond_prob_sampling_recipient, outfile.figures)
  
  # estimated probability of detecting a transmission event from source i to recipient j
  prob_thinning <- find_summary_output_by_round(samples, 'log_thinning_prob', c('INDEX_DIRECTION', 'INDEX_ROUND', 'INDEX_AGE'), transform = 'exp')
  plot_prob_thinning(prob_thinning, outfile.figures)
}


#
## force of infection
#

cat("\nPlot force of infection\n")

# 2D for all categories
force_infection <- find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_ROUND', 'INDEX_AGE'), transform = 'exp')
plot_force_infection(force_infection, outfile.figures)

# shift in sex-specific transmission dynamics by period
force_infection_sex_source <- find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_ROUND'), transform = 'exp')
plot_force_infection_sex_source(force_infection_sex_source, outfile.figures)

# shift in age source by period
force_infection_age_source <- find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'), transform = 'exp')
plot_force_infection_sex_age_source(force_infection_age_source, outfile.figures)

# shift in age source by round
force_infection_age_recipient <-  find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), transform = 'exp')
plot_force_infection_sex_age_recipient(force_infection_age_recipient, outfile.figures)


#
# median age of source
#

cat("\nPlot median age at transmission of the source by age at infection of recipient\n")

# by 1-year age band
median_age_source <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                  median_age_source = T)
plot_median_age_source(median_age_source, outfile.figures)

# by age groups
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
median_age_source_group <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                        transform = 'exp',
                                                        standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                        quantile_age_source = T, 
                                                        lab = '5yrs_age_band')
expected_contribution_age_group_source2 <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_ROUND'), 
                                                                        lab = '5yrs_age_band')
plot_median_age_source_group(median_age_source_group, expected_contribution_age_group_source2, reported_contact, outfile.figures)
plot_median_age_source_group_all_rounds(median_age_source_group, expected_contribution_age_group_source2, outfile.figures)

# by age group (other aggregation)
df_age_aggregated <- get.age.aggregated.map(c('15-29', '30-49'))
median_age_source_group2 <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                         transform = 'exp',
                                                         standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                         quantile_age_source = T)

# total
median_age_source <- find_summary_output_by_round(samples, 'log_lambda_latent', 
                                                  c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                                  quantile_age_source = T, 
                                                  lab = 'median_age_source_total')

# save
save_median_age_source(median_age_source_group, median_age_source_group2, median_age_source, outdir.table)


#
# median age difference male to female
#

cat("\nPlot median age difference at transmission\n")

# age group recipient
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
median_age_diff_group1 <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                        transform = 'exp',
                                                        standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                        quantile_age_diff = T)

df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))
median_age_diff_group2 <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                       transform = 'exp',
                                                       standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                       quantile_age_diff = T)

# total
median_age_diff <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'),
                                                      transform = 'exp',
                                                      standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                                      quantile_age_diff = T)
save_median_age_diff(median_age_diff_group1, median_age_diff_group2, outdir.table)
# Contribution to transmission
#

cat("\nPlot contribution\n")

# sex-specific contribution to transmission
contribution_sex_source <-  find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                                         standardised.vars = c('INDEX_ROUND'))
plot_contribution_sex_source(contribution_sex_source, unsuppressed_share_sex, prevalence_prop_sex, outfile.figures)

# age-specific contribution to transmission among all sources by sex
contribution_age_source <-  find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                         standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_source_unsuppressed(contribution_age_source, unsuppressed_share_sex_age, median_age_source, df_unsuppressed_median_age, outfile.figures)
plot_contribution_age_source(contribution_age_source, median_age_source, outfile.figures)

# contribution aggregated by age group of sources and recipients
contribution_age_group_source <- find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                               standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_group(contribution_age_group_source, outfile.figures)

# contribution aggregated by age group of recipients and classification of age of sources
contribution_age_classification_source <-  find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_classification(contribution_age_classification_source, outfile.figures)


#
# Expected Contribution to transmission
#

cat("\nPlot expected contribution\n")

# sex-specific contribution to transmission
expected_contribution_sex_source <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                                        transform = 'exp',
                                                        standardised.vars = c('INDEX_ROUND'))
plot_contribution_sex_source(expected_contribution_sex_source, unsuppressed_share_sex, prevalence_prop_sex, outfile.figures,'Expected_contribution')

# age-specific contribution to transmission among all age by sex
expected_contribution_age_source <- find_summary_output_by_round(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                 transform = 'exp',
                                                                 standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND'))

# age-specific contribution to transmission among all age source and reicpient by sex
expected_contribution_age_source_recipient <- find_summary_output_by_round(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'),
                                                                 transform = 'exp',
                                                                 standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND'))

# age-specific contribution to transmission among all age and  sex
expected_contribution_age_source2 <- find_summary_output_by_round(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                  transform = 'exp',
                                                                  standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_source_unsuppressed(expected_contribution_age_source2, unsuppressed_share_sex_age, median_age_source, df_unsuppressed_median_age, outfile.figures,'Expected_contribution')
plot_contribution_age_source(expected_contribution_age_source2, median_age_source, outfile.figures,'Expected_contribution_sex')

# age-specific sex ratio contribution to transmission
expected_contribution_age_source_sex_ratio <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                           c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                           transform = 'exp',
                                                                           standardised.vars = c('INDEX_ROUND'),
                                                                           sex_ratio= T)
plot_contribution_age_source_sex_ratio(expected_contribution_age_source_sex_ratio, outfile.figures,'Expected_contribution_sex_ratio')

# contribution aggregated by age group recipients and 1-year age band sources
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
expected_contribution_age_ungroup_source <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                       c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                       transform = 'exp',
                                                                       standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_ungroup(expected_contribution_age_ungroup_source, outfile.figures,'Expected_contribution')

# contribution aggregated by age group of sources and recipients
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))
expected_contribution_age_group_source <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_group(expected_contribution_age_group_source, outfile.figures,'Expected_contribution')

# contribution aggregated by age group of recipients and classification of age of sources
expected_contribution_age_classification_source <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_classification(expected_contribution_age_classification_source, outfile.figures,'Expected_contribution')

# contribution aggregated by age group of recipients and classification of age of male
expected_contribution_age_classification_male <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                                c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                                transform = 'exp',
                                                                                add_male_age_classification_nonsymmetric = T,
                                                                                standardised.vars = c('INDEX_ROUND'))
plot_contribution_age_classification_male(expected_contribution_age_classification_male, outfile.figures,'Expected_contribution')

#
# Make table 1
# contribution aggregated by age group source and sex
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
expected_contribution_age_group_source_total <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                             c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE'),
                                                                             transform = 'exp',
                                                                             standardised.vars = c('INDEX_ROUND', 'INDEX_DIRECTION'))
save_statistics_expected_contribution(expected_contribution_sex_source, expected_contribution_age_group_source_total, outdir.table)

#
# Make table 2
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))

# contribution aggregated by age group of recipient 
expected_contribution_sex_age_group_recipient <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                      c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                      transform = 'exp',
                                      standardised.vars = c('INDEX_ROUND'))
expected_contribution_age_group_recipient <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                      c('INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                      transform = 'exp',
                                      standardised.vars = c('INDEX_ROUND'))

# contribution aggregated by  classification of age of male
expected_contribution_sex_age_classification_male <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                      c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                      transform = 'exp',
                                      add_male_age_classification_nonsymmetric = T,
                                      standardised.vars = c('INDEX_ROUND'))
expected_contribution_age_classification_male_total <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                      c('INDEX_ROUND'),
                                      transform = 'exp',
                                      add_male_age_classification_nonsymmetric = T,
                                      standardised.vars = c('INDEX_ROUND'))

tt <- make_transmission_flows_table(expected_contribution_age_classification_male, expected_contribution_sex_age_group_recipient, 
                                          expected_contribution_sex_age_classification_male, expected_contribution_age_group_recipient, 
                                          expected_contribution_age_classification_male_total, expected_contribution_sex_source, outdir.table)

#
# Transmission flows (absolute)
#

cat("\nPlot tranmission flows (absolute)\n")

# find number of transmission by age group of sources and recipients
incidence_age_group_source  <- find_summary_output_by_round(samples, 'log_beta', 
                                                            c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                            transform = 'exp',
                                                            log_offset_round = log_offset_round)
plot_transmission_age_group(incidence_age_group_source, outfile.figures)

# find number of transmission by age group of recipients and classification of age of sources
incidence_age_classification_source  <- find_summary_output_by_round(samples, 'log_beta', 
                                                            c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                            transform = 'exp',
                                                            log_offset_round = log_offset_round)
plot_transmission_age_classification(incidence_age_classification_source, outfile.figures)

# find number of transmission by age group of recipients and classification of age of male
incidence_age_classification_male  <- find_summary_output_by_round(samples, 'log_beta', 
                                                                   c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                   transform = 'exp',
                                                                   add_male_age_classification_nonsymmetric = T,
                                                                   log_offset_round = log_offset_round)
plot_transmission_age_classification_male(incidence_age_classification_male, outfile.figures)


#
# Transmission risk per unsuppressed
#

cat("\nPlot transmission risk\n")


# sex-specific transmission risk
transmission_risk_sex_source <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                                             transform = 'exp',
                                                             log_offset_round = log_offset_round,
                                                             log_offset_formula = log_offset_formula,
                                                             per_unsuppressed = T)
plot_transmission_risk_sex_source(transmission_risk_sex_source, outfile.figures)

# sex and age-specific  transmission risk
transmission_risk_age_source <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                            transform = 'exp',
                                                            log_offset_round = log_offset_round,
                                                            log_offset_formula = log_offset_formula,
                                                            per_unsuppressed = T)
plot_transmission_risk_sex_age_source(transmission_risk_age_source, outfile.figures)


#
# Incidence infection
#

cat("\nPlot incidence infection and transmission\n")

# finer age bands
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))

#find incidence transmission
incidence_tranmission <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE'),
                                                      transform = 'exp',
                                                      log_offset_round = log_offset_round,
                                                      log_offset_formula = log_offset_formula,
                                                      relative_baseline = T,
                                                      per_susceptible = T)
plot_incidence_transmission(incidence_tranmission, outfile.figures)

#find incidence infection
incidence_infection <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                    transform = 'exp',
                                                    log_offset_round = log_offset_round,
                                                    log_offset_formula = log_offset_formula,
                                                    relative_baseline = T,
                                                    per_susceptible = T)
plot_incidence_infection(incidence_infection, outfile.figures)


#
# Counterfactual: comparison of targets
#

cat("\nPlot relative incidence infection if different groups of male are targeted\n")


# find incidence under the factual scenario by sex and age
incidence_factual <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  log_offset_formula = log_offset_formula_persusceptible,
                                                  log_offset_round = log_offset_round)

# identify the main spreaders, sources that contributes the most to incidence
spreaders <- find_spreaders(expected_contribution_age_source, outdir.table)

# generate counterfactual when only participant are treated 
# we compare treating main spreaders, male with the greatest diff in art coverage compared to female and random male
# in all scenario, the number of participants treated in the counterfactual are = number of participant sources that contribute to 1/3 incidence
counterfactuals_p_a <- make_counterfactual_target(samples, spreaders, log_offset_round, stan_data,
                                                  eligible_count_smooth, eligible_count_round, 
                                                  treatment_cascade, proportion_prevalence, participation,
                                                  only_participant = T, art_up_to_female = F, outdir.table)
plot_counterfactual_one(counterfactuals_p_a, incidence_factual, "Diagnosed unsuppressed", outfile.figures)

# generate counterfactual when participants and non-participant are treated
# we compare treating main spreaders, male with the greatest diff in art coverage compared to female and random male
# the number of male treated in the counterfactual are = number of male sources that contribute to 1/3 incidence
counterfactuals_a_a <- make_counterfactual_target(samples, spreaders, log_offset_round, stan_data,
                                                  eligible_count_smooth, eligible_count_round, 
                                                  treatment_cascade, proportion_prevalence,participation,
                                                  only_participant = F, art_up_to_female = F, outdir.table)
plot_counterfactual_one(counterfactuals_a_a, incidence_factual, "Unsuppressed", outfile.figures)
plot_counterfactual_strategy(counterfactuals_a_a, incidence_factual, 'Unsuppressed', outfile.figures)
  

#
# Counterfactual: comparison of strategies
#

cat("\nPlot relative incidence infection if different number of male are treated\n")

##
## generate counterfactual treating only men participant as much as female are diagnosed/treated/suppressed

if(exists("treatment_cascade_samples")){
  
  counterfactuals_p_f <- make_counterfactual(samples, log_offset_round, stan_data, 
                                             eligible_count_smooth, eligible_count_round, 
                                             treatment_cascade_samples, proportion_prevalence, participation,
                                             only_participant = T, art_up_to_female = 1, s959595 = NULL, s909090 = NULL, outdir.table)
  #  generate counterfactual treating only men participant half way to as much as female are diagnosed/treated/suppressed
  counterfactuals_p_f05 <- make_counterfactual(samples, log_offset_round, stan_data, 
                                               eligible_count_smooth, eligible_count_round, 
                                               treatment_cascade_samples, proportion_prevalence, participation,
                                               only_participant = T, art_up_to_female = 0.5, s959595 = NULL, s909090 = NULL, outdir.table)
  # generate counterfactual treating only men participant 95 95 95
  counterfactuals_p_959595 <- make_counterfactual(samples, log_offset_round, stan_data, 
                                                  eligible_count_smooth, eligible_count_round, 
                                                  treatment_cascade_samples, proportion_prevalence, participation,
                                                  only_participant = T, art_up_to_female = NULL, s959595 = 1, s909090 = NULL, outdir.table)
  # generate counterfactual treating only men participant 90 90 90
  counterfactuals_p_909090 <- make_counterfactual(samples, log_offset_round, stan_data, 
                                                  eligible_count_smooth, eligible_count_round, 
                                                  treatment_cascade_samples, proportion_prevalence, participation,
                                                  only_participant = T, art_up_to_female = NULL, s959595 = NULL, s909090 = 1, outdir.table)
  
  # plot
  plot_counterfactual(counterfactuals_p_f, counterfactuals_p_f05, counterfactuals_p_959595, counterfactuals_p_909090, 
                      incidence_factual, "Diagnosed unsuppressed", outfile.figures)
  
  # save
  save_counterfactual_results(counterfactuals_p_f, counterfactuals_p_f05, 
                              counterfactuals_p_959595, counterfactuals_p_909090, "Diagnosed unsuppressed", outdir.table)
  
}

##
## generate counterfactual treating all men as much as female are diagnosed/treated/suppressed

if(exists("treatment_cascade_samples")){
  
  counterfactuals_a_f <- make_counterfactual(samples, log_offset_round, stan_data, 
                                             eligible_count_smooth, eligible_count_round, 
                                             treatment_cascade_samples, proportion_prevalence, participation,
                                             only_participant = F, art_up_to_female = 1, s959595 = NULL, s909090 = NULL, outdir.table)
  #  generate counterfactual treating all men half way to as much as female are diagnosed/treated/suppressed
  counterfactuals_a_f05 <- make_counterfactual(samples, log_offset_round, stan_data, 
                                               eligible_count_smooth, eligible_count_round, 
                                               treatment_cascade_samples, proportion_prevalence, participation,
                                               only_participant = F, art_up_to_female = 0.5, s959595 = NULL, s909090 = NULL, outdir.table)
  # generate counterfactual treating all men 95 95 95
  counterfactuals_a_959595 <- make_counterfactual(samples, log_offset_round, stan_data, 
                                                  eligible_count_smooth, eligible_count_round, 
                                                  treatment_cascade_samples, proportion_prevalence, participation,
                                                  only_participant = F, art_up_to_female = NULL, s959595 = 1, s909090 = NULL,outdir.table)
  # generate counterfactual treating all men 90 90 90
  counterfactuals_a_909090 <- make_counterfactual(samples, log_offset_round, stan_data, 
                                                  eligible_count_smooth, eligible_count_round, 
                                                  treatment_cascade_samples, proportion_prevalence, participation,
                                                  only_participant = F, art_up_to_female = NULL, s959595 = NULL, s909090 = 1, outdir.table)
  # checks
  if(0){
    tmp <- counterfactuals_a_f$eligible_count_round.counterfactual[ROUND == 'R018' & SEX == 'F' & COMM == 'inland', .(ROUND, COMM, AGEYRS, PROP_UNSUPPRESSED_PARTICIPANTS_M, PROP_UNSUPPRESSED_NONPARTICIPANTS_M)]
    tmp1 <- counterfactuals_a_f$eligible_count_round.counterfactual[ROUND == 'R018' & SEX == 'M' & COMM == 'inland', .(ROUND, COMM, AGEYRS, INFECTED,  PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL)]
    tmp <- merge(tmp, tmp1, by = c('ROUND', 'COMM', 'AGEYRS'))
    tmp[PROP_UNSUPPRESSED_PARTICIPANTS_M != PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL]
    tmp[PROP_UNSUPPRESSED_NONPARTICIPANTS_M != PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL]
    
    tmp <- counterfactuals_a_f$eligible_count_round.counterfactual[ROUND == 'R018' & SEX == 'M' & COMM == 'inland', .(ROUND, COMM, AGEYRS, TREATED)]
    tmp[, TREATED_HALF := TREATED / 2]
    tmp1 <- counterfactuals_a_f05$eligible_count_round.counterfactual[ROUND == 'R018' & SEX == 'M' & COMM == 'inland', .(ROUND, COMM, AGEYRS, TREATED)]
    tmp <- merge(select(tmp, -'TREATED'), tmp1, by = c('ROUND', 'COMM', 'AGEYRS'))
    stopifnot(nrow(tmp[abs(TREATED - TREATED_HALF) > 1e-13]) == 0)
    
    tmp1 <- counterfactuals_a_959595$eligible_count_round.counterfactual[ROUND == 'R018' & SEX == 'M' & COMM == 'inland', .(ROUND, COMM, AGEYRS, INFECTED,  INFECTED_NON_SUPPRESSED)]
    tmp1[, RATIO := INFECTED_NON_SUPPRESSED / INFECTED]
    stopifnot(nrow(tmp1[abs(RATIO - (1-0.95*0.95*0.95)) > 1e-16]) == 0)
  }
  
  # plot
  plot_counterfactual(counterfactuals_a_f, counterfactuals_a_f05, counterfactuals_a_959595, counterfactuals_a_909090, 
                      incidence_factual, "Unsuppressed", outfile.figures)
  
  # save
  save_counterfactual_results(counterfactuals_a_f, counterfactuals_a_f05, 
                              counterfactuals_a_959595, counterfactuals_a_909090, "Unsuppressed", outdir.table)
}

#
# Find NNT
#

cat("\nPlot NNT\n")

# log offset formula (per year per unsuppressed)
log_offset_formula_perunsuppressed <- 'log_SUSCEPTIBLE'
if(!use_number_susceptible_offset)
  log_offset_formula_perunsuppressed = 'log_PROP_SUSCEPTIBLE'
if(use_contact_rates_prior)
  log_offset_formula_perunsuppressed = paste0(log_offset_formula_perunsuppressed, ' + log_CONTACT_RATES')

# NNT by 1 year age band
NNT <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'),
                                   transform = 'exp',
                                   log_offset_round = log_offset_round,
                                   log_offset_formula = log_offset_formula_perunsuppressed,
                                   invert = T)
plot_NNT(NNT, outfile.figures)

# NNT by age groups
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-29', '30-34', '35-39', '40-49'))
NNT_grouped <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'),
                                           transform = 'exp',
                                           log_offset_round = log_offset_round,
                                           log_offset_formula = log_offset_formula_perunsuppressed,
                                           invert = T)
plot_NNT_group(NNT_grouped, outfile.figures)


cat("End of postprocessing_figures.R")
