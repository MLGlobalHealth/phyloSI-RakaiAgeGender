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

jobname <- 'notsi'
stan_model <- 'gp_221002a'

indir <- "/rds/general/user/mm3218/home/git/phyloflows"
outdir <- paste0("/rds/general/user/mm3218/home/projects/2021/phyloflows/", stan_model, '-', jobname)

if(0){
  indir <- '~/git/phyloflows/'
  outdir <- paste0('~/Box\ Sync/2021/phyloflows/', stan_model, '-', jobname)
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

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_summary_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_plot_functions.R'))

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

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)

# temp
source(file.path(indir, 'functions', 'postprocessing_summary_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_plot_functions.R'))

#
# offset
#

log_offset_round <- find_log_offset_by_round(stan_data, eligible_count_round)



#
## PPC
#

cat("\nPlot PPC\n")

intensity_PP_sampled <- find_summary_output(samples, 'log_lambda', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'), transform = 'exp')
count_data <- prepare_count_data(stan_data)
plot_intensity_PP(intensity_PP_sampled, count_data, outfile.figures)

intensity_PP <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'INDEX_AGE'),
                                             transform = 'exp',
                                             log_offset_round = log_offset_round,
                                             log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')
plot_intensity_PP_by_round(intensity_PP, outfile.figures)

predict_y_source <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'))
predict_y_recipient <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'))
plot_PPC_observed_source(predict_y_source, count_data, outfile.figures)
plot_PPC_observed_recipient(predict_y_recipient, count_data, outfile.figures)

predict_z_source <- find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'))
predict_z_source_round <- find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'))
predict_z_recipient_round <- find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'))
predict_lambda_recipient_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                               transform = 'exp',
                                                               log_offset_round = log_offset_round,
                                                               log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')
incidence_cases_recipient_round <- prepare_incidence_cases(incidence_cases_round)
eligible_count_recipient <- prepare_eligible_count(eligible_count_round)
plot_PPC_augmented_recipient_round(predict_z_recipient_round, incidence_cases_recipient_round,
                                   eligible_count_recipient, outfile.figures)
plot_PPC_incidence_rate_round(predict_lambda_recipient_round, incidence_cases_recipient_round,
                                   eligible_count_recipient, outfile.figures)

unsuppressed_count <- prepare_unsuppressed(eligible_count)
plot_observed_to_augmented(predict_y_source, predict_z_source, unsuppressed_count, outfile.figures)


#
## force of infection
#

cat("\nPlot force of infection\n")

# 2D for all categories
force_infection <- find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'INDEX_AGE'), transform = 'exp')
plot_force_infection(force_infection, outfile.figures)

# shift in sex-specific transmission dynamics by period
force_infection_sex_source <- find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'), transform = 'exp')
plot_force_infection_sex_source(force_infection_sex_source, outfile.figures)

# shift in age source by period
force_infection_age_source <-  find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'), transform = 'exp')
plot_force_infection_age_source(force_infection_age_source, outfile.figures)

# shift in age source by round
force_infection_age_recipient <-  find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), transform = 'exp')
plot_force_infection_age_recipient(force_infection_age_recipient, outfile.figures)



#
# Contribution to transmission
#

cat("\nPlot contribution\n")

# sex-specific contribution to transmission
contribution_sex_source <-  find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                         standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
unsuppressed_prop_sex <- prepare_unsuppressed_proportion_by_round(file.unsuppressed.share, c('SEX'))
prevalence_prop_sex<- prepare_prevalence_proportion_by_round(file.prevalence.share, 'SEX')
plot_contribution_sex_source(contribution_sex_source, unsuppressed_prop_sex, prevalence_prop_sex, outfile.figures)

# age-specific contribution to transmission among all sources
contribution_age_source <-  find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                         standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
unsuppressed_prop_age <- prepare_unsuppressed_proportion_by_round(file.unsuppressed.share, c('SEX', 'AGEYRS'))
plot_contribution_age_source_unsuppressed(contribution_age_source, unsuppressed_prop_age, outfile.figures)
plot_contribution_age_source(contribution_age_source, outfile.figures)

# aggregated by agr group
contribution_age_group_source <-  find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                               standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_group(contribution_age_group_source, outfile.figures)

# aggregated by agr group and classified
contribution_age_classification_source <-  find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_classification(contribution_age_classification_source, outfile.figures)


#
# Expected Contribution to transmission
#

cat("\nPlot expected contribution\n")

# sex-specific contribution to transmission
expected_contribution_sex_source <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                        transform = 'exp',
                                                        standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_sex_source(expected_contribution_sex_source, unsuppressed_prop_sex, prevalence_prop_sex, outfile.figures,'Expected_contribution')

# age-specific contribution to transmission across sex
expected_contribution_age_source2 <-  find_summary_output_by_round(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                  transform = 'exp',
                                                                  standardised.vars = c( 'INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_source_unsuppressed(expected_contribution_age_source2, unsuppressed_prop_age, outfile.figures,'Expected_contribution')
plot_contribution_age_source(expected_contribution_age_source2, outfile.figures,'Expected_contribution_Sex')

# age-specific sex ration contribution to transmission
expected_contribution_age_source_sex_ratio <-  find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                   c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                   transform = 'exp',
                                                                   standardised.vars = c( 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                                   sex_ratio= T)
plot_contribution_age_source_sex_ratio(expected_contribution_age_source_sex_ratio, outfile.figures,'Expected_contribution_Sex_Ratio')

# aggregated by agr group
expected_contribution_age_group_source <-  find_summary_output_by_round(samples, 'log_lambda_latent',
                                                               c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                               transform = 'exp',
                                                               standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_group(expected_contribution_age_group_source, outfile.figures,'Expected_contribution')

# aggregated by agr group and classified
expected_contribution_age_classification_source <-  find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_classification(expected_contribution_age_classification_source, outfile.figures,'Expected_contribution')



#
# Transmission risk per unsuppressed
#

cat("\nPlot transmission risk\n")


# sex-specific transmission risk
transmission_risk_sex_source <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                             transform = 'exp',
                                                             log_offset_round = log_offset_round,
                                                             log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED',
                                                             per_unsuppressed = T)
plot_transmission_risk_sex_source(transmission_risk_sex_source, outfile.figures)

# age-specific  transmission risk
transmission_risk_age_source<- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                            transform = 'exp',
                                                            log_offset_round = log_offset_round,
                                                            log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED',
                                                            per_unsuppressed = T)
plot_transmission_risk_age_source(transmission_risk_age_source, outfile.figures)


#
# Incidence infection
#

cat("\nPlot incidence infection and transmission\n")

# finer age bands
# df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-29', '30-34', '35-39', '40-49'))
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))

#find incidence transmission
incidence_tranmission <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE'),
                                                      transform = 'exp',
                                                      log_offset_round = log_offset_round,
                                                      log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED',
                                                      relative_baseline = T,
                                                      per_eligible = T)
plot_incidence_transmission(incidence_tranmission, outfile.figures)

#find incidence infection
incidence_infection <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                    transform = 'exp',
                                                    log_offset_round = log_offset_round,
                                                    log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED',
                                                    relative_baseline = T,
                                                    per_eligible = T)
plot_incidence_infection(incidence_infection, outfile.figures)


#
# median age of source
#

cat("\nPlot median age at transmission of the source by age at infection of recipient\n")

# by 1-year age band
median_age_source <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                  median_age_source = T)
plot_median_age_source(median_age_source, outfile.figures)

# by ge groups
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
median_age_source_group <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                  quantile_age_source = T)
expected_contribution_age_group_source2 <-  find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
reported_contact <- clean_reported_contact(df_reported_contact)
plot_median_age_source_group(median_age_source_group, expected_contribution_age_group_source2, reported_contact, outfile.figures)

#
# Relative incidence infection if male had the same art uptake as female
#

cat("\nPlot relative incidence infection if male had the same art uptake as female\n")

# incidence actual 
incidence_factual <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  log_offset_round = log_offset_round,
                                                  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')
incidence_factual_all <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                  transform = 'exp',
                                                  log_offset_round = log_offset_round,
                                                  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')
# find age groups that contribute the most 
expected_contribution_age_source <- find_summary_output_by_round(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'), 
                                                                 transform = 'exp', 
                                                                 standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'))
spreaders <- find_spreaders(expected_contribution_age_source, outdir.table)

# counterfactual participant
## treated as much as female
counterfactuals_p_f <- make_counterfactual(samples, spreaders, log_offset_round, stan_data,
                                       eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                       only_participant = T, art_up_to_female = T, outdir.table)

## all treated
counterfactuals_p_a <- make_counterfactual(samples, spreaders, log_offset_round, stan_data,
                                           eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                           only_participant = T, art_up_to_female = F, outdir.table)

# plot
plot_counterfactual(counterfactuals_p_f, counterfactuals_p_a, eligible_count_round, incidence_factual, "Diagnosed unsuppressed", outfile.figures)
plot_counterfactual_all(counterfactuals_p_f, counterfactuals_p_a, eligible_count_round, incidence_factual_all, "Diagnosed unsuppressed", outfile.figures)


# counterfactual all males
## treated as much as female
counterfactuals_a_f <- make_counterfactual(samples, spreaders, log_offset_round, stan_data,
                                            eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                            only_participant = F, art_up_to_female = T, outdir.table)

## all treated
counterfactuals_a_a <- make_counterfactual(samples, spreaders, log_offset_round, stan_data,
                                           eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                           only_participant = F, art_up_to_female = F, outdir.table)

# plot
plot_counterfactual(counterfactuals_a_f, counterfactuals_a_a, eligible_count_round, incidence_factual, "Unsuppressed", outfile.figures)
plot_counterfactual_all(counterfactuals_a_f, counterfactuals_a_a, eligible_count_round, incidence_factual_all, "Unsuppressed", outfile.figures)


#
# Find NNT
#
cat("\nPlot NNT\n")
# NNT by 1 year age band
NNT <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'),
                                   transform = 'exp',
                                   log_offset_round = log_offset_round,
                                   log_offset_formula = 'log_PROP_SUSCEPTIBLE',
                                   invert = T)
plot_NNT(NNT, outfile.figures)

# NNT by age groups
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-29', '30-34', '35-39', '40-49'))
NNT_grouped <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'),
                                           transform = 'exp',
                                           log_offset_round = log_offset_round,
                                           log_offset_formula = 'log_PROP_SUSCEPTIBLE',
                                           invert = T)
plot_NNT_group(NNT_grouped, outfile.figures)


cat("End of postprocessing_figures.R")


