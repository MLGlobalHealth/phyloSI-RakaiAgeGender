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

jobname <- 'firstrun'
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
source(file.path(indir, 'functions', 'postprocessing_summary_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_plot_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_utils_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_statistics_functions.R'))

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
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_utils_functions.R'))
unsuppressed_share <- fread(file.unsuppressed.share)
infected_share <- fread(file.prevalence.share)
df_direction <- get.df.direction()
df_period[, PERIOD := factor(PERIOD, levels = PERIOD)]

#
# offset
#

log_offset_round <- find_log_offset_by_round(stan_data, eligible_count_round)


#
# Summarise data and merge to maps for figures
#

count_data <- prepare_count_data(stan_data)
incidence_cases_recipient_round <- prepare_incidence_cases(incidence_cases_round)
unsuppressed_share_sex <- prepare_unsuppressed_share(unsuppressed_share, c('SEX'))
unsuppressed_share_sex_age <- prepare_unsuppressed_share(unsuppressed_share, c('SEX', 'AGEYRS'))
prevalence_prop_sex<- prepare_infected_share(infected_share, 'SEX')
reported_contact <- clean_reported_contact(df_reported_contact)

#
## PPC
#

cat("\nPlot PPC\n")

# rate of observed transmission
intensity_PP_sampled <- find_summary_output(samples, 'log_lambda', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'), transform = 'exp')
plot_intensity_PP(intensity_PP_sampled, count_data, outfile.figures)

# rate of transmission
intensity_PP <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'INDEX_AGE'),
                                             transform = 'exp',
                                             log_offset_round = log_offset_round,
                                             log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')
plot_intensity_PP_by_round(intensity_PP, outfile.figures)

# predicted observed transmission
predict_y_source <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'))
predict_y_recipient <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'))
plot_PPC_observed_source(predict_y_source, count_data, outfile.figures)
plot_PPC_observed_recipient(predict_y_recipient, count_data, outfile.figures)

# predicted transmission
predict_z_source <- find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'))
predict_z_recipient_round <- find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'))
plot_PPC_augmented_recipient_round(predict_z_recipient_round, incidence_cases_recipient_round, outfile.figures)
predict_incidence_rate_round <- find_summary_output_by_round(samples, 'ir_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                             names = c('INDEX_AGE_INFECTION.RECIPIENT', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_PPC_incidence_rate_round(predict_incidence_rate_round, incidence_cases_recipient_round,outfile.figures)

# predicted observed transmission vs all transmission
plot_observed_to_augmented(predict_y_source, predict_z_source, outfile.figures)


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
force_infection_age_source <- find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'), transform = 'exp')
plot_force_infection_sex_age_source(force_infection_age_source, outfile.figures)

# shift in age source by round
force_infection_age_recipient <-  find_summary_output_by_round(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), transform = 'exp')
plot_force_infection_sex_age_recipient(force_infection_age_recipient, outfile.figures)



#
# Contribution to transmission
#

cat("\nPlot contribution\n")

# sex-specific contribution to transmission
contribution_sex_source <-  find_summary_output_by_round(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                         standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_sex_source(contribution_sex_source, unsuppressed_share_sex, prevalence_prop_sex, outfile.figures)

# age-specific contribution to transmission among all sources by sex
contribution_age_source <-  find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                         standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_source_unsuppressed(contribution_age_source, unsuppressed_share_sex_age, outfile.figures)
plot_contribution_age_source(contribution_age_source, outfile.figures)

# contribution aggregated by age group of sources and recipients
contribution_age_group_source <- find_summary_output_by_round(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                               standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_group(contribution_age_group_source, outfile.figures)

# contribution aggregated by age group of recipients and classification of age of sources
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
plot_contribution_sex_source(expected_contribution_sex_source, unsuppressed_share_sex, prevalence_prop_sex, outfile.figures,'Expected_contribution')

# age-specific contribution to transmission among all sources by sex
expected_contribution_age_source2 <- find_summary_output_by_round(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                  transform = 'exp',
                                                                  standardised.vars = c( 'INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_source_unsuppressed(expected_contribution_age_source2, unsuppressed_share_sex_age, outfile.figures,'Expected_contribution')
plot_contribution_age_source(expected_contribution_age_source2, outfile.figures,'Expected_contribution_sex')
save_statistics_expected_contribution(expected_contribution_sex_source, expected_contribution_age_source2, outdir.table)

# age-specific sex ratio contribution to transmission
expected_contribution_age_source_sex_ratio <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                           c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                           transform = 'exp',
                                                                           standardised.vars = c( 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                                           sex_ratio= T)
plot_contribution_age_source_sex_ratio(expected_contribution_age_source_sex_ratio, outfile.figures,'Expected_contribution_sex_ratio')

# contribution aggregated by age group recipients and 1-year age band sources
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
expected_contribution_age_ungroup_source <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                       c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                       transform = 'exp',
                                                                       standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_ungroup(expected_contribution_age_ungroup_source, outfile.figures,'Expected_contribution')

# contribution aggregated by age group of sources and recipients
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))
expected_contribution_age_group_source <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_contribution_age_group(expected_contribution_age_group_source, outfile.figures,'Expected_contribution')

# contribution aggregated by age group of recipients and classification of age of sources
expected_contribution_age_classification_source <- find_summary_output_by_round(samples, 'log_lambda_latent',
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

# sex and age-specific  transmission risk
transmission_risk_age_source <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                            transform = 'exp',
                                                            log_offset_round = log_offset_round,
                                                            log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED',
                                                            per_unsuppressed = T)
plot_transmission_risk_sex_age_source(transmission_risk_age_source, outfile.figures)


#
# Incidence infection
#

cat("\nPlot incidence infection and transmission\n")

# finer age bands
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))

#find incidence transmission
incidence_tranmission <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE'),
                                                      transform = 'exp',
                                                      log_offset_round = log_offset_round,
                                                      log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED',
                                                      relative_baseline = T,
                                                      per_susceptible = T)
plot_incidence_transmission(incidence_tranmission, outfile.figures)

#find incidence infection
incidence_infection <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                    transform = 'exp',
                                                    log_offset_round = log_offset_round,
                                                    log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED',
                                                    relative_baseline = T,
                                                    per_susceptible = T)
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

# by age groups
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
median_age_source_group <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                  quantile_age_source = T)
expected_contribution_age_group_source2 <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_COMMUNITY', 'INDEX_ROUND'))
plot_median_age_source_group(median_age_source_group, expected_contribution_age_group_source2, reported_contact, outfile.figures)


#
# Counterfactual: comparison of targeted male
#

cat("\nPlot relative incidence infection if different groups of male are targeted\n")

# in the counterfactual we assume that 95% of art user are unsuppressed (previously we assumed 100)
# we need to recalculate the number of infected unsuppressed with this new assumption 
# and the log offset with this new assumption 
# and define it as the factual scenario 
eligible_count_round_95suppression_given_ART <- find_eligible_count_round_95suppression_given_ART(eligible_count_smooth, proportion_prevalence, proportion_unsuppressed)
log_offset_round_95suppression_given_ART <- find_log_offset_by_round(stan_data, eligible_count_round_95suppression_given_ART)

# find incidence under the factual scenario by sex and age
incidence_factual <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                  transform = 'exp',
                                                  log_offset_round = log_offset_round_95suppression_given_ART,
                                                  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')

# find incidence under the factual scenario by sex
incidence_factual_all <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                  transform = 'exp',
                                                  log_offset_round = log_offset_round_95suppression_given_ART,
                                                  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')

# find age groups that contribute the most
expected_contribution_age_source <- find_summary_output_by_round(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
                                                                 transform = 'exp',
                                                                 standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'))

# identify the main spreaders, sources that contributes the most to incidence
spreaders <- find_spreaders(expected_contribution_age_source, outdir.table)

# generate counterfactual when only participant are treated 
# we compare treating main spreaders, male with the greatest diff in art coverage compared to female and random male
# in all scenario, the number of participants treated in the counterfactual are = number of participant sources that contribute to 1/3 incidence
counterfactuals_p_a <- make_counterfactual_target(samples, spreaders, log_offset_round_95suppression_given_ART, stan_data,
                                           eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                           only_participant = T, art_up_to_female = F, outdir.table)
plot_counterfactual_one(counterfactuals_p_a, eligible_count_round, incidence_factual, "Diagnosed unsuppressed", outfile.figures)

# generate counterfactual when participants and non-participant are treated
# we compare treating main spreaders, male with the greatest diff in art coverage compared to female and random male
# the number of male treated in the counterfactual are = number of male sources that contribute to 1/3 incidence
counterfactuals_a_a <- make_counterfactual_target(samples, spreaders, log_offset_round_95suppression_given_ART, stan_data,
                                           eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                           only_participant = F, art_up_to_female = F, outdir.table)
plot_counterfactual_one(counterfactuals_a_a, eligible_count_round, incidence_factual, "Unsuppressed", outfile.figures)


#
# Counterfactual: comparison of the number of male treated
#

cat("\nPlot relative incidence infection if different number of male are treated\n")

# in the previous counterfactual, we found that the most effective intervention is to
# target male with greatest art diff compared to female
# we now explore different number of age groups targeted in this scenario

# males with greatest art diff compared to female stratified by categories
# including more and more age groups until reaching all males
targeted.males <- find_male_with_greatest_art_difference_category(eligible_count_round_95suppression_given_ART, outdir.table)

# generate counterfactual treating only participant treated as much as female are treated
counterfactuals_p_f <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                           eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                           only_participant = T, art_up_to_female = 1, s959595 = NULL, outdir.table)
# generate counterfactual treating only participant treated half way to a much as female are treated
counterfactuals_p_f05 <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                             eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                             only_participant = T, art_up_to_female = 0.5, s959595 = NULL, outdir.table)
# generate counterfactual treating only participant 95 95 95
counterfactuals_p_959595 <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                             eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                             only_participant = T, art_up_to_female = NULL, s959595 = 1, outdir.table)
# generate counterfactual treating only participant treated half way to 95 95 95
counterfactuals_p_95959505 <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                             eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                             only_participant = T, art_up_to_female = NULL, s959595 = 0.5, outdir.table)

plot_counterfactual(counterfactuals_p_f, counterfactuals_p_f05, counterfactuals_p_959595, counterfactuals_p_95959505, 
                    eligible_count_round_95suppression_given_ART, incidence_factual, "Diagnosed unsuppressed", outfile.figures)

# generate counterfactual treating all men treated as much as female are treated
counterfactuals_a_f <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                           eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                           only_participant = F, art_up_to_female = 1, s959595 = NULL, outdir.table)
# generate counterfactual treating all men treated half way to a much as female are treated
counterfactuals_a_f05 <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                             eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                             only_participant = F, art_up_to_female = 0.5, s959595 = NULL, outdir.table)
# generate counterfactual treating all men 95 95 95
counterfactuals_a_959595 <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                                eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                                only_participant = F, art_up_to_female = NULL, s959595 = 1, outdir.table)
# generate counterfactual treating all men treated half way to 95 95 95
counterfactuals_a_95959505 <- make_counterfactual(samples, targeted.males, log_offset_round_95suppression_given_ART, stan_data,
                                                  eligible_count_smooth, proportion_unsuppressed, proportion_prevalence,
                                                  only_participant = F, art_up_to_female = NULL, s959595 = 0.5, outdir.table)

plot_counterfactual(counterfactuals_a_f, counterfactuals_a_f05, counterfactuals_a_959595, counterfactuals_a_95959505, 
                    eligible_count_round_95suppression_given_ART, incidence_factual, "Unsuppressed", outfile.figures)


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


