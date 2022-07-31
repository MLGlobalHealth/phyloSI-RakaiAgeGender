cat("Start of postprocessing_figures.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)

jobname <- 'firstrun'
stan_model <- 'gp_220721'

indir <- "/rds/general/user/mm3218/home/git/phyloflows"
outdir <- paste0("/rds/general/user/mm3218/home/projects/2021/phyloflows/", stan_model, '-', jobname)

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

# temporary 
source(file.path(indir, 'functions', 'summary_functions.R'))
df_round <- make.df.round(df_round, df_period)


#
## PPC
#

intensity_PP <- find_summary_output(samples, 'log_lambda', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'), df_direction, df_community, df_period, df_age, transform = 'exp')
count_data <- prepare_count_data(stan_data, df_direction, df_community, df_period, df_age)
plot_intensity_PP(intensity_PP, count_data, outfile.figures)

predict_y_source <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'), df_direction, df_community, df_period, df_age)
predict_y_recipient <- find_summary_output(samples, 'y_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age)
plot_PPC_observed_source(predict_y_source, count_data, outfile.figures)
plot_PPC_observed_recipient(predict_y_recipient, count_data, outfile.figures)

predict_z_source <- find_summary_output(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'), df_direction, df_community, df_period, df_age)
predict_z_recipient <- find_summary_output(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age)
incidence_cases_recipient <- prepare_incidence_cases(incidence_cases)
susceptible_recipient <- prepare_susceptible_count(eligible_count)
plot_PPC_augmented_recipient(predict_z_recipient, incidence_cases_recipient, susceptible_recipient, outfile.figures)


# find log offset by round 
log_offset_round <- find_log_offset_by_round(stan_data, eligible_count_round, df_age, df_direction, df_community, df_period)
predict_z_recipient_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'), 
                                                           df_direction, df_community, df_period, df_age, 
                                                           transform = 'exp', operation = function(x) rpois(1, sum(x)),
                                                           log_offset_round = log_offset_round)
incidence_cases_recipient_round <- prepare_incidence_cases(incidence_cases_round)
susceptible_recipient_count <- prepare_susceptible_count(eligible_count_round)
plot_PPC_augmented_recipient_round(predict_z_recipient_round, incidence_cases_recipient_round, susceptible_recipient_count, outfile.figures)

#
# Total infected
#

unsuppressed_count <- prepare_unsuppressed(eligible_count)
plot_observed_to_augmented(predict_y_source, predict_z_source, unsuppressed_count, outfile.figures)


#
## force of infection
#

force_infection <- find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'), df_direction, df_community, df_period, df_age, transform = 'exp')
plot_force_infection(force_infection, outfile.figures)

# shift in sex-specific transmission dynamics
force_infection_sex_source <- find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'), df_direction, df_community, df_period, df_age, transform = 'exp')
plot_force_infection_sex_source(force_infection_sex_source, outfile.figures)

# shift in age-specific transmission dynamics
force_infection_age_source <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'), df_direction, df_community, df_period, df_age, transform = 'exp')
plot_force_infection_age_source(force_infection_age_source, outfile.figures)

force_infection_age_source[, is_among_5 := M %in% sort(M, decreasing = T)[1:5], by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
tmp <- force_infection_age_source[is_among_5 == T, list(total_M = paste0(round(sum(M)*100, 2), '%'), age_group = paste0(min(AGE_TRANSMISSION.SOURCE), '-', max(AGE_TRANSMISSION.SOURCE))), by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
print(tmp[order(LABEL_DIRECTION, PERIOD)])

# aggregated by agr group
force_infection_aggregated_age_group <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age, transform = 'exp')
plot_force_infection_age_group(force_infection_aggregated_age_group, outfile.figures)
  
# aggregated by agr group and classified
force_infection_aggregated_age_classification <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age, transform = 'exp')
plot_force_infection_age_classification(force_infection_aggregated_age_classification, outfile.figures)

# compare to empirical
force_infection_age_recipient <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age, transform = 'exp')
force_infection_age_recipient_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'), 
                                                          df_direction, df_community, df_period, df_age, transform = 'exp')
plot_force_infection_age_recipient(force_infection_age_recipient, crude_force_infection_age_recipient, outfile.figures)
plot_force_infection_age_recipient_by_round(force_infection_age_recipient_round, crude_force_infection_age_recipient_round, outfile.figures)

crude_force_infection_age_recipient_round[, INDEX_TIME := 0]
crude_force_infection_age_recipient_round[ROUND == 'R015', INDEX_TIME := 1]
crude_force_infection_age_recipient_round[ROUND == 'R016', INDEX_TIME := 2]
crude_force_infection_age_recipient_round[, type := 'byround']
crude_force_infection_age_recipient[, type := 'byperiod']
tmp <- rbind(crude_force_infection_age_recipient_round, crude_force_infection_age_recipient, fill=TRUE)

ggplot(tmp[COMM == 'inland'], aes(x = AGE_INFECTION.RECIPIENT)) +
  geom_line(aes(y = INCIDENT_CASES/PERIOD_SPAN, col = type)) +
  labs(x = 'Age', y = 'Force of infection received', fill = '') +
  theme_bw() +
  facet_grid(INDEX_TIME~LABEL_DIRECTION) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)),
        legend.position = 'bottom') +
  ggsci::scale_fill_npg()


tmp[type == 'byround' & INDEX_TIME == 2]
#
# Contribution to transmission
#

# sex-specific contribution to transmission
contribution_sex_source <-  find_summary_output(samples, 'z_predict', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'), df_direction, df_community, df_period, df_age, standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME'))
eligible_prop_sex <- prepare_eligible_proportion(eligible_count, c('INDEX_TIME', 'COMM', 'SEX'), c('INDEX_TIME', 'COMM'))
plot_contribution_sex_source(contribution_sex_source, eligible_prop_sex, outfile.figures)

# age-specific contribution to transmission
contribution_age_source <-  find_summary_output(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'), df_direction, df_community, df_period, df_age, standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'))
eligible_prop_age <- prepare_eligible_proportion(eligible_count, c('INDEX_TIME', 'COMM', 'SEX', 'AGEYRS'), c('INDEX_TIME', 'COMM', 'SEX'))
plot_contribution_age_source(contribution_age_source, eligible_prop_age, outfile.figures)

contribution_age_source[, is_among_5 := M %in% sort(M, decreasing = T)[1:5], by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
tmp <- contribution_age_source[is_among_5 == T, list(total_M = paste0(round(sum(M)*100, 2), '%'), age_group = paste0(min(AGE_TRANSMISSION.SOURCE), '-', max(AGE_TRANSMISSION.SOURCE))), by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
print(tmp[order(LABEL_DIRECTION, PERIOD)])

# aggregated by agr group
contribution_age_group_source <-  find_summary_output(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age, standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME'))
plot_contribution_age_group(contribution_age_group_source, outfile.figures)

# aggregated by agr group and classified
contribution_age_classification_source <-  find_summary_output(samples, 'z_predict',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age, standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME'))
plot_contribution_age_classification(contribution_age_classification_source, outfile.figures)


#
# Expected Contribution to transmission
#

# sex-specific contribution to transmission
expected_contribution_sex_source <- find_summary_output(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'), df_direction, df_community, df_period, df_age,transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME'))
plot_contribution_sex_source(expected_contribution_sex_source, eligible_prop_sex, outfile.figures,'Expected contribution to infection')

# age-specific contribution to transmission
expected_contribution_age_source <-  find_summary_output(samples, 'log_lambda_latent',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'), df_direction, df_community, df_period, df_age, transform = 'exp', standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'))
plot_contribution_age_source(expected_contribution_age_source, eligible_prop_age, outfile.figures,'Expected contribution to infection')

# aggregated by agr group
expected_contribution_age_group_source <-  find_summary_output(samples, 'log_lambda_latent',
                                                               c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), 
                                                               df_direction, df_community, df_period, df_age,transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_DIRECTION', 'AGE_GROUP_INFECTION.RECIPIENT'))
plot_contribution_age_group(expected_contribution_age_group_source, outfile.figures,'Expected contribution to infection')

# aggregated by agr group and classified
expected_contribution_age_classification_source <-  find_summary_output(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), 
                                                                        df_direction, df_community, df_period, df_age, transform = 'exp',standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_DIRECTION', 'AGE_GROUP_INFECTION.RECIPIENT'))
plot_contribution_age_classification(expected_contribution_age_classification_source, outfile.figures,'Expected contribution to infection')


#
# Expected Contribution to transmission by round
#

# sex-specific contribution to transmission
expected_contribution_sex_source_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND'), 
                                                                       df_direction, df_community, df_period, df_age, 
                                                                       transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'ROUND'), 
                                                                       log_offset_round = log_offset_round)
unsuppressed_prop_sex <- prepare_unsuppressed_proportion_by_round(eligible_count_round, c('ROUND', 'COMM', 'SEX'), c('ROUND', 'COMM'))
plot_contribution_sex_source_by_round(expected_contribution_sex_source_round, unsuppressed_prop_sex, outfile.figures, 'Expected contribution to infection')

# age-specific contribution to transmission
expected_contribution_age_source_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_TRANSMISSION.SOURCE'), 
                             df_direction, df_community, df_period, df_age, 
                             transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'ROUND', 'INDEX_DIRECTION'), 
                             log_offset_round = log_offset_round)
unsuppressed_prop_age <- prepare_unsuppressed_proportion_by_round(eligible_count_round, c('ROUND', 'COMM', 'SEX', 'AGEYRS'), c('ROUND', 'COMM', 'SEX'))
plot_contribution_age_source_by_round(expected_contribution_age_source_round, unsuppressed_prop_age, outfile.figures,'Expected contribution to infection')


# aggregated by agr group
expected_contribution_age_group_source_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), 
                                                                             df_direction, df_community, df_period, df_age, 
                                                                             transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'ROUND'), 
                                                                             log_offset_round = log_offset_round)
plot_contribution_age_group_by_round(expected_contribution_age_group_source_round, outfile.figures,'Expected contribution to infection')

# aggregated by age group and classified
expected_contribution_age_classification_source_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), 
                             df_direction, df_community, df_period, df_age, 
                             transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'ROUND'), 
                             log_offset_round = log_offset_round)
plot_contribution_age_classification_by_round(expected_contribution_age_classification_source_round, outfile.figures,'Expected contribution to infection')


#
# Expected transmission risk
#

# by sex
transmission_risk_sex_source_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND'), 
                                                                       df_direction, df_community, df_period, df_age, 
                                                                       transform = 'exp', 
                                                                       log_offset_round = log_offset_round, log_offset_name = 'log_INFECTED_NON_SUPPRESSED')
plot_transmission_risk_sex_source_by_round(transmission_risk_sex_source_round, outfile.figures)

# by age
transmission_risk_round <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'INDEX_AGE'), 
                                                        df_direction, df_community, df_period, df_age, 
                                                        transform = 'exp', 
                                                        log_offset_round = log_offset_round, log_offset_name = 'log_INFECTED_NON_SUPPRESSED')
plot_transmission_risk_age_source_by_round(transmission_risk_round, outfile.figures)


#
# median age of source
#

median_age_source <- find_median_age_source(samples, 'log_beta', df_age, df_direction, df_community, df_period)
plot_median_age_source(median_age_source, outfile.figures)


cat("End of postprocessing_figures.R")