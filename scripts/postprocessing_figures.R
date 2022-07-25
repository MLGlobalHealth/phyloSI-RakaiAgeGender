cat("Start of postprocessing_figures.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)

jobname <- 'secondrun'
stan_model <- 'gp_220720'

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
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))


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
plot_PPC_augmented_recipient(predict_z_recipient, incidence_cases_recipient, outfile.figures)


#
# Total infected
#

plot_observed_to_augmented(predict_y_source, predict_z_source, outfile.figures)


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


#
# Standardised force of infection 
#

force_infection_standardised <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'), df_direction, df_community, df_period, df_age, transform = 'exp', standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'))
plot_force_infection(force_infection_standardised, outfile.figures, 'standardised')

# shift in sex-specific transmission dynamics
profile_sex_source <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'), df_direction, df_community, df_period, df_age, transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME'))
plot_force_infection_sex_source(profile_sex_source, outfile.figures, 'Standardised force of infection')

# shift in age-specific transmission dynamics
profile_age_source <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'), df_direction, df_community, df_period, df_age, transform = 'exp', standardised.vars = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'))
plot_force_infection_age_source(profile_age_source, outfile.figures, 'Standardised force of infection')

profile_age_source[, is_among_5 := M %in% sort(M, decreasing = T)[1:5], by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
tmp <- profile_age_source[is_among_5 == T, list(total_M = paste0(round(sum(M)*100, 2), '%'), age_group = paste0(min(AGE_TRANSMISSION.SOURCE), '-', max(AGE_TRANSMISSION.SOURCE))), by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
print(tmp[order(LABEL_DIRECTION, PERIOD)])

# aggregated by agr group
profile_age_group_source_recipient <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_GROUP_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age, transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME'))
plot_force_infection_age_group(profile_age_group_source_recipient, outfile.figures, 'Standardised force of infection')

# aggregated by agr group and classified
profile_age_classification_source_recipient <-  find_summary_output(samples, 'log_beta',c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_CLASSIFICATION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'), df_direction, df_community, df_period, df_age, transform = 'exp', standardised.vars = c('INDEX_COMMUNITY', 'INDEX_TIME'))
plot_force_infection_age_classification(profile_age_classification_source_recipient, outfile.figures, 'Standardised force of infection')


#
# median age of source
#
median_age_source <- find_median_age_source(samples, 'log_beta', df_age, df_direction, df_community, df_period)
plot_median_age_source(median_age_source, outfile.figures)


cat("End of postprocessing_figures.R")