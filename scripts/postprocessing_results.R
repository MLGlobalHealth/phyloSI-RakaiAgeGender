cat("Start of postprocessing_results.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)

jobname <- '2014_IpriorGP'
stan_model <- 'gp_220108'
DEBUG <- F

indir <- "/rds/general/user/mm3218/home/git/phyloflows"
outdir <- paste0("/rds/general/user/mm3218/home/projects/2021/phyloflows/", stan_model, '-', jobname)

if(0)
{
  indir <- '~/git/phyloflows'
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

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))

# paths
path.to.stan.output = paste0(outfile, "-stanout_", jobname, ".rds")
outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))

# load data
path.to.stan.data <- paste0(outfile, "-stanin_",jobname,".RData")
load(path.to.stan.data)

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)

## convergence diagnostics 
make_convergence_diagnostics_stats(fit, outdir.table)
range_age_observed <- find_range_age_observed(copy(pairs), df_group)
age_aggregated <- c('5-14', '15-29', '30-39', '40-49')
df_age_aggregated <- data.table(expand.grid(age_group_infection.RECIPIENT = age_aggregated, age_group_transmission.SOURCE = age_aggregated))
df_age_aggregated[, age_from.RECIPIENT := gsub('(.+)-.*', '\\1', age_group_infection.RECIPIENT)]
df_age_aggregated[, age_from.SOURCE := gsub('(.+)-.*', '\\1', age_group_transmission.SOURCE)]
df_age_aggregated[, age_to.RECIPIENT := gsub('.*-(.+)', '\\1', age_group_infection.RECIPIENT)]
df_age_aggregated[, age_to.SOURCE := gsub('.*-(.+)', '\\1', age_group_transmission.SOURCE)]
tmp <- df_age_aggregated[, list(age_infection.RECIPIENT = age_from.RECIPIENT:age_to.RECIPIENT), by = c('age_group_infection.RECIPIENT')]
tmp1 <- df_age_aggregated[, list(age_transmission.SOURCE = age_from.SOURCE:age_to.SOURCE), by = c('age_group_transmission.SOURCE')]
df_age_aggregated <- merge(df_age_aggregated, tmp, by = 'age_group_infection.RECIPIENT', allow.cartesian=TRUE)
df_age_aggregated <- merge(df_age_aggregated, tmp1, by = 'age_group_transmission.SOURCE', allow.cartesian=TRUE)

## intensity of the poisson process
cat("\nPlot transmission intensity\n")
intensity_PP <- summarise_var_by_age_group(samples, 'log_lambda', df_group, df_age, transform = 'exp')
count_data <- prepare_count_data(stan_data, df_age, df_group)
plot_intensity_PP(intensity_PP, count_data, outfile.figures)

## relative intensity of the poisson process
relative_intensity_PP <- find_relative_intensity_PP(samples, df_group, df_age)
plot_relative_intensity_PP(standardised_intensity_PP, outfile.figures)
relative_intensity_PP_aggregated <- find_relative_intensity_PP_aggregated(samples, df_group, df_age)

## number of incident cases
cat("\nPlot incident cases\n")
incident_cases <- find_incident_cases(samples, df_group, df_age, incidence)
plot_incident_cases(incident_cases, outfile.figures)
relative_incident_cases_aggregated <- find_relative_incident_cases_aggregated(samples, df_group, df_age, incidence)
plot_relative_intensity_PP_standardised(relative_intensity_PP_aggregated, relative_incident_cases_aggregated, df_age_aggregated, outfile.figures)


## shift in age-specific transmission dynamics
cat("\nPlot age-specific transmission dynamics\n")

# median age of source
age_source <- find_age_source_by_age_group(samples, df_group, df_age)
plot_median_age_source(age_source, outfile.figures)

age_source_difference <- find_age_source_difference_by_age_group(samples, df_group, df_age)
plot_median_age_source_difference(age_source_difference, outfile.figures)

age_source_overall <- find_age_source_by_group(samples, df_group, df_age, incidence, range_age_observed)
plot_median_age_source_overall(age_source_overall, outfile.figures)

# median age of recipient
age_recipient <- find_age_recipient_by_age_group(samples, df_group, df_age)
plot_median_age_recipient(age_recipient, outfile.figures)

age_recipient_difference <- find_age_recipient_difference_by_age_group(samples, df_group, df_age)
plot_median_age_recipient_difference(age_recipient_difference, outfile.figures)

age_recipient_overall <- find_age_recipient_by_group(samples, df_group, df_age, incidence, range_age_observed)
plot_median_age_recipient_overall(age_recipient_overall, outfile.figures)


cat("End of postprocessing_results.R")

