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

# load data
path.to.stan.data <- file.path(outfile, paste0("stanin_",jobname,".RData"))
load(path.to.stan.data)

# paths
path.to.stan.output = paste0(outfile, "-stanout_", jobname, ".rds")
outdir.fig <- file.path(outdir, 'figures', paste0(stan_model,'-', .JOBID))
outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', .JOBID))

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)

## convergence diagnostics 
make_convergence_diagnostics_stats(fit, outdir.table)

range_age_observed <- find_range_age_observed(copy(pairs), df_group)

## intensity of the poisson process
cat("\nPlot transmission intensity\n")
intensity_PP <- summarise_var_by_age_group(samples, 'log_lambda', df_group, df_age, transform = 'exp')
count_data <- prepare_count_data(stan_data, df_age, df_group)
plot_intensity_PP(intensity_PP, count_data, outdir.fig)

## number of incident cases
cat("\nPlot incident cases\n")
incident_cases <- find_incident_cases_by_group(samples, df_group, df_age, incidence, range_age_observed)
plot_incident_cases(incident_cases, outdir.fig)


## shift in age-specific transmission dynamics
cat("\nPlot age-specific transmission dynamics\n")

# median age of source
age_source <- find_age_source_by_age_group(samples, df_group, df_age)
plot_median_age_source(age_source, outdir.fig)

age_source_difference <- find_age_source_difference_by_age_group(samples, df_group, df_age)
plot_median_age_source_difference(age_source_difference, outdir.fig)

age_source_overall <- find_age_source_by_group(samples, df_group, df_age, incidence, range_age_observed)
plot_median_age_source_overall(age_source_overall, outdir.fig)

# median age of recipient
age_recipient <- find_age_recipient_by_age_group(samples, df_group, df_age)
plot_median_age_recipient(age_recipient, outdir.fig)

age_recipient_difference <- find_age_recipient_difference_by_age_group(samples, df_group, df_age)
plot_median_age_recipient_difference(age_recipient_difference, outdir.fig)

age_recipient_overall <- find_age_recipient_by_group(samples, df_group, df_age, incidence, range_age_observed)
plot_median_age_recipient_overall(age_recipient_overall, outdir.fig)


cat("End of postprocessing_results.R")

