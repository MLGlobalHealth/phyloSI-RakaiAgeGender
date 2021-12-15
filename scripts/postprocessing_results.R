cat("Start of postprocessing_results.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)

lab <- "MRC_FALSE_OnlyHTX_TRUE_threshold_0.5"
.stan_model <- 'gp_211207'
DEBUG <- F
.JOBID = 21114

.indir <- "/rds/general/user/mm3218/home/git/phyloflows"
datadir <- "/rds/general/user/mm3218/home/projects/2021/phyloflows"
.outdir <- file.path(datadir, lab, paste0(.stan_model,'-', .JOBID))

if(0){
  .indir <- '~/git/phyloflows'
  datadir <- '~/Box\ Sync/2021/phyloflows/'
  .outdir <- file.path(datadir, lab)
}

datadir <- file.path(datadir, lab)

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-datadir')
  stopifnot(args_line[[5]]=='-outdir')
  stopifnot(args_line[[7]]=='-stan_model')
  stopifnot(args_line[[9]]=='-JOBID')
  stopifnot(args_line[[11]]=='-lab')
  .indir <- args_line[[2]]
  datadir <- args_line[[4]]
  .outdir <- args_line[[6]]
  .stan_model <- args_line[[8]]
  .JOBID <- args_line[[10]]
  lab <- args_line[[12]]
}

# load functions
source(file.path(.indir, 'functions', 'postprocessing_summary_functions.R'))
source(file.path(.indir, 'functions', 'postprocessing_plot_functions.R'))

# load data
path.to.stan.data <- file.path(datadir, paste0("stanin_",lab,".RData"))
load(path.to.stan.data)

# paths
path.to.stan.output = file.path(datadir, paste0(.stan_model,'-', .JOBID), paste0(.stan_model,'-', .JOBID, '_', lab, '.rds'))
outdir.fig <- file.path(.outdir, 'figures', paste0(.stan_model,'-', .JOBID))
outdir.table <- file.path(.outdir, 'tables', paste0(.stan_model,'-', .JOBID))

# table
create.table.reference(stan_data, df_age_time)
range_age_observed <- find_range_age_observed(copy(pairs.all), df_direction)

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)

# convergence diagnostics 
make_convergence_diagnostics_stats(fit, outdir.table)

# intensity of the poisson process
intensity_PP <- summarise_var_by_agextime_direction(samples, 'log_lambda', df_direction, df_age_time, transform = 'exp')
plot_intensity_PP(intensity_PP, outdir.fig)

# intensity of the poisson process reduced
intensity_PP_reduced <- summarise_var_by_agextime_direction(samples, 'log_lambda_reduced', df_direction, df_age_time_reduced, transform = 'exp')
count_data <- prepare_count_data(stan_data, df_age_time_reduced)
plot_intensity_reduced_PP(intensity_PP_reduced, count_data, outdir.fig)

# median age of source
age_source <- find_age_source_by_agextime_direction(samples, df_direction, df_age_time)
age_source <- age_source[!date_infection_evaluated.RECIPIENT %in% range(date_infection_evaluated.RECIPIENT)]
plot_mean_age_source(age_source, range_age_observed, outdir.fig)



range(pairs.all$age_infection.RECIPIENT)
age_source_overall <- find_age_source_by_time_direction(samples, df_direction, df_age_time)
age_source_overall <- age_source_overall[!date_infection_evaluated.RECIPIENT %in% range(date_infection_evaluated.RECIPIENT)]
plot_mean_age_source_overall(age_source_overall, outdir.fig)

cat("End of postprocessing_results.R")

