library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)

lab <- "MRC_FALSE_OnlyHTX_TRUE_threshold_0.5"
.stan_model <- 'gp_211207'
DEBUG <- F
JOBID = 12

.indir <- "/rds/general/user/mm3218/home/git/phyloflows"
datadir <- "/rds/general/user/mm3218/home/projects/2021/phyloflows"

if(0){
  .indir <- '~/git/phyloflows'
  datadir <- '~/Box\ Sync/2021/phyloflows/'
}

.outdir <- file.path(datadir, lab)
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

# paths
path.to.stan.data <- file.path(datadir, paste0("stanin_",lab,".RData"))
path.to.stan.output = file.path(.outdir, paste0(.stan_model,'-', .JOBID), paste0(.stan_model,'-', .JOBID, '_', lab, '.rds'))
outdir.fig <- file.path(.outdir, 'figures', paste0(.stan_model,'-', .JOBID))
outdir.table <- file.path(.outdir, 'table', paste0(.stan_model,'-', .JOBID))

# load data
load(path.to.stan.data)

# table
df_direction <- data.table(index_direction = 1:2, is_mf = stan_data$is_mf)
df_direction[, label_direction := ifelse(is_mf == 1, 'Male -> Female', 'Female -> Male')]

df_age_time$index_age_time <- 1:nrow(df_age_time); 
df_age_time[, year_infection_reduced.RECIPIENT := as.numeric(format(date_infection_reduced.RECIPIENT, '%Y'))]
time_bands = unique(diff(unique(df_age_time$year_infection_reduced.RECIPIENT)))
df_age_time[, date_infection_reduced_name.RECIPIENT := paste0('Jan ', year_infection_reduced.RECIPIENT, '-Dec ', year_infection_reduced.RECIPIENT + time_bands - 1)]

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)

# convergence diagnostics 
make_convergence_diagnostics_stats(fit, outdir.table)

# intensity of the poisson process
intensity_PP <- summarise_var_by_agextime_direction(samples, 'log_lambda', df_direction, df_age_time, transform = 'exp')
count_data <- prepare_count_data(stan_data, df_age_time)
plot_intensity_PP(intensity_PP, count_data, outdir.fig)
  
# median age of source
# age_source <- find_age_source(samples, df_direction, df_age)
  



