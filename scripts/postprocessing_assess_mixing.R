cat("Start of postprocessing_assess_mixing.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)
library(bayesplot)

jobname <- 'informativeprior_incrate0524'
stan_model <- 'gp_220427'
DEBUG <- F

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
source(file.path(gitdir.functions, 'functions_transmission_flow', 'postprocessing_summary_functions.R'))
source(file.path(gitdir.functions, 'functions_transmission_flow', 'postprocessing_plot_functions.R'))
source(file.path(gitdir.functions, 'functions_transmission_flow', 'postprocessing_utils_functions.R'))

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


#
## convergence diagnostics 
#

make_convergence_diagnostics_stats(fit, samples, outdir.table)


#
# Trace plots
#
p <- bayesplot::mcmc_trace(fit, regex_pars = c('rho_gp', 'alpha_gp')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc-trace_plots.png'), w  = 8, h = 8)


#
# Interval plot
#

# baseline parameters
p <- bayesplot::mcmc_intervals(fit, pars = c('log_beta_baseline', 'log_beta_baseline_contrast_direction', 
                                             paste0('log_beta_baseline_contrast_round[', 1:(stan_data$N_ROUND-1), ']'),
                                             'log_beta_baseline_contrast_period')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc-intervals_plots-baseline.png'), w  = 8, h = 8)

# hyperparameters
p <- bayesplot::mcmc_intervals(fit, regex_pars = c('rho_gp', 'alpha_gp')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc-intervals_plots.png'), w  = 8, h = 8)


#
# Pairs plot
#

# baseline parameters
p <- bayesplot::mcmc_pairs(fit, pars = c('log_beta_baseline', 'log_beta_baseline_contrast_direction', 
                                         paste0('log_beta_baseline_contrast_round[', 1:(stan_data$N_ROUND-1), ']'),
                                         'log_beta_baseline_contrast_period')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc-pairs_plots-baseline.png'), w  = 8, h = 8)

# hyperparameters
p <- bayesplot::mcmc_pairs(fit, regex_pars = c('rho_gp', 'alpha_gp')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc-pairs_plots.png'), w  = 8, h = 8)


#
# Force of infection contrasts
#

#
# # period contrast
add.vars = NULL
if(length(dim(samples[['log_beta_period_contrast']])) == 3){
  add.vars = c('INDEX_DIRECTION', add.vars)
}

log_period_contrast <- find_summary_output(samples, 'log_beta_period_contrast', c(add.vars, 'INDEX_AGE'), names = c(add.vars, 'INDEX_AGE'))
plot_2D_contrast(log_period_contrast, outfile.figures, paste0(df_period[INDEX_TIME == 2, PERIOD], ' period contrast'), 'period')

log_period_contrast_source <- find_summary_output(samples, 'log_beta_period_contrast', c(add.vars, 'AGE_TRANSMISSION.SOURCE'), 
                                                  names = c(add.vars, 'INDEX_AGE'), operation = 'mean')
plot_source_contrast(log_period_contrast_source, outfile.figures, paste0(df_period[INDEX_TIME == 2, PERIOD], ' period contrast'), 'period')

log_period_contrast_recipient<- find_summary_output(samples, 'log_beta_period_contrast', c(add.vars, 'AGE_INFECTION.RECIPIENT'), 
                                                  names = c(add.vars, 'INDEX_AGE'), operation = 'mean')
plot_recipient_contrast(log_period_contrast_recipient, outfile.figures, paste0(df_period[INDEX_TIME == 2, PERIOD], ' period contrast'), 'period')

#
# # round contrast
add.vars = 'INDEX_ROUND'
if(length(dim(samples[['log_beta_round_contrast']])) == 4){
  add.vars = c('INDEX_DIRECTION', add.vars)
}

log_round_contrast <- find_summary_output_by_round(samples, 'log_beta_round_contrast', c(add.vars, 'INDEX_AGE'),  names = c(add.vars, 'INDEX_AGE'))
log_round_contrast <- remove_first_round(log_round_contrast)
plot_2D_contrast(log_round_contrast, outfile.figures, paste0('Round contrast'), 'round')

log_round_contrast_source <- find_summary_output_by_round(samples, 'log_beta_round_contrast', c(add.vars, 'AGE_TRANSMISSION.SOURCE'), 
                                                  names = c(add.vars, 'INDEX_AGE'), operation = 'mean')
log_round_contrast_source <- remove_first_round(log_round_contrast_source)
plot_source_contrast(log_round_contrast_source, outfile.figures, 'Round constrast', 'round')

log_round_contrast_recipient<- find_summary_output_by_round(samples, 'log_beta_round_contrast', c(add.vars, 'AGE_INFECTION.RECIPIENT'), 
                                                    names = c(add.vars, 'INDEX_AGE'), operation = 'mean')
log_round_contrast_recipient <- remove_first_round(log_round_contrast_recipient)
plot_recipient_contrast(log_round_contrast_recipient, outfile.figures, 'Round constrast', 'round')



cat("End of postprocessing_assess_mixing.R")

