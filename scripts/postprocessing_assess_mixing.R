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
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_summary_functions.R'))

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


#
## convergence diagnostics 
#

make_convergence_diagnostics_stats(fit, outdir.table)

#
# Trace plots
#
p <- bayesplot::mcmc_trace(fit, regex_pars = c('beta_baseline', 'beta_community', 'beta_period', 'rho_gp', 'alpha_gp')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc_trace_plots.png'), w  = 8, h = 8)

#
# Pairs plot
#
p <- bayesplot::mcmc_pairs(fit, regex_pars = c('beta_baseline', 'beta_community', 'beta_period', 'rho_gp', 'alpha_gp')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc_pairs_plots.png'), w  = 8, h = 8)

#
# Interval plot
#
p <- bayesplot::mcmc_intervals(fit, regex_pars = c('beta_baseline', 'beta_community', 'beta_period', 'rho_gp', 'alpha_gp')) + theme_bw()
ggsave(p, file = paste0(outfile.figures, '-mcmc_intervals_plots.png'), w  = 8, h = 8)



cat("End of postprocessing_assess_mixing.R")

