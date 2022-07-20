cat("Start of postprocessing_figures.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)

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

#
## PPC
#

intensity_PP <- summarise_var_by_age_group(samples, 'log_lambda',  df_direction, df_community, df_period, df_age, transform = 'exp')
count_data <- prepare_count_data(stan_data, df_direction, df_community, df_period, df_age)
plot_intensity_PP(intensity_PP, count_data, outfile.figures)

predict_y <- summarise_var_by_age_recipient_group(samples, 'y_predict',  df_direction, df_community, df_period, df_age)
plot_PPC(predict_y, count_data, outfile.figures)

#
# Total infected
#

predict_z <- summarise_var_by_age_recipient_group(samples, 'z_predict',  df_direction, df_community, df_period, df_age)
plot_observed_to_augmented(predict_y, predict_z, outfile.figures)

#
## force of transmission
#

force_infection <- summarise_var_by_age_group(samples, 'log_beta', df_direction, df_community, df_period, df_age, transform = 'exp')
plot_force_infection(force_infection, outfile.figures)

#
## shift in sex-specific transmission dynamics
#

force_infection_sex_source <- summarise_var_by_sex_source(samples, var, df_age, df_direction, df_community, df_period)
plot_force_infection_sex_source(force_infection_sex_source, outfile.figures)


#
## shift in age-specific transmission dynamics
#

cat("\nPlot age-specific transmission dynamics\n")

# aggregated by agr group
force_infection_aggregated_age_group <- summarise_var_by_aggregated_age_group(samples, 'log_beta', df_direction, df_community, df_age, df_period, df_age_aggregated)
force_infection_aggregated_age_classification <- summarise_var_by_aggregated_age_classification(samples, 'log_beta', df_direction, df_community, df_age, df_period, df_age_aggregated)

# aggregated by age of the source                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
force_infection_age_source <- summarise_var_by_age_source(samples, 'log_beta', df_age, df_direction, df_community, df_period)
plot_force_infection_age_source(force_infection_age_source, outfile.figures)

force_infection_age_source[, is_among_5 := M %in% sort(M, decreasing = T)[1:5], by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
tmp <- force_infection_age_source[is_among_5 == T, list(total_M = paste0(round(sum(M)*100, 2), '%'), 
                                                        age_group = paste0(min(AGE_TRANSMISSION.SOURCE), '-', max(AGE_TRANSMISSION.SOURCE))), by = c('LABEL_COMMUNITY', 'PERIOD', 'LABEL_DIRECTION')]
print(tmp[order(LABEL_DIRECTION, PERIOD)])

# median age of source
median_age_source <- find_median_age_source(samples, 'log_beta', df_age, df_direction, df_community, df_period)
plot_median_age_source(median_age_source, outfile.figures)


cat("End of postprocessing_figures.R")