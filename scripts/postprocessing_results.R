cat("Start of postprocessing_results.R")

library(rstan)
library(data.table)	
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)

jobname <- 'onlyinland_cutoff2014'
stan_model <- 'gp_220108'
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

# maps
range_age_observed_cat <- find_range_age_observed(copy(pairs), df_group)
range_age_observed <- pairs[, list(min_age = min(c(age_infection.RECIPIENT, age_transmission.SOURCE)), 
                                                 max_age = max(c(age_infection.RECIPIENT, age_transmission.SOURCE)))]
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'), incidence)

# samples 
fit <- readRDS(path.to.stan.output)
samples <- rstan::extract(fit)


#
## convergence diagnostics 
#

make_convergence_diagnostics_stats(fit, outdir.table)


#
## intensity of the poisson process
#

cat("\nPlot transmission intensity\n")

intensity_PP <- summarise_var_by_age_group(samples, 'log_lambda', df_group, df_age, transform = 'exp')
count_data <- prepare_count_data(stan_data, df_age, df_group)
plot_intensity_PP(intensity_PP, count_data, range_age_observed, outfile.figures)

transmission_flows <- find_transmission_flows(samples, df_group, df_age)
plot_transmission_flows(transmission_flows = transmission_flows, count_data = count_data, range_age_observed = range_age_observed, 
                        with_contour = T, outdir = outfile.figures)
plot_transmission_flows_vertical(transmission_flows = transmission_flows, count_data = count_data, range_age_observed = range_age_observed, 
                        with_contour = T, outdir = outfile.figures)

transmission_flows_aggregated <- find_transmission_flows_aggregated(samples, df_group, df_age, df_age_aggregated)
transmission_flows_aggregated2 <- find_transmission_flows_aggregated2(samples, df_group, df_age, df_age_aggregated)


#
## Standardised transmission flows startified by time period
#

cat("\nPlot Standardised transmission flows\n")

standardised_transmission_flows <- find_standardised_transmission_flows(samples, df_group, df_age, incidence)
plot_transmission_flows(transmission_flows = standardised_transmission_flows, lab = 'Standardised', 
                        range_age_observed=range_age_observed,outdir = outfile.figures, with_contour = T)

standardised_transmission_flows_aggregated <- find_standardised_transmission_flows_aggregated(samples, df_group, df_age, df_age_aggregated, incidence)
standardised_transmission_flows_aggregated2 <- find_standardised_transmission_flows_aggregated2(samples, df_group, df_age, df_age_aggregated, incidence)

plot_transmission_flows_aggregated(transmission_flows_aggregated, standardised_transmission_flows_aggregated, df_age_aggregated, outfile.figures)
plot_transmission_flows_aggregated2(transmission_flows_aggregated2, standardised_transmission_flows_aggregated2, df_age_aggregated, outfile.figures)


#
## Standardised transmission flows statified by round
#

cat("\nPlot Standardised transmission flows by round\n")

standardised_transmission_flows_by_round <- find_standardised_transmission_flows_by_round(samples, df_group, df_age, incidence)
plot_transmission_flows_by_round(transmission_flows = standardised_transmission_flows_by_round, lab = 'Standardised', 
                        range_age_observed=range_age_observed,outdir = outfile.figures, with_contour = T)

standardised_transmission_flows_aggregated_by_round <- find_standardised_transmission_flows_aggregated_by_round(samples, df_group, df_age, df_age_aggregated, incidence)
standardised_transmission_flows_aggregated2_by_round <- find_standardised_transmission_flows_aggregated2_by_round(samples, df_group, df_age, df_age_aggregated, incidence)

plot_transmission_flows_aggregated_by_round(standardised_transmission_flows_aggregated_by_round, df_age_aggregated, outfile.figures)
plot_transmission_flows_aggregated2_by_round(standardised_transmission_flows_aggregated2_by_round, df_age_aggregated, outfile.figures)

# transmission flow from age and sex 
standardised_transmission_flows_across_age_by_round <- find_standardised_transmission_flows_across_age_by_round(samples, df_group, df_age, incidence)
plot_standardised_transmission_flows_across_age_by_round(standardised_transmission_flows_across_age_by_round, outfile.figures)
  




#
## shift in sex-specific transmission dynamics
#

sex_source <- find_sex_source(samples, df_group, df_age, incidence)
sex_source_standardised <- find_sex_source_standardised(samples, df_group, df_age, incidence)
plot_sex_source_standardised(sex_source, sex_source_standardised, outfile.figures)

sex_source_standardised_by_round <- find_sex_source_standardised_by_round(samples, df_group, df_age, incidence)
plot_sex_source_standardised_by_round(sex_source_standardised_by_round, outfile.figures)

#
## shift in age-specific transmission dynamics
#

cat("\nPlot age-specific transmission dynamics\n")

# median age of source
age_source <- find_age_source_by_age_group(samples, df_group, df_age)
plot_median_age_source(age_source, outfile.figures)
plot_median_age_source_with_empirical_data(age_source, pairs, outfile.figures)


# median age of recipient
age_recipient <- find_age_recipient_by_age_group(samples, df_group, df_age)
plot_median_age_recipient_with_empirical_data(age_recipient, pairs, outfile.figures)

age_recipient_standardised <- find_age_recipient_by_age_group_standardised(samples, df_group, df_age, incidence)
plot_median_age_recipient(age_recipient, age_recipient_standardised, outfile.figures)


# 
# age_source_difference <- find_age_source_difference_by_age_group(samples, df_group, df_age)
# plot_median_age_source_difference(age_source_difference, outfile.figures)
# 
# age_source_overall <- find_age_source_by_group(samples, df_group, df_age, di, range_age_observed)
# plot_median_age_source_overall(age_source_overall, outfile.figures)

# median age of recipient



cat("End of postprocessing_results.R")

