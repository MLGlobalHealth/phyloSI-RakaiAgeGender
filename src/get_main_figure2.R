cat("Start of get_main_figure2.R")

library(rstan)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)
library(ggnewscale)

jobname <- 'central3'
stan_model <- 'gp_221201d'

indir <- "/rds/general/user/mm3218/home/git/phyloflows"
outdir <- paste0("/rds/general/user/mm3218/home/projects/2021/phyloflows/", stan_model, '-', jobname)

if(1){
    indir <- "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live/temporary_for_Andrea/phyloflows/"
    outdir <- file.path(indir,"results", paste0(stan_model, '-', jobname))
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
source(file.path(indir, 'functions', 'postprocessing_utils_functions.R'))
source(file.path(indir, 'functions', 'postprocessing_statistics_functions.R'))

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))

# paths
path.to.stan.output = paste0(outfile, "-stanout_", jobname, ".rds")
.outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
.outdir.table <-  file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
path.to.suboutput <- '~/Downloads/subsample.rds'

# load data
path.to.stan.data <- paste0(outfile, "-stanin_",jobname,".RData")
load(path.to.stan.data)
outfile.figures <- .outfile.figures
outdir.table <- .outdir.table

# samples 
if(file.exists(path.to.suboutput))
{
    .outfile.figures <- '~/Downloads' #file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
    .outdir.table <- '~/Downloads' # file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
    samples <- readRDS( file=path.to.suboutput)
}else{
    fit <- readRDS(path.to.stan.output)
    samples <- rstan::extract(fit)
}

if(0)
{
    index_array <- function(x, dim, value, drop = FALSE) { 
      # Create list representing arguments supplied to [
      # bquote() creates an object corresponding to a missing argument
      indices <- rep(list(bquote()), length(dim(x)))
      indices[[dim]] <- value

      # Generate the call to [
      call <- as.call(c(
        list(as.name("["), quote(x)),
        indices,
        list(drop = drop)))
      # Print it, just to make it easier to see what's going on
      print(call)

      # Finally, evaluate it
      eval(call)
    }

    .f <- function(x) 
    {
        index_array(x, 1, 1:1000)
    }

    samples <- lapply(samples, .f)
}


# need to be able to get to: 
# a Expected_contribution_sex_age_inland.pdf 
# b MedianAgeSource_ByAgeGroupRecipient_inland.pdf 
# c Expected_contribution_age_inland.pdf 
# d Expected_contribution_sex_barplot.png 


# temporary
source(file.path(indir, 'functions', 'summary_functions.R'))
df_direction <- get.df.direction()

#
# offset
#

log_offset_round <- find_log_offset_by_round(stan_data, eligible_count_round, df_estimated_contact_rates, 
                                             use_number_susceptible_offset, use_contact_rates_prior)

# log offset formula (per year)
log_offset_formula <- 'log_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED'
if(!use_number_susceptible_offset)
  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED'
if(use_contact_rates_prior)
  log_offset_formula = paste0(log_offset_formula, ' + log_CONTACT_RATES')

# log offset formula (per year per unsuppressed)
log_offset_formula_perunsuppressed <- 'log_SUSCEPTIBLE'
if(!use_number_susceptible_offset)
  log_offset_formula_perunsuppressed = 'log_PROP_SUSCEPTIBLE'
if(use_contact_rates_prior)
  log_offset_formula_perunsuppressed = paste0(log_offset_formula_perunsuppressed, ' + log_CONTACT_RATES')

# log offset formula (per year per susceptible)
log_offset_formula_persusceptible <- 'log_INFECTED_NON_SUPPRESSED'
if(use_contact_rates_prior)
  log_offset_formula_persusceptible = paste0(log_offset_formula_persusceptible, ' + log_CONTACT_RATES')


#
# Summarise data and merge to maps for figures
#

count_data <- prepare_count_data(stan_data)
incidence_cases_recipient_round <- prepare_incidence_cases(incidence_cases_round)
unsuppressed_share_sex <- prepare_unsuppressed_share(unsuppressed_share, c('SEX'))
unsuppressed_share_sex_age <- prepare_unsuppressed_share(unsuppressed_share, c('SEX', 'AGEYRS'))
prevalence_prop_sex<- prepare_infected_share(infected_share, 'SEX')
reported_contact <- clean_reported_contact(df_reported_contact)
df_unsuppressed_median_age<-prepare_unsuppressed_median_age(unsuppressed_median_age)

median_age_source <- find_summary_output_by_round(samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'),
    transform = 'exp',
    standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
    median_age_source = T)

# by age groups
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))
median_age_source_group <- find_summary_output_by_round(samples, 'log_lambda_latent', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                        transform = 'exp',
                                                        standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                        quantile_age_source = T)
expected_contribution_age_group_source2 <- find_summary_output_by_round(samples, 'log_lambda_latent',
                                                                        c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
                                                                        transform = 'exp',
                                                                        standardised.vars = c('INDEX_ROUND'))


# load plot requirements
naturemed_reqs()

p_b <- plot_median_age_source_group(median_age_source_group,
    expected_contribution_age_group_source2,
    reported_contact,
    outfile.figures,
    nm_reqs=TRUE) + reqs

# age-specific contribution to transmission among all sources by sex
contribution_age_source <-  find_summary_output_by_round(
    samples,
    'z_predict',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
    standardised.vars = c('INDEX_ROUND'))

# age-specific contribution to transmission among all age and  sex
expected_contribution_age_source2 <- find_summary_output_by_round(
    samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
    transform = 'exp',
                                                                  standardised.vars = c('INDEX_ROUND'))


expected_contribution_sex_source <- find_summary_output_by_round(
    samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND'),
    transform = 'exp',
    standardised.vars = c('INDEX_ROUND'))


p_c1 <- plot_contribution_age_source_unsuppressed(
    expected_contribution_age_source2,
    unsuppressed_share_sex_age,
    median_age_source,
    df_unsuppressed_median_age,
    outfile.figures,
    'Expected_contribution', 
    nm_reqs=TRUE) + reqs

p_c2 <- plot_contribution_sex_source(
    expected_contribution_sex_source,
    unsuppressed_share_sex,
    prevalence_prop_sex,
    outfile.figures,
    'Expected_contribution',
    nm_reqs=TRUE) + reqs

p_a <- plot_contribution_age_source(
    expected_contribution_age_source2,
    median_age_source,
    outfile.figures,
    'Expected_contribution_sex',
    nm_reqs = TRUE) + reqs

# use patchwork to manage everything 
library(patchwork)
p_c <-  (p_c1 / p_c2) + plot_layout(heights=c(6,1)) 
right_side <- ( (p_b + labs(subtitle = 'b')) / ( p_c  ) + plot_layout( heights = c(6, 7)) )
p_2 <- (p_a + labs(subtitle = 'a') |  right_side) + plot_layout( widths = c(2, 4)) 

ggsave_nature(p_2, filename=paste0(outfile.figures, 'MainFigure2.pdf'))
