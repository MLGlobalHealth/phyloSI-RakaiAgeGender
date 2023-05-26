cat("Start of get_main_figure2.R")

suppressPackageStartupMessages({
    library(rstan)
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(gridExtra)
    library(matrixStats)
    library(dplyr)
    library(lubridate)
    library(ggnewscale)
})

gitdir <- here::here()
source(file.path(gitdir, 'config.R'))

# load functions
source(file.path(gitdir.R.flow, 'postprocessing_summary_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_plot_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_utils_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_statistics_functions.R'))

jobname <- 'central3'
stan_model <- 'gp_221201d'
suffix <- paste(stan_model, jobname, sep='-')

# if on HPC, we know the username
if(dir.exists('/rds/')){
    if(usr == 'ab1820')
        outdir <- paste0("/rds/general/user/",usr,"/home/projects/2022/phyloflows/", suffix)
    if(usr == 'mm3218')
        outdir <- paste0("/rds/general/user/",usr,"/home/projects/2021/phyloflows/", suffix)
} else {
    if(usr == 'andrea'){
        outdir <- file.path(outdir.phyloflows, suffix) 
    }
}


outfile <- file.path(outdir, paste0(stan_model,'-', jobname))

# paths
path.to.stan.output = paste0(outfile, "-stanout_", jobname, ".rds")
.outfile.figures <- file.path(outdir, 'figures', suffix)
.outdir.table <-  file.path(outdir, 'tables', suffix)
if(0)
    path.to.suboutput <- '~/Downloads/subsample.rds'

# load data, then overwrite outdirs 
path.to.stan.data <- paste0(outfile, "-stanin_",jobname,".RData")
load(path.to.stan.data)

# samples: select the thinned sample if specified
if(exists('path.to.suboutput') )
{
    stopifnot("Thinned chains not found"=file.exists(path.to.suboutput))
    .outfile.figures <- '~/Downloads' 
    .outdir.table <- '~/Downloads' 
    samples <- readRDS( file=path.to.suboutput)
}else{
    fit <- readRDS(path.to.stan.output)
    samples <- rstan::extract(fit)
}
outfile.figures <- .outfile.figures
outdir.table <- .outdir.table


# temporary
# source(file.path(indir, 'functions', 'summary_functions.R'))
# df_direction <- get.df.direction()

#
# offset
#

log_offset_round <- find_log_offset_by_round(
    stan_data,
    eligible_count_round, 
    df_estimated_contact_rates,
    use_number_susceptible_offset, use_contact_rates_prior)

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

median_age_source <- find_summary_output_by_round(
    samples,
    'log_lambda_latent', 
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
    transform = 'exp',
    standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND'),
    quantile_age_source = T, 
    save_output = FALSE)

# by age groups
df_age_aggregated <- get.age.aggregated.map(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'))

median_age_source_group <- find_summary_output_by_round(
    samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE', 'AGE_GROUP_INFECTION.RECIPIENT'),
    transform = 'exp',
    standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
    quantile_age_source = T, save_output = FALSE ) 


expected_contribution_age_group_source2 <- find_summary_output_by_round(
    samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_GROUP_INFECTION.RECIPIENT'),
    transform = 'exp',
    standardised.vars = c('INDEX_ROUND'),
    save_output = FALSE)


# load plot requirements
naturemed_reqs()

p_b <- plot_median_age_source_group(median_age_source_group,
    expected_contribution_age_group_source2,
    reported_contact,
    outfile.figures,
    nm_reqs=TRUE) + reqs

# ggsave('~/Downloads/tmp.pdf', p_b, w=12.5, h=11, unit='cm')

# age-specific contribution to transmission among all sources by sex
contribution_age_source <-  find_summary_output_by_round(
    samples,
    'z_predict',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
    standardised.vars = c('INDEX_ROUND'))

expected_contribution_sex_source <- find_summary_output_by_round(
    samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND'),
    transform = 'exp',
    standardised.vars = c('INDEX_ROUND'))

# age-specific contribution to transmission among all age and  sex
expected_contribution_age_source2 <- find_summary_output_by_round(
    samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
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
