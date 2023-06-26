cat("Start of postprocessing_figure_time_trends_sources.R\n")

library(rstan)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(matrixStats)
library(dplyr)
library(lubridate)
library(ggnewscale)
library(patchwork)

usr <- Sys.info()[['user']]

if(usr == 'melodiemonod'){# if on laptop
  indir <- here::here()
  stan_model <- 'gp_230602'
  jobname <- 'newdetectionprob'
  outdir <- file.path('~/Box\ Sync/2021/phyloflows/', paste0(stan_model,'-', jobname))
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

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))

# paths
outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))

# load functions
gitdir <- copy(indir)
source(file.path(gitdir, 'config.R'))
source(file.path(gitdir.R.flow, 'plotting_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_summary_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_plot_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_utils_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_statistics_functions.R'))


#
# Load files
#

# name files data
file.reported_contact <- paste0(outdir.table, '-data-reported_contact.rds')
file.unsuppressed_share_sex_age <- paste0(outdir.table, '-data-unsuppressed_share_by_sex_ageyrs.rds')
file.unsuppressed_share_sex <- paste0(outdir.table, '-data-unsuppressed_share_by_sex.rds')
file.df_unsuppressed_median_age <- paste0(outdir.table, '-data-median_age_unsuppressed.rds')
file.prevalence_prop_sex <- paste0(outdir.table, '-data-share_prevalence_by_sex.rds')

# name files with estimated quantities
file.median_age_source <- paste0(outdir.table, '-output-log_lambda_latentby_direction_round_quantilestandardisedby_direction_round_median_age_source_total.rds')
file.median_age_source_group <- paste0(outdir.table, '-output-log_lambda_latentby_direction_round_age_group_infection.recipient_quantilestandardisedby_direction_round_age_group_infection.recipient_5yrs_age_band.rds')
file.expected_contribution_age_group_source2 <- paste0(outdir.table, '-output-log_lambda_latentby_direction_round_age_group_infection.recipientstandardisedby_round_5yrs_age_band.rds')
file.expected_contribution_sex_source <- paste0(outdir.table, '-output-log_lambda_latentby_direction_roundstandardisedby_round.rds')
file.expected_contribution_age_source2 <- paste0(outdir.table, '-output-log_lambda_latentby_direction_round_age_transmission.sourcestandardisedby_round.rds')

all.files.exist <- all(file.exists(c(
  file.median_age_source ,
  file.median_age_source_group,
  file.expected_contribution_age_group_source2,
  file.expected_contribution_sex_source,
  file.expected_contribution_age_source2, 
  file.reported_contact,
  file.unsuppressed_share_sex_age,
  file.unsuppressed_share_sex,
  file.df_unsuppressed_median_age,
  file.prevalence_prop_sex
  )) )

stopifnot(all.files.exist)

# median age source by 1-year age band
median_age_source <- readRDS(file.median_age_source)

# median age source by age group
median_age_source_group <- readRDS(file.median_age_source_group)

# expected contribution
expected_contribution_age_group_source2 <- readRDS(file.expected_contribution_age_group_source2)

# sex-specific contribution to transmission
expected_contribution_sex_source <- readRDS(file.expected_contribution_sex_source)

# age-specific contribution to transmission among all age and  sex
expected_contribution_age_source2 <- readRDS(file.expected_contribution_age_source2)
  
# estimated number of sexual contacts
reported_contact <- readRDS(file.reported_contact)

# unsuppressed share by sex and age
unsuppressed_share_sex_age <- readRDS(file.unsuppressed_share_sex_age)
unsuppressed_share_sex <- readRDS(file.unsuppressed_share_sex)

# median age unsuppressed
df_unsuppressed_median_age <- readRDS(file.df_unsuppressed_median_age)

# share prevalence by sex
prevalence_prop_sex <- readRDS(file.prevalence_prop_sex)


#
# Plot
#

# load plot requirements
naturemed_reqs()

# find color palette of rounds
find_palette_round()

# find range age
range_age_non_extended <- range(unsuppressed_share_sex_age$AGEYRS)

p_b <- plot_median_age_source_group(median_age_source_group,
                                    expected_contribution_age_group_source2,
                                    reported_contact,
                                    outfile.figures,
                                    nm_reqs=TRUE) + reqs

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
p_c <-  (p_c1 / p_c2) + plot_layout(heights=c(6,1)) 
right_side <- ( (p_b + labs(subtitle = 'b')) / ( p_c  ) + plot_layout( heights = c(6, 7)) )
p_2 <- (p_a + labs(subtitle = 'a') |  right_side) + plot_layout( widths = c(2, 4)) 

ggsave_nature(p_2, filename=paste0(outfile.figures, '-output-MainFigure2.pdf'))

cat("End of postprocessing_figure_time_trends_sources.R\n")

