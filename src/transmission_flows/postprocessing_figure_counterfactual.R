cat("Start of postprocessing_figure_counterfactual.R\n")

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
  outdir <- '~/Box\ Sync/2021/phyloflows/'
  stan_model <- 'gp_230602'
  jobname <- 'usedep'
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
source(file.path(gitdir.R.flow, 'postprocessing_summary_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_plot_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_utils_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_statistics_functions.R'))


#
# Load files
#

# name files 
file.counterfactuals_a_f <- paste0(outdir.table, '-output-counterfactuals_art_up_to_female1.rds')
file.counterfactuals_a_f05 <- paste0(outdir.table, '-output-counterfactuals_art_up_to_female0.5.rds')
file.counterfactuals_a_959595 <- paste0(outdir.table, '-output-counterfactuals_9595951.rds')
file.counterfactuals_a_909090 <- paste0(outdir.table, '-output-counterfactuals_9090901.rds')
file.incidence_factual <- paste0(outdir.table, '-output-log_betaby_direction_round_age_infection.recipient.rds')
  
all.files.exist <- all(file.exists(c(
  file.counterfactuals_a_f ,
  file.counterfactuals_a_f05 ,
  file.counterfactuals_a_959595,
  file.counterfactuals_a_909090, 
  file.incidence_factual)) )

stopifnot(all.files.exist)

# counterfactual treating all men as much as female are diagnosed/treated/suppressed
counterfactuals_a_f <- readRDS(file.counterfactuals_a_f)

#  counterfactual treating all men half way to as much as female are diagnosed/treated/suppressed
counterfactuals_a_f05 <- readRDS(file.counterfactuals_a_f05)

# counterfactual treating all men 95 95 95
counterfactuals_a_959595 <- readRDS(file.counterfactuals_a_959595)

# counterfactual treating all men 90 90 90
counterfactuals_a_909090 <- readRDS(file.counterfactuals_a_909090)

# incidence under the factual scenario 
incidence_factual <- readRDS(file.incidence_factual)

#
# plot
#

naturemed_reqs()
fig4_list <- plot_counterfactual(
  counterfactuals_a_f,
  counterfactuals_a_f05,
  counterfactuals_a_959595,
  counterfactuals_a_909090,
  incidence_factual,
  "Unsuppressed",
  outfile.figures, 
  nm_reqs = TRUE
)

# remember to add labels
p_a <-fig4_list[['a']] + reqs + theme(strip.text = element_blank())
p_b <-fig4_list[['b']] 
p_c <-fig4_list[['c']]
p_d <-fig4_list[['d']]

# Prepare legend grabbing it from 
cols <- c('#F1A661', '#C55300', '#749F82')
p_a_labs <- c(
    'Remaining virally unsuppressed in R18',
    'Closing the suppression gap\nin men relative to women',
    'Closing half the suppression gap\nin men relative to women',
    '95-95-95 in men',
    'Already virally suppressed in R18'
)

p_a_legend <- p_a  +
    labs(fill='Counterfactual interventions', color='Counterfactual interventions') +
    theme(legend.position='bottom') +
    scale_fill_manual(values = c('grey50', cols[c(2,1,3)], 'grey80'), labels=p_a_labs)

p_d_legend <- p_d +
    theme(legend.position='bottom') + 
    guides(fill='none', color='none')

legend_a <- get_legend(p_a_legend)
legend_d <- get_legend(p_d_legend)

legends <- ggarrange_nature(
    as_ggplot(legend_a),
    as_ggplot(legend_d),
    ncol=2)

# consider removing legend for a 

p_bcd <- ggarrange_nature(
    p_b , p_c , p_d,
    ncol=1,
    heights=c(3.3, 3.3, 4.6),
    labels = c('b', 'c', 'd'))

p_a_with_legend <- ggarrange_nature(
    p_a , legends,
    ncol = 1,
    heights = c(4, 2),
    labels = c('a', ''), 
    add_reqs = F
)

p_abcd <- ggarrange_nature(
    p_a_with_legend, p_bcd,
    ncol=2,
    widths=c(3.2, 3.4)
)

ggsave_nature(p_abcd,
    filename=paste0(outfile.figures, '-output-MainFigure4.pdf'),
    w=18, h=15)


#
# Prepare excel for Global UNAIDS report
#

save_counterfactual_results_for_UNAIDS(  counterfactuals_a_f,
                                         counterfactuals_a_f05,
                                         counterfactuals_a_959595,
                                         counterfactuals_a_909090,
                                         incidence_factual,
                                         "Unsuppressed", 
                                         outdir.table)

cat("End of postprocessing_figure_counterfactual.R\n")
