# copying from run_stan.R 

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
require(lubridate)
library(rstan)
library(gridExtra)
library(lognorm)
library(ggExtra)
library(Hmisc)

# laptop
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir <- '~/git/phyloflows'
  outdir <- '~/Box\ Sync/2021/phyloflows/'

  jobname <- 'new_treatment_cascade'
  stan_model <- 'gp_221201d'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/home/andrea'))
{
  indir <-'~/git/phyloflows'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'test'
  stan_model <- 'gp_221201d'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
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
outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
if(!dir.exists(dirname(outfile.figures))) dir.create(dirname(outfile.figures))
if(!dir.exists(dirname(outdir.table))) dir.create(dirname(outdir.table))

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'statistics_functions.R'))
source(file.path(indir, 'functions', 'stan_utils.R'))

# indicators -- fixed
only.transmission.after.start.observational.period <- T
only.transmission.before.stop.observational.period <- T
only.one.community <- 'inland'

# indicators -- sensitivity analyses
remove.pairs.from.rounds <- NULL
pairs_replicates.seed <- NULL
use_30com_pairs <- F
use_tsi_non_refined <- F
# load files
file.pairs <- file.path(indir, 'data', 'pairsdata_toshare_d1_w11_netfrompairs_postponessrem.rds')

#
# Define start time, end time and cutoff
#

file.path.round.timeline <- file.path(indir, 'data', 'RCCS_round_timeline_220905.RData')
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]

start_first_period_inland <- df_round_inland[round == 'R010', min_sample_date] # "2003-09-26"
stop_first_period_inland <- df_round_inland[round == 'R015', max_sample_date] # "2013-07-05"
start_second_period_inland <-df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"
stop_second_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

stopifnot(start_first_period_inland < stop_first_period_inland)
stopifnot(stop_first_period_inland < start_second_period_inland)
stopifnot(start_second_period_inland < stop_second_period_inland)

# 2nd line needed for ROUND_SPANYRS
df_round <- make.df.round(df_round_inland)
df_period <- make.df.period(start_first_period_inland, stop_first_period_inland, 
                            start_second_period_inland, stop_second_period_inland, 
                            df_round)


# load pairs
pairs.all <- read_pairs(file.pairs)


#
# Find phylo pairs 
#
pairs <- copy(pairs.all)
if(1)
{
  cat('Keep only heterosexual pairs\n')
  cat('Removing ', nrow(pairs[! ((SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M'))]), ' pairs\n')
  pairs <- pairs[(SEX.RECIPIENT == 'M' & SEX.SOURCE == 'F') | (SEX.RECIPIENT == 'F' & SEX.SOURCE == 'M')]
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
  
  cat('Keep only RCCS participants\n')
  cat('Removing ', nrow(pairs[(COMM.SOURCE == 'neuro' | COMM.RECIPIENT == "neuro")]), ' pairs\n')
  pairs <- pairs[COMM.SOURCE != 'neuro' & COMM.RECIPIENT != "neuro"]
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(!is.null(only.one.community)){
  cat('\nExcluding sources and recipients in ',   pairs[COMM.RECIPIENT != only.one.community, unique(COMM.RECIPIENT)] ,'\n')
  cat('Removing ', nrow(pairs[!(COMM.SOURCE == only.one.community & COMM.RECIPIENT == only.one.community)]), ' pairs\n')
  pairs <- pairs[COMM.SOURCE == only.one.community & COMM.RECIPIENT == only.one.community]
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(use_30com_pairs){
  cat('\nExcluding sources and recipients outside of the 30 continuously surveyed communities\n')
  comm_continuously_surveyed <- c(1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 
                                  89, 94, 106, 107, 108, 120, 391, 602, 754)
  cat('Removing ', nrow(pairs[!(COMM_NUM.SOURCE %in% comm_continuously_surveyed & COMM_NUM.RECIPIENT %in% comm_continuously_surveyed)]), ' pairs\n')
  pairs <- pairs[(COMM_NUM.SOURCE %in% comm_continuously_surveyed & COMM_NUM.RECIPIENT %in% comm_continuously_surveyed)]
}
if(only.transmission.after.start.observational.period){
  cat('\nFor inland excluding recipients infected before ', as.character(start_first_period_inland), '\n')
  cat('Removing ', nrow(pairs[DATE_INFECTION.RECIPIENT < start_first_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs <- pairs[!(DATE_INFECTION.RECIPIENT < start_first_period_inland & COMM.RECIPIENT == 'inland')]

  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(only.transmission.before.stop.observational.period){
  cat('\nFor inland excluding recipients infected after ', as.character(stop_second_period_inland), '\n')
  cat('Removing ', nrow(pairs[DATE_INFECTION.RECIPIENT > stop_second_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
  pairs <- pairs[!(DATE_INFECTION.RECIPIENT > stop_second_period_inland & COMM.RECIPIENT == 'inland')]
  
  cat('resulting in a total of ', nrow(pairs),' pairs\n\n')
}
if(!is.null(remove.pairs.from.rounds)){
  cat('\nExcluding pairs in inland community from round', remove.pairs.from.rounds, '\n')
  tmp <- df_round_inland[round %in% remove.pairs.from.rounds, list(min_exclusion = min(min_sample_date), 
                                                            max_exclusion = max(max_sample_date))]
  pairs[, DATE.COLLECTION.PAIR := max(c(DATE.COLLECTION.SOURCE, DATE.COLLECTION.RECIPIENT)), by = c('RECIPIENT', 'SOURCE')]
  cat('Removing ', nrow(pairs[COMM.RECIPIENT == 'inland' & DATE.COLLECTION.PAIR <= tmp[, max_exclusion] & DATE.COLLECTION.PAIR >= tmp[,min_exclusion ]]), ' pairs\n')
  pairs <- pairs[!(COMM.RECIPIENT == 'inland' & DATE.COLLECTION.PAIR <= tmp[, max_exclusion] & DATE.COLLECTION.PAIR >= tmp[,min_exclusion ])]
}

print.which.NA(pairs)
print.statements.about.pairs(copy(pairs))

# keep only pairs with source-recipient with a time of infection
pairs <- pairs[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
pairs[, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT := DATE_INFECTION.RECIPIENT < start_second_period_inland]
tab <- pairs[, list(count = .N), by = c('DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT', 'COMM.RECIPIENT', 'SEX.RECIPIENT')]
print_table(tab[order(DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT, COMM.RECIPIENT, SEX.RECIPIENT)])

# replace pairs with a bootstrap sample 
if(!is.null(pairs_replicates.seed)){
  cat('\nSeed for replicate ', pairs_replicates.seed)
  set.seed(pairs_replicates.seed)
  stopifnot(all(duplicated(pairs) == F))
  pairs <- pairs[sample(nrow(pairs), replace = T)]
  stopifnot(any(duplicated(pairs) == T))
}



outfile.figures <- '~/Downloads/EDF5/inkscape'

find_palette_round()
p <- plot_pairs(pairs, outfile.figures, nm_reqs = TRUE)
p1 <- plot_pairs_all(pairs.all, outfile.figures, nm_reqs=TRUE)
p2 <- plot_transmission_events_over_time(pairs, outfile.figures, nm_reqs=TRUE)

# extract legend
legend <- cowplot::get_legend(p2)
p2 <- p2 + theme(legend.position = 'none')

top <- (p1 + p2 + plot_layout(widths = c(1.2,3))) / 
    legend + 
    plot_layout( heights = c(100,1)) + 
    plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=8, face='bold', family='sans'))

all <- ggarrange_nature(top, p , heights=c(2,3), ncol=1)
all

ggsave_nature('~/Downloads/EDF5/extended_data_figure_pairs.pdf',all, w=18, h=21)
ggsave_nature
