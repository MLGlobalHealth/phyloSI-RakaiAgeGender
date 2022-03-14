library(data.table)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
library(grid)
library(ggtree)
library(ggnet) 
require(lubridate)

# change as appropriate
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
  outdir <- '~/Box\ Sync/2021/phyloflows/'
  
  jobname <- 'cutoff_2014'
  stan_model <- 'gp_220108'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/home/andrea'))
{
  indir <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata/'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'cutoff_2014'
  stan_model <- 'gp_220108'
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

# indicators 
cutoff_date <- as.Date('2014-01-01')
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.5
use.tsi.estimates <- F
remove.inconsistent.infection.dates <- F
remove.young.individuals <- T
only.inland <- T

# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.path.meta <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_20220308.csv')
file.path.tsiestimates <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', 'TSI_estimates_220119.csv')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.incidence <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_incident_cases_220311.csv')

# load functions
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'stan_utils.R'))
source(file.path(indir, 'functions', 'check_potential_TNet.R'))

# load chains
load(file.path.chains.data)
dchain <- as.data.table(dchain)

# load meta data
aik <- .read(file.anonymisation.keys); aik$X <- NULL
meta_data <- .read(file.path.meta)

# load Tanya's estimate time since infection using phylogenetic data
time.since.infection <- make.time.since.infection(as.data.table(read.csv(file.path.tsiestimates)))

# load incidence from adam
incidence <- read.csv(file.incidence)


#
# TRANSFORM AND MERGE DATA
#

# get time of infection (using Tanya's estimate if use.tsi.estimates == T)
meta_data <- find.time.of.infection(meta_data, time.since.infection, use.tsi.estimates)

# get likely transmission pairs
chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)

# merge meta data to source and recipienx
pairs.all <- pairs.get.meta.data(chain, meta_data, aik)

if(include.only.heterosexual.pairs){
  cat('Keep only heterosexual pairs\n')
  cat('Removing ', nrow(pairs.all[! ((sex.RECIPIENT == 'M' & sex.SOURCE == 'F') | (sex.RECIPIENT == 'F' & sex.SOURCE == 'M'))]), ' pairs\n')
  pairs.all <- pairs.all[(sex.RECIPIENT == 'M' & sex.SOURCE == 'F') | (sex.RECIPIENT == 'F' & sex.SOURCE == 'M')]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(remove.inconsistent.infection.dates){
  plot_pairs_infection_dates(pairs.all)
  cat('Remove infections for which estimated date at infection of source is more than one year after the estimated date at infection of the recipient.\n ')
  cat('Removing ', nrow(pairs.all[ date_infection.SOURCE >= date_infection.RECIPIENT + 365]), ' pairs\n')
  pairs.all <- pairs.all[! date_infection.SOURCE >= date_infection.RECIPIENT + 365 ]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(remove.young.individuals){
  # exclude young indivis
  cat('\nExcluding very young individuals')
  cat('Removing ', nrow(pairs.all[age_infection.SOURCE < 11 | age_infection.RECIPIENT < 11]), ' pairs\n')
  pairs.all <- pairs.all[age_infection.SOURCE >= 11 & age_infection.RECIPIENT >= 11]
  cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
}
if(only.inland){
  cat('\nInclude only inland recipients')
  pairs.all <- pairs.all[comm.RECIPIENT == 'inland']
}

print.which.NA(pairs.all)
print.statements.about.pairs(copy(pairs.all))

# which base frequency files we have on the HPC
# atm gives error: maybe TODO when I understand more about PHSC pipeline
# missing_bff <- print.statements.about.basefreq.files(pairs.all)

# keep only pairs with source-recipient with proxy for the time of infection
pairs <- pairs.all[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
pairs[, date_infection_before_cutoff.RECIPIENT := date_infection.RECIPIENT < cutoff_date]


#
# PREPARE MAPS
#

# prepare age map
df_age <- get.age.map(pairs, age_bands_reduced = 4)

# prepare group map
df_group <- get.group.map()


#
# MAKE EXPLANATORY PLOTS
#

if(1){
  plot_hist_age_infection(copy(pairs), outfile)
  plot_hist_time_infection(copy(pairs), cutoff_date, outfile)
  plot_age_infection_source_recipient(pairs[sex.SOURCE == 'M' & sex.RECIPIENT == 'F'], 'Male -> Female', 'MF', outfile)
  plot_age_infection_source_recipient(pairs[sex.SOURCE == 'F' & sex.RECIPIENT == 'M'], 'Female -> Male', 'FM', outfile)
  plot_CI_age_infection(pairs, outfile)
  plot_CI_age_transmission(pairs, outfile)
  # phsc.plot.transmission.network(copy(as.data.table(dchain)), copy(as.data.table(dc)), pairs,outdir=outfile, arrow=arrow(length=unit(0.02, "npc"), type="open"), edge.size = 0.1)
}


#
# PREPARE STAN DATA
#

# prepare stan data
stan_data <- prepare_stan_data(pairs, df_age, df_group)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3, 
                                      n_knots_rows = 8, n_knots_columns = 8, 
                                      X = unique(df_age$age_transmission.SOURCE),
                                      Y = unique(df_age$age_infection.RECIPIENT))
stan_data <- add_prior_gp_mean(stan_data, df_age, outfile)

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.jobname, paste0("stanin_",jobname,".RData")) )


