library(data.table)
library(phyloscannerR)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
# Needed for Network plot:
library(grid)
library(ggtree)
library(ggnet) 
  
# change as appropriate
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir.repository <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
  outdir <- '~/Box\ Sync/2021/phyloflows/'
}
if(dir.exists('~/Documents/ratmann_deepseq_analyses'))
{
  indir.repository <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '~/Documents/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequence_analyses_MRC   <- '~/Documents/PANGEA2_MRCUVRI'
  indir.deepsequencedata <- '~/Documents/ratmann_pangea_deepsequencedata/'
  outdir <- '~/Documents/2021/phyloflows'
}

# indicators 
include.mrc <- T
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.5
date_implementation_UTT <- as.Date('2016-12-01')
lab <- paste0('MRC_', include.mrc, '_OnlyHTX_', include.only.heterosexual.pairs, '_threshold_', threshold.likely.connected.pairs)

# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.path.meta.data.rccs.1 <- file.path(indir.deepsequence_analyses, 'RakaiPangeaMetaData_v2.rda')
file.path.meta.data.rccs.2 <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')
file.path.meta.data.mrc <- file.path(indir.deepsequencedata, 'PANGEA2_MRC','200319_pangea_db_sharing_extract_mrc.csv')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.time.first.positive <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '211111_pangea_db_sharing_extract_rakai_age_firstpos_lastneg.csv')
outdir.lab <- file.path(outdir, lab); dir.create(outdir.lab)
  
# load functions
source(file.path(indir.repository, 'functions', 'summary_functions.R'))
source(file.path(indir.repository, 'functions', 'plotting_functions.R'))
source(file.path(indir.repository, 'functions', 'stan_utils.R'))
source(file.path(indir.repository, 'functions', 'check_potential_TNet.R'))

# load files
load(file.path.chains.data)
load(file.path.meta.data.rccs.1)
meta.rccs.1 <- as.data.table(rccsData)
meta.rccs.2 <- as.data.table( read.csv(file.path.meta.data.rccs.2))
meta.mrc <- as.data.table( read.csv(file.path.meta.data.mrc))


# load keys
anonymisation.keys <- as.data.table(read.csv(file.anonymisation.keys))
community.keys <- as.data.table( read.csv(file.community.keys) )

time.first.positive <- as.data.table( read.csv(file.time.first.positive))
time.first.positive <- make.time.first.positive(time.first.positive)

# get meta data
meta_data <- get.meta.data(meta.rccs.1, meta.rccs.2, meta.mrc, time.first.positive, anonymisation.keys, community.keys)

# get likely transmission pairs
chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)

# study whether couple found in previous analyses are included in the same Potential Transmission Network cluster.
# commented by Melodie: this gives me an error
# print.statements.about.potential.TNet()

# merge meta data to source and recipient
pairs.all <- pairs.get.meta.data(chain, meta_data)

if(!include.mrc){
  cat('Keep only pairs in RCCS\n')
  pairs.all <- pairs.all[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS']
}
if(include.only.heterosexual.pairs){
  cat('Keep only heterosexual pairs\n')
  pairs.all <- pairs.all[(sex.RECIPIENT == 'M' & sex.SOURCE == 'F') | (sex.RECIPIENT == 'F' & sex.SOURCE == 'M')]
}

print.statements.about.pairs(copy(pairs.all), outdir.lab)

# keep only pairs with source-recipient with proxy for the time of infection
pairs <- pairs.all[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
pairs[, date_infection_before_UTT.RECIPIENT := date_infection.RECIPIENT <= date_implementation_UTT]

# prepare age map
get.age.map(pairs)

# make some explanatory plots
plot_hist_age_infection(copy(pairs), outdir.lab)
plot_hist_time_infection(copy(pairs), date_implementation_UTT, outdir.lab)
plot_age_infection_source_recipient(pairs[sex.SOURCE == 'M' & sex.RECIPIENT == 'F'], 'Male -> Female', 'MF', outdir.lab)
plot_age_infection_source_recipient(pairs[sex.SOURCE == 'F' & sex.RECIPIENT == 'M'], 'Female -> Male', 'FM', outdir.lab)
plot_CI_age_infection(pairs, outdir.lab)

# Unfortunately have to run 'manually' at the moment. Need to understand what s wrong with this one.
dchain <- as.data.table(dchain)
dc <- as.data.table(dc)
phsc.plot.transmission.network(copy(dchain), copy(dc),outdir=outdir, arrow=arrow(length=unit(0.02, "npc"), type="open"), edge.size = 0.1)

# prepare stan data
stan_data <- prepare_stan_data(pairs, df_age)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3, n_knots_rows = 15, n_knots_columns = 15, unique(df_age$age_infection.SOURCE))

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.lab, paste0("stanin_",lab,".RData")) )
