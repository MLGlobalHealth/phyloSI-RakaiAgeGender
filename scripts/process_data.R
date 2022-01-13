library(data.table)
library(phyloscannerR)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
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
include.mrc <- F
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.5
use.tsi.estimates <- T
cutoff_date <- as.Date('2015-01-01')
jobname <- 'cutoff2015priorgp'
lab <- paste0('MRC_', include.mrc, '_OnlyHTX_', include.only.heterosexual.pairs, '_threshold_', threshold.likely.connected.pairs, '_jobname_', jobname)

# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.path.meta.data.rccs.1 <- file.path(indir.deepsequence_analyses, 'RakaiPangeaMetaData_v2.rda')
file.path.meta.data.rccs.2 <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')
file.path.meta.data.mrc <- file.path(indir.deepsequencedata, 'PANGEA2_MRC','200319_pangea_db_sharing_extract_mrc.csv')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.time.first.positive <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '211111_pangea_db_sharing_extract_rakai_age_firstpos_lastneg.csv')
file.path.phscinput <- file.path(indir.deepsequence_analyses, '210120_RCCSUVRI_phscinput_runs.rds')
file.path.bflocs <- file.path(indir.repository, 'data/bfloc2hpc_20220103.rds')
file.path.tsiestimates <- file.path(indir.repository, 'data/TSI_estimates_20220601.csv')

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
anonymisation.keys <- as.data.table(read.csv(file.anonymisation.keys))[, X:=NULL]
print.anonkey.statements()
community.keys <- as.data.table( read.csv(file.community.keys) )

# get age at infection
time.first.positive <- as.data.table( read.csv(file.time.first.positive))
time.since.infection <- as.data.table(read.csv(file.path.tsiestimates))
time.since.infection <- merge(time.since.infection, anonymisation.keys, by='AID')
if(! use.tsi.estimates)
{
  cat("Do not use TSI estimates ")
  time.since.infection[, TSI_estimated_mean := 1]
}
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

if(file.exists(file.path.phscinput) & file.exists(file.path.bflocs))
{
  missing_bff <- print.statements.about.basefreq.files(pairs.all)
  missing_bff <- merge(anonymisation.keys, missing_bff[HPC_EXISTS == FALSE, ], by='AID')
  missing_bff[, HPC_EXISTS := NULL]
  
  # if want to send Tanya:
  tmp <- missing_bff[,  .(PT_ID, PREFIX)]
  name <- file.path(indir.repository, 'data/missing_bf_files_20220106.csv')
  if(!file.exists(name)){write.csv(tmp, name, row.names = F)}
  # commented out here, but all missing bf's have RCCS2 prefixes.
}


# keep only pairs with source-recipient with proxy for the time of infection
pairs <- pairs.all[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
pairs[, date_infection_before_cutoff.RECIPIENT := date_infection.RECIPIENT < cutoff_date]

# prepare age map
df_age <- get.age.map(pairs, age_bands_reduced = 4)

# prepare group map
df_group <- get.group.map()

# make some explanatory plots
plot_hist_age_infection(copy(pairs), outdir.lab)
plot_hist_time_infection(copy(pairs), cutoff_date, outdir.lab)
plot_age_infection_source_recipient(pairs[sex.SOURCE == 'M' & sex.RECIPIENT == 'F'], 'Male -> Female', 'MF', outdir.lab)
plot_age_infection_source_recipient(pairs[sex.SOURCE == 'F' & sex.RECIPIENT == 'M'], 'Female -> Male', 'FM', outdir.lab)
plot_CI_age_infection(pairs, outdir.lab)
phsc.plot.transmission.network(copy(as.data.table(dchain)), copy(as.data.table(dc)),outdir=outdir.lab, arrow=arrow(length=unit(0.02, "npc"), type="open"), edge.size = 0.1)

# extract incidence rate in rakai
incidence <- find_incidence_rate(range(df_age$age_infection.RECIPIENT))

# prepare stan data
stan_data <- prepare_stan_data(pairs, df_age, df_group)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3, 
                                      n_knots_rows = 8, n_knots_columns = 8, 
                                      AGES = unique(df_age$age_transmission.SOURCE))
stan_data <- add_prior_gp_mean(stan_data, df_age)

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.lab, paste0("stanin_",lab,".RData")) )
