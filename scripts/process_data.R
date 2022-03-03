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
require(lubridate)

# change as appropriate
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir.repository <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
  outdir <- '~/Box\ Sync/2021/phyloflows/'
}

if(dir.exists('/home/andrea'))
{
  indir.repository <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata/'
  outdir <- '~/Documents/Box/2021/phyloflows'
  file.path(indir.deepsequencedata,'RCCS_R15_18')
}

# indicators 
include.mrc <- F
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.5
use.tsi.estimates <- T
remove.inconsistent.infection.dates <- F
make.TSI.linear.adjustment <- F
remove.young.individuals <- T

cutoff_date <- as.Date('2014-01-01')
jobname <- '2014_IpriorGP'
lab <- paste0('OnlyHTX_', include.only.heterosexual.pairs, '_threshold_', threshold.likely.connected.pairs, '_jobname_', jobname)


# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
# file.path.meta.data.rccs.1 <- file.path(indir.deepsequence_analyses, 'RakaiPangeaMetaData_v2.rda')
# file.path.meta.data.rccs.2 <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.time.first.positive <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '211111_pangea_db_sharing_extract_rakai_age_firstpos_lastneg.csv')
file.path.phscinput <- file.path(indir.deepsequence_analyses, '210120_RCCSUVRI_phscinput_runs.rds')
file.path.bflocs <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', 'bfloc2hpc_20220103.rds')
file.path.tsiestimates <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', 'TSI_estimates_220119.csv')

outdir.lab <- file.path(outdir, lab); dir.create(outdir.lab)

# Latest data from Rakai's CCS (Joseph's data)
file.path.allhiv <- file.path(indir.deepsequencedata, 'RCCS_R15_18', 'All_HIVpcr_for_questR15_R18_220129.csv')
file.path.flow <- file.path(indir.deepsequencedata, 'RCCS_R15_18', 'FlowR15_R18_VoIs_220129.csv')
file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_R15_18', 'HIV_R15_R18_VOIs_220129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_R15_18', 'quest_R15_R18_VoIs_220129.csv')

# load functions
source(file.path(indir.repository, 'functions', 'summary_functions.R'))
source(file.path(indir.repository, 'functions', 'plotting_functions.R'))
source(file.path(indir.repository, 'functions', 'stan_utils.R'))
source(file.path(indir.repository, 'functions', 'check_potential_TNet.R'))

# load chains
load(file.path.chains.data)
dchain <- as.data.table(dchain)

# Load data
aik <- .read(file.anonymisation.keys); aik$X <- NULL
colnames(aik) <- tolower(colnames(aik))
community.keys <-.read(file.community.keys)

allhiv <- .read(file.path.allhiv) 
flow <- .read(file.path.flow)
hiv <- .read(file.path.hiv)
quest <- .read(file.path.quest)

# Add 'RK-' to all data.tables study_ids
invisible(lapply(list(hiv, allhiv, quest, flow),
                 function(dt) {
                   dt[!grepl('RK-', study_id), study_id := paste0('RK-', study_id)]
                 }))


# process quest and make date.birth
quest <- process.quest(quest)
date.birth <- quest[,list(birthdat=unique(birthdat)),by='study_id']

# process allhiv and make time.first.positive
allhiv <- process.allhiv(allhiv)

# process hiv and check time.first.positive obtained as such:
hiv <- process.hiv(hiv)

# what s the difference between hiv and allhiv?
# hiv probably contains all tests, while allhiv only those from people ever tested positive?
# There are 6 entries with different first positive diagnoses between hiv and allhivl
date.first.positive <- compare.hiv.allhiv.firstpositivedates(hiv, allhiv)
stopifnot(.vars.with.multiple.values(date.first.positive, 'study_id')[, .N == 0])
date.first.last.visit <- make.date.first.last.visit(hiv)

# load Tanya's estimate time since infection using phylogenetic data
time.since.infection <- as.data.table(read.csv(file.path.tsiestimates))
time.since.infection <- make.time.since.infection(time.since.infection)

# get meta data
meta_data <- get.meta.data(quest, aik, date.first.positive, time.since.infection)

# get likely transmission pairs
chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)

# study whether couple found in previous analyses are included in the same Potential Transmission Network cluster.
# commented by Melodie: this gives me an error
# print.statements.about.potential.TNet()

# merge meta data to source and recipienx
pairs.all <- pairs.get.meta.data(chain, meta_data)

# 99 study id NOT in Joseph Data. All appear in meta.rccs.2 and in 66% in meta.rccs.1
# But none of them appears in any new datesets (quest, hiv, allhiv and flow)
# tmp <- c("RK-A066261", "RK-J055504", "RK-C054888", "RK-J068023", "RK-H065970", "RK-H133273", "RK-H132927", "RK-G063880", "RK-J054622", "RK-B132191", "RK-B095525", "RK-F005239", "RK-D050431", "RK-F115356", "RK-K059296", "RK-A133666", "RK-K132741", "RK-F100447", "RK-J069565", "RK-A094422", "RK-D021259", "RK-C094814", "RK-B067364", "RK-K100417", "RK-K066707", "RK-B038263", "RK-C038434", "RK-J133656", "RK-C051281", "RK-H030219", "RK-H046775", "RK-C100842", "RK-A132124", "RK-H133061", "RK-D115354", "RK-C073172", "RK-E102438", "RK-J133698", "RK-D131863", "RK-K094560", "RK-C056943", "RK-F101418", "RK-G100383", "RK-C133359", "RK-F061554", "RK-E115410", "RK-B046282", "RK-K027734", "RK-G054078", "RK-D131588", "RK-H039230", "RK-B063380", "RK-K058835", "RK-D039919", "RK-H096194", "RK-C132388", "RK-G131743", "RK-K100823", "RK-H061522", "RK-E100111", "RK-J188436", "RK-K094633", "RK-D095072", "RK-J036456", "RK-K088460", "RK-J067882", "RK-B009138", "RK-K096293", "RK-F026612", "RK-F038893", "RK-A100696", "RK-H045013", "RK-H096209", "RK-F059276", "RK-E115389", "RK-H068933", "RK-B095818", "RK-G055934", "RK-K102776", "RK-J096771", "RK-C100834", "RK-D067324", "RK-E053649", "RK-E062736", "RK-F064949", "RK-D102788", "RK-B195356", "RK-E036460", "RK-B002793", "RK-D066117", "RK-E096442", "RK-J005965", "RK-C101596", "RK-E025226", "RK-C064035", "RK-K091620", "RK-K004814", "RK-A066774", "RK-E031889")

# MRC do not appear on meta_data so no need to add anything here
# if(!include.mrc){
#   cat('Keep only pairs in RCCS\n')
#   cat('removing ', nrow(pairs.all[cohort.RECIPIENT == 'MRC' | cohort.SOURCE == 'MRC']), ' pairs\n')
#   pairs.all <- pairs.all[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS']
#   cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
# }
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

print.which.NA(pairs.all,name='pairs.all')
print.statements.about.pairs(copy(pairs.all), outdir.lab)

# which base frequency files we have on the HPC
# atm gives error: maybe TODO when I understand more about PHSC pipeline
# missing_bff <- print.statements.about.basefreq.files(pairs.all)


# keep only pairs with source-recipient with proxy for the time of infection
pairs <- pairs.all[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
pairs[, date_infection_before_cutoff.RECIPIENT := date_infection.RECIPIENT < cutoff_date]



#
# MAKE EXPLANATORY PLOTS
#

# prepare age map
df_age <- get.age.map(pairs, age_bands_reduced = 4)
df_age

# prepare group map
df_group <- get.group.map()

# Purely for the sake of making plots:
pairs[, cohort_round.SOURCE := "R15-18"]
pairs[, cohort_round.RECIPIENT := "R15-18"]

# make some explanatory plots
plot_hist_age_infection(copy(pairs), outdir.lab)
plot_hist_time_infection(copy(pairs), cutoff_date, outdir.lab)
plot_age_infection_source_recipient(pairs[sex.SOURCE == 'M' & sex.RECIPIENT == 'F'], 'Male -> Female', 'MF', outdir.lab)
plot_age_infection_source_recipient(pairs[sex.SOURCE == 'F' & sex.RECIPIENT == 'M'], 'Female -> Male', 'FM', outdir.lab)
plot_CI_age_infection(pairs, outdir.lab)
plot_CI_age_transmission(pairs, outdir.lab)
phsc.plot.transmission.network(copy(as.data.table(dchain)), copy(as.data.table(dc)), pairs,outdir=outdir.lab, arrow=arrow(length=unit(0.02, "npc"), type="open"), edge.size = 0.1)



#
# PREPARE STAN DATA
#

# extract incidence rate in rakai from https://www.rhsp.org/research/rccs/explore-rccs-data
incidence <- find_incidence_rate(range(df_age$age_infection.RECIPIENT))

# prepare stan data
stan_data <- prepare_stan_data(pairs, df_age, df_group)
stan_data <- add_2D_splines_stan_data(stan_data, spline_degree = 3, 
                                      n_knots_rows = 8, n_knots_columns = 8, 
                                      X = unique(df_age$age_transmission.SOURCE),
                                      Y = unique(df_age$age_infection.RECIPIENT))
stan_data <- add_prior_gp_mean(stan_data, df_age, outdir.lab)

## save image before running Stan
tmp <- names(.GlobalEnv)
tmp <- tmp[!grepl('^.__|^\\.|^model$',tmp)]
save(list=tmp, file=file.path(outdir.lab, paste0("stanin_",lab,".RData")) )


