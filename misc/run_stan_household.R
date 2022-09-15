library(data.table)
library(ggplot2)

# laptop
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
  outdir <- '~/Box\ Sync/2021/phyloflows/'

  jobname <- 'test_new'
  stan_model <- 'gp_220905a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/home/andrea'))
{
  indir <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata/'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'test'
  stan_model <- 'gp_220901a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/rds/general'))
{
  indir.deepsequence_analyses   <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
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

# TODO: set an output directory you want
# outfile <- file.path(outdir, paste0(stan_model,'-', jobname))
# outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
# outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
# if(!dir.exists(dirname(outfile.figures))) dir.create(dirname(outfile.figures))
# if(!dir.exists(dirname(outdir.table))) dir.create(dirname(outdir.table))

# indicators
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.5
use.tsi.estimates <- F
remove.inconsistent.infection.dates <- F
remove.young.individuals <- T
remove.missing.community.recipient <- T

remove.neuro.individuals <- T
only.transmission.after.start.observational.period <- T
only.transmission.before.stop.observational.period <- T
only.transmission.same.community <- F

# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')


# from misc/
file.path.meta <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_20220329.RData')

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))

# load anonymous aid TODO: add file.anonymisation.keys path
# aik <- .read(file.anonymisation.keys); aik$X <- NULL

# load chains
load(file.path.chains.data)
dchain <- as.data.table(dchain)
dchain

load(file.path.meta)
meta_data



if(0)
{
        # NEED TO GET PAIRS pairs.all

        if(only.transmission.after.start.observational.period){
          cat('\nFor inland excluding recipients infected before ', as.character(start_observational_period_inland), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland')]
          
          cat('\nFor fishing excluding recipients infected before ', as.character(start_observational_period_fishing), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT < start_observational_period_fishing & COMM.RECIPIENT == 'fishing']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT < start_observational_period_fishing & COMM.RECIPIENT == 'fishing')]
          
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        }
        if(only.transmission.before.stop.observational.period){
          cat('\nFor inland excluding recipients infected after ', as.character(stop_observational_period_inland), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland')]
          
          cat('\nFor fishing excluding recipients infected after ', as.character(stop_observational_period_fishing), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT > stop_observational_period_fishing & COMM.RECIPIENT == 'fishing']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT > stop_observational_period_fishing & COMM.RECIPIENT == 'fishing')]
          
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        }
        if(remove.missing.community.recipient){
          cat('\nExcluding recipients without community \n')
          cat('Removing ', nrow(pairs.all[is.na(COMM.RECIPIENT)]), ' pairs\n')
          pairs.all <- pairs.all[!is.na(COMM.RECIPIENT)]
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        } 
        if(only.transmission.same.community ){
          cat('\nExcluding transmission events between communities (I->F or F->I) \n')
          cat('Removing ', nrow(pairs.all[COMM.SOURCE != COMM.RECIPIENT]), ' pairs\n')
          pairs.all <- pairs.all[COMM.SOURCE == COMM.RECIPIENT]
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        }
}


