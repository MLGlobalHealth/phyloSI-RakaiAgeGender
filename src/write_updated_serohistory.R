library(data.table)
library(lubridate)

if(dir.exists('/Users/melodiemonod'))
{
  indir.repository <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
}

if(dir.exists('/home/andrea'))
{
    base.hpc <- '/home/andrea/HPC'
    # out.dir <- '/home/andrea/HPC/Documents/Box/2022/genintervals/'
    out.dir <- '/home/andrea/HPC/ab1820/home/projects/2022/genintervals'
    indir.deepsequence_xiaoyue   <- file.path(base.hpc, 'project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI')
    indir.deepsequencedata <- file.path(base.hpc, 'project/ratmann_pangea_deepsequencedata/live')
    indir.deepsequence_analyses   <- file.path(base.hpc, 'project/ratmann_deepseq_analyses/live')
    indir <- '/home/andrea/git/phyloflows'
}

if(dir.exists('/rds/general/user/'))
{
  indir.repository <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
}

################
#    HELPERS   #
################

source(file.path(indir.repository, 'functions', 'utils.R'))
source(file.path(indir.repository, 'functions', 'summary_functions.R'))

if(0)
{   
    # Joseph claims that 3 people never tested positive and 2 are seroreverters, althoguh they all appear in the pangea samples
    file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
    aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))
    ddates <- get.sample.collection.dates()
    names(ddates) <- toupper(names(ddates))
    ddates <- merge(ddates, aik)

    # however they do not appear in the resulting pairs
    filename <-  file.path(indir.deepsequencedata, "RCCS_R15_R18/pairsdata_toshare_d1_w11_netfrompairs_seropairs.rds")
    dresults <- readRDS( filename)    
    idx <- ddates[ PT_ID %in% paste0('RK-', c('B077943', 'B106184', 'J010515', 'A008742', 'E118889')), AID]
    dresults[SOURCE %in% idx]
    dresults[RECIPIENT %in% idx]
}

new_dates <- data.table(study_id = paste0('RK-', c('K067249','C117824', 'C119303')),
                        lastnegdat=as.Date(c('2010-07-09', NA, NA)), 
                        firstposdat=as.Date(c('2014-10-02', '2014-09-18', '2015-01-08')))

filename <- file.path( indir.deepsequencedata, 'RCCS_R15_R18/221128_requested_updated_serohistory.csv')
fwrite(new_dates, filename)
