require(data.table)
require(here)

#########
# PATHS # 
#########

usr <- Sys.info()[['user']]
indir <- here::here()

if(usr == 'andrea')
{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepanalyses_xiaoyue <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
}

# TODO: better on HPC, as anyways do not push this
outdir.confidential <- file.path(indir.deepsequencedata, 'RCCS_R15_R18')
outdir.data <- file.path(indir, 'data')

# path from phylo analyses
path.anonymisation.keys <- file.path(indir.deepanalyses_xiaoyue, 'important_anonymisation_keys_210119.csv')
path.selected.samples <- file.path(indir.deepanalyses_xiaoyue,"210120_RCCSUVRI_phscinput_samples.rds" )

# path from PANGEA data
path.sdates.rccs <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')
path.sdates.mrc <- file.path(indir.deepsequencedata, 'PANGEA2_MRC','200319_pangea_db_sharing_extract_mrc.csv')

###########
# HELPERS #
###########

catn <- function(x) cat('\n----  ', x, '  ----\n')

perturb_dates <- function(date, .min=-1, .max=1)
{
    N <- length(date)
    noise <- as.integer(runif(n=N, min=.min, max=.max)*365)
    return(date + noise)
}

##############
# READ FILES #
##############

catn("Get Anonimization Keys")
aik <- fread(path.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))
names(aik) <- tolower(names(aik)) 

catn("Get selected samples for Xiaoyue's analysis")
dsamples <- readRDS(path.selected.samples) |>
    subset(select=c('PANGEA_ID', 'RENAME_ID'))
names(dsamples) <- tolower(names(dsamples))
dsamples[, aid := gsub( '-fq[0-9]$', '', rename_id) ]
dsamples[, pangea_id := gsub( '^(.*)_', '',pangea_id) ]
stopifnot(dsamples[, uniqueN(pangea_id) == .N])
dsamples <- merge(aik, dsamples, by='aid')

catn("Get Blood sample ids and visit dates")
files <- c(path.sdates.rccs, path.sdates.mrc)
cols <- c('pt_id', 'pangea_id', 'visit_dt')
ddates <- rbindlist(lapply(files, fread, select=cols))
ddates <- unique(ddates) |>
    merge(dsamples, by='pangea_id') |>
    subset(select=c('rename_id','aid' ,'visit_dt'))

catn("Save the original in encrypted data folder")
filename <- file.path(outdir.confidential, 'sequences_collection_dates.rds')
saveRDS(ddates, file=filename)

catn("Perturb dates in (-3, 3) months")
ddates[, visit_dt_perturbed := perturb_dates(visit_dt, -.25, .25) ]
ddates_perturbed <- ddates[, .(
    aid = aid, 
    rename_id = rename_id,
    visit_dt = visit_dt_perturbed
)]

catn("Save the perturbed dates in the encrypted data folder")
filename <- file.path(outdir.data, 'sequences_collection_dates_randomized.rds')
saveRDS(ddates_perturbed, file=filename)
