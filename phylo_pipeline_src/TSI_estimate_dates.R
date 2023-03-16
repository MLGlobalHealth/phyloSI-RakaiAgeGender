cat('\n\n=====  TSI_estimate_dates.R =====\n\n')

# Objective of this script is to load every TSI prediction produced in the previous step
# And then to attribute estimates of dates of infectios.
# To this end, we require a dataframe containing the dates of collection of each sample.
library(data.table)
library(lubridate)
library(here)

option_list <- list(
  optparse::make_option(
    "--tsi_out_dir",
    type = "character",
    default = "/home/andrea/HPC/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T",
    help = "Absolute file path to directory containing phylo-TSI results", 
    dest = 'rel.dir'
  ),
  optparse::make_option(
    "--confidential",
    type = "logical",
    default = TRUE, 
    help = "Whether using the confidential sample collection dates or not",
    dest = 'confidential'
  )
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))

# if dataset with dates of collection exists, also compute date of infection estimates!
user <- Sys.info()[['user']]
use.randomized <- ! args$confidential

if(user != 'andrea')
{
    indir.deepseqdata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
}else{
    indir.deepseqdata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    # indir.deepseqanalyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live'
}

gitdir <- here()
# db.sharing.path.rccs <- file.path(indir.deepseqdata, 'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
# db.sharing.path.mrc <- file.path( indir.deepseqdata, 'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')

path.collection.dates <- fifelse(args$confidential, 
    yes=file.path(indir.deepseqdata, 'RCCS_R15_R18', 'sequences_collection_dates.rds'),
    no=file.path(gitdir, 'data', 'sequences_collection_dates_randomized.rds'),
)

################
# Testing
################

# if(user == 'andrea')
# {
#     args <- list(
#         confidential=FALSE,
#         rel.dir=file.path(indir.deepseqanalyses, 'PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_00'),
#         ps.samples=file.path(indir.deepseqanalyses, "PANGEA2_RCCS_MRC_UVRI_TSI/220331_RCCSUVRI_phscinput_samples_with_bf.rds"),
#         date='2022-08-22'
#     )
# }


################
# MAIN
################

stopifnot("Path to TSI-output-directory does not exist"=dir.exists(args$rel.dir))
# stopifnot(dir.exists(args$out.dir))

# aggregate outputs from HIV_phylo_tsi and prettify
dfiles <- list.files(args$rel.dir, pattern='_tsi.csv$', full.names = TRUE)
dpreds <- lapply(dfiles, fread)
.gs <- function(x) gsub('[A-z]|_|\\.', '', basename(x) )
names(dpreds) <- .gs(dfiles)
cols <- grep('host.id|^RF', colnames(dpreds[[1]]), value=T)
dpreds <- lapply(dpreds, function(DT) subset(DT, select=cols) )
dpreds <- rbindlist(dpreds, idcol = 'PTY')
dpreds[, `:=` (AID = gsub('-fq[0-9]$','',host.id), PTY=as.integer(PTY))]
setnames(dpreds, 'host.id', 'RENAME_ID')
if(dpreds[1, RENAME_ID == AID])
        dpreds[, RENAME_ID := NULL]

# Take median of predictions and cc's in sqrt space, then transform to linear space
# Alternatively could pick prediciton with minimum MAE
if( dpreds[, .N,by='RENAME_ID'][, any(N>1)] ) 
{
    cols <- grep('RF',colnames(dpreds), value=T)
    cols <- grep('linear', cols, value=T, invert = TRUE)
    dpreds <- dpreds[, lapply(.SD, median) ,by='RENAME_ID', .SDcols=cols]
    cols = grep('pred_sqrt|cc',cols, value=T)
    cols1 <- gsub('sqrt', 'linear',cols)
    cols1 <- gsub('cc025', 'pred_min_linear',cols1)
    cols1 <- gsub('cc975', 'pred_max_linear',cols1)
    dpreds [, (cols1):=lapply(.SD, function(x) x^2), .SDcols=cols, by='RENAME_ID']
}

# get sample collection dates and merge
ddates <- readRDS(path.collection.dates) 
names(ddates) <- toupper(names(ddates)) 
ddates  <- subset(ddates, select=c('RENAME_ID', 'VISIT_DT')) 
dpreds <- merge(dpreds, ddates, by="RENAME_ID")
dpreds[, PTY := NULL]

# Transform TSI estimates in day and then estimate DateofInfection
cols <- grep('linear', colnames(dpreds), value=T)
cols1 <- gsub('linear', 'days', cols)
tmp <- dpreds[, lapply(.SD, function(x) as.integer(x*365)) , .SDcols=cols, by=c('RENAME_ID', 'VISIT_DT')]
setnames(tmp, cols, cols1)
tmp <- tmp[, lapply(.SD, function(x){VISIT_DT - x}) , .SDcols=cols1, by='RENAME_ID']
setnames(tmp, grep('pred_min', names(tmp),value=T), 'pred_doi_max')
setnames(tmp, grep('pred_max', names(tmp),value=T), 'pred_doi_min')
setnames(tmp, grep('RF', names(tmp),value=T), 'pred_doi_mid')
dpreds <- merge(dpreds, tmp, by='RENAME_ID')

# Extract only columns of interest for ....  
cols <- c( "AID", "RENAME_ID",
    "RF_pred_sqrt", "RF_std", "RF_cc025", "RF_cc975", "RF_pred_MAE",
    "RF_pred_linear", "RF_pred_min_linear", "RF_pred_max_linear")
dpreds  <- subset(dpreds, select=cols) 
dtsi <- fread(filename, select = cols)


cat("\n---- Save in gitdir.data ----\n")
suffix <- ''
if(is.randomized)
    suffix <- '_randomized'
filename <- paste0("TSI_estimates", suffix,".csv")
filename <- file.path(gitdir, 'data', filename)

if(! file.exists(filename))
{
    cat('\nSaving',filename,'...\n')
    fwrite(dtsi, file=filename)
}else{
    cat('\nFile',filename,'already exists.\n')
}

cat('\n\n=====  end of script =====\n\n')
