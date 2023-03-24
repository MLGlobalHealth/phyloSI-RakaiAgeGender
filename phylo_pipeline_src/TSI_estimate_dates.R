cat('\n\n=====  TSI_estimate_dates.R =====\n\n')

# Objective of this script is to load every TSI prediction produced in the previous step
# And then to attribute estimates of dates of infectios.
# To this end, we require a dataframe containing the dates of collection of each sample.
library(data.table) |> suppressPackageStartupMessages() 
library(lubridate) |> suppressPackageStartupMessages()
library(here) |> suppressPackageStartupMessages()

option_list <- list(
  optparse::make_option(
    "--tsi-indir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to directory containing phylo-TSI results", 
    dest = 'rel.dir'
  ),
  optparse::make_option(
    "--confidential",
    type = "logical",
    default = TRUE, 
    help = "Whether using the confidential sample collection dates or not",
    dest = 'confidential'
  ),
  optparse::make_option(
    "--outdir",
    type = "character",
    default = NULL,
    help = "Absolute file path to output directory to save the scripts outputs in", 
    dest = 'outdir'
  )
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))

# if dataset with dates of collection exists, also compute date of infection estimates!
user <- Sys.info()[['user']]
is.randomized <- ! args$confidential

if(user != 'andrea')
{
    indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
}else{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
}

gitdir <- here()
source(file.path(gitdir, 'paths.R'))

# if output directory is null, set it to gitdir.data
if( is.null(args$outdir) )
    args$outdir <- gitdir.data

path.collection.dates <- fifelse(args$confidential, 
    yes=path.collection.dates.confidential,
    no=path.collection.dates.randomized)


################
# MAIN
################

stopifnot("Path to TSI-output-directory does not exist"=dir.exists(args$rel.dir))

# aggregate outputs from HIV_phylo_tsi and prettify
dfiles.csv <- list.files(args$rel.dir, pattern='_tsi.csv$', full.names = TRUE)
dfiles.rds <- list.files(args$rel.dir, pattern='_tsi.rds$', full.names = TRUE)

.gs <- function(x) gsub('[A-z]|_|\\.', '', basename(x) )
if( length(dfiles.csv) )
{
    cat('extracting',length(dfiles.csv),'.csv files\n')
    dpreds <- lapply(dfiles.csv, fread)
    names(dpreds) <- .gs(dfiles.csv)
}else if( length(dfiles.rds) ){
    cat('extracting',length(dfiles.rds),'.rds files\n')
    dpreds <- lapply(dfiles.rds, readRDS)
    names(dpreds) <- .gs(dfiles.rds)
}else{
    stop('no `*tsi.csv` or `tsi.rds` files found in args$rel.dir .\n')
}
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


cat("\n---- Save in outdir (default:gitdir.data) ----\n")
suffix <- ''
if(is.randomized)
    suffix <- '_randomized'
filename <- paste0("TSI_estimates", suffix,".csv")
filename <- file.path(args$outdir, filename)

if(! file.exists(filename))
{
    cat('\nSaving',filename,'...\n')
    fwrite(dpreds, file=filename)
}else{
    cat('\nFile',filename,'already exists.\n')
}

cat('\n\n=====  end of script =====\n\n')
