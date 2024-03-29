require(data.table)
require(here)

gitdir <- here()
source(file.path(gitdir, 'config.R'))

file.exists(c(
    file.anonymisation.keys ,
    path.selected.samples ,
    path.sdates.rccs ,
    path.sdates.mrc ))  |> all() |> stopifnot()

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
aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))
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
filename <- path.collection.dates.confidential
if(! file.exists(filename))
{
    cat("Saving file:", filename, '\n')
    saveRDS(ddates, file=filename)
}else{
    cat("File:", filename, "already exists...\n")
}

catn("Perturb dates in (-3, 3) months")
ddates[, visit_dt_perturbed := perturb_dates(visit_dt, -.25, .25) ]
ddates_perturbed <- ddates[, .(
    aid = aid, 
    rename_id = rename_id,
    visit_dt = visit_dt_perturbed
)]

catn("Save the perturbed dates in the zenodo dir")
filename <- path.collection.dates.randomized 
if(! file.exists(filename))
{
    cat("Saving file:", filename, '\n')
    saveRDS(ddates_perturbed, file=filename)
}else{
    cat("File:", filename, "already exists...\n")
}
