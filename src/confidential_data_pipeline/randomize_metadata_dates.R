# load previous metadata and then substitute values for:
# date of birth
# date first postiive
# date last negative

library(data.table)
library(here)

gitdir <- here()
source(file.path(gitdir, 'config.R'))

file.exists( path.meta.confidential) |> stopifnot()


# load original meta data.

meta_env <- new.env()
load(path.meta.confidential, envir = meta_env)
ls(meta_env)

columns_of_interest <- c(
    "aid",
    "sex",
    "comm",
    "round",
    "sample_date",
    "date_birth",
    "date_first_positive",
    "date_last_negative"
)

meta <- subset(meta_env$meta_data, select = columns_of_interest)
meta <- unique(meta[!is.na(aid)])

# initialize new metadata data.table. 
# keep the round, sex and comm constant.
constant_cols <- c('aid', 'sex', 'round', 'comm')
meta_randomized <- subset(meta,select = constant_cols)
no_repeated_combinations <- meta_randomized[, .N == uniqueN(.SD)]
stopifnot(no_repeated_combinations)

cat("-- within rounds, shuffle the sample date. --\n")
cols <- c('aid', 'round', 'sample_date')
shuffled_visit_dates <- subset(meta, select=cols)[, 
    .(aid=aid, sample_date=sample_date[sample(1:.N)]) ,
by='round']
meta_randomized <- merge(meta_randomized, shuffled_visit_dates, by=c('aid','round'))
setkey(meta_randomized, aid, sex, round)

# check consistency
meta_randomized[round!='neuro', round_i := as.integer(gsub('[A-Z]', '', round))]
check.visits.in.different.rounds.are.successive <- meta_randomized[
    ! round == 'R015S', 
    {
        z <- sort(round_i, index.return=TRUE )$ix
        list(sorted=!is.unsorted(sample_date[z]))
    }, by='aid' ][, all(sorted)] 
stopifnot(check.visits.in.different.rounds.are.successive)
meta_randomized[, round_i := NULL]

# 
cat("-- shifts date of birth and serohistory uniformly by [-6, 6] months --\n")
variable_cols <- c('aid', 'date_birth', 'date_first_positive', 'date_last_negative')
ddates <- subset(meta,select=variable_cols) |> unique()
stopifnot(ddates[, uniqueN(aid) == .N ])
add_dob_noise <- function(N, min=-2.5, max=2.5)
    as.integer(runif( n=N, min=min*365, max=max*365))
ddates[, date_birth := date_birth + add_dob_noise(.N, min=.25, max=.25)]

cat("-- do the same for dates first positive and dates of last negative test --\n")
ddates[, date_first_positive  := date_first_positive + add_dob_noise(.N, min=.25)]
ddates[, date_last_negative  := date_last_negative + add_dob_noise(.N, max=.25)]

# check consistency
while(ddates[date_first_positive < date_last_negative, .N ])
{
    ddates[, date_first_positive := date_first_positive + add_dob_noise(.N, min=0)]
}

# merge back
meta_randomized <- merge(meta_randomized, ddates, by='aid')
# rename to 'meta'
meta_data <- copy(meta_randomized)

filename <- path.meta.randomized
if( !file.exists(file.name))
{
    cat('\n Careful: This data should already exist exist in ', file.name  )
    cat('\n check that your Zenodo path is correctly specified in config.R ' )
    cat('\nIf you wish to proceed, and save this file anyway run the commented line below')
    #     save(meta_data, file=filename)
}else{
    cat('\n Output file', file.name,'already exists.\n')
}

