# load previous metadata and then substitute values for:
# date of birth
# date first postiive
# date last negative

library(data.table)
library(here)

indir <- here()
indir.data <- file.path(indir, 'data')
usr <- Sys.info()[['user']]
if(usr=='andrea')
{
    indir.deepsequencedata <- "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live"
}

file.path.meta.original <- file.path(
    indir.deepsequencedata,
    "RCCS_R15_R18",
    "Rakai_Pangea2_RCCS_Metadata_20221128.RData"
)

file.exists(
    file.path.meta.original,
    indir,
    indir.data
) |> all() |> stopifnot()


# load original meta data.

meta_env <- new.env()
load(file.path.meta.original, envir = meta_env)
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
cat("-- shifts date of birth uniformly by [-2.5, 2.5] years --\n")
variable_cols <- c('aid', 'date_birth', 'date_first_positive', 'date_last_negative')
ddates <- subset(meta,select=variable_cols) |> unique()
stopifnot(ddates[, uniqueN(aid) == .N ])
add_dob_noise <- function(N, min=-2.5, max=2.5)
    as.integer(runif( n=N, min=min*365, max=max*365))
ddates[, date_birth := date_birth + add_dob_noise(.N)]

cat("-- do the same for dates first positive and dates of last negative test --\n")
ddates[, date_first_positive  := date_first_positive + add_dob_noise(.N)]
ddates[, date_last_negative  := date_last_negative + add_dob_noise(.N)]

# check consistency
while(ddates[date_first_positive < date_last_negative, .N ])
{
    ddates[, date_first_positive := date_first_positive + add_dob_noise(.N, min=0)]
}

# merge back
meta_randomized <- merge(meta_randomized, ddates, by='aid')
# rename to 'meta'
meta_data <- copy(meta_randomized)

filename <- file.path(indir.data, 'Rakai_Pangea2_RCCS_Metadata_randomized.RData')
cat("Saving",filename, '\n')
save(meta_data, file=filename)
