library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "paths.R"))

# files
file.exists(file.treatment.cascade) |> stopifnot()
file.exists(file.prevalence) |> stopifnot()
file.exists(file.eligible.count) |> stopifnot()

# load census eligible ount
eligible_count <- fread(file.eligible.count)

# load proportion prevalence
proportion_prevalence <- as.data.table(readRDS(file.prevalence))

# load unsuppressed proportion 
treatment_cascade <- as.data.table(readRDS(file.treatment.cascade))

ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')

#############################

# FIND INFECTED UNSUPPRESSED

############################

# define round
df <- merge(proportion_prevalence, treatment_cascade, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX', 'iterations'))
df[, ROUND := gsub('R0(.+)', '\\1', ROUND)]

# merge number of eligible and the prevalence
df <- merge(eligible_count[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)], df, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))

# find infected
df[, INFECTED := ELIGIBLE * PREVALENCE_POSTERIOR_SAMPLE]

# find infected unsuppressed
df[, PROP_UNSUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_POSTERIOR_SAMPLE]
df[, UNSUPPRESSED := INFECTED * PROP_UNSUPPRESSED_POSTERIOR_SAMPLE]


#####################################################

# FIND MEDIAN AGE OF INFECTED UNSUPPRESSED 

#####################################################

# find share of unsuppressed by sex across age
df[, TOTAL_UNSUPPRESSED := sum(UNSUPPRESSED), by = c('ROUND', 'COMM', 'SEX', 'iterations')]
df[, UNSUPPRESSED_SHARE := UNSUPPRESSED / TOTAL_UNSUPPRESSED, by = c('ROUND', 'COMM', 'AGEYRS','SEX', 'iterations')]

# find qunatile of unsuppressed
df <- df[, list(value = Hmisc::wtd.quantile(x = AGEYRS, weight = UNSUPPRESSED_SHARE,
                                                probs = c(0.1, 0.25, 0.5, 0.75, 0.9), normwt = TRUE),
                    quantile = c('C10', 'C25', 'C50', 'C75', 'C90')), by = c('iterations', 'ROUND', 'COMM', 'SEX')]
# summarise
sing = df[, list(q= quantile(value, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'quantile')]
sing = as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

#########################################

# SAVE

#########################################

# file.name <- file.path(gitdir.fit,'RCCS_unsuppressed_median_age_221208.csv')
file.name <- file.unsuppressed_median_age  
if(! file.exists(fil.name))
{
    cat("Saving file:", file.name, '\n')
    write.csv(sing, file = file.name, row.names = F)
}else{
    cat("File:", file.name, "already exists...\n")
}
