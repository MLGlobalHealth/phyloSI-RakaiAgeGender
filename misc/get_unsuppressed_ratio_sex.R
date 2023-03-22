library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "paths.R"))

# directory to repository
gitdir <- getwd()

# files
file.exists(file.treatment.cascade) |> stopifnot()
file.exists(file.prevalence) |> stopifnot()
file.exists(file.eligible.count ) |> stopifnot()

# load census eligible ount
eligible_count <- fread(file.eligible.count)

# load proportion prevalence
proportion_prevalence <- as.data.table(readRDS(file.prevalence))

# load unsuppressed proportion 
treatment_cascade <- as.data.table(readRDS(file.treatment.cascade))


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

# UNSUPPRESSED BY 3 AGE GROUP

#####################################################

# aggregated by age groups
df[, INDEX_AGE_GROUP := 3]
df[AGEYRS < 35, INDEX_AGE_GROUP := 2]
df[AGEYRS < 25, INDEX_AGE_GROUP := 1]
df[, AGE_GROUP := c('15-24', '25-34', '35-49')[INDEX_AGE_GROUP]]
df.agg <- df[, list(INFECTED = sum(INFECTED), 
                UNSUPPRESSED = sum(UNSUPPRESSED), 
                ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'COMM', 'AGE_GROUP', 'SEX', 'iterations')]

# find unsuppressed rate
df.agg[, SUPPRESSED := INFECTED - UNSUPPRESSED]
df.agg[, SUPPRESSION_RATE := SUPPRESSED / INFECTED]
df.agg[, UNSUPPRESSION_RATE := UNSUPPRESSED / INFECTED]

# find unsuppressed rate relative to round 10
df.agg[COMM == 'inland', UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == '10'], by = c('COMM', 'AGE_GROUP', 'SEX', 'iterations')]
df.agg[COMM == 'fishing', UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == '15'], by = c('COMM', 'AGE_GROUP', 'SEX', 'iterations')]
df.agg[, UNSUPPRESSION_RATE_REL := UNSUPPRESSION_RATE / UNSUPPRESSION_RATE.REF]

# find ratio of art uptake male to female
df.agg <- dcast(df.agg, ROUND + COMM + AGE_GROUP + iterations ~ SEX, value.var = 'UNSUPPRESSION_RATE_REL')
setnames(df.agg, c('M', 'F'), c('UNSUPPRESSION_RATE_REL_M', 'UNSUPPRESSION_RATE_REL_F'))
df.agg[, UNSUPPRESSION_RATE_RATIO := UNSUPPRESSION_RATE_REL_M / UNSUPPRESSION_RATE_REL_F ]

# summarise
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.age = df.agg[, list(q= quantile(UNSUPPRESSION_RATE_RATIO, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'AGE_GROUP')]
sing.age = as.data.table(reshape2::dcast(sing.age, ... ~ q_label, value.var = "q"))

# name
setnames(sing.age, qlab, paste0('UNSUPPRESSION_RATE_RATIO_BY_AGE_', qlab))

# plot
ggplot(sing.age, aes(x = ROUND)) + 
  geom_point(aes(y = UNSUPPRESSION_RATE_RATIO_BY_AGE_M)) + 
  geom_errorbar(aes(ymin = UNSUPPRESSION_RATE_RATIO_BY_AGE_CL, ymax = UNSUPPRESSION_RATE_RATIO_BY_AGE_CU), alpha = 0.5) +
  facet_grid(COMM~AGE_GROUP) + 
  theme_bw()


#####################################################

# FIND SEX RATIO OF INFECTED UNSUPPRESSED ACROSS AGE

#####################################################

df.agg <- df[, list(INFECTED = sum(INFECTED), 
                    UNSUPPRESSED = sum(UNSUPPRESSED), 
                    ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'COMM', 'SEX', 'iterations')]

# find unsuppressed rate
df.agg[, SUPPRESSED := INFECTED - UNSUPPRESSED]
df.agg[, SUPPRESSION_RATE := SUPPRESSED / INFECTED]
df.agg[, UNSUPPRESSION_RATE := UNSUPPRESSED / INFECTED]

# find unsuppressed rate relative to round 10
df.agg[COMM == 'inland', UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == '10'], by = c('COMM', 'SEX', 'iterations')]
df.agg[COMM == 'fishing', UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == '15'], by = c('COMM', 'SEX', 'iterations')]
df.agg[, UNSUPPRESSION_RATE_REL := UNSUPPRESSION_RATE / UNSUPPRESSION_RATE.REF]

# find ratio of art uptake male to female
df.agg <- dcast(df.agg, ROUND + COMM + iterations ~ SEX, value.var = 'UNSUPPRESSION_RATE_REL')
setnames(df.agg, c('M', 'F'), c('UNSUPPRESSION_RATE_REL_M', 'UNSUPPRESSION_RATE_REL_F'))
df.agg[, UNSUPPRESSION_RATE_RATIO := UNSUPPRESSION_RATE_REL_M / UNSUPPRESSION_RATE_REL_F ]

# summarise
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing = df.agg[, list(q= quantile(UNSUPPRESSION_RATE_RATIO, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM')]
sing = as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

# name
setnames(sing, qlab, paste0('UNSUPPRESSION_RATE_RATIO_RATIO_', qlab))

# plot
ggplot(sing, aes(x = ROUND)) + 
  geom_point(aes(y = UNSUPPRESSION_RATE_RATIO_RATIO_M)) + 
  geom_errorbar(aes(ymin = UNSUPPRESSION_RATE_RATIO_RATIO_CL, ymax = UNSUPPRESSION_RATE_RATIO_RATIO_CU), alpha = 0.5) + 
  facet_grid(COMM~.) + 
  theme_bw()

#########################################

# SAVE

#########################################

tmp <- merge(sing.age, sing, by=c('ROUND', 'COMM'))
# file.name <- file.path(gitdir.fit,'RCCS_unsuppressed_ratio_sex_221208.csv')
file.name <- file.unsuppressed_rate_ratio 
if(! file.exists(file.name))
{
    cat("Saving file:", file.name, '\n')
    write.csv(tmp, file = file.name, row.names = F)
}else{
    cat("File:", file.name, "already exists...\n")
}
