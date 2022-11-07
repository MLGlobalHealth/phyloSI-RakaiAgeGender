library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

file.treatment.cascade <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_treatment_cascade_population_posterior_samples_221101.csv'))
file.prevalence <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('RCCS_prevalence_posterior_sample_220818.csv'))
file.eligible.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220830.csv')

# load census eligible ount
eligible_count <- as.data.table(read.csv(file.eligible.count))

# load proportion prevalence
proportion_prevalence <- as.data.table(read.csv(file.prevalence))

# load unsuppressed proportion 
treatment_cascade <- as.data.table(read.csv(file.treatment.cascade))


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

# ART UPTAKE BY 3 AGE GROUP

#####################################################

# aggregated by age groups
df[, INDEX_AGE_GROUP := 3]
df[AGEYRS < 35, INDEX_AGE_GROUP := 2]
df[AGEYRS < 25, INDEX_AGE_GROUP := 1]
df[, AGE_GROUP := c('15-24', '25-34', '35-49')[INDEX_AGE_GROUP]]
df.agg <- df[, list(INFECTED = sum(INFECTED), 
                UNSUPPRESSED = sum(UNSUPPRESSED), 
                ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'COMM', 'AGE_GROUP', 'SEX', 'iterations')]

# find ART uptake
df.agg[, SUPPRESSED := INFECTED - UNSUPPRESSED]
df.agg[, SUPPRESSION_RATE := SUPPRESSED / INFECTED]
df.agg[, UNSUPPRESSION_RATE := UNSUPPRESSED / INFECTED]

# find ratio of art uptake male to female
df.agg <- dcast(df.agg, ROUND + COMM + AGE_GROUP + iterations ~ SEX, value.var = 'UNSUPPRESSION_RATE')
setnames(df.agg, c('M', 'F'), c('UNSUPPRESSION_RATE_M', 'UNSUPPRESSION_RATE_F'))
df.agg[, UNSUPPRESSION_RATE_RATIO := UNSUPPRESSION_RATE_M / UNSUPPRESSION_RATE_F ]

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

# find ART uptake
df.agg[, SUPPRESSED := INFECTED - UNSUPPRESSED]
df.agg[, SUPPRESSION_RATE := SUPPRESSED / INFECTED]
df.agg[, UNSUPPRESSION_RATE := UNSUPPRESSED / INFECTED]

# find ratio of art uptake male to female
df.agg <- dcast(df.agg, ROUND + COMM + iterations ~ SEX, value.var = 'UNSUPPRESSION_RATE')
setnames(df.agg, c('M', 'F'), c('UNSUPPRESSION_RATE_M', 'UNSUPPRESSION_RATE_F'))
df.agg[, UNSUPPRESSION_RATE_RATIO := UNSUPPRESSION_RATE_M / UNSUPPRESSION_RATE_F ]

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
file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_unsuppressed_ratio_sex_221101.csv'))
write.csv(tmp, file = file.name, row.names = F)



