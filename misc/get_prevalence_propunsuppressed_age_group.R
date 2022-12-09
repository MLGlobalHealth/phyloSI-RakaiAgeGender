library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

# directory repository
indir.repository <- '~/git/phyloflows'

# files
file.treatment.cascade <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_population_posterior_samples_221208.rds'))
file.prevalence <- file.path(indir.repository, 'fit', paste0('RCCS_prevalence_posterior_sample_221116.rds'))
file.eligible.count <- file.path(indir.repository, 'data', 'RCCS_census_eligible_individuals_221116.csv')

# load census eligible ount
eligible_count <- as.data.table(read.csv(file.eligible.count))

# load proportion prevalence
proportion_prevalence <- as.data.table(readRDS(file.prevalence))

# load unsuppressed proportion 
treatment_cascade <- as.data.table(readRDS(file.treatment.cascade))


#############################

# FIND INFECTED UNSUPPRESSED

############################

# define round
df <- merge(proportion_prevalence, treatment_cascade, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX', 'iterations'), all.x = T)
df[, ROUND := gsub('R0(.+)', '\\1', ROUND)]

# merge number of eligible and the prevalence
df <- merge(eligible_count[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)], df, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))

# find infected
df[, INFECTED := ELIGIBLE * PREVALENCE_POSTERIOR_SAMPLE]

# find infected unsuppressed
df[, PROP_UNSUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_POSTERIOR_SAMPLE]
df[, UNSUPPRESSED := INFECTED * PROP_UNSUPPRESSED_POSTERIOR_SAMPLE]


#####################################################

# FIND PROPORTION UNSUPPRESSED BY AGE GROUP

#####################################################

# age groups
df_age_aggregated <- data.table(AGEYRS = df[, sort(unique(AGEYRS))])
df_age_aggregated[, INDEX_AGE_GROUP := 3]
df_age_aggregated[AGEYRS < 35, INDEX_AGE_GROUP := 2]
df_age_aggregated[AGEYRS < 25, INDEX_AGE_GROUP := 1]
df_age_aggregated[, AGE_GROUP := c('15-24', '25-34', '35-49')[INDEX_AGE_GROUP]]

# aggregated by age groups
df <- merge(df, df_age_aggregated, by = 'AGEYRS')
df <- df[, list(INFECTED = sum(INFECTED), 
                UNSUPPRESSED = sum(UNSUPPRESSED), 
                ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'COMM', 'AGE_GROUP', 'SEX', 'iterations')]

# find prevalence by age groups across age
df[, PROP_UNSUPPRESSED := UNSUPPRESSED / ELIGIBLE]
df[, PREVALENCE := INFECTED / ELIGIBLE]

# summarise unsuppressed
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.age = df[, list(q= quantile(na.omit(PROP_UNSUPPRESSED), prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')]
sing.age = as.data.table(reshape2::dcast(sing.age, ... ~ q_label, value.var = "q"))
setnames(sing.age, qlab, paste0('PROP_UNSUPPRESSED_AGE_GROUP_', qlab))

# summarise prevalence
sinf.age = df[, list(q= quantile(PREVALENCE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')]
sinf.age = as.data.table(reshape2::dcast(sinf.age, ... ~ q_label, value.var = "q"))
setnames(sinf.age, qlab, paste0('PREVALENCE_AGE_GROUP_', qlab))

# plot
ggplot(sing.age, aes(x = AGE_GROUP, group = SEX)) + 
  geom_errorbar(aes(ymin = PROP_UNSUPPRESSED_AGE_GROUP_CL, ymax = PROP_UNSUPPRESSED_AGE_GROUP_CU), alpha = 0.5, width = 0.3, position=position_dodge(width = 0.3)) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_AGE_GROUP_M, col = SEX), position=position_dodge(width = 0.3)) + 
  facet_grid(ROUND~COMM) + 
  theme_bw()

ggplot(sinf.age, aes(x = AGE_GROUP, group = SEX)) + 
  geom_errorbar(aes(ymin = PREVALENCE_AGE_GROUP_CL, ymax = PREVALENCE_AGE_GROUP_CU), alpha = 0.5, width = 0.3, position=position_dodge(width = 0.3)) + 
  geom_point(aes(y = PREVALENCE_AGE_GROUP_M, col = SEX), position=position_dodge(width = 0.3)) + 
  facet_grid(ROUND~COMM) + 
  theme_bw()

#########################################

# SAVE

#########################################

file.name <- file.path(indir.repository, 'fit', paste0('RCCS_propunsuppressed_age_group_221208.csv'))
write.csv(sing.age, file = file.name, row.names = F)

file.name <- file.path(indir.repository, 'fit', paste0('RCCS_prevalence_age_group_221116.csv'))
write.csv(sinf.age, file = file.name, row.names = F)

