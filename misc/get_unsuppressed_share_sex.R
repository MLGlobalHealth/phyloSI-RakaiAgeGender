library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

# directory to repository
indir.repository <- '~/git/phyloflows'

# files
file.treatment.cascade <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_population_posterior_samples_221116.rds'))
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

# FIND SEX SHARE OF INFECTED UNSUPPRESSED ACROSS AGE

#####################################################

# find share of unsuppressed by sex across age
df[, TOTAL_UNSUPPRESSED := sum(UNSUPPRESSED), by = c('ROUND', 'COMM', 'iterations')]
df[, UNSUPPRESSED_SHARE := UNSUPPRESSED / TOTAL_UNSUPPRESSED, by = c('ROUND', 'COMM', 'AGEYRS', 'iterations')]

df[COMM == 'inland' & ROUND == '18' & iterations == 1 & UNSUPPRESSED_SHARE > 0.05]
ggplot(df[COMM == 'inland' & ROUND == '18' & iterations == 1, ]) + 
  geom_point(aes(x = UNSUPPRESSED_SHARE, y = PROP_UNSUPPRESSED_POSTERIOR_SAMPLE))


# summarise
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.age = df[, list(q= quantile(UNSUPPRESSED_SHARE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGEYRS')]
sing.age = as.data.table(reshape2::dcast(sing.age, ... ~ q_label, value.var = "q"))

# name
setnames(sing.age, qlab, paste0('UNSUPPRESSED_SHARE_AGE_AND_SEX_', qlab))

# plot
ggplot(sing.age[COMM == 'inland'], aes(x = AGEYRS)) + 
  geom_line(aes(y = UNSUPPRESSED_SHARE_AGE_AND_SEX_M)) + 
  geom_ribbon(aes(ymin = UNSUPPRESSED_SHARE_AGE_AND_SEX_CL, ymax = UNSUPPRESSED_SHARE_AGE_AND_SEX_CU), alpha = 0.5) + 
  facet_grid(ROUND~COMM+SEX) + 
  theme_bw()


#########################################

# FIND SEX SHARE OF INFECTED 

#########################################

# find share of unsuppressed by sex across age
tmp <- df[, list(UNSUPPRESSED = sum(UNSUPPRESSED)), by = c('ROUND', 'COMM', 'iterations', 'SEX')]
tmp[, UNSUPPRESSED_SHARE := UNSUPPRESSED / sum(UNSUPPRESSED), by = c('ROUND', 'COMM', 'iterations')]

# summarise
sing = tmp[, list(q= quantile(UNSUPPRESSED_SHARE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX')]
sing = as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

# name
setnames(sing, qlab, paste0('UNSUPPRESSED_SHARE_SEX_', qlab))

# plot
ggplot(sing, aes(x = ROUND)) + 
  geom_point(aes(y = UNSUPPRESSED_SHARE_SEX_M)) + 
  geom_errorbar(aes(ymin = UNSUPPRESSED_SHARE_SEX_CL, ymax = UNSUPPRESSED_SHARE_SEX_CU), alpha = 0.5) + 
  facet_grid(SEX~COMM) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,1))


#########################################

# SAVE

#########################################

tmp <- merge(sing.age, sing, by=c('ROUND', 'COMM', 'SEX'))
file.name <- file.path(indir.repository, 'fit', paste0('RCCS_unsuppressed_share_sex_221116.csv'))
write.csv(tmp, file = file.name, row.names = F)
