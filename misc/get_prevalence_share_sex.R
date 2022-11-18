library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

# directory to reposity
indir.repository <- '~/git/phyloflows'

# files
file.prevalence <- file.path(indir.repository, 'fit', paste0('RCCS_prevalence_posterior_sample_221116.rds'))
file.eligible.count <- file.path(indir.repository, 'data', 'RCCS_census_eligible_individuals_221116.csv')

# load census eligible ount
eligible_count <- as.data.table(read.csv(file.eligible.count))

# load proportion prevalence
proportion_prevalence <- as.data.table(readRDS(file.prevalence))


####################

# FIND INFECTED

###################

# define round
df <- copy(proportion_prevalence)
df[, ROUND := gsub('R0(.+)', '\\1', ROUND)]

# merge number of eligible and the prevalence
df <- merge(eligible_count[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)], df, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))

# find infected
df[, INFECTED := ELIGIBLE * PREVALENCE_POSTERIOR_SAMPLE]


#########################################

# FIND SEX SHARE OF INFECTED ACROSS AGE

#########################################

# find share of infected by sex across age
df[, TOTAL_INFECTED := sum(INFECTED), by = c('ROUND', 'COMM', 'iterations')]
df[, INFECTED_SHARE := INFECTED / TOTAL_INFECTED, by = c('ROUND', 'COMM', 'AGEYRS', 'iterations')]

# summarise
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.age = df[, list(q= quantile(INFECTED_SHARE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGEYRS')]
sing.age = as.data.table(reshape2::dcast(sing.age, ... ~ q_label, value.var = "q"))

# name
setnames(sing.age, qlab, paste0('PREVALENCE_SHARE_SEX_AND_AGE_', qlab))

# plot
ggplot(sing.age, aes(x = AGEYRS)) + 
  geom_line(aes(y = PREVALENCE_SHARE_SEX_AND_AGE_M)) + 
  geom_ribbon(aes(ymin = PREVALENCE_SHARE_SEX_AND_AGE_CL, ymax = PREVALENCE_SHARE_SEX_AND_AGE_CU), alpha = 0.5) + 
  facet_grid(ROUND~COMM+SEX) + 
  theme_bw()


#########################################

# FIND SEX SHARE OF INFECTED 

#########################################

# find share of infected by sex across age
tmp <- df[, list(INFECTED = sum(INFECTED)), by = c('ROUND', 'COMM', 'iterations', 'SEX')]
tmp[, INFECTED_SHARE_SEX := INFECTED / sum(INFECTED), by = c('ROUND', 'COMM', 'iterations')]

# summarise
sing = tmp[, list(q= quantile(INFECTED_SHARE_SEX, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX')]
sing = as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

# name
setnames(sing, qlab, paste0('PREVALENCE_SHARE_SEX_', qlab))

# plot
ggplot(sing, aes(x = ROUND)) + 
  geom_point(aes(y = PREVALENCE_SHARE_SEX_M)) + 
  geom_errorbar(aes(ymin = PREVALENCE_SHARE_SEX_CL, ymax = PREVALENCE_SHARE_SEX_CU), alpha = 0.5) + 
  facet_grid(SEX~COMM) + 
  theme_bw() + 
  scale_y_continuous(limits = c(0,1))


#########################################

# SAVE

#########################################

tmp <- merge(sing.age, sing, by=c('ROUND', 'COMM', 'SEX'))

file.name <- file.path(indir.repository, 'fit', paste0('RCCS_prevalence_share_sex_221116.csv'))
write.csv(tmp, file = file.name, row.names = F)

