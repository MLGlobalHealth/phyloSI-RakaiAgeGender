library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

file.art.coverage <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_artcoverage_posterior_samples_220906.csv'))
file.prevalence <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('RCCS_prevalence_posterior_sample_220818.csv'))
file.eligible.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220830.csv')

# load census eligible ount
eligible_count <- as.data.table(read.csv(file.eligible.count))

# load proportion prevalence
proportion_prevalence <- as.data.table(read.csv(file.prevalence))

# load unsuppressed proportion 
proportion_unsuppressed <- as.data.table(read.csv(file.art.coverage))


#############################

# FIND INFECTED UNSUPPRESSED

############################

# define round
df <- merge(proportion_prevalence, proportion_unsuppressed, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX', 'iterations'))
df[, ROUND := gsub('R0(.+)', '\\1', ROUND)]

# merge number of eligible and the prevalence
df <- merge(eligible_count[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)], df, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))

# find infected
df[, INFECTED := ELIGIBLE * PREVALENCE_POSTERIOR_SAMPLE]

# find infected unsuppressed
df[, UNSUPPRESSED := INFECTED * PROP_UNSUPPRESSED_POSTERIOR_SAMPLE]


#####################################################

# FIND SEX SHARE OF INFECTED UNSUPPRESSED ACROSS AGE

#####################################################

# find share of infected by sex across age
df[, TOTAL_UNSUPPRESSED := sum(UNSUPPRESSED), by = c('ROUND', 'COMM', 'iterations')]
df[, UNSUPPRESSED_SHARE := UNSUPPRESSED / TOTAL_UNSUPPRESSED, by = c('ROUND', 'COMM', 'AGEYRS', 'iterations')]

# summarise
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.age = df[, list(q= quantile(UNSUPPRESSED_SHARE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGEYRS')]
sing.age = as.data.table(reshape2::dcast(sing.age, ... ~ q_label, value.var = "q"))

# name
setnames(sing.age, qlab, paste0('UNSUPPRESSED_SHARE_AGE_AND_SEX_', qlab))

# plot
ggplot(sing.age, aes(x = AGEYRS)) + 
  geom_line(aes(y = UNSUPPRESSED_SHARE_AGE_AND_SEX_M)) + 
  geom_ribbon(aes(ymin = UNSUPPRESSED_SHARE_AGE_AND_SEX_CL, ymax = UNSUPPRESSED_SHARE_AGE_AND_SEX_CU), alpha = 0.5) + 
  facet_grid(ROUND~COMM+SEX) + 
  theme_bw()


#########################################

# FIND SEX SHARE OF INFECTED 

#########################################

# find share of infected by sex across age
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

file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_artcoverage_share_sex_220906.csv'))
write.csv(tmp, file = file.name, row.names = F)
