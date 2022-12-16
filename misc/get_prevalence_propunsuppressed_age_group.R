library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

# directory repository
indir.repository <- '~/git/phyloflows'

# outdir
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

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

# FIND  PREVALENCE & PROPORTION UNSUPPRESSED BY AGE GROUP '15-24', '25-34', '35-49'

#####################################################

# age groups
df_age_aggregated <- data.table(AGEYRS = df[, sort(unique(AGEYRS))])
df_age_aggregated[, INDEX_AGE_GROUP := 3]
df_age_aggregated[AGEYRS < 35, INDEX_AGE_GROUP := 2]
df_age_aggregated[AGEYRS < 25, INDEX_AGE_GROUP := 1]
df_age_aggregated[, AGE_GROUP := c('15-24', '25-34', '35-49')[INDEX_AGE_GROUP]]

# aggregated by age groups
df_1 <- merge(df, df_age_aggregated, by = 'AGEYRS')
df_1 <- df_1[, list(INFECTED = sum(INFECTED), 
                UNSUPPRESSED = sum(UNSUPPRESSED), 
                ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'COMM', 'AGE_GROUP', 'SEX', 'iterations')]

# find prevalence by age groups across age
df_1[, PROP_UNSUPPRESSED := UNSUPPRESSED / ELIGIBLE]
df_1[, PREVALENCE := INFECTED / ELIGIBLE]

# summarise unsuppressed
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.age = df_1[, list(q= quantile(na.omit(PROP_UNSUPPRESSED), prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')]
sing.age = as.data.table(reshape2::dcast(sing.age, ... ~ q_label, value.var = "q"))
setnames(sing.age, qlab, paste0('PROP_UNSUPPRESSED_AGE_GROUP_', qlab))

# summarise prevalence
sinf.age = df_1[, list(q= quantile(PREVALENCE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')]
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



#####################################################

# FIND  PREVALENCE &  PROPORTION UNSUPPRESSED BY AGE GROUP '15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49'

#####################################################

# age groups
df_age_aggregated <- data.table(AGEYRS = df[, sort(unique(AGEYRS))])
df_age_aggregated[, INDEX_AGE_GROUP := 7]
df_age_aggregated[AGEYRS < 45, INDEX_AGE_GROUP := 6]
df_age_aggregated[AGEYRS < 40, INDEX_AGE_GROUP := 5]
df_age_aggregated[AGEYRS < 35, INDEX_AGE_GROUP := 4]
df_age_aggregated[AGEYRS < 30, INDEX_AGE_GROUP := 3]
df_age_aggregated[AGEYRS < 25, INDEX_AGE_GROUP := 2]
df_age_aggregated[AGEYRS < 20, INDEX_AGE_GROUP := 1]
df_age_aggregated[, AGE_GROUP := c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49')[INDEX_AGE_GROUP]]

# aggregated by age groups
df_1 <- merge(df, df_age_aggregated, by = 'AGEYRS')
df_1 <- df_1[, list(INFECTED = sum(INFECTED), 
                UNSUPPRESSED = sum(UNSUPPRESSED), 
                ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'COMM', 'AGE_GROUP', 'SEX', 'iterations')]

# find prevalence by age groups across age
df_1[, PROP_UNSUPPRESSED := UNSUPPRESSED / INFECTED]
df_1[, PREVALENCE := INFECTED / ELIGIBLE]
df_1[, DIFF_PROP_UNSUPPRESSED := PROP_UNSUPPRESSED[SEX == 'M'] -  PROP_UNSUPPRESSED[SEX == 'F'], by = c('ROUND', 'COMM', 'AGE_GROUP', 'iterations')]

# summarise unsuppressed
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.age = df_1[, list(q= quantile(na.omit(PROP_UNSUPPRESSED), prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')]
sing.age = as.data.table(reshape2::dcast(sing.age, ... ~ q_label, value.var = "q"))
setnames(sing.age, qlab, paste0('PROP_UNSUPPRESSED_AGE_GROUP_', qlab))

# summarise prevalence
sinf.age = df_1[, list(q= quantile(PREVALENCE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')]
sinf.age = as.data.table(reshape2::dcast(sinf.age, ... ~ q_label, value.var = "q"))
setnames(sinf.age, qlab, paste0('PREVALENCE_AGE_GROUP_', qlab))

# summarise difference prop unsuppressed
sind.age <- unique(df_1[, .(ROUND, COMM, AGE_GROUP, iterations, DIFF_PROP_UNSUPPRESSED)])
sind.age = sind.age[, list(q= quantile(DIFF_PROP_UNSUPPRESSED, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'AGE_GROUP')]
sind.age = as.data.table(reshape2::dcast(sind.age, ... ~ q_label, value.var = "q"))
setnames(sind.age, qlab, paste0('DIFF_PROP_UNSUPPRESSED_AGE_GROUP_', qlab))


#########################################

# SAVE FOR ROUND 18

#########################################

n_digits <- 1

tmp <- sing.age[ROUND == '18' & COMM == 'inland']
tmp[, `:=` (PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CL*100, n_digits), nsmall = n_digits), 
            PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CU*100, n_digits), nsmall = n_digits), 
            PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(PROP_UNSUPPRESSED_AGE_GROUP_M*100, n_digits), nsmall = n_digits))]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_CL := gsub(' ', '', PROP_UNSUPPRESSED_AGE_GROUP_CL)]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_CU := gsub(' ', '', PROP_UNSUPPRESSED_AGE_GROUP_CU)]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_M := gsub(' ', '', PROP_UNSUPPRESSED_AGE_GROUP_M)]
file.name <- file.path(outdir, paste0('RCCS_propunsuppressed_age_group_5years_R18_221215.rds'))
saveRDS(tmp, file = file.name)

tmp <- sinf.age[ROUND == '18' & COMM == 'inland']
tmp[, `:=` (PREVALENCE_AGE_GROUP_CL = format(round(PREVALENCE_AGE_GROUP_CL*100, n_digits), nsmall = n_digits), 
            PREVALENCE_AGE_GROUP_CU = format(round(PREVALENCE_AGE_GROUP_CU*100, n_digits), nsmall = n_digits), 
            PREVALENCE_AGE_GROUP_M = format(round(PREVALENCE_AGE_GROUP_M*100, n_digits), nsmall = n_digits))]
tmp[, PREVALENCE_AGE_GROUP_CL := gsub(' ', '', PREVALENCE_AGE_GROUP_CL)]
tmp[, PREVALENCE_AGE_GROUP_CU := gsub(' ', '', PREVALENCE_AGE_GROUP_CU)]
tmp[, PREVALENCE_AGE_GROUP_M := gsub(' ', '', PREVALENCE_AGE_GROUP_M)]
file.name <- file.path(outdir, paste0('RCCS_prevalence_age_group_5years_R18_221215.rds'))
saveRDS(tmp, file = file.name)

tmp <- sind.age[ROUND == '18' & COMM == 'inland']
tmp[, `:=` (DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL*100, n_digits), nsmall = n_digits), 
            DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU*100, n_digits), nsmall = n_digits), 
            DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M*100, n_digits), nsmall = n_digits))]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL := gsub(' ', '', DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL)]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU := gsub(' ', '', DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU)]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M := gsub(' ', '', DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M)]
file.name <- file.path(outdir, paste0('RCCS_diffpropunsuppressed_age_group_5years_R18_221215.rds'))
saveRDS(tmp, file = file.name)



#####################################################

# FIND PREVALENCE & PROPORTION UNSUPPRESSED ACROSS AGE

#####################################################

# aggregated across ages
df_1 <- df[, list(INFECTED = sum(INFECTED), 
                    UNSUPPRESSED = sum(UNSUPPRESSED), 
                    ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'COMM', 'SEX', 'iterations')]

# find prevalence and unsuppressed
df_1[, PROP_UNSUPPRESSED := UNSUPPRESSED / INFECTED]
df_1[, PREVALENCE := INFECTED / ELIGIBLE]
df_1[, DIFF_PROP_UNSUPPRESSED := PROP_UNSUPPRESSED[SEX == 'M'] -  PROP_UNSUPPRESSED[SEX == 'F'], by = c('ROUND', 'COMM', 'iterations')]

# summarise unsuppressed
ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')
sing.t = df_1[, list(q= quantile(na.omit(PROP_UNSUPPRESSED), prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX')]
sing.t = as.data.table(reshape2::dcast(sing.t, ... ~ q_label, value.var = "q"))
setnames(sing.t, qlab, paste0('PROP_UNSUPPRESSED_AGE_GROUP_', qlab))

# summarise prevalence
sinf.t = df_1[, list(q= quantile(PREVALENCE, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM', 'SEX')]
sinf.t = as.data.table(reshape2::dcast(sinf.t, ... ~ q_label, value.var = "q"))
setnames(sinf.t, qlab, paste0('PREVALENCE_AGE_GROUP_', qlab))

# summarise difference prop unsuppressed
sind.t <- unique(df_1[, .(ROUND, COMM, iterations, DIFF_PROP_UNSUPPRESSED)])
sind.t = sind.t[, list(q= quantile(DIFF_PROP_UNSUPPRESSED, prob=ps, na.rm = T), q_label=qlab), by=c('ROUND', 'COMM')]
sind.t = as.data.table(reshape2::dcast(sind.t, ... ~ q_label, value.var = "q"))
setnames(sind.t, qlab, paste0('DIFF_PROP_UNSUPPRESSED_AGE_GROUP_', qlab))


#########################################

# SAVE FOR ROUND 18

#########################################

n_digits <- 1

tmp <- sing.t[ROUND == '18' & COMM == 'inland']
tmp[, `:=` (PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CL*100, n_digits), nsmall = n_digits), 
            PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CU*100, n_digits), nsmall = n_digits), 
            PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(PROP_UNSUPPRESSED_AGE_GROUP_M*100, n_digits), nsmall = n_digits))]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_CL := gsub(' ', '', PROP_UNSUPPRESSED_AGE_GROUP_CL)]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_CU := gsub(' ', '', PROP_UNSUPPRESSED_AGE_GROUP_CU)]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_M := gsub(' ', '', PROP_UNSUPPRESSED_AGE_GROUP_M)]
file.name <- file.path(outdir, paste0('RCCS_propunsuppressed_total_R18_221215.rds'))
saveRDS(tmp, file = file.name)

tmp <- sinf.t[ROUND == '18' & COMM == 'inland']
tmp[, `:=` (PREVALENCE_AGE_GROUP_CL = format(round(PREVALENCE_AGE_GROUP_CL*100, n_digits), nsmall = n_digits), 
            PREVALENCE_AGE_GROUP_CU = format(round(PREVALENCE_AGE_GROUP_CU*100, n_digits), nsmall = n_digits), 
            PREVALENCE_AGE_GROUP_M = format(round(PREVALENCE_AGE_GROUP_M*100, n_digits), nsmall = n_digits))]
tmp[, PREVALENCE_AGE_GROUP_CL := gsub(' ', '', PREVALENCE_AGE_GROUP_CL)]
tmp[, PREVALENCE_AGE_GROUP_CU := gsub(' ', '', PREVALENCE_AGE_GROUP_CU)]
tmp[, PREVALENCE_AGE_GROUP_M := gsub(' ', '', PREVALENCE_AGE_GROUP_M)]
file.name <- file.path(outdir, paste0('RCCS_prevalence_total_R18_221215.rds'))
saveRDS(tmp, file = file.name)

tmp <- sind.t[ROUND == '18' & COMM == 'inland']
tmp[, `:=` (DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL*100, n_digits), nsmall = n_digits), 
            DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU*100, n_digits), nsmall = n_digits), 
            DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M*100, n_digits), nsmall = n_digits))]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL := gsub(' ', '', DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL)]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU := gsub(' ', '', DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU)]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M := gsub(' ', '', DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M)]
file.name <- file.path(outdir, paste0('RCCS_diffpropunsuppressed_total_R18_221215.rds'))
saveRDS(tmp, file = file.name)

