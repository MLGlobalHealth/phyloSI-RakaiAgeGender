library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')

# posterior samples
file.unsuppressedviralload <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.csv'))
file.selfreportedart <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_art_posterior_samples_221101.csv'))

ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')


############################

# COMBINE POSTERIOR SAMPLE

############################

# load proportion unsuppressed viral load
uns <- as.data.table(read.csv(file.unsuppressedviralload))

# load proportion art
sre <- as.data.table(read.csv(file.selfreportedart))

# merge (fyi for round < 15, we do not have viral load info)
df <- merge(uns, sre, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'), all.y = T)

# for round 10, set art coverage and unsuppressed as round 11 as questions art related q were introduced from roudn 11
tmp <- df[ROUND == 'R011']
tmp[, ROUND := 'R010']
df <- rbind(tmp, df[ROUND != 'R010'])

# restrict age 
df <- df[AGEYRS > 14 & AGEYRS < 50]

# find proportion suppressed among HIV positive participants
df[, PROP_SUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_POSTERIOR_SAMPLE]

# find suppression rate among HIV positive participants
df[, SUPPRESSION_RATE_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE / PROP_ART_COVERAGE_POSTERIOR_SAMPLE]
df[, SUPPRESSION_RATE_EMPIRICAL := (1-PROP_UNSUPPRESSED_EMPIRICAL) / PROP_ART_COVERAGE_EMPIRICAL ]
df[SUPPRESSION_RATE_EMPIRICAL > 1 + 1e-6, table(ROUND)] # the suppression rate should be < 1


##################################################################

# FIND PROPORTION OF DIAGNOSED AMONG HIV POSITIVE PARTICIPANTS

##################################################################

# all participants are diagnosed
df[, PROP_DIAGNOSED_POSTERIOR_SAMPLE := 1]


######################################################################################

# FIND ART COVERERAGE GIVEN DIAGNOSED AMONG HIV POSITIVE PARTICIPANTS FOR ALL ROUND

######################################################################################

# for round 15 find % on art by using the same suppression rate as round 16
df[, SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[, SUPPRESSION_RATE_EMPIRICAL_R016 := SUPPRESSION_RATE_EMPIRICAL[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]

df[ROUND %in% c('R015'), SUPPRESSION_RATE_POSTERIOR_SAMPLE := SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), SUPPRESSION_RATE_EMPIRICAL := SUPPRESSION_RATE_EMPIRICAL_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[SUPPRESSION_RATE_EMPIRICAL > 1 + 1e-6, table(ROUND)] # the suppression rate should be < 1

# find art coverage with the same suppression rate as round 16
df[ROUND %in% c('R015'), PROP_ART_COVERAGE_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE / SUPPRESSION_RATE_POSTERIOR_SAMPLE, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
set(df, NULL, 'SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016', NULL)
set(df, NULL, 'SUPPRESSION_RATE_EMPIRICAL_R016', NULL)

# add constraint that suppression rate must be smaller than 1
df <- df[PROP_ART_COVERAGE_POSTERIOR_SAMPLE <= 1]


##########################################################################

# FIND SUPPRESSION RATE AMONG HIV POSITIVE PARTICIPANTS FOR ALL ROUND

##########################################################################

# for round <15 find % suppressed by using the same suppression rate as round 16
df[, SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), SUPPRESSION_RATE_POSTERIOR_SAMPLE := SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_SUPPRESSED_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_POSTERIOR_SAMPLE * SUPPRESSION_RATE_POSTERIOR_SAMPLE, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_UNSUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_POSTERIOR_SAMPLE]
set(df, NULL, 'SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016', NULL)

# add constraint that suppression rate must be smaller than 1
df <- df[SUPPRESSION_RATE_POSTERIOR_SAMPLE <= 1]

# find prop suppressed given diagnosed among hiv positive participants
df[, PROP_SUPPRESSED_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE / PROP_DIAGNOSED_POSTERIOR_SAMPLE]


####################################

# SUMMARISE

####################################

# melt
df <- melt.data.table(df, id.vars= c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations', 'PROP_UNSUPPRESSED_EMPIRICAL', 'PROP_ART_COVERAGE_EMPIRICAL', 'SUPPRESSION_RATE_EMPIRICAL'))
df[, variable := gsub('(.+)_POSTERIOR_SAMPLE', '\\1', variable)]

# round the empirical otherwise we get two entries 
df[, PROP_UNSUPPRESSED_EMPIRICAL := round(PROP_UNSUPPRESSED_EMPIRICAL, 7)] 
df[, PROP_ART_COVERAGE_EMPIRICAL := round(PROP_ART_COVERAGE_EMPIRICAL, 7)] 
df[, SUPPRESSION_RATE_EMPIRICAL := round(SUPPRESSION_RATE_EMPIRICAL, 7)] 

# summarise
ns = df[, list(q= quantile(value, prob=ps, na.rm = T), q_label=paste0(variable, '_', qlab)), by=c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'PROP_UNSUPPRESSED_EMPIRICAL', 'PROP_ART_COVERAGE_EMPIRICAL', 'SUPPRESSION_RATE_EMPIRICAL', 'variable')]
ns = as.data.table(reshape2::dcast(ns, AGEYRS + SEX + COMM + ROUND + PROP_UNSUPPRESSED_EMPIRICAL + PROP_ART_COVERAGE_EMPIRICAL + SUPPRESSION_RATE_EMPIRICAL ~ q_label, value.var = "q"))

# check all entries are complete
stopifnot(nrow(ns[COMM == 'inland']) == ns[, length(unique(AGEYRS))] * ns[, length(unique(SEX))] * ns[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(ns[COMM == 'fishing']) == ns[, length(unique(AGEYRS))] * ns[, length(unique(SEX))] * ns[COMM == 'fishing', length(unique(ROUND))])


####################################

# PLOT

####################################

tmp <- copy(ns)
tmp[, ROUND_LABEL := paste0('ROUND ', gsub('R0(.+)', '\\1', ROUND))]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

# SUPPRESSION RATE
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_line(aes(y = SUPPRESSION_RATE_M)) + 
  geom_ribbon(aes(ymin = SUPPRESSION_RATE_CL, ymax = SUPPRESSION_RATE_CU), alpha = 0.5) + 
  geom_point(aes(y = (1-PROP_UNSUPPRESSED_EMPIRICAL) / PROP_ART_COVERAGE_EMPIRICAL), alpha = 0.5, col = 'darkred') + 
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  labs(x = 'Age', y = 'Suppression rate among HIV-positive participants', fill = '') +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) +
  scale_y_continuous(limits= c(0,1))
ggsave(file=file.path(outdir, paste0('smooth_suppression_rate_participants_221101.png')), w=9, h=8)


# PROPORTION OF ON ART 
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_line(aes(y = PROP_ART_COVERAGE_M)) + 
  geom_ribbon(aes(ymin = PROP_ART_COVERAGE_CL, ymax = PROP_ART_COVERAGE_CU), alpha = 0.5) + 
  geom_point(aes(y = PROP_ART_COVERAGE_EMPIRICAL), alpha = 0.5, col = 'darkred') + 
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  labs(x = 'Age', y = 'Proportion of HIV-positive participants with ART use', fill = '') +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) + 
  scale_y_continuous(labels = scales::percent, limits= c(0,1))
ggsave(file=file.path(outdir, paste0('smooth_art_coverage_participants_221101.png')), w=9, h=8)


# PROPORTION OF UNSUPPRESSED
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU), alpha = 0.5) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL), alpha = 0.5, col = 'darkred') + 
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  labs(x = 'Age', y = 'Proportion of HIV-positive participants with unsuppressed viral load', fill = '') +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) + 
  scale_y_continuous(labels = scales::percent, limits= c(0,1))
ggsave(file=file.path(outdir, paste0('smooth_unsuppressed_proportion_participants_221101.png')), w=9, h=8)


# PROPORTION OF SUPPRESSED
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_line(aes(y = PROP_SUPPRESSED_M)) + 
  geom_ribbon(aes(ymin = PROP_SUPPRESSED_CL, ymax = PROP_SUPPRESSED_CU), alpha = 0.5) + 
  geom_point(aes(y = 1 - PROP_UNSUPPRESSED_EMPIRICAL), alpha = 0.5, col = 'darkred') + 
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  labs(x = 'Age', y = 'Proportion of HIV-positive participants with suppressed viral load', fill = '') +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) + 
  scale_y_continuous(labels = scales::percent, limits= c(0,1))
ggsave(file=file.path(outdir, paste0('smooth_suppressed_proportion_participants_221101.png')), w=9, h=8)


# save
file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_treatment_cascade_participants_posterior_samples_221101.csv'))
write.csv(df, file = file.name, row.names = F)

file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_treatment_cascade_participants_estimates_221101.csv'))
write.csv(ns, file = file.name, row.names = F)
