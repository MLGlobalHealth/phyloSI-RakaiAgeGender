library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

# directory to repository
indir.repository <- "~/git/phyloflows"

# outdir to save figures
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')

# posterior samples participants
file.unsuppressedviralload <- file.path(indir.repository, 'fit', paste0('RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.rds'))
file.selfreportedart <- file.path(indir.repository, 'fit', paste0('RCCS_art_posterior_samples_221101.rds'))

# posterior samples non-participants
file.unsuppressedviralload.np <- file.path(indir.repository, 'fit', paste0('RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_newlyregistered_221101.rds'))
file.selfreportedart.np <- file.path(indir.repository, 'fit', paste0('RCCS_art_posterior_samples_newlyregistered_221101.rds'))

# participation rate
file.participation <- file.path(indir.repository, 'data', 'RCCS_participation_220915.csv')

ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')


###########################

# COMBINE POSTERIOR SAMPLE

###########################

# load proportion unsuppressed viral load 
uns.p <- as.data.table(readRDS(file.unsuppressedviralload))
uns.np <- as.data.table(readRDS(file.unsuppressedviralload.np))

# load proportion art
sre.p <- as.data.table(readRDS(file.selfreportedart))
sre.np <- as.data.table(readRDS(file.selfreportedart.np))

# add label
setnames(uns.p, 'PROP_UNSUPPRESSED_POSTERIOR_SAMPLE', 'PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE')
setnames(uns.np, 'PROP_UNSUPPRESSED_POSTERIOR_SAMPLE', 'PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE')
setnames(sre.p, 'PROP_ART_COVERAGE_POSTERIOR_SAMPLE', 'PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE')
setnames(sre.np, 'PROP_ART_COVERAGE_POSTERIOR_SAMPLE', 'PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE')

# remove crude estimates
set(uns.p, NULL, 'PROP_UNSUPPRESSED_EMPIRICAL', NULL)
set(uns.np, NULL, 'PROP_UNSUPPRESSED_EMPIRICAL', NULL)
set(sre.p, NULL, 'PROP_ART_COVERAGE_EMPIRICAL', NULL)
set(sre.np, NULL, 'PROP_ART_COVERAGE_EMPIRICAL', NULL)

# merge participants (fyi for round < 15, we do not have viral load info)
df.p <- merge(uns.p, sre.p, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'), all.y = T)

# merge non-participants (fyi for round < 15, we do not have viral load info)
df.np <- merge(uns.np, sre.np, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'), all.y = T)

# merge population
df <- merge(df.p, df.np, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'))

# for round 10, set art coverage and unsuppressed as round 11 as questions art related q were introduced from roudn 11
tmp <- df[ROUND == 'R011']
tmp[, ROUND := 'R010']
df <- rbind(tmp, df[ROUND != 'R010'])

# restrict age 
df <- df[AGEYRS > 14 & AGEYRS < 50]

# add participation
participation <- as.data.table(read.csv(file.participation))
participation[, ROUND := paste0('R0', ROUND)]
df <- merge(df, participation, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND'))

# find proportion suppressed
df[, PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]
df[, PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# find suppression rate
df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE / PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE]
df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE / PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# add constraint that suppression rate must at maximum 1
df[SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE > 1, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := 1]
df[SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1]


####################################

# FIND PROPORTION OF DIAGNOSED

####################################

# all participants are diagnosed
df[, PROP_DIAGNOSED_PARTICIPANT_POSTERIOR_SAMPLE := 1]

# assume that proportion diagnosed in non -participants = proportion on art
df[, PROP_DIAGNOSED_NONPARTICIPANT_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# find proportion diagnosed in the population
df[, PROP_DIAGNOSED_POSTERIOR_SAMPLE := PROP_DIAGNOSED_PARTICIPANT_POSTERIOR_SAMPLE * PARTICIPATION + PROP_DIAGNOSED_NONPARTICIPANT_POSTERIOR_SAMPLE * (1-PARTICIPATION)]


####################################

# FIND PROPORTION OF ART USE GIVEN DIAGNOSED

####################################

# for round 15 find % on art by using the same suppression rate as round 16
df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE := min(1, PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE / SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE), by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations')]
set(df, NULL, 'SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)

df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE := min(1, PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE / SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE), by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations')]
set(df, NULL, 'SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)

# add constraint that art coverage must be smaller than prop suppression (can occur if prop suppression > art coverage)
df[PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE < PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE, PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]
df[PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE < PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE, PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]
stopifnot(nrow(df[PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)

# find art coverage in population
df[, PROP_ART_COVERAGE_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE * PARTICIPATION + PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE * (1-PARTICIPATION)]

# find art coverage given diagnosed in population
df[, PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE := min(1, PROP_ART_COVERAGE_POSTERIOR_SAMPLE / PROP_DIAGNOSED_POSTERIOR_SAMPLE), by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations')]

stopifnot(nrow(df[PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE>1]) == 0)


####################################

# FIND SUPPRESSION RATE FOR ALL ROUND

####################################

# for round <15 find % suppressed by using the same suppression rate as round 16
df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE * SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]
set(df, NULL, 'SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)

df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE * SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]
set(df, NULL, 'SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)

stopifnot(nrow(df[SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)
stopifnot(nrow(df[SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)


####################################

# FIND PROPORTION OF SUPPRESSION RATE 

####################################

# find suppression in the population
df[, PROP_SUPPRESSED_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE * PARTICIPATION + PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE * (1-PARTICIPATION)]

# find suppression given art use in the population (i.e., suppression rate)
df[, PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE / PROP_ART_COVERAGE_POSTERIOR_SAMPLE ]

stopifnot(nrow(df[PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE > 1]) == 0)


####################################

# SUMMARISE

####################################

tmp <- df[, .(ROUND, SEX, COMM, AGEYRS, iterations, PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE, 
             PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE, PROP_DIAGNOSED_POSTERIOR_SAMPLE)]
tmp <- melt.data.table(tmp, id.vars= c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'))
tmp[, variable := gsub('(.+)_POSTERIOR_SAMPLE', '\\1', variable)]
ns = tmp[, list(q= quantile(value, prob=ps, na.rm = T), q_label=qlab), by=c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'variable')]
ns = as.data.table(reshape2::dcast(ns, AGEYRS + SEX + COMM + ROUND + variable ~ q_label, value.var = "q"))

# check all entries are complete
stopifnot(nrow(ns[COMM == 'inland']) == ns[, length(unique(AGEYRS))] * ns[, length(unique(SEX))] * ns[, length(unique(variable))] * ns[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(ns[COMM == 'fishing']) == ns[, length(unique(AGEYRS))] * ns[, length(unique(SEX))] * ns[, length(unique(variable))] * ns[COMM == 'fishing', length(unique(ROUND))])
stopifnot(nrow(df[COMM == 'inland']) == df[, length(unique(AGEYRS))] * df[, length(unique(SEX))] * df[, length(unique(iterations))] * df[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(df[COMM == 'fishing']) == df[, length(unique(AGEYRS))] * df[, length(unique(SEX))] * df[, length(unique(iterations))] * df[COMM == 'fishing', length(unique(ROUND))])


####################################

# PLOT

####################################

# label
diagnosed.label = 'estimated proportion of diagnosed'
on.art.label = "proportion of infected who report ART use"
suppressed.label = "proportion of infected reporting ART use who had a viral load < 1,000 cps / mL blood"

# add label
tab <- copy(ns)
tab[, VARIABLE_LEVEL := diagnosed.label]
tab[variable == 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED', VARIABLE_LEVEL := on.art.label]
tab[variable == 'PROP_SUPPRESSION_GIVEN_ART', VARIABLE_LEVEL := suppressed.label]
tab[, VARIABLE_LEVEL := factor(VARIABLE_LEVEL, 
                               levels = c(diagnosed.label, on.art.label, suppressed.label))]
tab[,SEX_LABEL := 'Men']
tab[SEX == 'F',SEX_LABEL := 'Women']
tab[,COMM_LABEL := 'Inland communities']
tab[COMM == 'fishing',COMM_LABEL := 'Fishing communities']

tmp <- tab[ROUND == 'R018' & COMM == 'inland']
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_hline(aes(yintercept = 0.9), linetype='dashed', alpha = 0.5) +
  geom_hline(aes(yintercept = 0.95), linetype='dashed', alpha = 0.5) +
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = VARIABLE_LEVEL), alpha = 0.25) + 
  geom_line(aes(y = M, col = VARIABLE_LEVEL)) + 
  facet_grid(~SEX_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = '', col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  ggsci::scale_color_jama()+
  ggsci::scale_fill_jama()+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  guides(color = guide_legend(byrow = T, nrow = 3),
         fill = guide_legend(byrow = T, nrow = 3))
ggsave(file = file.path(outdir, 'treatment_cascade_inland_population_221108.pdf'), w = 6.5, h = 4.5)  

tmp <- tab[COMM == 'inland']
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_hline(aes(yintercept = 0.9), linetype='dashed', alpha = 0.5) +
  geom_hline(aes(yintercept = 0.95), linetype='dashed', alpha = 0.5) +
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = VARIABLE_LEVEL), alpha = 0.25) + 
  geom_line(aes(y = M, col = VARIABLE_LEVEL)) + 
  facet_grid(ROUND~SEX_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = '', col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  ggsci::scale_color_jama()+
  ggsci::scale_fill_jama()+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  guides(color = guide_legend(byrow = T, nrow = 3),
         fill = guide_legend(byrow = T, nrow = 3))


####################################

# SAVE

####################################

file.name <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_population_posterior_samples_221101.rds'))
saveRDS(df, file = file.name)

file.name <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_population_estimates_221101.csv'))
write.csv(ns, file = file.name, row.names = F)

