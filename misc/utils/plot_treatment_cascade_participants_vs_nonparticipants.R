library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(stringr)
library(ggsci)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')

# treatment cascade participants
file.cascade.participants <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_treatment_cascade_participants_estimates_221101.csv'))

# treatment cascade non-participants
file.cascade.nonparticipants <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_treatment_cascade_nonparticipants_estimates_221101.csv'))


##############

# PRE-PROCESS

##############

# load files
cascade.p <- as.data.table(read.csv(file.cascade.participants))
cascade.np <- as.data.table(read.csv(file.cascade.nonparticipants))

# combine
cascade.p[, type := 'Participants']
cascade.np[, type := 'Non participants']
cascade.all.var <- rbind(cascade.p, cascade.np)
cascade.all.var[, type := factor(type, c('Participants', 'Non participants'))]

# select variable of interst for cascade
cascade <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,
                       PROP_DIAGNOSED_M, PROP_DIAGNOSED_CL, PROP_DIAGNOSED_CU,
                       PROP_ART_COVERAGE_M, PROP_ART_COVERAGE_CL, PROP_ART_COVERAGE_CU,
                       SUPPRESSION_RATE_M, SUPPRESSION_RATE_CL, SUPPRESSION_RATE_CU)]
setnames(cascade, c('PROP_ART_COVERAGE_M', 'PROP_ART_COVERAGE_CL', 'PROP_ART_COVERAGE_CU'),
         c('PROP_ART_COVERAGE_GIVEN_DIAGNOSED_M', 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED_CL', 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED_CU'))
setnames(cascade, c('SUPPRESSION_RATE_M', 'SUPPRESSION_RATE_CL', 'SUPPRESSION_RATE_CU'),
         c('PROP_SUPPRESSION_GIVEN_ART_M', 'PROP_SUPPRESSION_GIVEN_ART_CL', 'PROP_SUPPRESSION_GIVEN_ART_CU'))

# melt
cascade <- melt.data.table(cascade, id.vars = c('ROUND', 'SEX', 'AGEYRS', 'COMM', 'type'))

# find quantiles index
cascade[, q := gsub('_', '', str_sub(variable, start= -2))]
cascade[, variable := gsub('(.+)_.*', '\\1',variable)]

# cast
cascade <- dcast.data.table(cascade, ROUND + SEX + AGEYRS + COMM + type + variable ~ q, value.var = 'value')


###########################

# PLOT TREATMENT CASCADE

###########################

# label
diagnosed.label = 'estimated proportion of diagnosed'
on.art.label = "proportion of infected who report ART use"
suppressed.label = "proportion of infected reporting ART use who had a viral load < 1,000 cps / mL blood"

# add label
tab <- copy(cascade)
tab[, VARIABLE_LEVEL := diagnosed.label]
tab[variable == 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED', VARIABLE_LEVEL := on.art.label]
tab[variable == 'PROP_SUPPRESSION_GIVEN_ART', VARIABLE_LEVEL := suppressed.label]
tab[, VARIABLE_LEVEL := factor(VARIABLE_LEVEL, 
                               levels = c(diagnosed.label, on.art.label, suppressed.label))]
tab[,SEX_LABEL := 'Men']
tab[SEX == 'F',SEX_LABEL := 'Women']
tab[,COMM_LABEL := 'Inland communities']
tab[COMM == 'fishing',COMM_LABEL := 'Fishing communities']

tab <- tab[VARIABLE_LEVEL != diagnosed.label]
tmp <- tab[ROUND == 'R018' & COMM == 'inland']

colors <- pal_jama("default")(3)[2:3]

ggplot(tmp, aes(x = AGEYRS)) + 
  geom_hline(aes(yintercept = 0.9), linetype='dashed', alpha = 0.5) +
  geom_hline(aes(yintercept = 0.95), linetype='dashed', alpha = 0.5) +
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = VARIABLE_LEVEL), alpha = 0.25) + 
  geom_line(aes(y = M, col = VARIABLE_LEVEL)) + 
  facet_grid(type~SEX_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = '', col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  guides(color = guide_legend(byrow = T, nrow = 3),
         fill = guide_legend(byrow = T, nrow = 3))
ggsave(file = file.path(outdir, 'treatment_cascade_inland_participants_vs_nonparticipants.pdf'), w = 6.5, h = 6.5)  


###########################

# PLOT PROPORTION OF UNSUPPRESSED

###########################

# select variable of interst for cascade
uns <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,
                               PROP_UNSUPPRESSED_M, PROP_UNSUPPRESSED_CL, PROP_UNSUPPRESSED_CU)]

# for participants
tab <- uns[type == 'Participants']
tab[,SEX_LABEL := 'Men']
tab[SEX == 'F',SEX_LABEL := 'Women']
tab[,COMM_LABEL := 'Inland communities']
tab[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
tab[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU, fill = SEX_LABEL), alpha = 0.25) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M, col = SEX_LABEL)) + 
  facet_grid(ROUND_LABEL~COMM_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of unsuppressed among HIV-positive participants', 
       col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) 
ggsave(file = file.path(outdir, 'treatment_cascade_inland_participants_vs_nonparticipants.pdf'), w = 6.5, h = 6.5)  

