library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(stringr)
library(ggsci)
library(gridExtra)
library(ggpubr)

# directory repository
indir.repository <- '~/git/phyloflows'

# directory to save the figures
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')

# treatment cascade participants
file.cascade.participants <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_participants_estimates_221116.csv'))

# treatment cascade non-participants
file.cascade.nonparticipants <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_nonparticipants_estimates_221116.csv'))


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
ggsave(file = file.path(outdir, 'treatment_cascade_inland_participants_vs_nonparticipants_221108.pdf'), w = 6.5, h = 6.5)  


###########################

# PLOT PROPORTION OF UNSUPPRESSED

###########################

# select variable of interst for cascade
uns <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,PROP_UNSUPPRESSED_EMPIRICAL,
                               PROP_UNSUPPRESSED_M, PROP_UNSUPPRESSED_CL, PROP_UNSUPPRESSED_CU)]
uns[,SEX_LABEL := 'Men']
uns[SEX == 'F',SEX_LABEL := 'Women']
uns[,COMM_LABEL := 'Inland communities']
uns[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
uns[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

# for participants
tab <- uns[type == 'Participants' & COMM == 'inland' & ROUND != 'R015S']
ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL, col = SEX_LABEL, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M, col = SEX_LABEL, linetype= 'Fit')) + 
  facet_wrap(~ROUND_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of HIV-positive participants with unsuppressed viral load', 
       col = '', fill = '', linetype = '', shape = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))
ggsave(file = file.path(outdir, 'participants_smooth_unsuppressed_221108.pdf'), w = 7, h = 7)  

# for non-participants
tab <- uns[type == 'Non participants' & COMM == 'inland' & ROUND != 'R015S']
ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL, col = SEX_LABEL, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M, col = SEX_LABEL, linetype= 'Fit')) + 
  facet_wrap(~ROUND_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of newly registered HIV-positive participants\nwith unsuppressed viral load', 
       col = '', fill = '', linetype = '', shape = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))
ggsave(file = file.path(outdir, 'nonparticipants_smooth_unsuppressed_221108.pdf'), w = 7, h = 7)  


###########################

# PLOT PROPORTION OF ON ART

###########################

# select variable of interst for cascade
art <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,PROP_ART_COVERAGE_EMPIRICAL,
                           PROP_ART_COVERAGE_M, PROP_ART_COVERAGE_CL, PROP_ART_COVERAGE_CU)]
art[,SEX_LABEL := 'Men']
art[SEX == 'F',SEX_LABEL := 'Women']
art[,COMM_LABEL := 'Inland communities']
art[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
art[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]
art[ROUND == 'R010',PROP_ART_COVERAGE_EMPIRICAL := NA ]# round 10 was fixed to be the same as round 11

# for participants
tab <- art[type == 'Participants' & COMM == 'inland' & ROUND != 'R015S']
ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_ART_COVERAGE_CL, ymax = PROP_ART_COVERAGE_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_point(aes(y = PROP_ART_COVERAGE_EMPIRICAL, col = SEX_LABEL, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_ART_COVERAGE_M, col = SEX_LABEL, linetype = 'Fit')) + 
  facet_wrap(~ROUND_LABEL, ncol = 3) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of HIV-positive RCCS participants on ART', 
       col = '', fill = '', shape = '', linetype= '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))
ggsave(file = file.path(outdir, 'participants_smooth_art_221108.pdf'), w = 7, h = 7)  

# for non-participants
tab <- art[type == 'Non participants' & COMM == 'inland' & ROUND != 'R015S']
ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_ART_COVERAGE_CL, ymax = PROP_ART_COVERAGE_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_point(aes(y = PROP_ART_COVERAGE_EMPIRICAL, col = SEX_LABEL, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_ART_COVERAGE_M, col = SEX_LABEL, linetype = 'Fit')) + 
  facet_wrap(~ROUND_LABEL, ncol = 3) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of newly registered HIV-positive RCCS participants on ART', 
       col = '', fill = '', shape = '', linetype= '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))
ggsave(file = file.path(outdir, 'nonparticipants_smooth_art_221108.pdf'), w = 7, h = 7)  


###########################

# PLOT SUPPRESSION RATE

###########################

# select variable of interst for cascade
art <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,
                           SUPPRESSION_RATE_M, SUPPRESSION_RATE_CL, SUPPRESSION_RATE_CU)]
art[,SEX_LABEL := 'Men']
art[SEX == 'F',SEX_LABEL := 'Women']
art[,COMM_LABEL := 'Inland communities']
art[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
art[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

#
# for participants
#

# suppression rate
tab <- art[type == 'Participants' & COMM == 'inland' & !ROUND %in% c('R010', 'R011', 'R012', 'R013',  'R014', 'R015S', 'R015')]
p1 <- ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = SUPPRESSION_RATE_CL, ymax = SUPPRESSION_RATE_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_line(aes(y = SUPPRESSION_RATE_M, col = SEX_LABEL)) + 
  facet_wrap(~ROUND_LABEL, ncol = 1) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Viral suppression rate in HIV-infected participants', 
       col = '', fill = '', shape = '', linetype= '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))
p1 <- ggarrange(p1, labels = c('a'))

# suppression rate ratio
tab <- dcast.data.table(tab, ROUND_LABEL + AGEYRS + COMM ~ SEX, value.var = 'SUPPRESSION_RATE_M')
tab[, SUPPRESSION_RATE_RATIO := `F`/`M`]
cols <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(4)
p2 <- ggplot(tab, aes(x = AGEYRS)) + 
  geom_hline(aes(yintercept = 1), linetype = 'dashed', col = 'grey50') + 
  geom_line(aes(y = SUPPRESSION_RATE_RATIO, col = ROUND_LABEL)) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Female-to-male viral suppression rate\nratio in HIV-infected participants', 
       col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = cols) + 
  scale_y_log10() + 
  scale_x_continuous(expand = c(0,0))
p2 <- ggarrange(p2, labels = c('b'))

p <- grid.arrange(p1,p2,layout_matrix = rbind(c(1, 2), c(1, NA)), heights =c(0.6, 0.4), widths = c(0.45, 0.55))
ggsave(p, file = file.path(outdir, 'participants_suppression_rate_221108.pdf'), w = 8, h = 7)


#
# for non-participants
#

# suppression rate
tab <- art[type == 'Non participants' & COMM == 'inland' & !ROUND %in% c('R010', 'R011', 'R012', 'R013',  'R014', 'R015S', 'R015')]
p1 <- ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = SUPPRESSION_RATE_CL, ymax = SUPPRESSION_RATE_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_line(aes(y = SUPPRESSION_RATE_M, col = SEX_LABEL)) + 
  facet_wrap(~ROUND_LABEL, ncol = 1) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Viral suppression rate in HIV-infected participants', 
       col = '', fill = '', shape = '', linetype= '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))
p1 <- ggarrange(p1, labels = c('a'))

# suppression rate ratio
tab <- dcast.data.table(tab, ROUND_LABEL + AGEYRS + COMM ~ SEX, value.var = 'SUPPRESSION_RATE_M')
tab[, SUPPRESSION_RATE_RATIO := `F`/`M`]
cols <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(4)
p2 <- ggplot(tab, aes(x = AGEYRS)) + 
  geom_hline(aes(yintercept = 1), linetype = 'dashed', col = 'grey50') + 
  geom_line(aes(y = SUPPRESSION_RATE_RATIO, col = ROUND_LABEL)) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Female-to-male viral suppression rate\nratio in HIV-infected participants', 
       col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = cols) + 
  scale_y_log10() + 
  scale_x_continuous(expand = c(0,0))
p2 <- ggarrange(p2, labels = c('b'))

p <- grid.arrange(p1,p2,layout_matrix = rbind(c(1, 2), c(1, NA)), heights =c(0.6, 0.4), widths = c(0.45, 0.55))
ggsave(p, file = file.path(outdir, 'nonparticipants_suppression_rate_221108.pdf'), w = 8, h = 7)

