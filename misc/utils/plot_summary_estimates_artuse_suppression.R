library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(gridExtra)
library("haven")
library(ggpubr)

# directory of the repository
indir.repository <- '~/git/phyloflows'

# outdir directory for stan fit
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')

# files
artuse.path.np <- file.path(indir.repository, 'data', 'aggregated_newlyregistered_count_art_coverage.csv')
artuse.path.p <- file.path(indir.repository, 'data', 'aggregated_participants_count_art_coverage.csv')
file.cascade.participants <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_participants_estimates_221101.csv'))
file.cascade.nonparticipants <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_nonparticipants_estimates_221101.csv'))
file.cascade.population <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_population_estimates_221101.csv'))

# find count of newly registered participants who reported art use
rart.nonpart <- as.data.table( read.csv(artuse.path.np) )

# find count of participants who reported art use
rart.part <- as.data.table( read.csv(artuse.path.p) )

# treatment cascade participants 
cascade.participants  <- as.data.table( read.csv(file.cascade.participants) )

# treatment cascade non-participants
cascade.nonparticipants  <- as.data.table( read.csv(file.cascade.nonparticipants) )

# treatment cascade population
cascade.population  <- as.data.table( read.csv(file.cascade.population) )


#################################

# PLOT PARTICIPATION COUNT BY ART USE #

#################################

# participants
tmp <- copy(rart.part)
tmp[, `Do not use` := TOTAL_COUNT - COUNT] 
setnames(tmp, 'COUNT', 'Use')
tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'COMM', 'SEX', 'AGEYRS', 'TOTAL_COUNT'))
tmp[, variable := factor(variable, levels = c('Use', 'Do not use'))]
tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp[, SEX_LABEL := 'Women']
tmp[SEX== 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
tmp <- tmp[AGEYRS > 14 & AGEYRS < 50]

tmp1 <- tmp[!ROUND %in% c("06", "07", "08", "09", '10') & COMM == 'inland' & SEX == 'M']
pP.M <- ggplot(tmp1, aes(x = AGEYRS, y = value)) +
  geom_bar(aes(fill = variable), stat = 'identity') + 
  labs(x = 'Age', y = 'Count infected male participants', fill = '') +
  facet_grid(ROUND_LABEL~.) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        legend.justification = "left",
        strip.background = element_rect(colour="white", fill="white"),
        strip.text= element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 40)) + 
  scale_x_continuous(expand = c(0,0))+ 
  scale_fill_manual(values = c('#90B77D', '#425F57'), 
                    labels = c('Reported ART use', 'Did not report ART use')) 

# non participants
tmp <- copy(rart.nonpart)
tmp[, `Do not use` := TOTAL_COUNT - COUNT] 
setnames(tmp, 'COUNT', 'Use')
tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'COMM', 'SEX', 'AGEYRS', 'TOTAL_COUNT'))
tmp[, variable := factor(variable, levels = c('Use', 'Do not use'))]
tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp[, SEX_LABEL := 'Women']
tmp[SEX== 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
tmp <- tmp[AGEYRS > 14 & AGEYRS < 50]

tmp1 <- tmp[!ROUND %in% c("06", "07", "08", "09", '10') & COMM == 'inland' & SEX == 'M']
pNP.M <- ggplot(tmp1, aes(x = AGEYRS, y = value)) +
  geom_bar(aes(fill = variable), stat = 'identity') + 
  labs(x = 'Age', y = 'Count infected first-time male participants', fill = '') +
  facet_grid(ROUND_LABEL~.) +
  theme_bw() +
  theme(legend.position = 'none', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text= element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 40)) + 
  scale_x_continuous(expand = c(0,0))+ 
  scale_fill_manual(values = c('#90B77D', '#425F57'), 
                    labels = c('Reported ART use', 'Did not report ART use')) 


#################################

# PLOT ART USE ESTIMATES #

#################################

# process population
cap <- cascade.population[variable == 'PROP_ART_COVERAGE']
cap <- select(cap, -c('CL', 'CU'))
setnames(cap, c('M'), c('PROP_ART_COVERAGE_M'))
set(cap, NULL, 'variable', NULL)
cap[,SEX_LABEL := 'Men']
cap[SEX == 'F',SEX_LABEL := 'Women']
cap[,COMM_LABEL := 'Inland communities']
cap[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
cap[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

#process among participants
capa <- cascade.participants[, .(AGEYRS, SEX, COMM, ROUND, PROP_ART_COVERAGE_M, 
                                                 PROP_ART_COVERAGE_CL, PROP_ART_COVERAGE_CU, PROP_ART_COVERAGE_EMPIRICAL)]
capa[,SEX_LABEL := 'Men']
capa[SEX == 'F',SEX_LABEL := 'Women']
capa[,COMM_LABEL := 'Inland communities']
capa[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
capa[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]
capa[ROUND == 'R010',PROP_ART_COVERAGE_EMPIRICAL := NA ]# round 10 was fixed to be the same as round 11

# combine
label.population <- 'Estimates in\ncensus eligible population'
label.participant <- 'Estimates in\nparticipants'

# plot
tab <- capa[COMM == 'inland' & ROUND != 'R015S']
tabpop <- cap[COMM == 'inland' & ROUND != 'R015S']
p.ART <-ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_ART_COVERAGE_CL, ymax = PROP_ART_COVERAGE_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_point(aes(y = PROP_ART_COVERAGE_EMPIRICAL, col = SEX_LABEL, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_ART_COVERAGE_M, col = SEX_LABEL, linetype = label.participant)) + 
  geom_line(data = tabpop, aes(y = PROP_ART_COVERAGE_M, col = SEX_LABEL, size = label.population), linetype = 'dashed') + 
  facet_grid(ROUND_LABEL~.) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of infected individuals who report ART use', 
       col = '', fill = '', shape = '', linetype= '', size = '') +
  theme(legend.position = 'none', 
        # legend.box = 'vertical',
        panel.spacing.y = unit(0.6, "lines"),
        strip.background = element_rect(colour="white", fill="white"), 
        strip.text= element_blank()) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_size_manual(values = 0.5) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         size = guide_legend(order = 3), 
         color = guide_legend(order = 4),fill = guide_legend(order = 4))

#################################

# PLOT SUPPRESSED ESTIMATES #

#################################

# process population
cap <- cascade.population[variable == 'PROP_SUPPRESSED']
cap <- select(cap, -c('CL', 'CU'))
setnames(cap, 'M', 'PROP_SUPPRESSED_M')
set(cap, NULL, 'variable', NULL)
cap[,SEX_LABEL := 'Men']
cap[SEX == 'F',SEX_LABEL := 'Women']
cap[,COMM_LABEL := 'Inland communities']
cap[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
cap[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

#process among participants
capa <- cascade.participants[, .(AGEYRS, SEX, COMM, ROUND, PROP_SUPPRESSED_M, 
                                 PROP_SUPPRESSED_CL, PROP_SUPPRESSED_CU, PROP_UNSUPPRESSED_EMPIRICAL)]
capa[, PROP_SUPPRESSED_EMPIRICAL := 1- PROP_UNSUPPRESSED_EMPIRICAL]
capa[,SEX_LABEL := 'Men']
capa[SEX == 'F',SEX_LABEL := 'Women']
capa[,COMM_LABEL := 'Inland communities']
capa[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
capa[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]
capa[ROUND == 'R010',PROP_ART_COVERAGE_EMPIRICAL := NA ]# round 10 was fixed to be the same as round 11

# combine
label.population <- 'Estimates in\ncensus eligible population'
label.participant <- 'Estimates in\nparticipants'

# plot
tab <- capa[COMM == 'inland' & ROUND != 'R015S']
tabpop <- cap[COMM == 'inland' & ROUND != 'R015S']
p.US <-ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_SUPPRESSED_CL, ymax = PROP_SUPPRESSED_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_point(aes(y = PROP_SUPPRESSED_EMPIRICAL, col = SEX_LABEL, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_SUPPRESSED_M, col = SEX_LABEL, linetype = label.participant)) + 
  geom_line(data = tabpop, aes(y = PROP_SUPPRESSED_M, col = SEX_LABEL, size = label.population), linetype = 'dashed') + 
  facet_grid(ROUND_LABEL~.) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of infected individuals with suppressed virus', 
       col = '', fill = '', shape = '', linetype= '', size = '') +
  theme(legend.position = 'bottom', 
        legend.direction = 'vertical',
        legend.justification = "right",
        legend.spacing.x = unit(1, "pt"),
        legend.margin=margin(4,4,4,4),
        legend.box.spacing = margin(0.5),
        strip.text = element_text(size = rel(1)),
        panel.spacing.y = unit(0.6, "lines"),
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_size_manual(values = 0.5) + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         size = guide_legend(order = 3), 
         color = guide_legend(order = 4),fill = guide_legend(order = 4))

#################################

# COMBINE #

#################################

pP.M <- ggarrange(pP.M, labels = c('a'), vjust = 1)
pNP.M <- ggarrange(pNP.M, labels = c('b'), vjust = 1)
p.ART <- ggarrange(p.ART, labels = c('c'), vjust = 1)
p.US <- ggarrange(p.US, labels = c('d'), vjust = 1)

p <- grid.arrange(pP.M, pNP.M, p.ART, p.US, layout_matrix = rbind(c(1, 2, 3, 4), 
                                                                  c(1, NA, NA, 4), 
                                                                  c(NA, NA, NA, 4)), 
                  heights = c(0.95, 0.05, 0.02), 
                  widths =  c(rep(0.22, 3),0.25))
ggsave(p, file = file.path(outdir, 'summary_estimates_art_suppressed.pdf'), w = 9, h = 12)



