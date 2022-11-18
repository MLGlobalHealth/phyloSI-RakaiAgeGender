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
file.cascade.population <- file.path(indir.repository, 'fit', paste0('RCCS_treatment_cascade_population_estimates_221116.csv'))

# treatment cascade population
ns  <- as.data.table( read.csv(file.cascade.population) )



####################################

# PLOT TREATMENT CASCADE

####################################

# label
diagnosed.label = 'estimated proportion of diagnosed'
on.art.label = "proportion of infected who report ART use"
suppressed.label = "proportion of infected reporting ART use who had a viral load < 1,000 cps / mL blood"

# add label
tab <- copy(ns[variable %in% c('PROP_SUPPRESSION_GIVEN_ART',
                               'PROP_ART_COVERAGE_GIVEN_DIAGNOSED',
                               'PROP_DIAGNOSED')])
tab[, VARIABLE_LEVEL := diagnosed.label]
tab[variable == 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED', VARIABLE_LEVEL := on.art.label]
tab[variable == 'PROP_SUPPRESSION_GIVEN_ART', VARIABLE_LEVEL := suppressed.label]
tab[, VARIABLE_LEVEL := factor(VARIABLE_LEVEL, 
                               levels = c(diagnosed.label, on.art.label, suppressed.label))]
tab[,SEX_LABEL := 'Men']
tab[SEX == 'F',SEX_LABEL := 'Women']
tab[,COMM_LABEL := 'Inland communities']
tab[COMM == 'fishing',COMM_LABEL := 'Fishing communities']

# ROUND 18
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

# ALL ROUNDS
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

# PLOT SUPPRESSION RATE

####################################

tab <- copy(ns[variable %in% c('PROP_SUPPRESSION_GIVEN_DIAGNOSED')])
tab[,SEX_LABEL := 'Men']
tab[SEX == 'F',SEX_LABEL := 'Women']
tab[,COMM_LABEL := 'Inland communities']
tab[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
tab[, ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

tmp <- tab[ COMM == 'inland' & !ROUND %in% c( 'R015S')]
p1 <- ggplot(tmp, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_line(aes(y = M, col = SEX_LABEL)) + 
  facet_wrap(~ROUND_LABEL, ncol = 2) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Viral suppression rate in infected population', 
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
tab <- dcast.data.table(tmp, ROUND_LABEL + AGEYRS + COMM ~ SEX, value.var = 'M')
tab[, SUPPRESSION_RATE_RATIO := `F`/`M`]
tmp <- tab[!grepl('11|13|15|17', ROUND_LABEL)]
cols <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(9)[c(1,3,5,7,9)]
p2 <- ggplot(tmp, aes(x = AGEYRS)) + 
  geom_hline(aes(yintercept = 1), linetype = 'dashed', col = 'grey50') + 
  geom_line(aes(y = SUPPRESSION_RATE_RATIO, col = ROUND_LABEL)) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Female-to-male viral suppression rate\nratio in infected population', 
       col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = cols) + 
  scale_y_log10() + 
  scale_x_continuous(expand = c(0,0)) +
  guides(color = guide_legend(byrow = T, nrow = 2))
p2 <- ggarrange(p2, labels = c('b'))

p <- grid.arrange(p1,p2,layout_matrix = rbind(c(1, 2), c(1, NA)), heights =c(0.6, 0.4), widths = c(0.45, 0.55))
ggsave(p, file = file.path(outdir, 'population_suppression_rate_221118.pdf'), w = 8, h = 8)


