library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(gridExtra)
library("haven")
library(ggpubr)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

# directory to save the figure
outdir <- file.path('/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live','PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')
if(usr == 'melodiemonod'){
  outdir <- file.path('~/Box\ Sync/2021/ratmann_deepseq_analyses/live/', 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# check files exist
file.exists(c(file.treatment.cascade.population))  |> all() |> stopifnot()

# treatment cascade population
ns  <- as.data.table( read.csv(file.treatment.cascade.population) )



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

# save
file.name <- file.path(outdir, 'treatment_cascade_inland_population_221208.pdf')
if (! file.exists(file.name) || config$overwrite.existing.files) {
  cat("Saving file:", file.name, "\n")
  ggsave(file = file.name, w = 6.5, h = 4.5)  
} else {
  cat("File:", file.name, "already exists...\n")
}


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

# save
file.name <- file.path(outdir, 'population_suppression_rate_221208.pdf')
if (! file.exists(file.name) || config$overwrite.existing.files) {
  cat("Saving file:", file.name, "\n")
  ggsave(p, file = file.name, w = 8, h = 8)
} else {
  cat("File:", file.name, "already exists...\n")
}


####################################

# PLOT SUPPRESSION RATE ROUND 18

####################################

tab <- copy(ns[variable %in% c('PROP_SUPPRESSED')])
tab[,SEX_LABEL := 'Men']
tab[SEX == 'F',SEX_LABEL := 'Women']
tab[,COMM_LABEL := 'Inland communities']
tab[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
tab[, ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

tmp <- tab[ COMM == 'inland' & ROUND == c( 'R018')]
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_line(aes(y = M, col = SEX_LABEL)) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of individuals with HIV\nwho have suppressed virus in round 18', 
       col = '', fill = '', shape = '', linetype= '') +
  theme(legend.position = c(0.85, 0.17), 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) + 
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))

# save
file.name <- file.path(outdir, 'population_suppression_rate_R18_221213.pdf')
if (! file.exists(file.name) || config$overwrite.existing.files) {
  cat("Saving file:", file.name, "\n")
  ggsave(file = file.name, w = 4, h = 3.2)
} else {
  cat("File:", file.name, "already exists...\n")
}


####################################

# PLOT SUPPRESSION RATE GIVEN ART 

####################################


# add label
tab <- copy(ns[variable %in% c('PROP_SUPPRESSION_GIVEN_ART')])
tab[, VARIABLE_LEVEL := suppressed.label]

tab[,SEX_LABEL := 'Men']
tab[SEX == 'F',SEX_LABEL := 'Women']
tab[,COMM_LABEL := 'Inland communities']
tab[COMM == 'fishing',COMM_LABEL := 'Fishing communities']

tmp <- tab[COMM == 'inland' & ROUND %in% c('R015', 'R016', 'R017', 'R018')]
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_hline(aes(yintercept = 0.9), linetype='dashed', alpha = 0.5) +
  geom_hline(aes(yintercept = 0.95), linetype='dashed', alpha = 0.5) +
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = ROUND), alpha = 0.25) + 
  geom_line(aes(y = M, col = ROUND)) + 
  facet_grid(~SEX_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = "proportion of infected reporting ART use\nwho had a viral load < 1,000 cps / mL blood", 
       col = '', fill = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  ggsci::scale_color_jama()+
  ggsci::scale_fill_jama()+
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  guides(color = guide_legend(byrow = T, nrow = 3),
         fill = guide_legend(byrow = T, nrow = 3))

# save
file.name <- file.path(outdir, 'population_suppression_rate_overtime_221208.pdf')
if (! file.exists(file.name) || config$overwrite.existing.files) {
  cat("Saving file:", file.name, "\n")
  ggsave(file = file.name, w = 6, h = 5)
} else {
  cat("File:", file.name, "already exists...\n")
}



####################################

# SAVE PROPORTION OF SUPPRESSED

####################################

n_digit <- 1

tab <- ns[variable == 'PROP_SUPPRESSED']
tab <- tab[COMM == 'inland' & ROUND == 'R018']
tab <- tab[, .(ROUND, SEX, AGEYRS, M, CL, CU)]
tab <- tab[, `:=` (M = format(round(M*100, n_digit), n_digit), 
                   CL  = format(round(CL*100, n_digit), n_digit), 
                   CU = format(round(CU*100, n_digit), n_digit))]
tab <- tab[order(ROUND, SEX, AGEYRS)]

# save
file.name <- file.path(outdir, 'prop_suppressed.rds')
if (! file.exists(file.name) || config$overwrite.existing.files) {
  cat("Saving file:", file.name, "\n")
  saveRDS(tab, file = file.name)
} else {
  cat("File:", file.name, "already exists...\n")
}



