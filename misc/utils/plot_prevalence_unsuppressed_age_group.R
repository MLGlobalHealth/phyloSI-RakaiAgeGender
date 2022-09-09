library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(ggpubr)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

file.path.round.timeline <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'RCCS_round_timeline_220905.RData')
infile.unsuppressed <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_art_coverage_age_group_220906.csv'))
infile.prevalence <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('RCCS_prevalence_age_group_220830.csv'))

load(file.path.round.timeline)
prevalence <- as.data.table(read.csv(infile.prevalence))
unsuppressed <- as.data.table(read.csv(infile.unsuppressed))

# format
tmp <- data.table(reshape2::melt(prevalence, id.vars = c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')))
tmp[, TYPE := 'HIV-1 prevalence']
tmp1 <- data.table(reshape2::melt(unsuppressed, id.vars = c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')))
tmp1[, TYPE := 'HIV-positive with unsuppressed viral load']
tmp <- rbind(tmp, tmp1)
tmp[, variable := stringi::stri_sub(variable,-2,-1)]
tmp[variable == '_M', variable := 'M']
tmp <- dcast.data.table(tmp, ROUND + COMM + SEX + AGE_GROUP + TYPE ~ variable, value.var = 'value')
tmp[, ROUND := paste0('R0', ROUND)]
tmp[, TYPE := factor(TYPE, levels = c('HIV-1 prevalence', 'HIV-positive with unsuppressed viral load'))]

# round labels
df_round <- rbind(df_round_inland, df_round_fishing)
colnames(df_round) <- toupper(colnames(df_round))
df_round[, MIN_SAMPLE_DATE_LABEL := format(MIN_SAMPLE_DATE, '%b %Y')]
df_round[, MAX_SAMPLE_DATE_LABEL := format(MAX_SAMPLE_DATE - 31, '%b %Y')]
# df_round[, LABEL_ROUND := paste0('Round ', gsub('R0', '', ROUND), '\n', MIN_SAMPLE_DATE_LABEL, '-', MAX_SAMPLE_DATE_LABEL)]
df_round[, LABEL_ROUND := paste0('Round ', gsub('R0', '', ROUND))]
df_round[, LABEL_ROUND := factor(LABEL_ROUND, levels = df_round[order(ROUND), LABEL_ROUND])]
df_round[, INDEX_ROUND := 1:length(ROUND), by = 'COMM']

tmp <- merge(tmp, df_round, by = c('ROUND', 'COMM'))

# label sex and communities
tmp[, SEX_LABEL := 'Female']
tmp[SEX == 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Inland communities']
tmp[COMM == 'fishing', COMM_LABEL := 'Fishing communities']
tmp[, COMM_LABEL := factor(COMM_LABEL, levels = c('Inland communities', 'Fishing communities'))]
tmp[, AGE_LABEL := paste0('Age: ', AGE_GROUP)]

# plot
communities <- tmp[, unique(COMM)]

for(i in seq_along(communities)){
  tmp1 <- tmp[COMM == communities[i]]
  
  # prevalence
  p1 <- ggplot(tmp1[TYPE == 'HIV-1 prevalence'], aes(x = INDEX_ROUND, group = SEX_LABEL)) + 
    geom_line(aes(y = M), linetype = 'dashed', position=position_dodge(width = 0.5), alpha = 0.5) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0.5, position=position_dodge(width = 0.5)) + 
    geom_point(aes(y = M, col = SEX_LABEL), shape = 17, size = 2, position=position_dodge(width = 0.5)) + 
    facet_grid(.~AGE_LABEL) + 
    scale_color_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)), 
          axis.text.x = element_blank(),, 
          axis.title.x = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          legend.position = 'bottom', 
          legend.box = 'vertical', 
          legend.title = element_blank(), 
          legend.spacing.y= unit(0.00001, 'cm')) + 
    scale_y_continuous(labels = scales::percent, limits = c(0,  tmp1[TYPE == 'HIV-1 prevalence', max(CU)]), expand = expansion(mult = c(0, 0.1))) + 
    labs(y = 'HIV-1 prevalence\nin census eligible population', col= '', shape = '') + 
    scale_x_continuous(labels = df_round[COMM == communities[i], LABEL_ROUND], breaks =  df_round[COMM == communities[i], INDEX_ROUND]) + 
    guides(color = guide_legend(override.aes = list(shape = 16)))
  
  # unsuppressed
  p2 <- ggplot(tmp1[TYPE == 'HIV-positive with unsuppressed viral load'], aes(x = INDEX_ROUND, group = SEX_LABEL)) + 
    geom_line(aes(y = M), linetype = 'dotted', position=position_dodge(width = 0.5), alpha = 0.5) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0.5, position=position_dodge(width = 0.5)) + 
    geom_point(aes(y = M, col = SEX_LABEL), shape = 15, size = 2, position=position_dodge(width = 0.5)) + 
    facet_grid(.~AGE_LABEL) + 
    scale_color_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_blank(), 
          axis.text.x = element_text(angle = 50, hjust = 1), 
          axis.title.x = element_blank(), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          legend.position = 'bottom', 
          legend.box = 'vertical', 
          legend.title = element_blank(), 
          legend.spacing.y= unit(0.00001, 'cm')) + 
    scale_y_continuous(labels = scales::percent, limits = c(0,  tmp1[TYPE == 'HIV-1 prevalence', max(CU)]), expand = expansion(mult = c(0, 0.1))) + 
    labs(y = 'HIV+ with unsuppressed viral load\nin census eligible population', col= '', shape = '') + 
    scale_x_continuous(labels = df_round[COMM == communities[i], LABEL_ROUND], breaks =  df_round[COMM == communities[i], INDEX_ROUND])
  
  
  p <- ggarrange(p1, p2, ncol = 1, heights = c(0.45, 0.55), common.legend = T, legend = 'bottom')
  ggsave(p, file = file.path(outdir, paste0('prevalence_unsuppressed_among_census_eligible_', communities[i], '2.png')), w = 7.5,h = 7)
  
  
  ### try to combine them
  p <- ggplot(tmp1, aes(x = INDEX_ROUND, group = interaction(TYPE, SEX_LABEL))) + 
    geom_line(aes(y = M, linetype = TYPE), position=position_dodge(width = 0.5), alpha = 0.5) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5, width = 0.5, position=position_dodge(width = 0.5)) + 
    geom_point(aes(y = M, col = SEX_LABEL, shape = TYPE), size = 2, position=position_dodge(width = 0.5), stroke = 1) + 
    facet_grid(.~AGE_LABEL) + 
    scale_color_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
    scale_linetype_manual(values = c('HIV-1 prevalence' = 'dotted', 'HIV-positive with unsuppressed viral load' = 'solid')) + 
    scale_shape_manual(values = c('HIV-1 prevalence' = 1, 'HIV-positive with unsuppressed viral load' = 16)) + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          axis.text.x = element_text(angle = 30, hjust = 1), 
          axis.title.x = element_blank(), 
          legend.position = c(0.165, 0.85), 
          legend.key.size = unit(0.33, 'cm'),
          # legend.box = 'vertical', 
          legend.title = element_blank(), 
          legend.spacing.y= unit(0.00001, 'cm')) + 
    scale_y_continuous(labels = scales::percent, limits = c(0,  tmp1[, max(CU)]), expand = expansion(mult = c(0, 0.1))) + 
    labs(y = paste0('Percent among census eligible population\nin ',communities[i], ' communities'), col= '', shape = '', linetype = '') + 
    scale_x_continuous(labels = df_round[COMM == communities[i], LABEL_ROUND], breaks =  df_round[COMM == communities[i], INDEX_ROUND]) + 
    guides(color = guide_legend(override.aes = list(shape = 16), order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 2))
  ggsave(p, file = file.path(outdir, paste0('prevalence_unsuppressed_among_census_eligible_', communities[i], '.png')), w = 9,h = 4)
}


