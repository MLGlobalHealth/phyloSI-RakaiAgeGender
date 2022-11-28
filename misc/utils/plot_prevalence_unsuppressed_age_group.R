library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(ggpubr)

# directory repository
indir.repository <- '~/git/phyloflows'

# directory to save the figures
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

# files
file.path.round.timeline <- file.path(indir.repository, 'data', 'RCCS_round_timeline_220905.RData')
infile.unsuppressed <- file.path(indir.repository, 'fit', paste0('RCCS_propunsuppressed_age_group_221116.csv'))
infile.prevalence <- file.path(indir.repository, 'fit', paste0('RCCS_prevalence_age_group_221116.csv'))

load(file.path.round.timeline)
prevalence <- as.data.table(read.csv(infile.prevalence))
unsuppressed <- as.data.table(read.csv(infile.unsuppressed))

# format
tmp <- data.table(reshape2::melt(prevalence, id.vars = c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')))
tmp[, TYPE := 'HIV prevalence']
tmp1 <- data.table(reshape2::melt(unsuppressed, id.vars = c('ROUND', 'COMM', 'SEX', 'AGE_GROUP')))
tmp1[, TYPE := 'HIV-positive with\nunsuppressed viral load']
tmp <- rbind(tmp, tmp1)
tmp[, variable := stringi::stri_sub(variable,-2,-1)]
tmp[variable == '_M', variable := 'M']
tmp <- dcast.data.table(tmp, ROUND + COMM + SEX + AGE_GROUP + TYPE ~ variable, value.var = 'value')
tmp[, ROUND := paste0('R0', ROUND)]
tmp[, TYPE := factor(TYPE, levels = c('HIV prevalence', 'HIV-positive with\nunsuppressed viral load'))]

# round labels
df_round <- rbind(df_round_inland, df_round_fishing)
colnames(df_round) <- toupper(colnames(df_round))
df_round[, MIN_SAMPLE_DATE_LABEL := format(MIN_SAMPLE_DATE, '%b %Y')]
df_round[, MAX_SAMPLE_DATE_LABEL := format(MAX_SAMPLE_DATE - 31, '%b %Y')]
# df_round[, LABEL_ROUND := paste0('Round ', gsub('R0', '', ROUND), '\n', MIN_SAMPLE_DATE_LABEL, '-', MAX_SAMPLE_DATE_LABEL)]
df_round[, LABEL_ROUND := paste0('Round ', gsub('R0', '', ROUND))]
df_round[, LABEL_ROUND := factor(LABEL_ROUND, levels = unique(df_round[order(ROUND), LABEL_ROUND]))]
df_round[, INDEX_ROUND := 1:length(ROUND), by = 'COMM']
df_round[, MIDPOINT_DATE := MIN_SAMPLE_DATE + as.numeric(MAX_SAMPLE_DATE - MIN_SAMPLE_DATE)/2]

tmp <- merge(tmp, df_round, by = c('ROUND', 'COMM'))

# label sex and communities
tmp[, SEX_LABEL := 'Women']
tmp[SEX == 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Inland communities']
tmp[COMM == 'fishing', COMM_LABEL := 'Fishing communities']
tmp[, COMM_LABEL := factor(COMM_LABEL, levels = c('Inland communities', 'Fishing communities'))]
tmp[, AGE_LABEL := paste0('Age: ', AGE_GROUP)]

if(0){
  tmp <- tmp[!(COMM == 'inland' & ROUND %in% c('R010', 'R011'))]
}

# plot
communities <- tmp[, unique(COMM)]

for(i in seq_along(communities)){
  tmp1 <- tmp[COMM == communities[i]]
  
  p <- ggplot(tmp1, aes(x = MIDPOINT_DATE, group = interaction(TYPE, SEX_LABEL))) + 
    geom_line(aes(y = M, linetype = TYPE), position=position_dodge(width = 500), alpha = 1, col = 'grey50', size = 0.4) + 
    geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 1, width = 500, position=position_dodge(width = 500), col = 'grey30', size = 0.4) + 
    geom_point(aes(y = M, shape = TYPE, fill = SEX_LABEL), size = 1.5, position=position_dodge(width = 500), stroke = 0.1, col = 'black') + 
    facet_grid(.~AGE_LABEL) + 
    scale_color_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    scale_linetype_manual(values = c('HIV prevalence' = 'dotted', 'HIV-positive with\nunsuppressed viral load' = 'solid')) + 
    scale_shape_manual(values = c('HIV prevalence' = 24, 'HIV-positive with\nunsuppressed viral load' = 21)) + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)), 
          # plot.title = element_text(hjust = 0.5), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          axis.text.x = element_text(angle = 30, hjust = 1),
          # axis.title.x = element_blank(), 
          legend.position = c(0.12, 0.8),
          # legend.position = 'bottom', 
          legend.key.size = unit(0.33, 'cm'),
          # legend.box = 'vertical',
          legend.direction = 'vertical',
          legend.title = element_blank(), 
          legend.spacing.y= unit(0.00001, 'cm')) + 
    scale_y_continuous(labels = scales::percent, limits = c(0,  tmp1[, max(CU)]), expand = expansion(mult = c(0, 0.1))) + 
    labs(y = paste0('Proportion in the corresponding age bracket\nof the census-eligible population'), col= '', shape = '', linetype = '', 
         x = 'Date (midpoint of survey interval)') + 
    scale_x_date(expand = c(0.03,0.03)) +
    guides(fill = guide_legend(override.aes = list(shape = 21), order = 1), 
           linetype = guide_legend(order = 2), 
           shape = guide_legend(order = 2, override.aes = list(fill = 'black')))
  ggsave(p, file = file.path(outdir, paste0('prevalence_unsuppressed_among_census_eligible_', communities[i], '_221116.png')), w = 8,h = 3.5)
  ggsave(p, file = file.path(outdir, paste0('prevalence_unsuppressed_among_census_eligible_', communities[i], '_221116.pdf')), w = 8,h = 3.5)
}

# find statistics
tmp1 <- tmp[COMM == 'inland' & ROUND == 'R018']
tmp1 <- dcast.data.table(tmp1, SEX + AGE_GROUP ~ TYPE, value.var = 'M')
tmp1[, UNSUPPRESSED_TO_PREVALENCE_RATIO := `HIV-positive with\nunsuppressed viral load` / `HIV prevalence` ]
tmp2 <- tmp1[order(AGE_GROUP, SEX), .(AGE_GROUP, SEX, UNSUPPRESSED_TO_PREVALENCE_RATIO)]
tmp2

tmp1[, PREVALENCE_TO_UNSUPPRESSED_RATIO := `HIV prevalence` / `HIV-positive with\nunsuppressed viral load`  ]
tmp2 <- tmp1[order(AGE_GROUP, SEX), .(AGE_GROUP, SEX, PREVALENCE_TO_UNSUPPRESSED_RATIO)]
tmp2

