library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(ggpubr)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'artcoverage_by_gender_loc_age')

file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_unsuppressed_ratio_sex_221101.csv'))
file.path.round.timeline <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'RCCS_round_timeline_220905.RData')

unsuppressed_ratio <- as.data.table(read.csv(file.name))
load(file.path.round.timeline)

# label sex and communities
unsuppressed_ratio[, COMM_LABEL := 'Inland communities']
unsuppressed_ratio[COMM == 'fishing', COMM_LABEL := 'Fishing communities']
unsuppressed_ratio[, COMM_LABEL := factor(COMM_LABEL, levels = c('Inland communities', 'Fishing communities'))]
unsuppressed_ratio[, AGE_LABEL := paste0('Age: ', AGE_GROUP)]
unsuppressed_ratio[, ROUND_LABEL := paste0('Round ', ROUND)]
unsuppressed_ratio <- unsuppressed_ratio[!(COMM == 'inland' & ROUND == '15S')]

# add midpoint date round
df_round <- rbind(df_round_inland, df_round_fishing)
unsuppressed_ratio[, round := paste0('R0', ROUND)]
unsuppressed_ratio <- merge(unsuppressed_ratio, df_round, by = c('round', 'COMM'))
unsuppressed_ratio[, MIDPOINT_DATE := min_sample_date + as.numeric(max_sample_date - min_sample_date)/2]

if(0){
  unsuppressed_ratio <- unsuppressed_ratio[!(COMM == 'inland' & round %in% c('R010', 'R011'))]
}

# plot
communities <- unsuppressed_ratio[, unique(COMM)]

for(i in seq_along(communities)){
  tmp1 <- unsuppressed_ratio[COMM == communities[i]]
  
  p<-ggplot(tmp1, aes(x = MIDPOINT_DATE, group = AGE_GROUP)) + 
    geom_hline(yintercept = 1, linetype = 'dashed', alpha= 0.5) + 
    geom_line(aes(y = UNSUPPRESSION_RATE_RATIO_BY_AGE_M),position=position_dodge(width = 300), alpha = 0.5) + 
    geom_errorbar(aes(ymin = UNSUPPRESSION_RATE_RATIO_BY_AGE_CL, ymax = UNSUPPRESSION_RATE_RATIO_BY_AGE_CU), alpha = 0.35, width = 300, position=position_dodge(width = 300)) + 
    geom_point(aes(y = UNSUPPRESSION_RATE_RATIO_BY_AGE_M, col = AGE_GROUP, shape = AGE_GROUP),  size = 2, position=position_dodge(width = 300)) + 
    scale_color_viridis_d(option = 'D') + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)), 
          # axis.title.x = element_blank(), 
          # axis.text.x = element_text(angle = 30, hjust = 1), 
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          legend.position = 'none', 
          # legend.title = element_blank(), 
          legend.spacing.y= unit(0.00001, 'cm')) + 
    scale_y_continuous( expand = expansion(mult = c(0.02, 0.1))) + 
    labs(y = 'Male-female ratio of the rate of infected\nindividuals that remain unsuppressed', col= 'Age', shape= 'Age', 
         x = 'Date (midpoint of survey interval)') 
    ggsave(p, file = file.path(outdir, paste0('unsuppressed_rate_ratio_', communities[i], '_221101.png')), w = 3.5,h = 3.15)
    ggsave(p, file = file.path(outdir, paste0('unsuppressed_rate_ratio_', communities[i], '_221101.pdf')), w = 3.5,h = 3.2)
    

}


