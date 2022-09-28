library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(ggpubr)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'artcoverage_by_gender_loc_age')

file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_artcoverage_ratio_sex_220926.csv'))

unsuppressed_ratio <- as.data.table(read.csv(file.name))

# label sex and communities
unsuppressed_ratio[, COMM_LABEL := 'Inland communities']
unsuppressed_ratio[COMM == 'fishing', COMM_LABEL := 'Fishing communities']
unsuppressed_ratio[, COMM_LABEL := factor(COMM_LABEL, levels = c('Inland communities', 'Fishing communities'))]
unsuppressed_ratio[, AGE_LABEL := paste0('Age: ', AGE_GROUP)]
unsuppressed_ratio[, ROUND_LABEL := paste0('Round ', ROUND)]
unsuppressed_ratio <- unsuppressed_ratio[!(COMM == 'inland' & ROUND == '15S')]

# plot
communities <- unsuppressed_ratio[, unique(COMM)]

for(i in seq_along(communities)){
  tmp1 <- unsuppressed_ratio[COMM == communities[i]]
  
  p<-ggplot(tmp1, aes(x = ROUND_LABEL, group = AGE_GROUP)) + 
    geom_hline(yintercept = 1, linetype = 'dashed', alpha= 0.5) + 
    geom_line(aes(y = UNSUPPRESSION_RATE_RATIO_BY_AGE_M),position=position_dodge(width = 0.5), alpha = 0.5) + 
    geom_errorbar(aes(ymin = UNSUPPRESSION_RATE_RATIO_BY_AGE_CL, ymax = UNSUPPRESSION_RATE_RATIO_BY_AGE_CU), alpha = 0.35, width = 0.5, position=position_dodge(width = 0.5)) + 
    geom_point(aes(y = UNSUPPRESSION_RATE_RATIO_BY_AGE_M, col = AGE_GROUP, shape = AGE_GROUP),  size = 2, position=position_dodge(width = 0.5)) + 
    scale_color_viridis_d(option = 'D') + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 30, hjust = 1), 
          plot.title = element_text(hjust = 0.5), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          legend.position = 'bottom', 
          legend.box = 'vertical', 
          # legend.title = element_blank(), 
          legend.spacing.y= unit(0.00001, 'cm')) + 
    scale_y_continuous( limits = c(0,  NA), expand = expansion(mult = c(0, 0.1))) + 
    labs(y = 'ART naive (unsuppressed) rate male to female ratio', col= 'Age', shape= 'Age') 
    ggsave(p, file = file.path(outdir, paste0('unsuppressed_rate_ratio_', communities[i], '.png')), w = 5,h = 5)
  

}


