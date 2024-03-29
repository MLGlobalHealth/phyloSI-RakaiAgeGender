library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(ggpubr)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

# directory to save the figure
outdir <- file.path('/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live','PANGEA2_RCCS', 'artcoverage_by_gender_loc_age')
if(usr == 'melodiemonod'){
  outdir <- file.path('~/Box\ Sync/2021/ratmann_deepseq_analyses/live/', 'PANGEA2_RCCS', 'artcoverage_by_gender_loc_age')
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# check files exist
file.exists(c(
  file.unsuppressed_rate_ratio ,
  file.path.round.timeline))  |> all() |> stopifnot()

# load files
unsuppressed_ratio <- fread(file.unsuppressed_rate_ratio)
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
    geom_line(aes(y = UNSUPPRESSION_RATE_RATIO_BY_AGE_M),position=position_dodge(width = 300),col = 'grey50') + 
    geom_errorbar(aes(ymin = UNSUPPRESSION_RATE_RATIO_BY_AGE_CL, ymax = UNSUPPRESSION_RATE_RATIO_BY_AGE_CU), col = 'grey30', width = 300, position=position_dodge(width = 300)) + 
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
    labs(y = 'Male-female ratio of the proportion of individuals\nwith HIV who have unsuppressed virus\nrelative to round 10', col= 'Age', shape= 'Age', 
         x = 'Date (midpoint of survey interval)') 
  
  # save
  file.name <- file.path(outdir, paste0('unsuppressed_rate_ratio_', communities[i], '_221208.png'))
  if (! file.exists(file.name) || config$overwrite.existing.files) {
    cat("Saving file:", file.name, "\n")
    ggsave(p, file = file.name, w = 3.8,h = 3.15)
  } else {
    cat("File:", file.name, "already exists...\n")
  }
  
  file.name <- file.path(outdir, paste0('unsuppressed_rate_ratio_', communities[i], '_221208.pdf'))
  if (! file.exists(file.name) || config$overwrite.existing.files) {
    cat("Saving file:", file.name, "\n")
    ggsave(p, file = file.name, w = 4.1,h = 3.4)
  } else {
    cat("File:", file.name, "already exists...\n")
  }

}


