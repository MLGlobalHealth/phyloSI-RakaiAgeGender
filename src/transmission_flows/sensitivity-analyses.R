library(data.table)
library(ggplot2)
library(gridExtra)

outdir.all <- "/rds/general/user/mm3218/home/projects/2021/phyloSI-RakaiAgeGender/"
# outdir.all <- "~/Box\ Sync/2021/phyloflows/"

# jobames
jobname <- c('usedep',  # central analysis
             'incloess', # incidence rate estimated with loess regression
             '30com', # incidence rate estimated in 30 communities
             'tsinonrefined', # non refined TSI
             'woR18', # remove pairs after round 17
             'woR1718', # remove pairs after round 16
             'woR161718', # remove pairs after round 15
             'seed12', # bootstrap sample of the pairs 1
             'seed13', # bootstrap sample of the pairs 2
             'seed14', # bootstrap sample of the pairs 3
             'usedep', # alternative form transmission rate 1
             'usedep', # alternative form transmission rate 2
             'firstrun', # alternative form transmission rate 3
             'firstrun', # alternative form transmission rate 4
             'firstrun', # alternative form transmission rate 5
             'firstrun', # alternative form transmission rate 6
             'nonparttreatedaspart', # non participants treated liked participants
             'nonpartnottreated', # non participants all unsupressed
             'nonpartfemalemale125moreinfectious', # non participants 25% higher prevalence
             'nonpartmale125moreinfectious',  # non participants men 25% higher prevalence
             'nonpartfemale125moreinfectious', # non participants women 25% higher prevalence
             'vl200', # viral load threshold 200ml
             'central' # not adjusting for sampling probability of source
             )

# stan models
stan_model <- c('gp_230602', 
                'gp_230602', 
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230614c',
                'gp_230614d',
                'gp_230614e',
                'gp_230614f',
                'gp_230614a',
                'gp_230614b',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_230602',
                'gp_221201d'
                )

# labels
label_models <- c('Central model', 
                  'Alternative model: incidence rates estimated with LOESS regression',
                  'Alternative model: incidence rates estimated on data subset to 28 continuously surveyed communities',
                  'Alternative model: non refined TSI estimates',
                  'Alternative model: without source-recipients pairs for which the source or recipient was sequenced after round 17',
                  'Alternative model: without source-recipients pairs for which the source or recipient was sequenced after round 16',
                  'Alternative model: without source-recipients pairs for which the source or recipient was sequenced after round 15',
                  'Alternative model: source-recipients pairs bootstrap sample 1',
                  'Alternative model: source-recipients pairs bootstrap sample 2',
                  'Alternative model: source-recipients pairs bootstrap sample 3',
                  'Alternative model: transmission rate form (11e)',
                  'Alternative model: transmission rate form (11f)',
                  'Alternative model: transmission rate form (11a)',
                  'Alternative model: transmission rate form (11b)',
                  'Alternative model: transmission rate form (11c)',
                  'Alternative model: transmission rate form (11d)',
                  'Alternative model: infected non-participants are assumed to be suppressed as the participants',
                  'Alternative model: infected non-participants are assumed to undiagnosed and not suppressed',
                  'Alternative model: 25% higher prevalence in non-participants',
                  'Alternative model: 25% higher prevalence in men non-participants',
                  'Alternative model: 25% higher prevalence in women non-participants',
                  'Alternative model: viral load threshold of 200ml',
                  'Alternative model: without adjustments for potentially unequal sampling of sources'
                  )


# log_lambda_latentby_direction_round_quantilestandardisedby_direction_round.rds
# z_predictby_direction_roundstandardisedby_round.rds
# counterfactuals_*
#log_betaby_direction_round.rds

outdir.here <- file.path(outdir.all, paste0(stan_model, '-', jobname))
N <- length(jobname)
# data_in <- file.path(outdir.here, paste0(stan_model, '-', jobname, '-stanin_', jobname, '.RData'))

prettyround <- function(x,n) format(round(x, n), nsmall = n)
results_sensitivity <- list()
cols <- c('black', viridis::viridis_pal()(N-1))

#
# load contribution from male
#

paths <- file.path(outdir.here, 'tables',
                   paste0(stan_model, '-', jobname, '-output-log_lambda_latentby_direction_roundstandardisedby_round.rds'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  output[[i]] <- as.data.table(readRDS(paths[i]))
  output[[i]][, type := label_models[i]]
}
output <- rbindlist(output)
output[, type := factor(type, levels= label_models)]

# keep only male contributions
output <- output[LABEL_SOURCE == 'Male sources' & COMM == 'inland']

# y axis label
type_cont <- 'Contribution from male sources\nto HIV incidence'

# plot
pdf(file.path(outdir.here[1], 'figures', paste0(stan_model[1], '-', jobname[1], '-sensitivity-contribution_male_sources.pdf')), onefile = TRUE, w = 8, h = 6)
for(i in 2:N){
  tmp <- output[type %in% label_models[c(1, i)]]
  p <- ggplot(tmp, aes(x = INDEX_ROUND, group = type)) +
    geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0, size = 0.5, position = position_dodge(width = 0.4))  +
    geom_point(aes(y = M, col = type), size= 1.7, position = position_dodge(width = 0.4)) + 
    labs(x = '', y = type_cont, col = '') + 
    theme_bw() +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title.x = element_blank(), 
          legend.position = 'bottom',
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          legend.margin = margin(),
          legend.title = element_blank())  + 
    scale_x_continuous(labels = output[order(INDEX_ROUND), unique(LABEL_ROUND)], 
                       breaks = output[order(INDEX_ROUND), unique(INDEX_ROUND)]) + 
    scale_color_manual(values = cols[c(1, i)]) + 
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0)), limits = c(0, 1)) + 
    guides(color = guide_legend(order = 1)) 
  do.call("grid.arrange", list(p)  )
}
dev.off()

# save for table
tmp <- output[ROUND %in% c('R010', 'R014', 'R018')]
tmp[, M := paste0(prettyround(M*100, 1), '\\%')]
tmp[, CI := paste0( '[', prettyround(CL*100, 1), '-', prettyround(CU*100, 1), ']')]
tmp1 <- dcast(tmp, type ~ ROUND, value.var = 'M')
tmp1 <- merge(tmp1, dcast(tmp, type ~ ROUND, value.var = 'CI'), by = 'type')
results_sensitivity[['contribution_male']] <- tmp1


#
# load median age of transmitter by age of the recipient
#

paths <- file.path(outdir.here, 'tables',
                   paste0(stan_model, '-', jobname, '-output-log_lambda_latentby_direction_round_quantilestandardisedby_direction_round_median_age_source_total.rds'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  output[[i]] <- as.data.table(readRDS(paths[i]))
  output[[i]][, type := label_models[i]]
}
output <- rbindlist(output)
output[, type := factor(type, levels= label_models)]

# keep only male contributions
output <- output[COMM == 'inland']

# dcast
mag <- dcast.data.table(output, ROUND + LABEL_ROUND + COMM + LABEL_RECIPIENT + type ~ quantile, value.var = 'M')
magm <- output[quantile == 'C50']

# plot
pdf(file.path(outdir.here[1], 'figures', paste0(stan_model[1], '-', jobname[1], '-sensitivity-MedianAgeSource.pdf')), onefile = TRUE, w = 8, h = 6)
for(i in 2:N){
  tmp <- mag[type %in% label_models[c(1, i)]]
  p <- ggplot(tmp) + 
    geom_boxplot(stat = "identity",
                 aes(x = LABEL_ROUND, fill = type,
                     lower  = C25,upper = C75, middle = C50, ymin = C10, ymax = C90, 
                     group= interaction(type, LABEL_ROUND)),
                 col = 'grey70', outlier.shape = NA, 
                 size = 0.2, position = position_dodge2(width = 0.4)) +
    facet_grid(.~LABEL_RECIPIENT) + 
    theme_bw() +
    labs(x = '', y = 'Age source') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom', 
          legend.direction = 'vertical',
          axis.text.x = element_text(angle= 40, hjust =1),
          panel.grid.minor.x = element_blank(), 
          legend.spacing.x = unit(0.1, 'cm'),
          legend.title = element_blank()) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = cols[c(1, i)]) + 
    guides(color = guide_legend(order = 1))
  do.call("grid.arrange", list(p)  )
}
dev.off()

# save for table
tmp <- mag[ROUND %in% c('R010', 'R014', 'R018')]
tmp[, RESULTS := paste0(prettyround(C50, 1))]
tmp[, CI := paste0( ' [', prettyround(C10, 1), '-', prettyround(C90, 1), ']')]
tmp1 <- dcast(tmp[LABEL_RECIPIENT == 'Female recipients'], type ~ ROUND, value.var = 'RESULTS')
tmp1 <- merge(tmp1, dcast(tmp[LABEL_RECIPIENT == 'Female recipients'], type ~ ROUND, value.var = 'CI'), by = 'type')
results_sensitivity[['age_male']] <- tmp1
tmp1 <- dcast(tmp[LABEL_RECIPIENT == 'Male recipients'], type ~ ROUND, value.var = 'RESULTS')
tmp1 <- merge(tmp1, dcast(tmp[LABEL_RECIPIENT == 'Male recipients'], type ~ ROUND, value.var = 'CI'), by = 'type')
results_sensitivity[['age_female']] <- tmp1


#
# load median age of transmitter by 1-year age of the recipient
#

paths <- file.path(outdir.here, 'tables',
                   paste0(stan_model, '-', jobname, '-output-log_lambda_latentby_direction_round_age_infection.recipient_quantilestandardisedby_direction_round_age_infection.recipient.rds'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  output[[i]] <- as.data.table(readRDS(paths[i]))
  output[[i]][, type := label_models[i]]
}
output <- rbindlist(output)
output[, type := factor(type, levels= label_models)]

# keep only male contributions
output <- output[COMM == 'inland' & LABEL_RECIPIENT == 'Female recipients']
magm <- magm[COMM == 'inland' & LABEL_RECIPIENT == 'Female recipients']

# keep up to round 15
output <- output[!ROUND %in% c('R016', 'R017', 'R018')]
magm <- magm[!ROUND %in% c('R016', 'R017', 'R018')]

# keep only counterfactual with right censoring
output <- output[grepl('Central|without source-recipients',type)]
magm <- magm[grepl('Central|without source-recipients',type)]
cols_res <- c('black', ggsci::pal_npg()(magm[, length(unique(type)) - 1]))

# plot
ggplot(output) + 
  geom_abline(slope = 1, intercept= 0, linetype = 'dashed', col = 'grey50') + 
  geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin = CL, ymax = CU, fill = type), alpha = 0.2) + 
  geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = type)) + 
  geom_errorbarh(data = magm, aes(xmin = CL, xmax = CU, y = 17, col= type, fill = type), size = 1.5) + 
  geom_point(data = magm, aes(x = M, y = 19, col= type, fill = type), shape = 25, size = 2) + 
  facet_grid(LABEL_ROUND~type) + 
  theme_bw() + 
  theme(legend.position = 'bottom', 
        legend.direction = 'vertical', 
        strip.background = element_rect(colour="white", fill="white"), 
        strip.text.x = element_blank()) + 
  labs(x = 'Age recipient', y = 'Age source', col = '', fill ='') + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = cols_res) +
  scale_fill_manual(values = cols_res)
ggsave(file = file.path(outdir.here[1], 'figures', paste0(stan_model[1], '-', jobname[1], '-sensitivity-MedianAgeSourceBy1yearAgeRecipient.png')), 
       w = 8, h = 10)
    


#
# load counterfactual budget
#

scenario.half.female = 'Closing half the gap in men\nrelative to women participants/new participants'
scenario.female = 'Closing the gap in men\nrelative to women participants/new participants'
scenario.959595 = '95-95-95 in men'

paths.1 <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-output-counterfactuals_9595951.rds'))
paths.2 <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-output-counterfactuals_art_up_to_female0.5.rds'))
paths.3 <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-output-counterfactuals_art_up_to_female1.rds'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  tmp <- readRDS(paths.1[i])
  output[[i]] <- cbind(tmp$budget[ROUND == 'R018' & SEX == 'M'], scenario = scenario.959595)
  
  tmp <- readRDS(paths.2[i])
  tmp <- cbind(tmp$budget[ROUND == 'R018' & SEX == 'M'], scenario = scenario.half.female)
  output[[i]] <- rbind(output[[i]], tmp)
  
  tmp <- readRDS(paths.3[i])
  tmp <- cbind(tmp$budget[ROUND == 'R018' & SEX == 'M'], scenario = scenario.female)
  output[[i]] <- rbind(output[[i]], tmp)
  
  output[[i]][, type := label_models[i]]
}
output <- rbindlist(output)
output[, type := factor(type, levels= label_models)]
output[, scenario := factor(scenario, levels = c(scenario.half.female, scenario.female, scenario.959595))]

# plot
pdf(file.path(outdir.here[1], 'figures', paste0(stan_model[1], '-', jobname[1], '-sensitivity-counterfactual_budget.pdf')), onefile = TRUE, w = 8, h = 6)
for(i in 2:N){
  tmp <- output[type %in% label_models[c(1, i)]]
  p <- ggplot(tmp, aes(x = scenario, group = type)) +
    geom_bar(aes(y = TREATED, fill = type), width= 0.7, position = position_dodge(width = 0.8), stat = 'identity') + 
    geom_errorbar(aes(ymin = TREATED_CL, ymax = TREATED_CU), width= 0.1, position = position_dodge(width = 0.8), col = 'grey50') + 
    labs(x = '', y = 'Additional number of men suppressed', col = '') + 
    theme_bw() +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title.x = element_blank(), 
          legend.position = 'bottom',
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          legend.margin = margin(),
          legend.title = element_blank())   + 
    scale_fill_manual(values = cols[c(1, i)]) + 
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, 0.01)))
  do.call("grid.arrange", list(p)  )
}
dev.off()

# save for table
tmp <- output[ROUND %in% c('R010', 'R014', 'R018')]
tmp[, RESULTS := paste0(prettyround(TREATED, 1))]
tmp[, CI := paste0( ' [', prettyround(TREATED_CL, 1), '-', prettyround(TREATED_CU, 1), ']')]
tmp[, CI := gsub(' ', '', CI)]
tmp1 <- dcast(tmp, type ~ scenario, value.var = 'RESULTS')
tmp1 <- merge(tmp1, dcast(tmp, type ~ scenario, value.var = 'CI'), by = 'type')
results_sensitivity[['budget']] <- tmp1


#
# load counterfactual results
#

paths.1 <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-output-counterfactuals_9595951.rds'))
paths.2 <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-output-counterfactuals_art_up_to_female0.5.rds'))
paths.3 <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-output-counterfactuals_art_up_to_female1.rds'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  tmp <- readRDS(paths.1[i])
  output[[i]] <- cbind(tmp$relative_incidence_counterfactual_all[ROUND == 'R018' & LABEL_RECIPIENT == 'Female recipients'], scenario = scenario.959595)
  
  tmp <- readRDS(paths.2[i])
  tmp <- cbind(tmp$relative_incidence_counterfactual_all[ROUND == 'R018' & LABEL_RECIPIENT == 'Female recipients'], scenario = scenario.half.female)
  output[[i]] <- rbind(output[[i]], tmp)
  
  tmp <- readRDS(paths.3[i])
  tmp <- cbind(tmp$relative_incidence_counterfactual_all[ROUND == 'R018' & LABEL_RECIPIENT == 'Female recipients'], scenario = scenario.female)
  output[[i]] <- rbind(output[[i]], tmp)
  
  output[[i]][, type := label_models[i]]
}
output <- rbindlist(output)
output[, type := factor(type, levels= label_models)]
output[, scenario := factor(scenario, levels = c(scenario.half.female, scenario.female, scenario.959595))]

# plot
pdf(file.path(outdir.here[1], 'figures', paste0(stan_model[1], '-', jobname[1], '-sensitivity-counterfactual_reduction_incidence.pdf')),, onefile = TRUE, w = 8, h = 6)
for(i in 2:N){
  tmp <- output[type %in% label_models[c(1, i)]]
  p <- ggplot(tmp, aes(x = scenario, group = type)) +
    geom_bar(aes(y = M, fill = type), width= 0.7, position = position_dodge(width = 0.8), stat = 'identity') + 
    geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0.2, size = 0.5, position = position_dodge(width = 0.8))  +
    labs(x = '', y = '% reduction in incidence in women', col = '') + 
    theme_bw() +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title.x = element_blank(), 
          legend.position = 'bottom',
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          legend.margin = margin(),
          legend.title = element_blank())   + 
    scale_fill_manual(values = cols[c(1, i)]) + 
    scale_y_continuous(labels = scales::percent, limits = c(0,NA), expand = expansion(mult = c(0, 0.01)))
  do.call("grid.arrange", list(p)  )
}
dev.off()


# save for table
output[, RESULTS := paste0(prettyround(M*100, 1), '\\%')]
output[, CI := paste0('[', prettyround(CL*100, 1), '-', prettyround(CU*100, 1), ']')]
tmp <- dcast(output, type ~ scenario, value.var = 'RESULTS')
tmp <- merge(tmp, dcast(output, type ~ scenario, value.var = 'CI'), by = 'type')
results_sensitivity[['counterfactual_redinc']] <- tmp



#
# load transmission risk from male
#

paths <- file.path(outdir.here, 'tables',
                   paste0(stan_model, '-', jobname, '-output-log_betaby_direction_round.rds'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  output[[i]] <- as.data.table(readRDS(paths[i]))
  output[[i]][, type := label_models[i]]
}
output <- rbindlist(output)
output[, type := factor(type, levels= label_models)]

# keep only male contributions
output <- output[LABEL_SOURCE == 'Male sources' & COMM == 'inland']

# y axis label
type_cont <- 'Transmission risk from male sources'

# plot
pdf(file.path(outdir.here[1], 'figures', paste0(stan_model[1], '-', jobname[1], '-sensitivity-transmission_risk_male_sources.pdf')), onefile = TRUE, w = 8, h = 6)
for(i in 2:N){
  tmp <- output[type %in% label_models[c(1, i)]]
  p <- ggplot(tmp, aes(x = INDEX_ROUND, group = type)) +
    geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0, size = 0.5, position = position_dodge(width = 0.8))  +
    geom_point(aes(y = M, col = type), size= 1.7, position = position_dodge(width = 0.8)) +
    labs(x = '', y = type_cont, col = '') +
    theme_bw() +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = 'bottom',
          legend.direction = 'vertical',
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.margin = margin(),
          legend.title = element_blank())  +
    scale_color_manual(values = cols[c(1, i)]) +
    scale_x_continuous(labels = output[order(INDEX_ROUND), unique(LABEL_ROUND)],
                       breaks = output[order(INDEX_ROUND), unique(INDEX_ROUND)]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    guides(color = guide_legend(order = 1))
  do.call("grid.arrange", list(p)  )
}
dev.off()

# save for table
tmp <- output[ROUND %in% c('R010', 'R014', 'R018')]
tmp[, RESULTS := paste0(prettyround(M*100, 2), '\\% [', prettyround(CL*100, 2), '-', prettyround(CU*100, 2), ']')]
tmp <- dcast(tmp, type ~ ROUND, value.var = 'RESULTS')
results_sensitivity[['transmission_risk_male']] <- tmp


###
# save table
###

saveRDS(results_sensitivity, file = file.path(outdir.here[1], 'tables', paste0(stan_model[1], '-', jobname[1], '-sensitivity-summary.rds')))

