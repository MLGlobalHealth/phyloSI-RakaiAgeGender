plot_age_infection_source_recipient <- function(data, title, lab, outdir){
  
  data <- data[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
  data[, `Cohort round recipient` := cohort_round.RECIPIENT]
  data[, `Cohort round source` := cohort_round.SOURCE]
  data[, `Community recipient` := comm.RECIPIENT]
  data[, `Community source` := comm.SOURCE]
  
  # all pairs
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) 
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_AllPairs_', lab, '.png')), w = 4, h = 4)
  
  # by cohort round
  p1 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round source`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p2 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round recipient`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p <- ggarrange(p1, p2, ncol = 2)
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CohortRound_', lab, '.png')), w = 9, h = 7)
  
  # by age infection round
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = date_infection_before_UTT.RECIPIENT)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient',
         col = 'Date infection recipient before UTT') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_DateInfectionRecipient_', lab, '.png')), w = 5, h = 5)
  
  
  # by community
  p1 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Community source`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p2 <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point(aes(col =`Community recipient`)) + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')
  
  p <- ggarrange(p1, p2, ncol = 2)
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CommunityRecipient_', lab, '.png')), w = 9, h = 7)
  
}

plot_hist_age_infection <- function(pairs, outdir){
  
  pairs[, Sex := sex.SOURCE]
  p1 <- ggplot(pairs, aes(x = age_infection.SOURCE)) + 
    geom_histogram(bins = 30) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection source') +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  pairs[, Sex := sex.RECIPIENT]
  p2 <- ggplot(pairs, aes(x = age_infection.RECIPIENT)) + 
    geom_histogram(bins = 30) +     
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  
  file = file.path(outdir, paste0('hist_age_infection_', lab, '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
  
  return(p)
}

plot_hist_time_infection <- function(pairs, outdir){
  
  pairs[, `Round source` := cohort_round.SOURCE]
  p1 <- ggplot(pairs, aes(x = date_infection.SOURCE)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round source`, nrow = length(unique(pairs$cohort_round.SOURCE))) +
    theme_bw() + 
    labs(x = 'Date of infection source') + 
    geom_vline(xintercept = date_implementation_UTT, linetype = 'dashed')
  
  pairs[, `Round recipient` := cohort_round.RECIPIENT]
  p2 <- ggplot(pairs, aes(x = date_infection.RECIPIENT)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round recipient`, nrow = length(unique(pairs$cohort_round.RECIPIENT))) +
    theme_bw() + 
    labs(x = 'Date of infection recipient') + 
    geom_vline(xintercept = date_implementation_UTT, linetype = 'dashed')

  p <- ggarrange(p1, p2, ncol = 2)
  
  file = file.path(outdir, paste0('hist_date_infection_', lab, '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
}

plot_hist_age_infection_diff_threshold <- function(pairs, outdir){
  
  chain <- keep.likely.transmission.pairs(as.data.table(dchain), 0.5)
  pairs <- pairs.get.meta.data(chain, meta_data)
  pairs$threshold = '0.5'
  
  chain <- keep.likely.transmission.pairs(as.data.table(dchain), 0.6)
  pairs2 <- pairs.get.meta.data(chain, meta_data)
  pairs2$threshold = '0.6'
  
  pairs = rbind(pairs, pairs2)
  
  if(!include.mrc){
    cat('Keep only pairs in RCCS')
    pairs <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS']
  }
  if(include.only.heterosexual.pairs){
    cat('Keep only heterosexual pairs')
    pairs <- pairs[(sex.RECIPIENT == 'M' & sex.SOURCE == 'F') | (sex.RECIPIENT == 'F' & sex.SOURCE == 'M')]
  }
  
  pairs[, Sex := sex.SOURCE]
  p1 <- ggplot(pairs, aes(x = age_infection.SOURCE)) + 
    geom_density(aes( group = threshold, fill = threshold), alpha = 0.5) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection source') +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  pairs[, Sex := sex.SOURCE]
  p2 <- ggplot(pairs, aes(x = age_infection.RECIPIENT)) + 
    geom_density(aes( group = threshold, fill = threshold), alpha = 0.5) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  
  file = file.path(outdir, paste0('hist_age_infection_source_', gsub('(.+)_threshold.*', '\\1', lab), '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
  
  return(p)
}

plot_CI_age_infection <- function(pairs, outdir){
  
  tmp <- copy(pairs)
  tmp[, age_infection.SOURCE := floor(age_infection.SOURCE)]
  tmp[, age_infection.RECIPIENT := floor(age_infection.RECIPIENT)]
  tmp <- merge(tmp, df_age, by = c('age_infection.RECIPIENT', 'age_infection.SOURCE'))
  
  ps <- c(0.5, 0.2, 0.8)
  p_labs <- c('M','CL','CU')
  
  quantile(subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M' & age_infection_reduced.RECIPIENT == 46)$age_infection.SOURCE, prob=ps)
  
  tmp = tmp[, list(q= quantile(age_infection.SOURCE, prob=ps, na.rm = T), q_label=p_labs), 
            by=c('sex.SOURCE', 'sex.RECIPIENT', 'date_infection_before_UTT.RECIPIENT', 'age_infection_reduced.RECIPIENT')]	
  tmp = dcast(tmp, sex.SOURCE + sex.RECIPIENT + date_infection_before_UTT.RECIPIENT + age_infection_reduced.RECIPIENT ~ q_label, value.var = "q")
  
  subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M' & age_infection_reduced.RECIPIENT == 46)
  
  # FM
  tmp1 <- subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_UTT.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_UTT.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection male recipient', y = 'Age at infection female source', col = 'Date infection of recipient before 2017') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 

  # MF
  tmp1 <- subset(tmp, sex.SOURCE == 'M' & sex.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_UTT.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_UTT.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection female recipient', y = 'Age at infection male source', col = 'Date infection of recipient before 2017') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 

  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = file.path(outdir, paste0('AgeInfectionSource_', lab, '.png')), w = 8, h = 5)
  
}

