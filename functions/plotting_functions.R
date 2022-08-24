plot_age_infection_source_recipient <- function(data, title, plotlab, outdir = NULL){
  
  plots = list()
  
  data <- data[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
  
  data[, `Community source` := COMM.SOURCE]
  data[, `Community recipient` := COMM.RECIPIENT]
  
  data[, COHORT_ROUND.SOURCE := substr(ROUND.SOURCE, 1, 4)]
  data[, COHORT_ROUND.RECIPIENT := substr(ROUND.RECIPIENT, 1, 4)]
  data[, `Cohort round recipient` := COHORT_ROUND.RECIPIENT]
  data[, `Cohort round source` := COHORT_ROUND.SOURCE]
  
  # all pairs
  p <- ggplot(data, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_point() + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) 
    if(!is.null(outdir))
      ggsave(p, filename = paste0(outdir, '-data-AgeInfection_AllPairs_', plotlab, '.png'), w = 4, h = 4)
  plots = c(plots, list(p))
  
  # by cohort round
  p1 <- ggplot(data, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round source`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')+
    guides(col=guide_legend(nrow=2,byrow=TRUE))

  p2 <- ggplot(data, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round recipient`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom') +
    guides(col=guide_legend(nrow=2,byrow=TRUE))

  p <- ggarrange(p1, p2, ncol = 2)
  if(!is.null(outdir))
    ggsave(p, filename = paste0(outdir, '-data-AgeInfection_CohortRound_', plotlab, '.png'), w = 9, h = 7)
  plots = c(plots, list(p))
  
  # by age infection round
  p <- ggplot(data, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_point(aes(col = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient',
         col = paste0('Date infection recipient before ', format(cutoff_date, '%Y'))) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')
  if(!is.null(outdir))
    ggsave(p, filename = paste0(outdir, '-data-AgeInfection_DateInfectionRecipient_', plotlab, '.png'), w = 4, h = 4)
  plots = c(plots, list(p))
  
  # by community
  p1 <- ggplot(data, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_point(aes(col = `Community source`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p2 <- ggplot(data, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_point(aes(col =`Community recipient`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p <- ggarrange(p1, p2, ncol = 2)
  if(!is.null(outdir))
    ggsave(p, filename = paste0(outdir, '-data-AgeInfection_CommunitySourceRecipient_', plotlab, '.png'), w = 9, h = 7)
  plots = c(plots, list(p))
  
  data[, date_infection_before_cutoff_name.RECIPIENT := 'After 2014']
  data[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == T, date_infection_before_cutoff_name.RECIPIENT := 'Before 2014']
  p2 <- ggplot(data, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_point(aes(col =`Community recipient`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$AGE_TRANSMISSION.SOURCE, data$AGE_INFECTION.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom') + 
    facet_grid(.~date_infection_before_cutoff_name.RECIPIENT)
  ggsave(p2, filename = paste0(outdir, '-data-AgeInfection_CommunityRecipient_', plotlab, '.png'), w = 5, h = 5)
  
  return(plots)
}

plot_hist_age_infection <- function(pairs, outdir = NULL){
  
  pairs[, Sex := SEX.SOURCE]
  p1 <- ggplot(pairs, aes(x = AGE_TRANSMISSION.SOURCE)) + 
    geom_histogram(bins = 30) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at transmission source') +
    scale_x_continuous(limits = range(c(pairs$AGE_TRANSMISSION.SOURCE, pairs$AGE_INFECTION.RECIPIENT)))
  
  pairs[, Sex := SEX.RECIPIENT]
  p2 <- ggplot(pairs, aes(x = AGE_INFECTION.RECIPIENT)) + 
    geom_histogram(bins = 30) +     
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$AGE_TRANSMISSION.SOURCE, pairs$AGE_INFECTION.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  
  if(!is.null(outdir)){
    file = paste0(outdir, '-data-hist_age_infection.png')
    cat('saving', file, '\n')
    ggsave(p, file = file, w = 6, h = 6)
  }
  
  return(p)
}

plot_hist_time_infection <- function(pairs, cutoff_date, outdir = NULL){
  
  pairs <- copy(pairs.all)
  
  # pairs[, COHORT_ROUND.SOURCE := substr(ROUND.SOURCE, 1, 4)]
  pairs[, COHORT_ROUND.RECIPIENT := substr(ROUND.RECIPIENT, 1, 4)]
  # pairs[, `Round source` := COHORT_ROUND.SOURCE]
  pairs[, Round_recipient := paste0('First round recipient:\n', COHORT_ROUND.RECIPIENT)]
  pairs <- pairs[COHORT_ROUND.RECIPIENT != 'neur']
 
  if(use.tsi.estimates){
    pairs[, type := 'With TSI estimates']
  }else{
    pairs[, type := 'Date infection = Date first visit - 1 year']
    pairs[!is.na(AGE_FIRST_POSITIVE.RECIPIENT) , type := 'Date infection = Date first positive - 1 year']
    pairs[AGE_FIRST_POSITIVE.RECIPIENT == AGE_FIRST_VISIT.RECIPIENT, type := 'Date infection = Date first visit(=first positive) - 1 year']
  }

  # inland
  tmp <- pairs[COMM.RECIPIENT == 'inland']
  tmp_round <- df_round[COMM == 'inland']
  p <- ggplot(tmp) +
    geom_rect(data = tmp_round, aes(ymin = -Inf, ymax = Inf, xmin = MIN_SAMPLE_DATE, 
                                    xmax = MAX_SAMPLE_DATE, fill = ROUND), alpha = 0.5) + 
    geom_histogram(aes(x = DATE_INFECTION.RECIPIENT), bins = 100) +
    facet_grid(Round_recipient~type) +
    theme_bw() + 
    labs(x = 'Date of infection recipient') + 
    geom_vline(xintercept = cutoff_date, linetype = 'dashed') + 
    scale_y_continuous(expand = expansion(mult = c(0, .05)))+ 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1))) +
    ggtitle('Inland communities')
  file = paste0(outdir, '-data-hist_date_infection_inland.png')
  ggsave(p, file = file, w = 10.5, h = 9.5)
  
  # fishing
  tmp <- pairs[COMM.RECIPIENT == 'fishing']
  tmp_round <- df_round[COMM == 'fishing']
  p <- ggplot(tmp) +
    geom_rect(data = tmp_round, aes(ymin = -Inf, ymax = Inf, xmin = MIN_SAMPLE_DATE, 
                                    xmax = MAX_SAMPLE_DATE, fill = ROUND), alpha = 0.5) + 
    geom_histogram(aes(x = DATE_INFECTION.RECIPIENT), bins = 100) +
    facet_grid(Round_recipient~type) +
    theme_bw() + 
    labs(x = 'Date of infection recipient', y = 'Count of phylo pairs') + 
    geom_vline(xintercept = cutoff_date, linetype = 'dashed') + 
    scale_y_continuous(expand = expansion(mult = c(0, .05)))+ 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)))  +
    ggtitle('Fishing communities')
  file = paste0(outdir, '-data-hist_date_infection_fishing.png')
  ggsave(p, file = file, w = 10.5, h = 9.5)

  return(p)
}

plot_CI_age_infection <- function(pairs, outdir = NULL){
  
  data <- copy(pairs)
  data[, AGE_TRANSMISSION.SOURCE := floor(AGE_TRANSMISSION.SOURCE)]
  data[, AGE_INFECTION.RECIPIENT := floor(AGE_INFECTION.RECIPIENT)]
  data <- merge(data, df_age, by = c('AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'))
  
  ps <- c(0.5, 0.2, 0.8)
  p_labs <- c('M','CL','CU')
  
  ## stratified by age of recipient
  tmp = data[, list(q= quantile(AGE_TRANSMISSION.SOURCE, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('SEX.SOURCE', 'SEX.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT', 'AGE_INFECTION_REDUCED.RECIPIENT')]	
  tmp = dcast(tmp, SEX.SOURCE + SEX.RECIPIENT + DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT + AGE_INFECTION_REDUCED.RECIPIENT ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = AGE_INFECTION_REDUCED.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection male recipient', y = 'Age at transmission female source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  # MF
  tmp1 <- subset(tmp, SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = AGE_INFECTION_REDUCED.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection female recipient', y = 'Age at transmission male source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = paste0(outdir, '-data-AgeTransmissionSource.png'), w = 8, h = 5)
  
  
  ## not stratified by age of recipient
  tmp = data[, list(q= quantile(AGE_TRANSMISSION.SOURCE, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('SEX.SOURCE', 'SEX.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT')]	
  tmp = dcast(tmp, SEX.SOURCE + SEX.RECIPIENT + DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT  ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Date infection recipient', y = 'Age at transmission female source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  # MF
  tmp1 <- subset(tmp, SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT)) + 
    geom_point(aes(y = M, col = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Date infection recipient', y = 'Age at transmission male source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = paste0(outdir, '-data-AgeTransmissionSource_Unstratified.png'), w = 8, h = 5)
  
}

plot_CI_age_transmission <- function(pairs, outdir = NULL){
  
  data <- copy(pairs)
  data[, AGE_TRANSMISSION.SOURCE := floor(AGE_TRANSMISSION.SOURCE)]
  data[, AGE_INFECTION.RECIPIENT := floor(AGE_INFECTION.RECIPIENT)]
  data <- merge(data, df_age, by = c('AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'))
  
  ps <- c(0.5, 0.2, 0.8)
  p_labs <- c('M','CL','CU')
  
  ## stratified by age of recipient
  tmp = data[, list(q= quantile(AGE_INFECTION.RECIPIENT, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('SEX.SOURCE', 'SEX.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT', 'AGE_TRANSMISSION_REDUCED.SOURCE')]	
  tmp = dcast(tmp, SEX.SOURCE + SEX.RECIPIENT + DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT + AGE_TRANSMISSION_REDUCED.SOURCE ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = AGE_TRANSMISSION_REDUCED.SOURCE)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at transmission female source', y = 'Age at infection male recipient', 
         col ='Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  # MF
  tmp1 <- subset(tmp, SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = AGE_TRANSMISSION_REDUCED.SOURCE)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at transmission male source', y = 'Age at infection female recipient', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = paste0(outdir, '-data-AgeInfectionRecipient.png'), w = 8, h = 5)
  
  
  ## not stratified by age of recipient
  tmp = data[, list(q= quantile(AGE_INFECTION.RECIPIENT, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('SEX.SOURCE', 'SEX.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT')]	
  tmp = dcast(tmp, SEX.SOURCE + SEX.RECIPIENT + DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT  ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = date_infection_before_cutoff_name.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Date infection recipient', y = 'Age at infection male recipient', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  # MF
  tmp1 <- subset(tmp, SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = date_infection_before_cutoff_name.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Date infection recipient', y = 'Age at infection female recipient', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = paste0(outdir, '-data-AgeInfectionRecipient_Unstratified.png'), w = 8, h = 5)
}

plot_pairs_infection_dates <- function(pairs.all, outdir)
{
  
  tmp <- copy(pairs.all)
  tmp[, DUMMY:= DATE_INFECTION.SOURCE - DATE_INFECTION.RECIPIENT]
  setkey(tmp, date_first_positive.SOURCE)
  tmp[, CONS:=ifelse(DUMMY < 0, 'consistent', 'inconsistent')]
  setkey(tmp, DUMMY)
  # tmp[, DUMMY:=1:.N]
  tmp[date_first_positive.SOURCE < date_first_positive.RECIPIENT ] # (316 pairs for which )
  tmp[DATE_INFECTION.SOURCE < DATE_INFECTION.RECIPIENT  ]
  
  p <- ggplot(tmp, aes(y=date_first_positive.SOURCE, col=CONS)) + 
    geom_point(aes(x=DATE_INFECTION.SOURCE)) + 
    geom_errorbarh(aes(xmin=DATE_INFECTION.SOURCE, xmax=DATE_INFECTION.RECIPIENT)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(data=tmp[date_first_positive.SOURCE > date_first_positive.RECIPIENT],aes(x=date_first_positive.RECIPIENT), color='black', pch=18) +
    scale_x_date(date_labels = '%Y', breaks = '12 months') +
    scale_y_date(date_labels = '%Y', breaks = '12 months') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position = 'bottom') +
    labs(x='Estimated infection dates', y=" Source's date of first positive test", col='source-date relationship')
  
  ggsave(p, filename = paste0(outdir, paste0('-data-PairsInfectionDates.png')), w = 8, h = 10)
}

# Tranmission Network plot from: https://github.com/olli0601/Phyloscanner.R.utilities/blob/7c58edac4812e53b0b67e69770f70e9d3826881d/R/phyloscan.fun.plotting.R
phsc.plot.transmission.network<- function(dchain, dc, pairs, outdir=NULL, point.size=10, point.size.couple=point.size*1.4, edge.gap=0.04, edge.size=0.4, curvature= -0.2, arrow=arrow(length=unit(0.02, "npc"), type="open"), curv.shift=0.08, label.size=3, node.label='ID', node.shape=NA_character_, node.fill=NA_character_, node.shape.values=c('M' = 15, 'F' = 17), node.fill.values=c('F'='hotpink2', 'M'='steelblue2') , threshold.linked=NA_real_, N_min = 4)
{	
  library(phyloscannerR)
  # translate dchain and dc in terms of df and di as required
  tmp.dchain <- copy(dchain)
  tmp.dc <- copy(dc)
  tmp.dchain[LINK_21 == 1, `:=` (SOURCE=H2, RECIPIENT=H1)]
  tmp.dchain[LINK_12 == 1, `:=` (SOURCE=H1, RECIPIENT=H2)]

  # Only include the pairs we want
  tmp.dchain <- merge(pairs[,.(RECIPIENT,SOURCE, SEX.RECIPIENT)], tmp.dchain, by=c('SOURCE', 'RECIPIENT'))
  # hist(tmp.dchain[, SCORE_LINKED]);hist(tmp.dchain[, SCORE_DIR_12])

  # make it conformable with the plot function arguments (forget about sex etc... atm)
  tmp.dchain[, TYPE := ifelse(LINK_12 ==1, '12', '21')]
  tmp.dchain[,.(IDCLU, H1, H2, TYPE, CATEGORISATION, SCORE_LINKED)]

  # get K_EFF from dc 
  tmp <- tmp.dc[CATEGORISATION == 'close.and.adjacent.cat' & !grepl('not',TYPE),]
  tmp[, `:=` (TYPE=NULL, CNTRL1=NULL, CNTRL2=NULL, CATEGORISATION=NULL, CATEGORICAL_DISTANCE=NULL, SCORE=NULL, PTY_RUN=NULL)]
  tmp.dchain <- merge(tmp.dchain, tmp, by=c('H1', 'H2'))

  setnames(tmp.dchain, c('H1', 'H2', 'K_EFF', 'SCORE_LINKED'), c('ID1', 'ID2', 'KEFF', 'POSTERIOR_SCORE'))
  df <- unique(tmp.dchain[, .(IDCLU, ID1, ID2, TYPE, KEFF, POSTERIOR_SCORE)])

  di <- data.table(ID = df[, unique( c(ID1,ID2))])
  di <- merge(di, pairs[, .(RECIPIENT, SEX.RECIPIENT)], by.x='ID' , by.y='RECIPIENT', all.x=TRUE)
  di <- unique(merge(di, pairs[, .(SOURCE, SEX.SOURCE)], by.x='ID' , by.y='SOURCE', all.x=TRUE))
  di[, NODE_FILL:=na.omit(c(SEX.RECIPIENT, SEX.SOURCE))[1], by=ID]
  di[, `:=`(SEX.RECIPIENT=NULL, SEX.SOURCE=NULL)]
  di[, `:=` (NODE_LABEL=ID, NODE_SHAPE=NODE_FILL)]

  # plot transmission Network:
  tmp <- df[, .N, by = IDCLU][N >= N_min,IDCLU]
  df <- df[IDCLU %in% tmp]
  di <- di[ID %in% df[, unique(c(ID1,ID2))] ]

  # Function from Phyloscanner
  if(is.na(node.label))
  {
    node.label<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
    set(di, NULL, node.label, NA_character_)
  }
  if(is.na(node.shape))
  {
    node.shape<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
    set(di, NULL, node.shape, 'NA')
  }
  if(is.na(node.fill))
  {
    node.fill<- paste0('DUMMY',1+length(which(grepl('DUMMY',colnames(di)))))
    set(di, NULL, node.fill, 'NA')
  }
  if(any(is.na(node.fill.values)))
  {
    z						<- unique(di[[node.fill]])
    node.fill.values		<- heat.colors(length(z))
    names(node.fill.values)	<- z
  }
  if(any(is.na(node.shape.values)))
  {
    z						<- unique(di[[node.shape]])
    node.shape.values		<- seq_along(z)
    names(node.shape.values)<- z
  }
  setnames(di, c(node.label, node.shape, node.fill), c('NODE_LABEL','NODE_SHAPE','NODE_FILL'))
  tmp	<- c('NODE_LABEL','NODE_SHAPE','NODE_FILL')[which(c(node.label, node.shape, node.fill)=='ID')]
  if(length(tmp))
    set(di, NULL, 'ID', di[[tmp]])
  di	<- subset(di, select=c(ID, NODE_LABEL, NODE_SHAPE, NODE_FILL))
  
  layout	<- as.data.table(ggnet2(network(unique(subset(df, select=c(ID1,ID2))), directed=FALSE, matrix.type="edgelist"))$data[,c("label", "x", "y")])
  setnames(layout, c('label','x','y'), c('ID1','ID1_X','ID1_Y'))
  df		<- merge(df, layout, by='ID1')
  setnames(layout, c('ID1','ID1_X','ID1_Y'), c('ID2','ID2_X','ID2_Y'))
  df		<- merge(df, layout, by='ID2')
  setnames(layout, c('ID2','ID2_X','ID2_Y'),  c('ID','X','Y'))	
  layout	<- merge(layout,di, by='ID')	
  
  df[, EDGETEXT_X:= (ID1_X+ID2_X)/2]
  df[, EDGETEXT_Y:= (ID1_Y+ID2_Y)/2]
  #
  #	calculate score for linked
  if(is.na(threshold.linked))
  {
    df	<- merge(df,df[, 	{
      z<- rep('edge_col_1', length(TYPE))
      z[which.max(POSTERIOR_SCORE)]	<- 'edge_col_2'
      list(EDGE_COL=z, TYPE=TYPE)	
    }, by=c('ID1','ID2')], by=c('ID1','ID2','TYPE'))		
  }
  if(!is.na(threshold.linked))
  {
    tmp	<- subset(df, TYPE!='not close/disconnected')[, list( EDGE_COL=as.character(factor(sum(POSTERIOR_SCORE)>=threshold.linked, levels=c(TRUE, FALSE), labels=c('edge_col_2','edge_col_1'))) ), by=c('ID1','ID2')]
    df	<- merge(df, tmp, by=c('ID1','ID2'))		
  }	
  #	for edges, move the start and end points on the line between X and Y
  #	define unit gradient
  df[, MX:= (ID2_X - ID1_X)]	
  df[, MY:= (ID2_Y - ID1_Y)]
  tmp		<- df[, sqrt(MX*MX+MY*MY)]
  set(df, NULL, 'MX', df[, MX/tmp])
  set(df, NULL, 'MY', df[, MY/tmp])	
  set(df, NULL, 'ID1_X', df[, ID1_X + MX*edge.gap])
  set(df, NULL, 'ID1_Y', df[, ID1_Y + MY*edge.gap])
  set(df, NULL, 'ID2_X', df[, ID2_X - MX*edge.gap])
  set(df, NULL, 'ID2_Y', df[, ID2_Y - MY*edge.gap])	
  #	label could just be move on the tangent vector to the line
  #	define unit tangent
  df[, TX:= -MY]
  df[, TY:= MX]
  tmp		<- df[, which(TYPE=='12')]
  set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X + TX*curv.shift])
  set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y + TY*curv.shift])
  tmp		<- df[, which(TYPE=='21')]
  set(df, tmp, 'EDGETEXT_X', df[tmp, EDGETEXT_X - TX*curv.shift])
  set(df, tmp, 'EDGETEXT_Y', df[tmp, EDGETEXT_Y - TY*curv.shift])
  #	
  
  p		<- ggplot() +			
    geom_point(data=layout, aes(x=X, y=Y, colour=NODE_FILL, pch=NODE_SHAPE), size=point.size) +
    geom_segment(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
    geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
    geom_curve(data=subset(df, EDGE_COL=='edge_col_1' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
    geom_segment(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='ambiguous' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), lineend="butt") +
    geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='12' & KEFF>0), aes(x=ID1_X, xend=ID2_X, y=ID1_Y, yend=ID2_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +
    geom_curve(data=subset(df, EDGE_COL=='edge_col_2' & TYPE=='21' & KEFF>0), aes(x=ID2_X, xend=ID1_X, y=ID2_Y, yend=ID1_Y, size=edge.size*KEFF, colour=EDGE_COL), curvature=curvature, arrow=arrow, lineend="butt") +									
    scale_colour_manual(values=c(node.fill.values, 'edge_col_1'='red', 'edge_col_2'='blue','NA'='green')) +
    scale_shape_manual(values=c(node.shape.values, 'NA'=21)) +
    scale_fill_manual(values=c(node.fill.values, 'NA'='grey50')) +
    scale_size_identity() +
    geom_text(data=subset(df, TYPE!='not close/disconnected' & KEFF>0), aes(x=EDGETEXT_X, y=EDGETEXT_Y, label=paste0(round(100*POSTERIOR_SCORE,d=1),'%')), size=label.size) +
    geom_text(data=layout, aes(x=X, y=Y, label=NODE_LABEL)) +
    theme_void() +
    guides(colour='none', fill='none',size='none', pch='none') 
  layout		<- subset(layout, select=c(ID,X,Y))
  setnames(layout, c('ID','X','Y'), c('label','x','y'))	
  # p$layout	<- layout
  
  if(!is.null(outdir))
    ggsave(p, filename = paste0(outdir, '-TransmissionNetworkClusters.png'), w = 10, h = 10)

  return(p)
  
}


plot_data_by_round <- function(eligible_count_round, proportion_unsuppressed, proportion_prevalence, incidence_cases_round, outdir){
  
  level_rounds <- c('R014', 'R015', 'R016', 'R017', 'R018')
  
  # round periods
  tmp <- rbind(df_round_inland, df_round_fishing)
  tmp1 <- copy(df_period)
  tmp1[, INDEX_TIME2 := paste0('Period: ', INDEX_TIME)]
  ggplot(tmp) + 
    geom_rect(data = tmp1, aes(ymin = -Inf, ymax = Inf, xmin = MIN_PERIOD_DATE, xmax = MAX_PERIOD_DATE, fill = INDEX_TIME2), alpha = 0.5) + 
    geom_errorbarh(aes(xmin =min_sample_date , xmax =  max_sample_date, col = as.factor(round), y = as.factor(round))) + 
    facet_grid(COMM~.) + 
    theme_bw() + 
    labs(y = 'Round', col = '', fill = '') 
  ggsave(paste0(outdir, '-data-period-round.png'), w = 7, h = 8)
  
  # Census eligible count 
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = ELIGIBLE, col = SEX)) +
    labs(y = 'Census eligible count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_count_round_sex.png'), w = 7, h = 8)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = ELIGIBLE, col = ROUND)) +
    labs(y = 'Census eligible count', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_count_round.png'), w = 7, h = 6)
  
  # Prevalence proportion
  ggplot(subset(proportion_prevalence, ROUND !='15S'), aes(x = AGEYRS)) +
    geom_point(aes(y =EMPIRICAL_PREVALENCE, col = SEX)) +
    geom_line(aes(y =PREVALENCE_M, col = SEX)) +
    geom_ribbon(aes(ymin =PREVALENCE_CL, ymax = PREVALENCE_CU, fill = SEX), alpha = 0.5) +
    labs(y = 'Prevlence among participant', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_infected_round_sex.png'), w = 7, h = 8)
  
  tmp <- as.data.table(proportion_prevalence)
  tmp[, ROUND := factor( ROUND, levels = level_rounds)]
  ggplot(subset(tmp, ROUND !='15S'), aes(x = AGEYRS)) +
    geom_point(aes(y =EMPIRICAL_PREVALENCE, col = ROUND), alpha = 0.5) +
    geom_line(aes(y =PREVALENCE_M, col = ROUND)) +
    geom_ribbon(aes(ymin =PREVALENCE_CL, ymax = PREVALENCE_CU, fill = ROUND), alpha = 0.1) +
    labs(y = 'Prevlence among participant', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')+ 
    scale_color_discrete(drop = FALSE) + 
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent)
  ggsave(paste0(outdir, '-data-census_eligible_prop_infected_round.png'), w = 7, h = 6)
  
  ggplot(subset(tmp, ROUND !='R015S'), aes(x = AGEYRS)) +
    geom_point(aes(y =1-EMPIRICAL_PREVALENCE, col = ROUND), alpha = 0.5) +
    geom_line(aes(y =1-PREVALENCE_M, col = ROUND)) +
    geom_ribbon(aes(ymin =1-PREVALENCE_CL, ymax = 1-PREVALENCE_CU, fill = ROUND), alpha = 0.1) +
    labs(y = 'Proportion of susceptibile', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')+ 
    scale_color_discrete(drop = FALSE) + 
    scale_fill_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent)
  ggsave(paste0(outdir, '-data-census_eligible_prop_susceptible_round.png'), w = 7, h = 6)
  
  tmp[, PROP_SUSCEPTIBLE:= 1-PREVALENCE_M]
  tmp[, TOTAL_PROP_SUSCEPTIBLE := sum(PROP_SUSCEPTIBLE), by = c('ROUND', 'SEX', 'COMM')]
  tmp[, distr := PROP_SUSCEPTIBLE / TOTAL_PROP_SUSCEPTIBLE]
  ggplot(tmp, aes(x = AGEYRS, y = distr)) + 
    geom_line(aes(col = ROUND)) + 
    facet_grid(SEX~COMM) + 
    theme_bw() + 
    labs(x = 'Age', y = 'Probability distribution of the age composition of proportion of susceptible')
  ggsave(paste0(outdir, '-data-distribution_function_age_composition_susceptible.png'), w = 7, h = 7)
  
  
  # HIV+ census eligible count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED , col = SEX)) +
    labs(y = 'Census eligible infected count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_infected_round_sex.png'), w = 7, h = 8)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED , col = ROUND)) +
    labs(y = 'Census eligible infected count', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_infected_round.png'), w = 7, h = 6)
  
  # Susceptible census eligible count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = SUSCEPTIBLE , col = SEX)) +
    labs(y = 'Census eligible susceptible count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_susceptible_round_sex.png'), w = 7, h = 8)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = SUSCEPTIBLE , col = ROUND)) +
    labs(y = 'Census eligible susceptible count', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_susceptible_round.png'), w = 7, h = 6)
  
  # proportion of unsuppressed
  ggplot(proportion_unsuppressed, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_UNSUPPRESSED_M , col = SEX)) +
    geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL , ymax = PROP_UNSUPPRESSED_CU , fill = SEX), alpha = 0.5) +
    geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL , col = SEX), alpha = 0.5) +
    labs(y = 'Proportion of unsupressed', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_unsuppressed_round_sex.png'), w = 7, h = 8)
  
  ggplot(proportion_unsuppressed, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_UNSUPPRESSED_M , col = ROUND)) +
    geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL , ymax = PROP_UNSUPPRESSED_CU , fill = ROUND), alpha = 0.5) +
    geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL , col = ROUND), alpha = 0.5) +
    labs(y = 'Proportion of unsupressed', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_unsuppressed_round.png'), w = 7, h = 6)
  
  #HIV+ unsupressed census eligible count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED_NON_SUPPRESSED , col = SEX)) +
    # geom_ribbon(aes(ymin = INFECTED_NON_SUPPRESSED_CL, ymax = INFECTED_NON_SUPPRESSED_CU , fill = SEX), alpha = 0.5) +
    labs(y = 'Number of HIV+ unsupressed census eligible', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_unsuppressed_round_sex.png'), w = 7, h = 8)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED_NON_SUPPRESSED , col = ROUND)) +
    # geom_ribbon(aes(ymin = INFECTED_NON_SUPPRESSED_CL, ymax = INFECTED_NON_SUPPRESSED_CU , fill = ROUND), alpha = 0.5) +
    labs(y = 'Number of HIV+ unsupressed census eligible', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_unsuppressed_round.png'), w = 7, h = 6)
  
  # incidence rate per person per round
  ggplot(incidence_cases_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENCE, col = ROUND)) +
    # geom_ribbon(aes(ymin = LB * ROUND_SPANYRS, ymax = UB* ROUND_SPANYRS, fill = SEX),  alpha = 0.5) +
    labs(y = 'Incidence rate per 1 person per year', x = 'Age') +
    facet_grid(COMM~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-incidence_rate_round.png'), w = 7, h = 6)
  
  # incidence cases
  ggplot(incidence_cases_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENT_CASES / ROUND_SPANYRS, col = ROUND)) +
    # geom_ribbon(aes(ymin = INCIDENT_CASES_LB , ymax = INCIDENT_CASES_UB , fill = SEX), alpha = 0.5) +
    labs(y = 'Number of incident cases per year', x = 'Age') +
    facet_grid(COMM~SEX, label = 'label_both', scale = 'free_y') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-incidence_case_round.png'), w = 7, h = 6)
  

  # empirical contribution and transmisison risk 
  tmp <- copy(incidence_cases_round)
  tmp <- tmp[, list(INCIDENT_CASES = sum(INCIDENT_CASES), 
                    INFECTED_NON_SUPPRESSED = sum(INFECTED_NON_SUPPRESSED)), by = c('COMM', 'SEX', 'ROUND', 'ROUND_SPANYRS')]
  tmp[, INFECTED_NON_SUPPRESSED_OPPOSITE_SEX := ifelse(SEX == 'F', INFECTED_NON_SUPPRESSED[SEX == 'M'], 
                                                       INFECTED_NON_SUPPRESSED[SEX == 'F']), by = c('COMM', 'ROUND', 'ROUND_SPANYRS')]
  tmp[, ADJ_INCIDENT_CASES := INCIDENT_CASES / INFECTED_NON_SUPPRESSED_OPPOSITE_SEX]
  tmp[,TOTAL_ADJ_CASES := sum(ADJ_INCIDENT_CASES), by = c('COMM', 'ROUND')]
  tmp[, PROP_ADJ_CASES := ADJ_INCIDENT_CASES / TOTAL_ADJ_CASES]
  tmp[,TOTAL_CASES := sum(INCIDENT_CASES), by = c('COMM', 'ROUND')]
  tmp[, PROP_CASES := INCIDENT_CASES / TOTAL_CASES]
  
  tmp <- melt.data.table(tmp, id.vars = c("COMM", 'SEX', 'ROUND', 'ROUND_SPANYRS'))
  tmp1 <- tmp[variable %in% c('PROP_ADJ_CASES', 'PROP_CASES')]
  tmp1[, type := 'not adjusted by the number of HIV+ unsuppressed']
  tmp1[grepl("ADJ", variable), type := 'adjusted by the number of HIV+ unsuppressed']
  
  tmp1[, SEX_SOURCE := 'Male']
  tmp1[SEX == 'M', SEX_SOURCE := 'Female']
  
  ggplot(tmp1, aes(x = ROUND, y = value, col = type)) + 
    geom_point() + 
    facet_grid(COMM~SEX_SOURCE) + 
    theme_bw() + 
    labs(y = 'Empirical contribution to infection', col = '')  +
    theme(legend.position = 'bottom') + 
    guides(col = guide_legend(byrow = T, nrow = 2)) + 
    scale_y_continuous(labels = scales::percent_format())
  ggsave(paste0(outdir, '-data-empirical_contribution_infection.png'), w = 7, h = 6)
  
  # empirical transmission risk
  tmp1 <- tmp[variable == 'ADJ_INCIDENT_CASES']
  tmp1[, SEX_SOURCE := 'Male']
  tmp1[SEX == 'M', SEX_SOURCE := 'Female']
  
  ggplot(tmp1, aes(x = ROUND, y = value / ROUND_SPANYRS, fill = SEX_SOURCE)) + 
    geom_bar(stat = 'identity') + 
    facet_grid(COMM~SEX_SOURCE, scale = 'free') + 
    theme_bw() + 
    scale_fill_manual(values = c('Male'='royalblue3','Female'='deeppink2')) + 
    labs(y = 'Empirical transmission risk per year', col = '')  +
    theme(legend.position = 'none') + 
    guides(col = guide_legend(byrow = T, nrow = 2))  + 
    scale_y_continuous(expand = expansion(mult = c(0, .05)))
  ggsave(paste0(outdir, '-data-empirical_transmission_risk_sex.png'), w = 7, h = 6)
  
}

plot_data_by_period <- function(incidence_cases, outdir){
  
  # incidence cases
  ggplot(incidence_cases, aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENT_CASES , col = PERIOD)) +
    # geom_ribbon(aes(ymin = INCIDENT_CASES_LB , ymax = INCIDENT_CASES_UB , fill = SEX), alpha = 0.5) +
    labs(y = 'Number of incident cases', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-incidence_case_period.png'), w = 7, h = 6)
}

plot_offset <- function(stan_data, outdir){
  
  tmp <- as.data.table(reshape2::melt(stan_data[['log_offset']]))
  setnames(tmp, 1:4, c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'INDEX_AGE'))
  tmp <- merge(tmp, df_direction, by = 'INDEX_DIRECTION')
  tmp <- merge(tmp, df_community, by = 'INDEX_COMMUNITY')
  tmp <- merge(tmp, df_round, by = c('INDEX_ROUND', 'COMM'))
  tmp <- merge(tmp, df_age, by = 'INDEX_AGE')

  tmp1 <- tmp[, list(value = sum(exp(value))), by = c('AGE_INFECTION.RECIPIENT', 'LABEL_DIRECTION', 'LABEL_COMMUNITY', 'ROUND')]
  
  ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_raster(aes(fill = exp(value))) +
    facet_grid( ROUND~LABEL_DIRECTION + LABEL_COMMUNITY) + 
    theme_bw() +
    scale_y_continuous(expand= c(0,0))+
    scale_x_continuous(expand= c(0,0)) +
    scale_fill_viridis_c()  + 
    theme(strip.background = element_rect(colour="black", fill="white"),
          strip.text = element_text(size = rel(1))) +
    labs(x= 'Age at infection recipient', y = 'Age at transmission source', fill = 'offset') 
  ggsave(paste0(outdir, '-offset-value.png'), w = 10, h = 12 )
}

plot_crude_force_infection <- function(crude_force_infection, outdir){
  
  communities <- crude_force_infection[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- crude_force_infection[ COMM == communities[i]]

    p <- ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = CRUDE_FOI)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = 'Crude force of\ninfection', 
           y= 'Age at transmission source') +
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) + 
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    
    ggsave(p, file = paste0(outdir, '-data-crude_FOI_',  communities[i], '.png'), w = 7, h = 7)
  }
  
  # by age recipient
  tmp <- crude_force_infection[, list(CRUDE_FOI = sum(CRUDE_FOI)), by= c('AGE_INFECTION.RECIPIENT', 'PERIOD', 'LABEL_DIRECTION', 'COMM')]
  ggplot(tmp, aes(x = AGE_INFECTION.RECIPIENT)) +
    geom_line(aes(y = CRUDE_FOI, col = PERIOD)) +
    labs(x = 'Age of recipient', y = 'Crude force of infection received', fill = '') +
    theme_bw() +
    facet_grid(COMM~LABEL_DIRECTION) +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg()
  ggsave(paste0(outdir, '-data-crude_FOI_age_recipient.png'), w = 7, h = 6)
  
  # by age source
  tmp <- crude_force_infection[, list(CRUDE_FOI = sum(CRUDE_FOI)), by= c('AGE_TRANSMISSION.SOURCE', 'PERIOD', 'LABEL_DIRECTION', 'COMM')]
  ggplot(tmp, aes(x = AGE_TRANSMISSION.SOURCE)) +
    geom_line(aes(y = CRUDE_FOI, col = PERIOD)) +
    labs(x = 'Age of recipient', y = 'Crude force of infection exerted', fill = '') +
    theme_bw() +
    facet_grid(COMM~LABEL_DIRECTION) +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg()
  ggsave(paste0(outdir, '-data-crude_FOI_age_source.png'), w = 7, h = 6)
  
}


plot_transmission_events_over_time <- function(eligible_count_round, incidence_cases_round, pairs, outdir){
  # timeline
  df_timeline <- copy(df_round)
  df_timeline[, MIDPOINT := as.Date(mean(c(MIN_SAMPLE_DATE, MAX_SAMPLE_DATE))), by = c('ROUND', 'COMM')]
  df_timeline <- df_timeline[, .(ROUND, MIDPOINT, COMM, INDEX_ROUND)]
  
  # age groups
  age_groups <- c('15-24', '25-34', '35-49')
  df_age_group <- data.table(AGEYRS = 15:49)
  df_age_group[, index_age_group := 3]
  df_age_group[AGEYRS < 35, index_age_group := 2]
  df_age_group[AGEYRS < 25, index_age_group := 1]
  df_age_group[, age_group := age_groups[index_age_group]]
  df_age_group[, AGE_GROUP_LABEL := paste0('Age: ', age_group)]
  
  # sex
  df_sex <- data.table(SEX = c('M', 'F'), SEX_LABEL = c('Male', 'Female'))
  
  # grid
  df_grid <- data.table(expand.grid(SEX_LABEL = df_sex[, unique(SEX)], 
                                    age_group = df_age_group[, unique(age_group)], 
                                    INDEX_ROUND = df_timeline[, unique(INDEX_ROUND)], 
                                    COMM=c('inland', 'fishing')))
  df_grid <- merge(df_grid, unique(df_age_group[, .(age_group, AGE_GROUP_LABEL)]), by = 'age_group')
  df_grid <- merge(df_grid, (df_timeline), by = c('INDEX_ROUND', 'COMM'))
  
  #
  # Prepare census eligible
  
  # sum across age group
  ecr <- merge(eligible_count_round, df_age_group, by = 'AGEYRS')
  ecr <- ecr[, list(ELIGIBLE = sum(ELIGIBLE)), by = c('SEX', 'ROUND', 'COMM', 'AGE_GROUP_LABEL')] 
  
  # merge labels
  ecr <- merge(ecr, df_timeline, by = c('ROUND', 'COMM'))
  ecr <- merge(ecr, df_sex, by = 'SEX')
  
  #
  # Prepare incidence cases
  icr <- merge(incidence_cases_round, df_age_group, by = 'AGEYRS')
  icr <- icr[, list(INCIDENT_CASES = sum(INCIDENT_CASES), 
                    INCIDENT_CASES_UB = sum(INCIDENT_CASES_UB),
                    INCIDENT_CASES_LB = sum(INCIDENT_CASES_LB)), by = c('SEX', 'ROUND', 'COMM', 'AGE_GROUP_LABEL')] 
  
  # merge labels
  icr <- merge(icr, df_timeline, by = c('ROUND', 'COMM'))
  icr <- merge(icr, df_sex, by = 'SEX')
  
  #
  # Prepare phylo pairs
  dp <- copy(pairs)
  setnames(dp, c('SEX.RECIPIENT', 'COMM.RECIPIENT', 'AGE_INFECTION.RECIPIENT', 'DATE_INFECTION.RECIPIENT'), 
           c('SEX', 'COMM', 'AGEYRS', 'DATE'))
  dp[, AGEYRS := floor(AGEYRS)]
  dp <- merge(dp, df_age_group, by = 'AGEYRS')
  dp[, DIRECTION := 'Male -> Female']
  dp[SEX.SOURCE == 'F', DIRECTION := 'Female -> Male']
  # aggregated  age groups
  # dp <- dp[, list(count = .N), by =  c('SEX', 'DATE', 'AGE_GROUP_LABEL', 'COMM')]
  
  
  
  #
  # Plot
  
  communities <- hivinc[, unique(COMM)]
  male_color <- 'lightblue3'
  female_color <- 'lightpink2'
  
  i = 1
  comm <- communities[i]
  df_round_comm <- df_round[COMM == comm]
  df_round_comm <- full_join(df_round_comm, unique(df_age_group[, .(AGE_GROUP_LABEL)]), by = character())
  
  # plot person year
  p1 <- ggplot(ecr[COMM == comm]) + 
    geom_bar(aes(x = MIDPOINT, y = ELIGIBLE, fill = SEX_LABEL), stat = 'identity', position = position_dodge()) + 
    facet_grid(.~AGE_GROUP_LABEL) + 
    labs(y = 'Census eligible count', x= 'Date (midpoint of survey interval)') + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          strip.background = element_rect(colour="white", fill="white"),
          # axis.text.x = element_text(angle= 70, hjust = 1),
          strip.text = element_text(size = rel(1.2)), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = c(0.93, 0.85), 
          legend.title = element_blank()) + 
    scale_fill_manual(values = c('Male'=male_color,'Female'=female_color)) +
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(limits = c(df_period[, min(MIN_PERIOD_DATE)], df_period[, max(MAX_PERIOD_DATE)]), expand = c(0,0)) 
  
  # plot incidence cases
  p2 <- ggplot(icr[COMM == comm]) + 
    geom_bar(aes(x = MIDPOINT, y = INCIDENT_CASES, fill = SEX_LABEL), stat = 'identity', position = position_dodge(width = 550)) +
    geom_errorbar(aes(x = MIDPOINT, ymin = INCIDENT_CASES_LB, ymax = INCIDENT_CASES_UB, group = SEX_LABEL), position = position_dodge(width = 550), width = 200, col = 'grey50') +
    facet_grid(.~AGE_GROUP_LABEL) + 
    labs(y = 'HIV incident cases\namong census eligible population', x= 'Date (midpoint of survey interval)') + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          strip.background = element_rect(colour="white", fill="white"),
          # axis.text.x = element_text(angle= 70, hjust = 1),
          strip.text = element_blank(), 
          legend.position = c(0.93, 0.85), 
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.title = element_blank()) + 
    scale_fill_manual(values = c('Male'=male_color,'Female'=female_color)) +
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(limits = c(df_period[, min(MIN_PERIOD_DATE)], df_period[, max(MAX_PERIOD_DATE)]), expand = c(0,0))
  
  # plot pairs 
  p3 <- ggplot(dp[COMM == comm]) + 
    geom_histogram(aes(x = DATE, fill = DIRECTION), bins = 30) + 
    facet_grid(.~AGE_GROUP_LABEL) + 
    labs(y = '\nDetected HIV-1 transmission events', x = 'Date at transmission') + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          # axis.text.x = element_text(angle= 70, hjust = 1),
          strip.text =  element_blank(), 
          legend.position = c(0.90, 0.85), 
          legend.title = element_blank()) + 
    scale_fill_manual(values = c('Male -> Female'=female_color,'Female -> Male'=male_color)) +
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, 0.05))) + 
    scale_x_date(limits = c(df_period[, min(MIN_PERIOD_DATE)], df_period[, max(MAX_PERIOD_DATE)]), expand = c(0,0))  
  
  # arrange
  p <- grid.arrange(p1, p2, p3, layout_matrix = rbind(c(NA,1,1), 
                                                      c(2, 2, 2), 
                                                      c(NA,NA,3)), 
                                                      widths = c(0.01, 0.008, 0.98), 
                                                      heights = c(0.35, 0.32,0.32))
  ggsave(p, file =  paste0(outdir, '-data-panel.png'), w = 8, h = 10)

}







