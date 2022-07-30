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
  
  pairs[, COHORT_ROUND.SOURCE := substr(ROUND.SOURCE, 1, 4)]
  pairs[, COHORT_ROUND.RECIPIENT := substr(ROUND.RECIPIENT, 1, 4)]
  pairs[, `Round source` := COHORT_ROUND.SOURCE]
  pairs[, `Round recipient` := COHORT_ROUND.RECIPIENT]
  
  p1 <- ggplot(pairs, aes(x = DATE_INFECTION.SOURCE)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round source`, nrow = length(unique(pairs$COHORT_ROUND.SOURCE))) +
    theme_bw() + 
    labs(x = 'Date of infection source') + 
    geom_vline(xintercept = cutoff_date, linetype = 'dashed')
  
  p2 <- ggplot(pairs, aes(x = DATE_INFECTION.RECIPIENT)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round recipient`, nrow = length(unique(pairs$COHORT_ROUND.RECIPIENT))) +
    theme_bw() + 
    labs(x = 'Date of infection recipient') + 
    geom_vline(xintercept = cutoff_date, linetype = 'dashed')

  p <- ggarrange(p1, p2, ncol = 2)
  
  if(!is.null(outdir)){
    file = paste0(outdir, '-data-hist_date_infection.png')
    cat('saving', file, '\n')
    ggsave(p, file = file, w = 6, h = 6)
  }

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


plot_data_by_round <- function(eligible_susceptible_count, proportion_unsuppressed, eligible_count_round, incidence_cases_round, outdir){
  
  # Census eligible count 
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = ELIGIBLE, col = SEX)) +
    labs(y = 'Census eligible count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_count_round.png'), w = 7, h = 8)
  
  # Prevalence proportion
  ggplot(eligible_susceptible_count, aes(x = AGEYRS)) +
    geom_point(aes(y =PREVALENCE_EMPIRICAL, col = SEX)) +
    geom_line(aes(y =PREVALENCE_PROPORTION, col = SEX)) +
    geom_ribbon(aes(ymin =PREVALENCE_PROPORTION_CL, ymax = PREVALENCE_PROPORTION_CU, fill = SEX), alpha = 0.5) +
    labs(y = 'Prevlence among participant', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-prevalence_proportion_round.png'), w = 7, h = 8)
  
  # HIV+ census eligible count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED , col = SEX)) +
    labs(y = 'Census eligible infected count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_infected_round.png'), w = 7, h = 8)
  
  # Susceptible census eligible count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = SUSCEPTIBLE , col = SEX)) +
    labs(y = 'Census eligible susceptible count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_susceptible_round.png'), w = 7, h = 8)
  
  eligible_count_round[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_SUSCEPTIBLE , col = ROUND)) +
    labs(y = 'Proportion of susceptible among census eligible individuals', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_susceptible_round.png'), w = 7, h = 8)
  
  # proportion of unsuppressed
  tmp1 <- as.data.table(proportion_unsuppressed)[ROUND != 'R015' & !(COMM == 'fishing' & ROUND == 'R016')]
  ggplot(tmp1, aes(x = AGEYRS)) +
    geom_line(aes(y = M , col = SEX)) +
    geom_ribbon(aes(ymin = CL , ymax = CU , fill = SEX), alpha = 0.5) +
    geom_point(aes(y = PROP_NON_SUPPRESSED_EMPIRICAL , col = SEX), alpha = 0.5) +
    labs(y = 'Proportion of unsupressed', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-unsuppressed_proportion_round.png'), w = 7, h = 8)
  
  #HIV+ unsupressed census eligible count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED_NON_SUPPRESSED , col = SEX)) +
    geom_ribbon(aes(ymin = INFECTED_NON_SUPPRESSED_CL, ymax = INFECTED_NON_SUPPRESSED_CU , fill = SEX), alpha = 0.5) +
    labs(y = 'Number of HIV+ unsupressed census eligible', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_unsuppressed_round.png'), w = 7, h = 8)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED_NON_SUPPRESSED , col = ROUND)) +
    geom_ribbon(aes(ymin = INFECTED_NON_SUPPRESSED_CL, ymax = INFECTED_NON_SUPPRESSED_CU , fill = ROUND), alpha = 0.5) +
    labs(y = 'Number of HIV+ unsupressed census eligible', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_unsuppressed_round2.png'), w = 7, h = 5)
  
  # incidence rate per person per round
  ggplot(incidence_cases_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENCE* ROUND_SPANYRS, col = SEX)) +
    geom_ribbon(aes(ymin = LB* ROUND_SPANYRS, ymax = UB* ROUND_SPANYRS, fill = SEX),  alpha = 0.5) +
    labs(y = 'Incidence rate per 1 person per round', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-incidence_rate_round.png'), w = 7, h = 8)
  
  # incidence cases
  ggplot(incidence_cases_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENT_CASES , col = SEX)) +
    geom_ribbon(aes(ymin = INCIDENT_CASES_LB , ymax = INCIDENT_CASES_UB , fill = SEX), alpha = 0.5) +
    labs(y = 'Number of incident cases', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-incidence_case_round.png'), w = 7, h = 8)
}

plot_data_by_period <- function(eligible_count, incidence_cases, proportion_sampling, outdir){
  
  eligible_count_wide <- dcast.data.table(eligible_count, SEX + COMM + AGEYRS + BEFORE_CUTOFF + PERIOD ~ variable, value.var = 'count')
  
  # Census eligible count 
  ggplot(eligible_count_wide, aes(x = AGEYRS)) +
    geom_line(aes(y = ELIGIBLE, col = SEX)) +
    labs(y = 'Census eligible count', x = 'Age') +
    facet_grid(PERIOD~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_count_period.png'), w = 7, h = 6)
  
  # Susceptible census eligible count
  ggplot(eligible_count_wide, aes(x = AGEYRS)) +
    geom_line(aes(y = SUSCEPTIBLE , col = SEX)) +
    labs(y = 'Census eligible susceptible count', x = 'Age') +
    facet_grid(PERIOD~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_susceptible_period.png'), w = 7, h = 6)
  
  eligible_count_wide[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  ggplot(eligible_count_wide, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_SUSCEPTIBLE , col = PERIOD)) +
    labs(y = 'Proportion of susceptible among census eligible individuals', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_susceptible_period.png'), w = 7, h = 6)
  
  #HIV+ unsupressed census eligible count
  ggplot(eligible_count_wide, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED_NON_SUPPRESSED , col = SEX)) +
    geom_ribbon(aes(ymin = INFECTED_NON_SUPPRESSED_CL, ymax = INFECTED_NON_SUPPRESSED_CU , fill = SEX), alpha = 0.5) +
    labs(y = 'Number of HIV+ unsupressed census eligible', x = 'Age') +
    facet_grid(PERIOD~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_unsuppressed_period.png'), w = 7, h = 6)
  
  # incidence cases
  ggplot(incidence_cases, aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENT_CASES , col = SEX)) +
    geom_ribbon(aes(ymin = INCIDENT_CASES_LB , ymax = INCIDENT_CASES_UB , fill = SEX), alpha = 0.5) +
    labs(y = 'Number of incident cases', x = 'Age') +
    facet_grid(PERIOD~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-incidence_case_period.png'), w = 7, h = 6)
  
  ggplot(proportion_sampling, aes(x = AGEYRS)) +
    geom_line(aes(y = prop_sampling , col = SEX)) +
    labs(y = 'Probability of observing a transmission event', x = 'Age') +
    facet_grid(PERIOD~COMM, label = 'label_both', scale = 'free_y') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-proportion_sampling_period.png'), w = 7, h = 6)
}

plot_offset <- function(stan_data, outdir){
  tmp <- as.data.table(reshape2::melt(stan_data[['log_offset']]))
  setnames(tmp, 1:4, c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'))
  tmp <- merge(tmp, df_direction, by = 'INDEX_DIRECTION')
  tmp <- merge(tmp, df_community, by = 'INDEX_COMMUNITY')
  tmp <- merge(tmp, df_period, by = 'INDEX_TIME')
  tmp <- merge(tmp, df_age, by = 'INDEX_AGE')
  tmp <- tmp[COMM == 'inland']
  tmp1 <- tmp[, list(value = sum(exp(value))), by = c('AGE_INFECTION.RECIPIENT', 'LABEL_DIRECTION', 'LABEL_COMMUNITY', 'PERIOD')]
  
  ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
    geom_raster(aes(fill = exp(value))) +
    facet_grid( PERIOD~LABEL_DIRECTION + LABEL_COMMUNITY) + 
    theme_bw() +
    labs(x= 'Age at infection recipient', y = 'Age at transmission source', fill = 'offset') 
  ggsave(paste0(outdir, '-offset-value.png'), w = 7, h = 8)
  
  tmp <- as.data.table(reshape2::melt(stan_data[['y']]))
  setnames(tmp, 1:4, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'))
  tmp <- merge(tmp, df_direction, by = 'INDEX_DIRECTION')
  tmp <- merge(tmp, df_community, by = 'INDEX_COMMUNITY')
  tmp <- merge(tmp, df_period, by = 'INDEX_TIME')
  tmp <- merge(tmp, df_age, by = 'INDEX_AGE')
  tmp <- tmp[COMM == 'inland']
  tmp2 <- tmp[, list(value = sum(value)), by = c('AGE_INFECTION.RECIPIENT', 'LABEL_DIRECTION', 'LABEL_COMMUNITY', 'PERIOD')]
  
  ggplot(tmp1, aes(x = AGE_INFECTION.RECIPIENT, y = (value), col = PERIOD)) + 
    geom_line() +
    facet_grid( LABEL_COMMUNITY~LABEL_DIRECTION ) + 
    theme_bw() +
    labs(x=  'Age at infection recipient', y = 'total offset') 
  ggsave(paste0(outdir, '-offset-value_aggregated.png'), w = 7, h = 5)
  
  ggplot(tmp2, aes(x = AGE_INFECTION.RECIPIENT, y = (value), col = PERIOD)) + 
    geom_line() +
    facet_grid( LABEL_COMMUNITY~LABEL_DIRECTION ) + 
    theme_bw() +
    labs(x=  'Age at infection recipient', y = 'total transmission events') 
  ggsave(paste0(outdir, '-offset-total_count.png'), w = 7, h = 5)
}
