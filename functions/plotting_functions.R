plot_age_infection_source_recipient <- function(data, title, plotlab, outdir = NULL){
  
  plots = list()
  
  data <- data[!is.na(age_transmission.SOURCE) & !is.na(age_infection.RECIPIENT)]
  
  data[, `Community source` := comm.SOURCE]
  data[, `Community recipient` := comm.RECIPIENT]
  
  data[, cohort_round.SOURCE := substr(round.SOURCE, 1, 4)]
  data[, cohort_round.RECIPIENT := substr(round.RECIPIENT, 1, 4)]
  data[, `Cohort round recipient` := cohort_round.RECIPIENT]
  data[, `Cohort round source` := cohort_round.SOURCE]
  
  # all pairs
  p <- ggplot(data, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) 
    if(!is.null(outdir))
      ggsave(p, filename = paste0(outdir, '-AgeInfection_AllPairs_', plotlab, '.png'), w = 4, h = 4)
  plots = c(plots, list(p))
  
  # by cohort round
  p1 <- ggplot(data, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round source`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')+
    guides(col=guide_legend(nrow=2,byrow=TRUE))

  p2 <- ggplot(data, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Cohort round recipient`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom') +
    guides(col=guide_legend(nrow=2,byrow=TRUE))

  p <- ggarrange(p1, p2, ncol = 2)
  if(!is.null(outdir))
    ggsave(p, filename = paste0(outdir, '-AgeInfection_CohortRound_', plotlab, '.png'), w = 9, h = 7)
  plots = c(plots, list(p))
  
  # by age infection round
  p <- ggplot(data, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_point(aes(col = date_infection_before_cutoff.RECIPIENT)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient',
         col = paste0('Date infection recipient before ', format(cutoff_date, '%Y'))) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')
  if(!is.null(outdir))
    ggsave(p, filename = paste0(outdir, '-AgeInfection_DateInfectionRecipient_', plotlab, '.png'), w = 4, h = 4)
  plots = c(plots, list(p))
  
  # by community
  p1 <- ggplot(data, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_point(aes(col = `Community source`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')

  p2 <- ggplot(data, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_point(aes(col =`Community recipient`)) + 
    labs(y = 'Age at transmission source', x = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_transmission.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    theme(legend.position = 'bottom')
  
  p <- ggarrange(p1, p2, ncol = 2)
  if(!is.null(outdir))
    ggsave(p, filename = paste0(outdir, '-AgeInfection_CommunityRecipient_', plotlab, '.png'), w = 9, h = 7)
  plots = c(plots, list(p))
  
  return(plots)
}

plot_hist_age_infection <- function(pairs, outdir = NULL){
  
  pairs[, Sex := sex.SOURCE]
  p1 <- ggplot(pairs, aes(x = age_transmission.SOURCE)) + 
    geom_histogram(bins = 30) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at transmission source') +
    scale_x_continuous(limits = range(c(pairs$age_transmission.SOURCE, pairs$age_infection.RECIPIENT)))
  
  pairs[, Sex := sex.RECIPIENT]
  p2 <- ggplot(pairs, aes(x = age_infection.RECIPIENT)) + 
    geom_histogram(bins = 30) +     
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$age_transmission.SOURCE, pairs$age_infection.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  
  if(!is.null(outdir)){
    file = paste0(outdir, '-hist_age_infection.png')
    cat('saving', file)
    ggsave(p, file = file, w = 6, h = 6)
  }
  
  return(p)
}

plot_hist_time_infection <- function(pairs, cutoff_date, outdir = NULL){
  
  pairs[, cohort_round.SOURCE := substr(round.SOURCE, 1, 4)]
  pairs[, cohort_round.RECIPIENT := substr(round.RECIPIENT, 1, 4)]
  pairs[, `Round source` := cohort_round.SOURCE]
  pairs[, `Round recipient` := cohort_round.RECIPIENT]
  
  p1 <- ggplot(pairs, aes(x = date_infection.SOURCE)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round source`, nrow = length(unique(pairs$cohort_round.SOURCE))) +
    theme_bw() + 
    labs(x = 'Date of infection source') + 
    geom_vline(xintercept = cutoff_date, linetype = 'dashed')
  
  p2 <- ggplot(pairs, aes(x = date_infection.RECIPIENT)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~`Round recipient`, nrow = length(unique(pairs$cohort_round.RECIPIENT))) +
    theme_bw() + 
    labs(x = 'Date of infection recipient') + 
    geom_vline(xintercept = cutoff_date, linetype = 'dashed')

  p <- ggarrange(p1, p2, ncol = 2)
  
  if(!is.null(outdir)){
    file = paste0(outdir, '-hist_date_infection.png')
    cat('saving', file)
    ggsave(p, file = file, w = 6, h = 6)
  }

  return(p)
}

plot_CI_age_infection <- function(pairs, outdir = NULL){
  
  data <- copy(pairs)
  data[, age_transmission.SOURCE := floor(age_transmission.SOURCE)]
  data[, age_infection.RECIPIENT := floor(age_infection.RECIPIENT)]
  data <- merge(data, df_age, by = c('age_infection.RECIPIENT', 'age_transmission.SOURCE'))
  
  ps <- c(0.5, 0.2, 0.8)
  p_labs <- c('M','CL','CU')
  
  ## stratified by age of recipient
  tmp = data[, list(q= quantile(age_transmission.SOURCE, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('sex.SOURCE', 'sex.RECIPIENT', 'date_infection_before_cutoff.RECIPIENT', 'age_infection_reduced.RECIPIENT')]	
  tmp = dcast(tmp, sex.SOURCE + sex.RECIPIENT + date_infection_before_cutoff.RECIPIENT + age_infection_reduced.RECIPIENT ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[date_infection_before_cutoff.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection male recipient', y = 'Age at transmission female source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  # MF
  tmp1 <- subset(tmp, sex.SOURCE == 'M' & sex.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at infection female recipient', y = 'Age at transmission male source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = paste0(outdir, '-AgeTransmissionSource.png'), w = 8, h = 5)
  
  
  ## not stratified by age of recipient
  tmp = data[, list(q= quantile(age_transmission.SOURCE, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('sex.SOURCE', 'sex.RECIPIENT', 'date_infection_before_cutoff.RECIPIENT')]	
  tmp = dcast(tmp, sex.SOURCE + sex.RECIPIENT + date_infection_before_cutoff.RECIPIENT  ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[date_infection_before_cutoff.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = date_infection_before_cutoff.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Date infection recipient', y = 'Age at transmission female source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  # MF
  tmp1 <- subset(tmp, sex.SOURCE == 'M' & sex.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = date_infection_before_cutoff.RECIPIENT)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Date infection recipient', y = 'Age at transmission male source', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = paste0(outdir, '-AgeTransmissionSource_Unstratified.png'), w = 8, h = 5)
  
}

plot_CI_age_transmission <- function(pairs, outdir = NULL){
  
  data <- copy(pairs)
  data[, age_transmission.SOURCE := floor(age_transmission.SOURCE)]
  data[, age_infection.RECIPIENT := floor(age_infection.RECIPIENT)]
  data <- merge(data, df_age, by = c('age_infection.RECIPIENT', 'age_transmission.SOURCE'))
  
  ps <- c(0.5, 0.2, 0.8)
  p_labs <- c('M','CL','CU')
  
  ## stratified by age of recipient
  tmp = data[, list(q= quantile(age_infection.RECIPIENT, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('sex.SOURCE', 'sex.RECIPIENT', 'date_infection_before_cutoff.RECIPIENT', 'age_transmission_reduced.SOURCE')]	
  tmp = dcast(tmp, sex.SOURCE + sex.RECIPIENT + date_infection_before_cutoff.RECIPIENT + age_transmission_reduced.SOURCE ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[date_infection_before_cutoff.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M') 
  p1 <- ggplot(tmp1, aes(x = age_transmission_reduced.SOURCE)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at transmission female source', y = 'Age at infection male recipient', 
         col ='Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  # MF
  tmp1 <- subset(tmp, sex.SOURCE == 'M' & sex.RECIPIENT == 'F') 
  p2 <- ggplot(tmp1, aes(x = age_transmission_reduced.SOURCE)) + 
    geom_point(aes(y = M, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5)) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2) +
    labs(x = 'Age at transmission male source', y = 'Age at infection female recipient', 
         col = 'Date infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    # coord_fixed() + 
    theme(legend.position = 'bottom') 
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  ggsave(p, filename = paste0(outdir, '-AgeInfectionRecipient.png'), w = 8, h = 5)
  
  
  ## not stratified by age of recipient
  tmp = data[, list(q= quantile(age_infection.RECIPIENT, prob=ps, na.rm = T), q_label=p_labs), 
             by=c('sex.SOURCE', 'sex.RECIPIENT', 'date_infection_before_cutoff.RECIPIENT')]	
  tmp = dcast(tmp, sex.SOURCE + sex.RECIPIENT + date_infection_before_cutoff.RECIPIENT  ~ q_label, value.var = "q")
  tmp[, date_infection_before_cutoff_name.RECIPIENT := paste0('Before ', format(cutoff_date, '%Y'))]
  tmp[date_infection_before_cutoff.RECIPIENT == 0, date_infection_before_cutoff_name.RECIPIENT := paste0('After ', format(cutoff_date, '%Y'))]
  tmp[, date_infection_before_cutoff_name.RECIPIENT := factor(date_infection_before_cutoff_name.RECIPIENT, levels = paste0(c('Before ', 'After '), format(cutoff_date, '%Y')))]
  
  # FM
  tmp1 <- subset(tmp, sex.SOURCE == 'F' & sex.RECIPIENT == 'M') 
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
  tmp1 <- subset(tmp, sex.SOURCE == 'M' & sex.RECIPIENT == 'F') 
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
  ggsave(p, filename = paste0(outdir, '-AgeInfectionRecipient_Unstratified.png'), w = 8, h = 5)
}

plot_pairs_infection_dates <- function(pairs.all)
{
  
  tmp <- copy(pairs.all)
  tmp[, DUMMY:= date_infection.SOURCE - date_infection.RECIPIENT]
  setkey(tmp, date_first_positive.SOURCE)
  tmp[, CONS:=ifelse(DUMMY < 0, 'consistent', 'inconsistent')]
  setkey(tmp, DUMMY)
  # tmp[, DUMMY:=1:.N]
  tmp[date_first_positive.SOURCE < date_first_positive.RECIPIENT ] # (316 pairs for which )
  tmp[date_infection.SOURCE < date_infection.RECIPIENT  ]
  
  p <- ggplot(tmp, aes(y=date_first_positive.SOURCE, col=CONS)) + 
    geom_point(aes(x=date_infection.SOURCE)) + 
    geom_errorbarh(aes(xmin=date_infection.SOURCE, xmax=date_infection.RECIPIENT)) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(data=tmp[date_first_positive.SOURCE > date_first_positive.RECIPIENT],aes(x=date_first_positive.RECIPIENT), color='black', pch=18) +
    scale_x_date(date_labels = '%Y', breaks = '12 months') +
    scale_y_date(date_labels = '%Y', breaks = '12 months') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position = 'bottom') +
    labs(x='Estimated infection dates', y=" Source's date of first positive test", col='source-date relationship')
  
  ggsave(p, filename = file.path(outdir, paste0('PairsInfectionDates.png')), w = 8, h = 10)
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
  tmp.dchain <- merge(pairs[,.(RECIPIENT,SOURCE, sex.RECIPIENT)], tmp.dchain, by=c('SOURCE', 'RECIPIENT'))
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
  di <- merge(di, pairs[, .(RECIPIENT, sex.RECIPIENT)], by.x='ID' , by.y='RECIPIENT', all.x=TRUE)
  di <- unique(merge(di, pairs[, .(SOURCE, sex.SOURCE)], by.x='ID' , by.y='SOURCE', all.x=TRUE))
  di[, NODE_FILL:=na.omit(c(sex.RECIPIENT, sex.SOURCE))[1], by=ID]
  di[, `:=`(sex.RECIPIENT=NULL, sex.SOURCE=NULL)]
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
