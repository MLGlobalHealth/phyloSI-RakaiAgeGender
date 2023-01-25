naturemed_reqs <- function() 
{
    # call this before doing your plots
    reqs <<- theme(axis.text = element_text(size=5, family='sans'), text=element_text(size=7,family='sans'), legend.text=element_text(size=7, family='sans'))
}

ggarrange_nature <- function(
  ...,
  plotlist = NULL,
  ncol = NULL,
  nrow = NULL,
  labels = NULL,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  align = c("none", "h", "v", "hv"),
  widths = 1,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL, 
  add_reqs=TRUE
){
    if(add_reqs)
        reqs <- theme(axis.text = element_text(size=5, family='sans'), text=element_text(size=7,family='sans'), legend.text=element_text(size=7, family='sans'))

    plots <- c(list(...), plotlist)

    if(add_reqs)
        plots <- lapply(plots, function(p){p + reqs})

    out <- ggarrange(plotlist = plots,
                     ncol = ncol,
                     nrow = nrow,
                     labels = labels,
                     label.x = label.x,
                     label.y = label.y,
                     hjust = hjust,
                     vjust = vjust,
                     font.label = list(size = 8, color = "black", face = "bold", family = 'sans'),
                     align = align,
                     widths = widths,
                     heights = heights,
                     legend = legend,
                     common.legend = common.legend,
                     legend.grob = legend.grob)
    return(out)
}

ggsave_nature <- function(filename, p, w=18,h=24, add_reqs=TRUE)
{
    # check size
    tmp <- sort(c(w,h))
    if(tmp[1] > 18 | tmp[2] > 24)
        warning('Plot is bigger than allowed for EDFs. Maximum size is 18cm x 24cm\n')
    if( tmp[1] < 10)
        warning('w and h represent cm units, not inches. Are you sure you want to save such a small plot?\n')

    # apply changes: for the moment only works for simple plots, but not for ggarrange.
    # Let's see if anyone answers this:
    # https://stackoverflow.com/questions/74379207/simultaneously-applying-same-modification-to-all-ggarrange-subplots
    if(add_reqs)
    {
        reqs <- theme(axis.text = element_text(size=5, family='sans'), text=element_text(size=7,family='sans'), legend.text=element_text(size=7, family='sans'))
        p <- p + reqs
    }

    # save
    ggsave(filename=filename, plot=p, width=w, height=h, units='cm', dpi=310)
}



find_palette_round <- function()
{
  # palette_round <- scales::viridis_pal(option = 'A', end= 0.9)(8)
  palette_round <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(10)
  palette_round_inland <<- palette_round[c(1:6, 8:10)]
  palette_round_fishing <<- palette_round[c(6:10)]
}


plot_age_infection_source_recipient <- function(data, title, plotlab, cutoff_date, outdir = NULL)
{
  
  plots = list()
  
  data <- data[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
  
  data[, `Community source` := COMM.SOURCE]
  data[, `Community recipient` := COMM.RECIPIENT]
  
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

plot_hist_age_infection <- function(pairs, outdir = NULL)
{
  
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

plot_hist_time_infection <- function(pairs, cutoff_date, outdir = NULL)
{
  
  # find round of infection
  tmp <- merge(pairs, df_round,by.x = 'COMM.RECIPIENT', by.y = 'COMM', allow.cartesian = T)
  tmp <- tmp[DATE_INFECTION.RECIPIENT >= MIN_SAMPLE_DATE & DATE_INFECTION.RECIPIENT <= MAX_SAMPLE_DATE]
  setnames(tmp, 'LABEL_ROUND', 'LABEL_ROUND_BIS') # to avoid rect to be only in one facet
  
  # plot
  p <- ggplot(tmp) +
    geom_rect(data = df_round, aes(ymin = -Inf, ymax = Inf, xmin = MIN_SAMPLE_DATE, 
                                    xmax = MAX_SAMPLE_DATE, fill = ROUND), alpha = 0.5) + 
    geom_histogram(aes(x = DATE_INFECTION.RECIPIENT), bins = 100) +
    facet_grid(LABEL_ROUND_BIS~.) +
    theme_bw() + 
    labs(x = 'Date of infection recipient') + 
    geom_vline(xintercept = cutoff_date, linetype = 'dashed') + 
    scale_y_continuous(expand = expansion(mult = c(0, .05)))+ 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1))) +
    ggtitle('Inland communities')
  file = paste0(outdir, '-data-hist_date_infection_inland.png')
  ggsave(p, file = file, w = 10.5, h = 12)
  
  return(p)
}

plot_CI_age_infection <- function(pairs, cutoff_date, outdir = NULL)
{
  
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
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2, col = 'grey30') +
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
    geom_errorbar(aes(ymin = CL, ymax = CU, col = date_infection_before_cutoff_name.RECIPIENT), position = position_dodge(1.5), width = 0.2, col = 'grey30') +
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

plot_CI_age_transmission <- function(pairs, cutoff_date, outdir = NULL)
{
  
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


plot_data_by_round <- function(eligible_count_round, treatment_cascade, proportion_prevalence, outdir)
{
  
  level_rounds <- c('R010', 'R011', 'R012', 'R013', 'R014', 'R015', 'R015S', 'R016', 'R017', 'R018')
  
  #
  # Calendar time of round survey 
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
  
  #
  # Census eligible count 
  tmp <- copy(eligible_count_round)
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp <- merge(tmp, df_round, by = c('COMM', 'ROUND'))
  tmp[, LABEL_ROUND := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  tmp[, SEX_LABEL := 'Female']
  tmp[SEX== 'M', SEX_LABEL := 'Male']
  ggplot(tmp, aes(x = AGEYRS)) +
    geom_line(aes(y = ELIGIBLE, col = SEX_LABEL)) +
    labs(y = 'Census eligible population count', x = 'Age', col = '') +
    facet_grid(LABEL_COMMUNITY~LABEL_ROUND) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1))) +
    scale_color_manual(values = c('Male'='royalblue3','Female'='deeppink')) 
  ggsave(paste0(outdir, '-data-census_eligible_count_round_sex.png'), w = 12, h = 7)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = ELIGIBLE, col = ROUND)) +
    labs(y = 'Census eligible count', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_count_round.png'), w = 7, h = 6)
  
  
  #
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
  ggplot(tmp, aes(x = AGEYRS)) +
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
  
  
  #
  # Susceptible proportion
  ggplot(tmp, aes(x = AGEYRS)) +
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
  
  #
  # Age distribution of susceptible
  tmp[, PROP_SUSCEPTIBLE:= 1-PREVALENCE_M]
  tmp[, TOTAL_PROP_SUSCEPTIBLE := sum(PROP_SUSCEPTIBLE), by = c('ROUND', 'SEX', 'COMM')]
  tmp[, distr := PROP_SUSCEPTIBLE / TOTAL_PROP_SUSCEPTIBLE]
  ggplot(tmp, aes(x = AGEYRS, y = distr)) + 
    geom_line(aes(col = ROUND)) + 
    facet_grid(SEX~COMM) + 
    theme_bw() + 
    labs(x = 'Age', y = 'Probability distribution of the age composition of proportion of susceptible')
  ggsave(paste0(outdir, '-data-distribution_function_age_composition_susceptible.png'), w = 7, h = 7)
  
  
  #
  # Infected count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED , col = SEX)) +
    labs(y = 'Census eligible infected count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_infected_round_sex.png'), w = 7, h = 9)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED , col = ROUND)) +
    labs(y = 'Census eligible infected count', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_infected_round.png'), w = 7, h = 6)
  
  
  #
  # Susceptible count
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = SUSCEPTIBLE , col = SEX)) +
    labs(y = 'Census eligible susceptible count', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_susceptible_round_sex.png'), w = 7, h = 9)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = SUSCEPTIBLE , col = ROUND)) +
    labs(y = 'Census eligible susceptible count', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_susceptible_round.png'), w = 7, h = 6)
  
  #
  # proportion of unsuppressed among participants
  ggplot(treatment_cascade, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_UNSUPPRESSED_PARTICIPANTS_M , col = SEX)) +
    geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_PARTICIPANTS_CL , ymax = PROP_UNSUPPRESSED_PARTICIPANTS_CU , fill = SEX), alpha = 0.5) +
    labs(y = 'Proportion of unsupressed', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_unsuppressed_participants_round_sex.png'), w = 7, h = 8)
  
  ggplot(treatment_cascade, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_UNSUPPRESSED_PARTICIPANTS_M , col = ROUND)) +
    geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_PARTICIPANTS_CL , ymax = PROP_UNSUPPRESSED_PARTICIPANTS_CU , fill = ROUND), alpha = 0.5) +
    labs(y = 'Proportion of unsupressed', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_unsuppressed_participants_round.png'), w = 7, h = 6)
  
  
  #
  # proportion of unsuppressed among non-participants
  ggplot(treatment_cascade, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_UNSUPPRESSED_NONPARTICIPANTS_M , col = SEX)) +
    geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_NONPARTICIPANTS_CL , ymax = PROP_UNSUPPRESSED_NONPARTICIPANTS_CU , fill = SEX), alpha = 0.5) +
    labs(y = 'Proportion of unsupressed', x = 'Age') +
    facet_grid(ROUND~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_unsuppressed_nonparticipants_round_sex.png'), w = 7, h = 8)
  
  ggplot(treatment_cascade, aes(x = AGEYRS)) +
    geom_line(aes(y = PROP_UNSUPPRESSED_NONPARTICIPANTS_M , col = ROUND)) +
    geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_NONPARTICIPANTS_CL , ymax = PROP_UNSUPPRESSED_NONPARTICIPANTS_CU , fill = ROUND), alpha = 0.5) +
    labs(y = 'Proportion of unsupressed', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_prop_unsuppressed_nonparticipants_round.png'), w = 7, h = 6)
  
  
  #
  # infected unsupressed count
  tmp <- copy(eligible_count_round)
  tmp[, SEX_LABEL := 'Women']
  tmp[SEX == 'M', SEX_LABEL := 'Men']
  tmp[, ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]
  ggplot(tmp, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED_NON_SUPPRESSED , col = SEX_LABEL)) +
    geom_ribbon(aes(ymin = INFECTED_NON_SUPPRESSED_CL, ymax = INFECTED_NON_SUPPRESSED_CU , fill = SEX_LABEL), alpha = 0.5) +
    labs(y = 'Number of HIV+ unsupressed census eligible', x = 'Age') +
    facet_wrap(~ROUND_LABEL, ncol = 2) +
    scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    scale_color_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    theme_bw() +
    theme(legend.position = 'bottom',    
          strip.background = element_rect(colour="white", fill="white"),
          legend.title = element_blank()) 
  ggsave(paste0(outdir, '-data-census_eligible_unsuppressed_round_sex.png'), w = 4, h = 7)
  
  ggplot(eligible_count_round, aes(x = AGEYRS)) +
    geom_line(aes(y = INFECTED_NON_SUPPRESSED , col = ROUND)) +
    # geom_ribbon(aes(ymin = INFECTED_NON_SUPPRESSED_CL, ymax = INFECTED_NON_SUPPRESSED_CU , fill = ROUND), alpha = 0.5) +
    labs(y = 'Number of HIV+ unsupressed census eligible', x = 'Age') +
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-data-census_eligible_unsuppressed_round.png'), w = 7, h = 6)
  
}

plot_data_by_period <- function(incidence_cases, outdir)
{
  
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

plot_offset <- function(stan_data, outdir)
{
  
  tmp <- as.data.table(reshape2::melt(stan_data[['log_offset']] + stan_data[['log_offset_time']]))
  setnames(tmp, 1:3, c('INDEX_DIRECTION', 'INDEX_ROUND', 'INDEX_AGE'))
  tmp <- merge(tmp, df_direction, by = 'INDEX_DIRECTION')
  tmp <- merge(tmp, df_round, by = c('INDEX_ROUND'))
  tmp <- merge(tmp, df_age, by = 'INDEX_AGE')
  tmp <- merge(tmp, df_community, by = 'COMM')

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


plot_transmission_events_over_time <- function(pairs, outdir, nm_reqs=FALSE){
  
  # timeline
  df_timeline <- copy(df_round)
  df_timeline[, MIDPOINT := as.Date(mean(c(MIN_SAMPLE_DATE, MAX_SAMPLE_DATE))), by = c('ROUND', 'COMM')]
  df_timeline <- df_timeline[, .(ROUND, MIDPOINT, COMM, INDEX_ROUND, ROUND_SPANYRS)]
  
  # age groups
  age_groups <- c('15-24', '25-34', '35-49')
  df_age_group <- data.table(AGEYRS = 15:49)
  df_age_group[, index_age_group := 3]
  df_age_group[AGEYRS < 35, index_age_group := 2]
  df_age_group[AGEYRS < 25, index_age_group := 1]
  df_age_group[, age_group := age_groups[index_age_group]]
  df_age_group[, AGE_GROUP_LABEL := paste0('Age recipient: ', age_group)]
  
  # sex
  df_sex <- data.table(SEX = c('M', 'F'), SEX_LABEL = c('Male', 'Female'))
  
  # grid
  df_grid <- data.table(expand.grid(SEX_LABEL = df_sex[, unique(SEX)], 
                                    age_group = df_age_group[, unique(age_group)], 
                                    INDEX_ROUND = df_timeline[, unique(INDEX_ROUND)], 
                                    COMM=c('inland', 'fishing')))
  df_grid <- merge(df_grid, unique(df_age_group[, .(age_group, AGE_GROUP_LABEL)]), by = 'age_group')
  df_grid <- merge(df_grid, (df_timeline), by = c('INDEX_ROUND', 'COMM'))
  
  # Prepare phylo pairs
  dp <- copy(pairs)
  setnames(dp, c('SEX.RECIPIENT', 'COMM.RECIPIENT', 'AGE_INFECTION.RECIPIENT', 'DATE_INFECTION.RECIPIENT'), 
           c('SEX', 'COMM', 'AGEYRS', 'DATE'))
  dp[, AGEYRS := floor(AGEYRS)]
  dp <- merge(dp, df_age_group, by = 'AGEYRS')
  dp[, DIRECTION := 'Male to Female']
  dp[SEX.SOURCE == 'F', DIRECTION := 'Female to Male']

  
  #
  # Plot
  
  communities <- df_round[, unique(COMM)]
  male_color <- 'lightblue3'
  female_color <- 'lightpink2'
  
  for(i in seq_along(communities)){
    
    comm <- communities[i]

    # plot pairs 
    p3 <- ggplot(dp[COMM == comm]) + 
      geom_histogram(aes(x = DATE, fill = DIRECTION), bins = 30) + 
      # geom_vline(aes(xintercept = pairs[, min(SAMPLE_DATE.RECIPIENT, SAMPLE_DATE.SOURCE)], linetype='Start deep-sequencing'), color = 'grey50') + 
      facet_grid(.~AGE_GROUP_LABEL) + 
      labs(y = paste0('Detected transmissions\nfrom deep-sequence data'), x = 'Date at transmission') + 
      theme_bw() + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            # axis.text.x = element_text(angle= 70, hjust = 1),
            strip.text = element_text(size = rel(1)), 
            legend.position = 'bottom',
            legend.title = element_blank()) + 
      scale_linetype_manual(values = 'dashed') + 
      scale_fill_manual(values = c('Male to Female'=male_color,'Female to Male'=female_color)) +
      scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, 0.05)),
                         breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) + 
      scale_x_date(limits = c(df_period[, min(MIN_PERIOD_DATE)], df_period[, max(MAX_PERIOD_DATE)]), expand = c(0,0))  

    
    ggsave(p3, file =  paste0(outdir, '-data-detected_transmission_events_', communities[i], '.pdf'), w = 5.2, h = 3.5)
    
  }
  
    if(nm_reqs){ p3 <- p3 + reqs }
    return(p3)
}

plot_incident_cases_over_time <- function(incidence_cases_round, participation, outdir){
  
  #
  # incidence cases among participants and non participants
  tmp <- copy(incidence_cases_round)
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp <- merge(tmp, df_round, by = c('COMM', 'ROUND', 'ROUND_SPANYRS'))
  tmp[, LABEL_ROUND := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  tmp[, SEX_LABEL := 'Female']
  tmp[SEX== 'M', SEX_LABEL := 'Male']
  tmp[, type := 'Census eligible population']
  
  pa <- copy(participation)
  pa[, ROUND := paste0('R0', ROUND)]
  tmp1 <- merge(tmp,pa, by = c('ROUND', 'SEX', 'AGEYRS', 'COMM'))
  tmp1[, `:=` (INCIDENT_CASES = INCIDENT_CASES * PARTICIPATION, 
               INCIDENT_CASES_LB = INCIDENT_CASES_LB * PARTICIPATION, 
               INCIDENT_CASES_UB = INCIDENT_CASES_UB * PARTICIPATION)]
  tmp1[, type := 'Participants']
  tmp <- rbind(tmp, tmp1, fill=TRUE)
  ggplot(tmp, aes(x = AGEYRS, group = interaction(SEX_LABEL, type))) +
    geom_line(aes(y = INCIDENT_CASES/ROUND_SPANYRS, col = SEX_LABEL, linetype = type)) +
    geom_ribbon(aes(ymin = INCIDENT_CASES_LB /ROUND_SPANYRS, ymax = INCIDENT_CASES_UB /ROUND_SPANYRS, 
                    fill = SEX_LABEL, group = interaction(SEX_LABEL, type)), alpha = 0.25)  +
    labs(y = 'Number of incident cases per year', x = 'Age', col= '', fill = '', linetype = 'Among') +
    facet_grid(LABEL_COMMUNITY~LABEL_ROUND) +
    theme_bw() +
    scale_color_manual(values = c('Male'='royalblue3','Female'='deeppink')) + 
    scale_fill_manual(values = c('Male'='lightblue3','Female'='lightpink1')) +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"))
  ggsave(paste0(outdir, '-data-incidence_case_round.png'), w = 12, h = 7)
  
  #
  # empirical contribution 
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
    scale_y_continuous(labels = scales::percent_format(), limits = c(0,1))
  ggsave(paste0(outdir, '-data-empirical_contribution_infection.png'), w = 7, h = 6)
  
  
  #
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

plot_date_collection_pairs <- function(pairs, df_round_inland, outdir){
  
  pairs[, DATE.COLLECTION.PAIR := max(c(DATE.COLLECTION.SOURCE, DATE.COLLECTION.RECIPIENT)), by = c('RECIPIENT', 'SOURCE')]
  
  ggplot(pairs) + 
    geom_rect(data = df_round_inland, aes(ymin = -Inf, ymax = Inf, xmin = min_sample_date, xmax = max_sample_date, fill = round), alpha = 0.5) + 
    geom_histogram(aes(x = DATE.COLLECTION.PAIR), bins = 100) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult= c(0, 0.01))) + 
    labs(x = 'Date collection pairs\nmax(date collection source, date collection recipient)') + 
    theme_bw()
  ggsave(paste0(outdir, '-data-pairs_date_collection_distribution.png'), w = 7, h = 6)
  
}


plot_incident_rates_over_time <- function(incidence_cases_round, 
                                          incidence_rates_round.samples,
                                          eligible_count_round,
                                          outdir, outdir.table){
  
  #
  # median age at infection
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  icr <- copy(incidence_rates_round.samples)
  rm(incidence_rates_round.samples)
  icr[COMM == 'fishing', REF.ROUND := 'R015']
  icr[COMM == 'inland', REF.ROUND := 'R010']
  
  medage <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
  medage[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
  medage[, WEIGHTED_INCIDENCE := INCIDENT_CASES / sum(INCIDENT_CASES), by = c('COMM', 'ROUND', 'SEX', 'iterations')]
  medage <- medage[, list(value = matrixStats::weightedMedian(AGEYRS, WEIGHTED_INCIDENCE )), by = c('iterations', 'COMM', 'ROUND', 'SEX')]
  medage = medage[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'SEX')]	
  medage = dcast(medage, ... ~ q_label, value.var = "q")
  
  
  #
  # incidence rate per person per round
  
  tmp <- copy(incidence_cases_round)
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp <- merge(tmp, df_round, by = c('COMM', 'ROUND'))
  tmp[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  tmp[, SEX_LABEL := 'Women']
  tmp[SEX== 'M', SEX_LABEL := 'Men']
  ggplot(tmp, aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENCE*100, col = LABEL_ROUND2)) +
    geom_ribbon(aes(ymin = LB *100, ymax = UB* 100, fill = LABEL_ROUND2),  alpha = 0.1) +
    labs(y = 'Incidence rate per 100 person-years', x = 'Age', fill = '', col = '') +
    facet_grid(LABEL_COMMUNITY~SEX_LABEL,  scale = 'free_y') +
    theme_bw() +
    scale_color_manual(values = palette_round) + 
    scale_fill_manual(values = palette_round) + 
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"))
  ggsave(paste0(outdir, '-data-incidence_rate_round.png'), w = 7, h = 6)
  
  # prepare median age
  max_y_limits <- ifelse(tmp[, max(INCIDENCE)] > 0.021, 2.7, 2.1)
  median_age <-  copy(medage)
  median_age[, SEX_LABEL := 'Women']
  median_age[SEX== 'M', SEX_LABEL := 'Men']
  set.seed(12)
  median_age[ROUND == 'R018' & SEX_LABEL == 'Women', M := M + runif(length(M), 0, 1)]
  median_age[ROUND == 'R018' & SEX_LABEL == 'Men', M := M + runif(length(M), -1, 0)]
  tmp1 <- median_age[COMM == 'inland' & ROUND %in% c('R010', 'R012', 'R014', 'R016', 'R018')]
  tmp1 <- merge(tmp1, df_round, by = c('COMM', 'ROUND'))
  
  ggplot(tmp[COMM == 'inland' & round %in% c(10, 12, 14, 16, 18)]) +
    geom_line(aes(x = AGEYRS, y = INCIDENCE*100, col = SEX_LABEL)) +
    geom_ribbon(aes(x = AGEYRS, ymin = LB *100, ymax = UB* 100, fill = SEX_LABEL),  alpha = 0.5) +
    # geom_errorbarh(data = tmp1[SEX_LABEL=='Men' ], aes(y = 0.01, xmin = CL, xmax = CU, col = SEX_LABEL), size =1.5) +
    # geom_errorbarh(data = tmp1[SEX_LABEL=='Men' & ROUND == 'R018'], aes(y = 0.025, xmin = CL, xmax = CU, col = SEX_LABEL), size =1.5) +
    geom_point(data = tmp1[SEX_LABEL=='Men' ], aes(y = 0.08, x = M, fill = SEX_LABEL,  col = SEX_LABEL), shape = 25, size =3,  alpha = 0.5) +
    geom_point(data = tmp1[SEX_LABEL=='Men' ], aes(y = 0.08, x = M, col = SEX_LABEL), shape = 6, size =3) +
    # geom_errorbarh(data = tmp1[SEX_LABEL=='Women' ], aes(y = 0.01, xmin = CL, xmax = CU, col = SEX_LABEL), size =1.5) +
    geom_point(data = tmp1[SEX_LABEL=='Women'], aes(y = 0.08, x = M,  fill = SEX_LABEL,col = SEX_LABEL), shape = 25, size =3,  alpha = 0.5) +
    geom_point(data = tmp1[SEX_LABEL=='Women'], aes(y = 0.08, x = M,  col = SEX_LABEL), shape = 6, size =3) +
    labs(y = 'Incidence rates\nper 100 person-years', x = 'Age') +
    facet_grid(.~LABEL_ROUND, scale = 'free_y') +
    theme_bw() +
    scale_color_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    theme(legend.position = 'none', 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.title = element_blank(), 
          strip.text = element_text(size = 9.3), 
          axis.title = element_text(size = 12)) + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(15, 49, 5))) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) + 
    coord_cartesian(ylim= c(0, max_y_limits))
  ggsave(paste0(outdir, '-data-incidence_rate_round_sex_inland_short.pdf'), w = 8, h = 3.5)
  
  tmp[, WEIGHTED_INCIDENCE := INCIDENCE / sum(INCIDENCE), by = c('COMM', 'ROUND', 'SEX')]
  median_age <- tmp[, list(MEDIAN_AGEYRS =matrixStats::weightedMedian(AGEYRS, WEIGHTED_INCIDENCE ) ), by = c('COMM', 'LABEL_ROUND', 'SEX_LABEL', 'ROUND', 'round')]
  ggplot(tmp[COMM == 'inland' & round %in% c(10, 12, 14, 16, 18)], aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENCE*100, col = SEX_LABEL)) +
    geom_ribbon(aes(ymin = LB *100, ymax = UB* 100, fill = SEX_LABEL),  alpha = 0.5) +
    geom_point(data = median_age[COMM == 'inland' & SEX_LABEL=='Men' & round %in% c(10, 12, 14, 16, 18)], aes(y = 0.08, x = MEDIAN_AGEYRS, fill = SEX_LABEL, col = SEX_LABEL), shape = 25, size =3) +
    geom_point(data = median_age[COMM == 'inland' & SEX_LABEL=='Women' & round %in% c(10, 12, 14, 16, 18)], aes(y = 0.08, x = MEDIAN_AGEYRS, fill = SEX_LABEL, col = SEX_LABEL), shape = 25, size =3) +
    labs(y = 'Incidence rates\nper 100 person-years', x = 'Age') +
    facet_grid(.~LABEL_ROUND, scale = 'free_y') +
    theme_bw() +
    scale_color_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    theme(legend.position = 'none', 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.title = element_blank(), 
          strip.text = element_text(size = 9.3), 
          axis.title = element_text(size = 12)) + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(15, 49, 5))) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) + 
    coord_cartesian(ylim= c(0, max_y_limits))
  ggsave(paste0(outdir, '-data-incidence_rate_round_sex_inland_short2.pdf'), w = 8, h = 3.5)
  
  ggplot(tmp[COMM == 'inland'], aes(x = AGEYRS)) +
    geom_line(aes(y = INCIDENCE*100, col = SEX_LABEL)) +
    geom_ribbon(aes(ymin = LB *100, ymax = UB* 100, fill = SEX_LABEL),  alpha = 0.5) +
    # geom_point(data = median_age[COMM == 'inland'], aes(y = 0.08, x = MEDIAN_AGEYRS, fill = SEX_LABEL, col = SEX_LABEL), shape = 25, size =3) +
    labs(y = 'Incidence rate per 100 person-years\nin inland communities', x = 'Age') +
    facet_wrap(.~LABEL_ROUND, nrow = 2) +
    theme_bw() +
    scale_color_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.title = element_blank()) + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(15, 49, 5))) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) + 
    coord_cartesian(ylim= c(0, 2.1))
  ggsave(paste0(outdir, '-data-incidence_rate_round_sex_inland.pdf'), w = 8, h = 6)
  
  #
  # female-male incidence rate ratio
  fmr <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
  fmr[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
  fmr <- fmr[, list(INCIDENT_RATE_SUSCEPTIBLE = sum(INCIDENT_CASES) / sum(SUSCEPTIBLE)), by = c('SEX', 'ROUND', 'COMM', 'iterations')] 
  fmr <- dcast.data.table(fmr, COMM + ROUND + iterations ~ SEX, value.var = 'INCIDENT_RATE_SUSCEPTIBLE')
  fmr[, INCIDENCE_RATIO := `F` / `M`]
  fmr = fmr[, list(q= quantile(INCIDENCE_RATIO, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND')]	
  fmr = dcast(fmr, ... ~ q_label, value.var = "q")
  fmr <- merge(fmr, df_community, by = 'COMM')
  fmr <- merge(fmr, df_round, by = c('COMM', 'ROUND'))
  fmr[, MIDPOINT_DATE := MIN_SAMPLE_DATE + (MAX_SAMPLE_DATE - MIN_SAMPLE_DATE)/2]
  
  #
  # incidence rate relative to first round per person per round by 1-year age group

  icr[, INCIDENCE_REL := INCIDENCE.DRAW / INCIDENCE.DRAW[ROUND == REF.ROUND], by = c('COMM', 'AGEYRS', 'SEX', 'iterations')]
  tmp1 = icr[, list(q= quantile(INCIDENCE_REL, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'AGEYRS', 'SEX')]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  tmp1 <- merge(tmp1, df_round, by = c('COMM', 'ROUND'))
  tmp1[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  tmp1[, SEX_LABEL := 'Female']
  tmp1[SEX== 'M', SEX_LABEL := 'Male']
  
  ggplot(tmp1[COMM == 'inland' & round %in% c(12, 14, 16, 18)], aes(x = AGEYRS)) +
    geom_hline(yintercept = 1, alpha = 0.2) + 
    geom_line(aes(y = M, col = SEX_LABEL)) +
    geom_ribbon(aes(ymin =CL, ymax = CU, fill = SEX_LABEL),  alpha = 0.3) +
    # geom_point(data = select(tmp1[COMM == 'inland' & ROUND == 'R012'& SEX_LABEL=='Male'], -'LABEL_ROUND'), aes(y = 0.12, x = MEDIAN_AGEYRS, col = SEX_LABEL), shape = 25, fill = 'grey50',  size =3) + 
    # geom_point(data = select(tmp1[COMM == 'inland' & ROUND == 'R012'& SEX_LABEL=='Female'], -'LABEL_ROUND'), aes(y = 0.06, x = MEDIAN_AGEYRS, col = SEX_LABEL), shape = 25, fill = 'grey50',  size =3) + 
    # geom_point(data = tmp1[COMM == 'inland' & SEX_LABEL=='Male'], aes(y = 0.12, x = MEDIAN_AGEYRS, fill = SEX_LABEL, col = SEX_LABEL), shape = 25, size =3) +
    # geom_point(data = tmp1[COMM == 'inland' & SEX_LABEL=='Female'], aes(y = 0.06, x = MEDIAN_AGEYRS, fill = SEX_LABEL, col = SEX_LABEL), shape = 25, size =3) +
    labs(y = 'Incidence rate relative to round 10', x = 'Age') +
    facet_grid(.~LABEL_ROUND, scale = 'free_y') +
    theme_bw() +
    scale_color_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
    scale_fill_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
    theme(legend.position = c(0.9, 0.83), 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.title = element_blank()) + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(15, 49, 5))) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) + 
    coord_cartesian(ylim= c(0, 2.1))
  ggsave(paste0(outdir, '-data-incidence_rate_relative_round_sex.pdf'), w = 7.5, h = 4)
  
  ggplot(tmp1[COMM == 'inland' & round %in% c(12, 14,16,18)], aes(x = AGEYRS)) +
    geom_hline(yintercept = 1, alpha = 0.2) + 
    geom_line(aes(y = M, col = LABEL_ROUND, linetype = SEX_LABEL)) +
    geom_ribbon(data = tmp1[COMM == 'inland' & round %in% c(18) & SEX == 'F'], aes(ymin = CL, ymax = CU, fill = LABEL_ROUND, group = interaction(LABEL_ROUND, SEX_LABEL)),  alpha = 0.3) +
    labs(y = 'Incidence rate relative to round 10', x = 'Age', col = '', fill = '') +
    theme_bw() +
    scale_color_manual(values = palette_round_inland[c(3, 5,7,9)]) + 
    scale_fill_manual(values =  palette_round_inland[c(9)]) + 
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"), 
          legend.title = element_blank(), 
          legend.direction = 'vertical') + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(15, 49, 5))) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) + 
    guides(fill = 'none')
  ggsave(paste0(outdir, '-data-incidence_rate_relative_round_sex2.png'), w = 5, h = 5)
  ggsave(paste0(outdir, '-data-incidence_rate_relative_round_sex2.pdf'), w = 3.5, h = 4.5)
  
  #
  # incidence rate relative to first round by 3 age groups
  
  # age groups
  age_groups <- c('15-24', '25-34', '35-49')
  df_age_group <- data.table(AGEYRS = 15:49)
  df_age_group[, index_age_group := 3]
  df_age_group[AGEYRS < 35, index_age_group := 2]
  df_age_group[AGEYRS < 25, index_age_group := 1]
  df_age_group[, age_group := age_groups[index_age_group]]
  df_age_group[, AGE_GROUP_LABEL := paste0('Age: ', age_group)]
  
  icrg <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
  icrg[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
  icrg <- merge(icrg, df_age_group, by = 'AGEYRS')
  icrg <- icrg[, list(INCIDENT_RATE_SUSCEPTIBLE = sum(INCIDENT_CASES) / sum(SUSCEPTIBLE)), by = c('SEX', 'ROUND', 'COMM', 'age_group', 'REF.ROUND', 'iterations')] 
  icrg[, INCIDENCE_REL := INCIDENT_RATE_SUSCEPTIBLE / INCIDENT_RATE_SUSCEPTIBLE[ROUND == REF.ROUND], by = c('COMM', 'age_group', 'SEX', 'iterations')]
  tmp1 = icrg[, list(q= quantile(INCIDENCE_REL, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'age_group', 'SEX')]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  tmp1 <- merge(tmp1, df_round, by = c('COMM', 'ROUND'))
  tmp1[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  tmp1[, SEX_LABEL := 'Women']
  tmp1[SEX== 'M', SEX_LABEL := 'Men']  
  tmp1[, MIDPOINT_DATE := MIN_SAMPLE_DATE + (MAX_SAMPLE_DATE - MIN_SAMPLE_DATE)/2]

  ggplot(tmp1[COMM == 'inland'], aes(x = MIDPOINT_DATE, group = interaction(age_group, SEX_LABEL))) + 
    geom_hline(yintercept = 1, linetype = 'dashed', alpha= 0.5) + 
    geom_line(aes(y = M, linetype = SEX_LABEL, group = interaction(SEX_LABEL, age_group)),position=position_dodge(width = 300), col = 'grey50') + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = interaction(SEX_LABEL, age_group)), width = 300, position=position_dodge(width = 300), col = 'grey30') + 
    geom_point(aes(y = M, col  = age_group, shape = SEX_LABEL, group = interaction(SEX_LABEL, age_group)),  
               size = 2, position=position_dodge(width = 300)) + 
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
          legend.direction = 'vertical',
          legend.position = 'bottom') + 
    # scale_y_continuous( limits = c(NA,  3), expand = expansion(mult = c(0.02, 0.02))) + 
    labs(y = 'Incidence rate relative to round 10', col= 'Age', shape= 'Gender', 
         x = 'Date (midpoint of survey interval)', linetype = 'Gender') + 
    scale_y_log10()
  ggsave(paste0(outdir, '-data-incidence_rate_relative_round_age_group_sex.pdf'), w = 4.2,h = 5.1)
  
  
  #
  # incidence rate relative to first round female to male ratio by 3 age groups
  
  # age groups
  age_groups <- c('15-24', '25-34', '35-49')
  df_age_group <- data.table(AGEYRS = 15:49)
  df_age_group[, index_age_group := 3]
  df_age_group[AGEYRS < 35, index_age_group := 2]
  df_age_group[AGEYRS < 25, index_age_group := 1]
  df_age_group[, age_group := age_groups[index_age_group]]
  df_age_group[, AGE_GROUP_LABEL := paste0('Age: ', age_group)]
  
  # incidence rate ration male to female
  icrr <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
  icrr[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
  icrr <- merge(icrr, df_age_group, by = 'AGEYRS')
  icrr <- icrr[, list(INCIDENT_RATE_SUSCEPTIBLE = sum(INCIDENT_CASES) / sum(SUSCEPTIBLE)), by = c('SEX', 'ROUND', 'COMM', 'age_group', 'REF.ROUND', 'iterations')] 
  icrr[, INCIDENCE_REL := INCIDENT_RATE_SUSCEPTIBLE / INCIDENT_RATE_SUSCEPTIBLE[ROUND == REF.ROUND], by = c('COMM', 'age_group', 'SEX', 'iterations')]
  icrrs <- dcast.data.table(icrr, COMM + ROUND + age_group + iterations ~ SEX, value.var = 'INCIDENCE_REL')
  icrrs[, INCIDENCE_REL_RATIO := `F` / `M`]
  icrrs = icrrs[, list(q= quantile(INCIDENCE_REL_RATIO, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'age_group')]	
  icrrs = dcast(icrrs, ... ~ q_label, value.var = "q")
  icrrs <- merge(icrrs, df_community, by = 'COMM')
  icrrs <- merge(icrrs, df_round, by = c('COMM', 'ROUND'))
  icrrs[, MIDPOINT_DATE := MIN_SAMPLE_DATE + (MAX_SAMPLE_DATE - MIN_SAMPLE_DATE)/2]
  # icrrs <- icrrs[!ROUND %in% c('R010', 'R011')]
  
  ggplot(icrrs[COMM == 'inland'], aes(x = MIDPOINT_DATE, group = age_group)) + 
    geom_hline(yintercept = 1, linetype = 'dashed', alpha= 0.5) + 
    geom_line(aes(y = M),position=position_dodge(width = 300), col = 'grey50') + 
    geom_errorbar(aes(ymin = CL, ymax = CU), width = 300, position=position_dodge(width = 300), col = 'grey30') + 
    geom_point(aes(y = M, col = age_group, shape = age_group),  size = 2, position=position_dodge(width = 300)) + 
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
          legend.position = 'bottom') + 
    scale_y_continuous( limits = c(NA,  3), expand = expansion(mult = c(0.02, 0.02))) + 
    labs(y = 'Female-male ratio of the incidence\nrate relative to round 10', col= 'Age', shape= 'Age', 
         x = 'Date (midpoint of survey interval)') 
  ggsave(paste0(outdir, '-data-incidence_rate_ratio_relative_round_sex.pdf'), w = 4.2,h = 4.2)
  
  #
  # incidence rate relative to first round female to male ratio total
  
  # incidence rate ration male to female
  icrrt <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
  icrrt[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
  icrrt <- icrrt[, list(INCIDENT_RATE_SUSCEPTIBLE = sum(INCIDENT_CASES) / sum(SUSCEPTIBLE)), by = c('SEX', 'ROUND', 'COMM', 'REF.ROUND', 'iterations')] 
  icrrt[, INCIDENCE_REL := INCIDENT_RATE_SUSCEPTIBLE / INCIDENT_RATE_SUSCEPTIBLE[ROUND == REF.ROUND], by = c('COMM', 'SEX', 'iterations')]
  icrrt <- dcast.data.table(icrrt, COMM + ROUND  + iterations ~ SEX, value.var = 'INCIDENCE_REL')
  icrrt[, INCIDENCE_REL_RATIO := `F` / `M`]
  icrrt = icrrt[, list(q= quantile(INCIDENCE_REL_RATIO, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND')]	
  icrrt = dcast(icrrt, ... ~ q_label, value.var = "q")
  icrrt <- merge(icrrt, df_community, by = 'COMM')
  icrrt <- merge(icrrt, df_round, by = c('COMM', 'ROUND'))
  
  
  #
  # save statistics
  
  save_statistics_incidence_rate_trends(icrr, icr, icrrs, icrrt, medage, fmr)
    
}


plot_incident_cases_to_unsuppressed_rate_ratio <- function(incidence_cases_round, unsuppressed_rate_ratio , outdir, outdir.table){
  
  # timeline
  df_timeline <- copy(df_round)
  df_timeline[, MIDPOINT := as.Date(mean(c(MIN_SAMPLE_DATE, MAX_SAMPLE_DATE))), by = c('ROUND', 'COMM')]
  df_timeline <- df_timeline[, .(ROUND, MIDPOINT, COMM, INDEX_ROUND, ROUND_SPANYRS)]
  
  # age groups
  age_groups <- c('15-24', '25-34', '35-49')
  df_age_group <- data.table(AGEYRS = 15:49)
  df_age_group[, index_age_group := 3]
  df_age_group[AGEYRS < 35, index_age_group := 2]
  df_age_group[AGEYRS < 25, index_age_group := 1]
  df_age_group[, age_group := age_groups[index_age_group]]
  df_age_group[, AGE_GROUP_LABEL := paste0('Age: ', age_group)]
  
  # Prepare incidence rates
  icr <- merge(incidence_cases_round, df_age_group, by = 'AGEYRS')
  icr <- icr[, list(INCIDENT_RATE_SUSCEPTIBLE = sum(INCIDENT_CASES) / sum(SUSCEPTIBLE), 
                    INCIDENT_RATE_UB_SUSCEPTIBLE = sum(INCIDENT_CASES_UB)  / sum(SUSCEPTIBLE),
                    INCIDENT_RATE_LB_SUSCEPTIBLE = sum(INCIDENT_CASES_LB) / sum(SUSCEPTIBLE)), by = c('SEX', 'ROUND', 'COMM', 'age_group')] 
  
  # merge labels
  icr <- merge(icr, df_timeline, by = c('ROUND', 'COMM'))
  
  # reference round
  icr[COMM == 'fishing', REFERENCE_ROUND := 'R015']
  icr[COMM == 'inland', REFERENCE_ROUND := 'R010']
  
  # take decline since first round
  icr[, INCIDENT_RATE_REF := INCIDENT_RATE_SUSCEPTIBLE / INCIDENT_RATE_SUSCEPTIBLE[ROUND == REFERENCE_ROUND], by = 'COMM']
  # icr <- icr[!(COMM == 'inland' & ROUND %in% c('R010', 'R011'))]
  
  # find incident cases ratio
  icr <- dcast(icr, ROUND + INDEX_ROUND + COMM + age_group ~ SEX, value.var = 'INCIDENT_RATE_REF')
  setnames(icr, c("M", 'F'), c('INCIDENT_RATE_REF_M', 'INCIDENT_RATE_REF_F'))
  icr[, INCIDENT_RATE_RATIO_REF := INCIDENT_RATE_REF_F / INCIDENT_RATE_REF_M]


  #
  # merge to unsuppressed rate ratio
  # urr <- unique(unsuppressed_rate_ratio[, .(ROUND, COMM, UNSUPPRESSION_RATE_RATIO_BY_AGE_M, UNSUPPRESSION_RATE_RATIO_BY_AGE_CL, UNSUPPRESSION_RATE_RATIO_BY_AGE_CU, AGE_GROUP)])
  urr <- unique(unsuppressed_rate_ratio[, .(ROUND, COMM, UNSUPPRESSION_RATE_RATIO_RATIO_M, UNSUPPRESSION_RATE_RATIO_RATIO_CL, UNSUPPRESSION_RATE_RATIO_RATIO_CU)])
  
  urr[, ROUND := paste0('R0', ROUND)]
  ic <- merge(icr, urr, by = c('ROUND', 'COMM'), allow.cartesian=TRUE)
  ic[, ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1',  ROUND))]
  
  #
  # Plot
  
  communities <- df_round[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    comm <- communities[i]
    tmp <- ic[COMM== comm]

    if(communities[i] == 'inland'){
      colors <- palette_round_inland[-1]
    }else{
      colors <- palette_round_fishing[-1]
    }
    
    p<-ggplot(tmp[INDEX_ROUND != min(INDEX_ROUND)]) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', alpha = 0.5) + 
      # geom_hline(yintercept = 1, linetype = 'dashed', alpha = 0.5) + 
      # geom_vline(xintercept = 1, linetype = 'dashed', alpha = 0.5) + 
      # geom_point(aes(col = ROUND_LABEL, shape = AGE_GROUP, x = UNSUPPRESSION_RATE_RATIO_BY_AGE_M, y = INCIDENT_RATE_RATIO_REF), size = 2) + 
      geom_point(aes(col = ROUND_LABEL, shape = age_group, x = UNSUPPRESSION_RATE_RATIO_RATIO_M, y = INCIDENT_RATE_RATIO_REF), size = 2) + 
      scale_y_continuous(limits = c(NA,2)) + 
      scale_color_manual(values = colors) + 
      labs(x = paste0('Male-female ratio of the rate of individuals\nwith HIV who have unsuppressed virus\nrelative to round ', tmp[INDEX_ROUND== min(INDEX_ROUND), gsub('R0(.+)', '\\1',unique(ROUND))]),
           y = paste0('Female-male ratio of the\nincidence rate relative\nto round ', tmp[INDEX_ROUND== min(INDEX_ROUND), gsub('R0(.+)', '\\1',unique(ROUND))])) + 
      theme_bw()  +
      theme(legend.position = 'none', 
            legend.title = element_blank(), 
            legend.text=element_text(size=rel(0.5))) +
      guides(shape=guide_legend(override.aes = list(size = 1), 
                                keywidth = unit(0.01, 'pt'), keyheight = unit(0.001, 'pt')), 
             color = guide_legend(byrow = T, nrow = 4,override.aes = list(size = 1), 
                                                keywidth = unit(0.01, 'pt'), keyheight = unit(0.001, 'pt'))) 
    
      # scale_y_log10()
    ggsave(p, file =  paste0(outdir, '-data-incidence_cases_rate_ratio_unsuppressed_', communities[i], '.pdf'), w = 3.7, h = 2.4)
    
  }
  
  #
  # save statistics
  save_statistics_incidence_rate_ratio_trends(ic, outdir.table)
}

plot_pairs_all <- function(pairs.all, outdir, nm_reqs=FALSE){
  
  tmp <- pairs.all[BOTH_PARTICIPATED == TRUE ]
  tmp[, {cat(sum(SEX.RECIPIENT != SEX.SOURCE)); cat(.N); table(SEX.RECIPIENT, SEX.SOURCE)/.N}]

  # find direction label
  tmp[, DIRECTION := 'Male to Female' ]
  tmp[SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M', DIRECTION := 'Female to Male' ]
  tmp[SEX.SOURCE == 'F' & SEX.RECIPIENT == 'F', DIRECTION := 'Female to Female' ]
  tmp[SEX.SOURCE == 'M' & SEX.RECIPIENT == 'M', DIRECTION := 'Male to Male' ]
  tmp[, DIRECTION := factor(DIRECTION, levels = c('Male to Male', 
                                                  'Female to Female', 
                                                  'Female to Male', 
                                                  'Male to Female'))]
  
  # find count and percentage
  tmp <- tmp[, list(COUNT = .N), by = c('DIRECTION')]
  tmp[, TOTAL_COUNT := sum(COUNT)]
  tmp[, PROPORTION := paste0(round(COUNT / TOTAL_COUNT*100, 1), '%')]
  
  # look at confidence interval
  # tmp1 <- tmp[, {
  #   intervals=Hmisc::binconf(COUNT, TOTAL_COUNT)
  #   list(M = intervals[1], CL = intervals[2], CU = intervals[3])}, by = c('COMM.RECIPIENT', 'DIRECTION')]
  # tmp1 <- tmp1[order(COMM.RECIPIENT ,DIRECTION)]
  # 
  # plot
  male_to_female_color <- 'lightblue3'
  female_to_male_color <- 'lightpink2'
  male_to_male_color <- 'grey70'
  female_to_female_color <- 'grey50'
  
  # tmp <- tmp[COMM.RECIPIENT == 'inland']
  labsize <- fifelse(nm_reqs, yes=3, no=4)
  p <- ggplot(tmp, aes(x = DIRECTION, y = COUNT, fill=DIRECTION, label=PROPORTION))+
    geom_col(width=0.6)+
    theme_bw() +
    geom_text(nudge_y= 7,color="black",size = labsize,fontface="bold") + 
    labs(y="Number of identified source-\nrecipient pairs in RCCS")+
    theme(legend.position='none', 
          axis.text.x = element_blank(), 
          axis.ticks.x= element_blank(), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), 
          axis.title.x = element_blank())+
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(values = c('Male to Female'=male_to_female_color,'Female to Male'=female_to_male_color, 
                                 'Female to Female'=female_to_female_color,'Male to Male'=male_to_male_color)) 
  ggsave(p, file = paste0(outdir, '-data-PairsAll.pdf'), w = 3.3, h = 2.5)

    if(nm_reqs){ p <- p + reqs }
    return(p)
  
}

plot_pairs <- function(pairs, outdir, nm_reqs=FALSE)
{
  
  # extend round periods to the beginning of next one.
  df_round_extended <- copy(df_round)
  df_round_extended[, MAX_SAMPLE_DATE := fcoalesce( shift(MIN_SAMPLE_DATE, -1), MAX_SAMPLE_DATE )]

    inf_date_range <- range(pairs$DATE_INFECTION.RECIPIENT)
    min_round_date <- min(df_round_extended$MIN_SAMPLE_DATE)
    if(inf_date_range[1] <= min_round_date & nm_reqs)
    {
        tmp <- data.table(
            ROUND='R000', 
            COMM='inland',
            MIN_SAMPLE_DATE = as.Date(-Inf),
            MAX_SAMPLE_DATE = min_round_date - 1,
            LABEL_ROUND = paste0('Before R10\n', 'Before ', format(min_round_date, '%b %Y'))
        )
        df_round_extended <- rbind( tmp,df_round_extended, fill=TRUE)
        palette_round_inland <- c( '#000000', palette_round_inland)
    }

    max_round_date <- max(df_round_extended$MAX_SAMPLE_DATE)
    if(inf_date_range[2] >= max_round_date & nm_reqs)
    {
        tmp <- data.table(
            ROUND='R111', 
            COMM='inland',
            MIN_SAMPLE_DATE = inf_date_range[2] + 1,
            MAX_SAMPLE_DATE = as.Date(Inf),
            LABEL_ROUND = paste0('After R18\n', 'After ', format(max_round_date, '%b %Y'))
        )
        df_round_extended <- rbind( df_round_extended, tmp, fill=TRUE)
        palette_round_inland <- c(palette_round_inland, '#FF0000')

    }

  # find round of infection
  tmp <- merge(pairs, df_round_extended, by.x = 'COMM.RECIPIENT', by.y = 'COMM', allow.cartesian = T)
  tmp <- tmp[DATE_INFECTION.RECIPIENT >= MIN_SAMPLE_DATE & DATE_INFECTION.RECIPIENT <= MAX_SAMPLE_DATE]
  
  # find direction label
  tmp[, DIRECTION := 'Male to female' ]
  tmp[SEX.SOURCE == 'F', DIRECTION := 'Female to male' ]
  
  COMMS <- pairs[, unique(COMM.RECIPIENT)]
  SEX <- c('M', 'F')
  for(i in seq_along(COMMS)){
    
    p <- list();index = 1
    for(j in seq_along(SEX)){
      
      comm <- COMMS[i]
      sex <- SEX[j]
      
      tmp1 <- tmp[COMM.RECIPIENT == comm & SEX.SOURCE == sex ]
      p[[index]] <- ggplot(tmp1, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
        geom_point(aes(col = LABEL_ROUND)) + 
        labs(y = 'Age at transmission source', x = 'Age at infection recipient', col = '') +
        geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
        theme_bw() + 
        coord_fixed() +
        ggtitle(tmp1[, unique(DIRECTION)]) + 
        scale_color_manual(values = palette_round_inland) + 
        scale_x_continuous(limits = c(15, 49))+
        scale_y_continuous(limits = c(15, 49)) +
        geom_label(x = 18, y = 49, label = paste0(paste0(nrow(tmp1), ' pairs')), 
                size=fifelse(nm_reqs, yes=3, no=NA_integer_),
                label.size = NA) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      
      if(j == 2){
        p[[index]] <- p[[index]] +
          theme(legend.position = 'none') 
      }else{
        p[[index]] <- p[[index]] +
          theme(legend.position = 'bottom', 
                legend.justification = "left") + 
          guides(color = guide_legend(byrow = T, nrow =2))
      }
      
      if(nm_reqs){ p[[index]] <- p[[index]] + reqs }

      p[[index]] <- ggExtra::ggMarginal(p[[index]], type = "histogram")
      
      index=index + 1
      
    }
    
    pp <- grid.arrange(grobs = p, layout_matrix = rbind(c(1,2), c(1,NA)), heights= c(0.83, 0.17))
    if(nm_reqs)
    {
        return(pp)
    }else{
        ggsave(pp, file = paste0(outdir, '-data-Pairs_', comm, '.pdf'), w = 8.2, h = 5)
    }
    
  }
}

plot_sources_histogram <- function(pairs, outdir)
{
  # prepare counts
  tmp <- as.data.table(reshape2::melt(stan_data[['y']]))
  setnames(tmp, 1:3, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_TIME'))
  tmp <- merge(tmp, df_direction, by = 'INDEX_DIRECTION')
  tmp <- merge(tmp, df_period, by = c('INDEX_TIME'))
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp <- merge(tmp, df_age, by = 'INDEX_AGE')
  setnames(tmp, 'value', 'count')

  # prepare percentage of sampling
  tmp1 <- as.data.table(reshape2::melt(exp(stan_data[['log_prop_sampling']])))
  setnames(tmp1, 1:3, c('INDEX_DIRECTION', 'INDEX_TIME', 'INDEX_AGE'))
  setnames(tmp1, 'value', 'prop')
  tmp <- merge(tmp, tmp1, by = c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_TIME'))
  
  # augment count with prop
  tmp[, count_augment := count / prop]
  
  # sum across recipients
  tmp <- tmp[, list(count_augment = sum(count_augment), 
                    count = sum(count)), by = c('LABEL_DIRECTION', 'COMM', 'PERIOD', 'AGE_TRANSMISSION.SOURCE')]
  
  # count density
  tmp[, count_augment_density := count_augment / sum(count_augment), by = 'PERIOD']
  tmp[, count_density := count / sum(count), by = 'PERIOD']
  
  # plot
  tmp1 <- tmp[count_density != 0.0000]
  ggplot(tmp1, aes(x = AGE_TRANSMISSION.SOURCE)) + 
    geom_line(aes(y = count_density, col = PERIOD)) + 
    geom_line(aes(y = count_augment_density, col = PERIOD), linetype = 'dashed') + 
    facet_grid(~LABEL_DIRECTION) + 
    theme_bw()


  # df <- copy(pairs)
  # df[, AGE_TRANSMISSION.SOURCE := floor(AGE_TRANSMISSION.SOURCE)]
  # df <- df[, list(count = .N), by = c('AGE_TRANSMISSION.SOURCE', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT', 'SEX.SOURCE')]
  # df[, count_density := count / sum(count), by = 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT']
  # 
  # ggplot(df, aes(x = AGE_TRANSMISSION.SOURCE)) + 
  #   geom_line(aes(y = count_density, col = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT)) + 
  #   facet_grid(~SEX.SOURCE) + 
  #   theme_bw()
  
  

} 

plot.coherent.tsi.estimates.with.seroconversion <- function(outdir)
{
        tmp <- meta_data[!is.na(date_infection)]
        cols <- grep('study_id|age|date', names(tmp), value=TRUE)
        tmp <- tmp[, ..cols]

        # Age collection does not strongly contradict age_first_positive, good
        tmp[ age_collection - age_first_positive < - 0.1]

        # Get visit date used to obtain estimate
        tmp <- merge(tmp, time.since.infection[, .(study_id, visit_dt)], by='study_id')
        cols <- c('visit_dt', 'date_infection')
        tmp[, (cols):=lapply(.SD,as.Date) , .SDcols=cols]
        tmp <- melt(tmp, measure.vars=c('visit_dt', 'date_infection'), value.name = "date")
        tmp[, variable:=fifelse(variable %like% 'visit', 'Collection date', 'Infection date')]

        idx <- tmp[!is.na(date_first_positive) & !is.na(date_last_negative) & !is.na(date), study_id]

        # study_id as.factor to set an order 
        tmp[, validity := fcase(
                           variable == 'Collection date', 'NA',
                           date <= date_first_positive & date >= date_last_negative, 'Coherent', 
                           rep(TRUE, .N), 'Incoherent'
        )] 
        
        tmp1 <- tmp[variable == 'Collection date', ] 
        setkey(tmp1, date)
        tmp1 <- tmp[, study_id := factor(study_id, levels=unique(tmp1$study_id), ordered=TRUE)]
        tmp1 <- tmp1[study_id %in% idx]
        .p <- function(x) paste0(round(x*100, 2), '%')
        lab <- tmp1[validity != 'NA', mean(validity=='Coherent')]
        lab <- paste0('\n\tCoherent predictions:\n\t',.p(lab))

        g <- ggplot(tmp1, aes(y=study_id, group=study_id)) +
                geom_segment(aes(x=date_last_negative, xend=date_first_positive, yend=study_id, linetype='seroconversion interval'), color='grey80') +
                geom_point(aes(x=date, pch=variable, color=validity)) +
                geom_text(aes(x=as.Date(-Inf), y=+Inf, label=lab), vjust=1, hjust=0) +
                scale_color_manual(values = c('green', 'red', 'blue'))+
                scale_x_date(breaks='6 months',expand=c(0,0), labels=scales:::date_format("%b %Y")) + 
                theme_bw() +
                theme(legend.position='bottom',
                      axis.text.x=element_text(angle=45, vjust=1, hjust=1),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank()) +
                labs(y='participants', x='', linetype='', pch='') 
        
        if(!is.null(outdir))
        {
                tmp <- '-TSI_all_seroconverters.png'
                if( file.path.tsiestimates %like% 'adjusted')
                        tmp <- gsub('TSI','TSI_adj',tmp)
                ggsave(g, filename = paste0(outdir, tmp), w = 10, h = 15)
        }

        # Coherence vs time difference between first positive and collection
        # get time difference
        tmp2 <- tmp1[, {
                z <- date[variable == 'Collection date']
                z1 <- z - date_first_positive 
                list(validity, deltaT=as.integer(z1))
                }, by='study_id'][validity != 'NA'] 
        tmp2 <- unique(tmp2)

        tmp2[, as.data.table(table(validity)), by=deltaT==0]

        # transform variable to wanted formats
        tmp2[, validity2 := fifelse(validity == 'Coherent', 1, 0)]

        # group together by time difference
        setkey(tmp2, deltaT)
        idx <- which(tmp2[, diff(deltaT) > 200])
        idx <- c(-1, tmp2$deltaT[idx], max(tmp2$deltaT))

        # compute binomial p estimates 
        tmp2[, deltaT2 := cut(tmp2$deltaT, idx) ]
        cols <- c('M', 'CL', 'CU')
        tmp2[, (cols) := 
             Hmisc::binconf(sum(validity2), .N, return.df=TRUE),
             by=deltaT2]
        tmp2[, midpoint2 := median(deltaT), by=deltaT2]

        # plot: Coherency drops with increasing time distance
        p <- ggplot(tmp2, aes(y=validity2, x=deltaT/365, col=deltaT2)) + 
                geom_jitter(pch='X', height=.01, width=0.02) +
                geom_point(data=tmp2[midpoint2 < 13 * 365],
                           aes(x=midpoint2/365, y=M), pch='triangle', size=3) +
                geom_linerange(data=tmp2[midpoint2 < 13 * 365],
                               aes(x=midpoint2/365, ymin=CL, ymax=CU)) +
                theme_bw() +
                theme(legend.position='none') +
                labs(
                     x= 'years between date of first positive test and sample collection (for TSI purposes)',
                     y='Proportion of TSI estimates in seroconversion interval',
                     title='Estimate Coherency drops as sample collection is delayed'
                )

        if(!is.null(outdir))
        {
                tmp <- '-TSI_vs_difftime_since_collection.png'
                if( file.path.tsiestimates %like% 'adjusted')
                        tmp <- gsub('TSI','TSI_adj',tmp)
                ggsave(p, filename = paste0(outdir, tmp), w = 10, h = 10)
        }

}

plot.tsi.relationships.among.source.recipient.pairs <- function(outdir)
{
        # Helpers #
        # _________

        .plot.pairs <- function(DT)
        {
                lab <- DT[, round(mean(DATE_INFECTION.SOURCE < DATE_INFECTION.RECIPIENT) * 100, 2)]
                lab <- data.table(lab = paste0('\n\tCoherent dates of infection:\n\t', lab, '%'))

                p <- ggplot(DT, aes(x=DATE_INFECTION.SOURCE, y=DATE_INFECTION.RECIPIENT)) +
                        geom_abline(slope=1, color='red', linetype='dashed') +
                        geom_point() +
                        geom_text(data=lab, aes(label=lab, x=as.Date(-Inf), y=as.Date(Inf), hjust=0, vjust=1)) +
                        scale_x_date(breaks='6 months',expand=c(0,0), labels=scales:::date_format("%b %Y")) + 
                        scale_y_date(breaks='6 months',expand=c(0,0), labels=scales:::date_format("%b %Y")) + 
                        theme_bw() +
                        theme(legend.position='bottom',
                              axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
                        labs(x='Source',
                             y='Recipient',
                             linetype='', pch='',
                             title='Transmission pairs: estimated dates of infection') 
        }

        adjust.filename <- function(x)
        {
                if( file.path.tsiestimates %like% 'adjusted')
                        x <- gsub('tsi','tsi-adj', x )
        }


        # Start
        # _____

        cols <- grep('AID|DATE_INFECTION|PRED_DOI|BIRTH',names(pairs.all), value=TRUE)
        pairs.tmp <- pairs.all[, ..cols] 
        # tmp[, uniqueN(.SD) == .N]

        # plot
        p <- .plot.pairs(pairs.tmp)
        file <- paste0(outdir, '-tsi-source_vs_recipient_from_HIVphyloTSI.png')
        file <- adjust.filename(file)
        ggsave(p, file = file, w = 10, h = 10)

        # Add 'rectangles'
        # ________________

        by_cols <- grep('^AID', names(pairs.tmp), value=TRUE)

        tmp <- pairs.tmp[,  find.center.of.mass.plausible.region(
                                xmin=PRED_DOI_MIN.SOURCE,
                                xmax=PRED_DOI_MAX.SOURCE,
                                ymin=PRED_DOI_MIN.RECIPIENT,
                                ymax=PRED_DOI_MAX.RECIPIENT                                                            ) , by=by_cols]

        stopifnot(tmp[, all(x < y, na.rm=TRUE)])

        cat('Adjust transmission pair-incoherent DOI\n')
        pairs.tmp <- merge(pairs.tmp, tmp, by=by_cols)

        idx <- pairs.tmp[, !is.na(x) & !is.na(y) & DATE_INFECTION.RECIPIENT <= DATE_INFECTION.SOURCE]
        pairs.tmp[idx, DATE_INFECTION.SOURCE := x]
        pairs.tmp[idx, DATE_INFECTION.RECIPIENT := y]
        pairs.tmp[, `:=` (x=NULL, y=NULL)]

        # This seems too optimistic
        g <- .plot.pairs(pairs.tmp)
        file = paste0(outdir, '-tsi-source_vs_recipient_from_HIVphyloTSI_centerofmassadj.png')
        file <- adjust.filename(file)
        ggsave(g, file = file, w = 10, h = 10)

        view.rectangle <- function(xmin2, xmax2, ymin2, ymax2)
        {
                dt <- data.table(xmin2=xmin2, xmax2=xmax2, ymin2=ymin2, ymax2=ymax2)
                ggplot(dt, aes(xmin=xmin2, xmax=xmax2, ymin=ymin2, ymax=ymax2)) +
                        geom_rect(fill=NA, color='black') +
                        geom_abline(slope=1, linetype='dashed', color='red')
        }

        # split DF by source and recipient and compute ages
        if(0)
        {
        
                .get.age <- function(x, birth) round(as.integer(x - birth)/365.25, 1)

                .get.CI.ages <- function(type)
                {
                        cols <- grep(type, names(tmp), value=TRUE)
                        tmp1 <- tmp[, ..cols]
                        names(tmp1) <- gsub(paste0('.', type), '', names(tmp1) )
                        cols <- grep('PRED|INFECTION', names(tmp1), value=TRUE)
                        newcols <- gsub('PRED|DATE', 'AGE', cols) 
                        tmp1[, (newcols) := lapply(.SD, .get.age , birth=DATE_BIRTH) , .SDcols=cols]

                        cols <- grep('AGE|AID', names(tmp1), value=TRUE)
                        newcols <- paste0(cols, '.', type)
                        tmp1[, ..cols]
                }
                
                tmp1 <- lapply(c('SOURCE', 'RECIPIENT'), .get.CI.ages)
                tmp1 <- rbindlist(tmp1)

                # Merge back to pairs.all data
                tmp2 <- copy(tmp1)
                setnames(tmp1, names(tmp1), paste0(names(tmp1), '.SOURCE'))
                setnames(tmp2, names(tmp2), paste0(names(tmp2), '.RECIPIENT'))
                tmp <- merge(tmp, tmp1, by='AID.SOURCE')
                tmp <- merge(tmp, tmp2, by='AID.RECIPIENT')

                # rename
                names(tmp) <- gsub('AGE_DOI', 'AGE_INFECTION', names(tmp))

        }


        # However, is this partially due to very unlikely transmission pairs? 
        .p <- function(x) round(100*x, 2)
        .f <- function(x) paste0( sum(x), "/", length(x), " (",  .p(mean(x)), "%)")

        cols <- grep('AID|DATE_FIRST_POSITIVE.RECIPIENT|DATE_LAST_NEGATIVE.SOURCE',
                     names(pairs.all), value=TRUE)
        tmp <- pairs.all[, ..cols]
        lab0 <- tmp[, paste0(.f(!is.na(DATE_LAST_NEGATIVE.SOURCE)), 'of date-last-negative for sources are defined\n')] 
        lab1 <- tmp[ !is.na(DATE_LAST_NEGATIVE.SOURCE),
            paste0('Of these:\n',
                .f(DATE_FIRST_POSITIVE.RECIPIENT < DATE_LAST_NEGATIVE.SOURCE),
                ', are incoherent in the sense that the source supposedely got infected after the recipient.\n')]
        cat(lab0)
        lab <- data.table(lab=paste0(lab0, gsub('source', 'source\n',lab1)))
        .f <- function(x) as.Date(x)
        p <- ggplot(tmp, aes(x=DATE_LAST_NEGATIVE.SOURCE, y=DATE_FIRST_POSITIVE.RECIPIENT)) +
                geom_abline(slope=1, color='red', linetype='dashed') +
                geom_point() +
                geom_text(data=lab, aes(label=lab, x=.f('2003-01-01'), y=.f(-Inf), hjust='left', vjust=0)) +
                scale_x_date(breaks='6 months',expand=c(0,0), labels=scales:::date_format("%b %Y")) + 
                scale_y_date(breaks='6 months',expand=c(0,0), labels=scales:::date_format("%b %Y")) + 
                theme_bw() +
                theme(legend.position='bottom',
                      axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
                labs(x='Source: date of last negative test',
                     y='Recipient: date of first postive test',
                     linetype='', pch='',
                     title='Transmission pairs: seroconversion date from interview data') 
        p
        file = paste0(outdir, '-tsi-source_vs_recipient_from_interviewdata.png')
        if( file.path.tsiestimates %like% 'adjusted')
                file <- gsub('tsi','tsi-adj',file)

        ggsave(p, file = file, w = 8, h = 10)
}

plot.phylopair.dates.scores <- function(DT, 
                                        add.rects=TRUE, add.dots=TRUE,
                                        only.rect.pairs=TRUE, only.contradict=FALSE, only.crossing=FALSE, only.coherent=FALSE, title=NULL, 
                                        daterange.vars=c('date_last_negative','date_first_positive'),
                                        doi.center.var='date_infection',
                                        return.DT=FALSE
)
{
        # DT <- copy(chain2)
        rgx <- paste0(daterange.vars, collapse='|')

        if(only.rect.pairs)
        {
                cols <- grep(rgx, names(DT), value=TRUE)
                tmp <- DT[, lapply(.SD, function(x) !is.na(x) ) , .SDcols=cols]
                idx <- which(apply(tmp, 1, all))
                DT <- DT[idx,]
                rm(idx)
        }

        # rename var
        vxmin <- paste0( daterange.vars[1], '.SOURCE')
        vxmax <- paste0( daterange.vars[2], '.SOURCE')
        vymin <- paste0( daterange.vars[1], '.RECIPIENT')
        vymax <- paste0( daterange.vars[2], '.RECIPIENT')

        idx <- 1:nrow(DT)
        if(only.contradict)
        {
                # the source cannot have been infected following the recipient
                idx <- DT[, ..vxmin][[1]] >= DT[, ..vymax][[1]]
        }else if(only.crossing){
                # 
                idx.1 <- DT[, ..vymax][[1]] >= DT[, ..vxmin][[1]]
                idx.2 <- DT[, ..vxmax][[1]] >= DT[, ..vymin][[1]]
                idx <- idx.1 & idx.2
        }else if(only.coherent){
                # 
                idx <- DT[, ..vxmax][[1]] <= DT[, ..vymin][[1]]
        }
        DT <- DT[idx, ]

        if(return.DT)
                return(DT)

        g <- ggplot(DT, aes(col=SCORE_DIR_SR, fill=SCORE_DIR_SR)) +
                theme_bw() + 
                #                 scale_fill_viridis_c( begin=0.5, end=1 )  + 
                #                 scale_color_viridis_c( begin=0.5, end=1 )  + 
                scale_fill_viridis_c( limits=c(.5, 1)) +
                scale_color_viridis_c( limits=c(.5, 1)) +
                theme(legend.position='bottom') +
                labs(x='date of infection(source)', y='date of infection(recipient)')

        if(add.rects)
        {
                g <- g + geom_rect(aes_string(ymin=vymin,
                                              ymax=vymax,
                                              xmin=vxmin,
                                              xmax=vxmax), 
                                   alpha=.2, color=NA)
        }
        
        if(add.dots)
        {
                vars <- paste0(doi.center.var, c('.SOURCE', '.RECIPIENT'))
                g <- g + 
                        geom_point(aes_string(x=vars[1],
                                              y=vars[2])) 
        }

        if(!is.null(title))
                g <- g + labs(title=title)

        g + 
                geom_abline(slope=1,linetype='dotted', color='red') 
}
