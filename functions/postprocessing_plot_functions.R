getLevel <- function(x,y,z,prob=0.95) {
  dx <- diff(unique(x)[1:2])
  dy <- diff(unique(y)[1:2])
  sz <- sort(z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

plot_intensity_PP <- function(intensity_PP, count_data, outdir){

  communities <- intensity_PP[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- intensity_PP[ COMM == communities[i]]
    tmp1 <- count_data[ COMM == communities[i]]
    
    p <- ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      geom_point(data = tmp1[count > 0], aes(size = count), col = 'grey50') +
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = 'Estimated median rate\nobserved transmission\nevents', 
           y= 'Age at transmission source',size='Pairs\ncount') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_size_continuous(range = c(1, 3), breaks = sort(unique(tmp1[count > 0]$count))) +
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) + 
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    
    ggsave(p, file = paste0(outdir, '-output-intensity_transmission_',  communities[i], '.png'), w = 7, h = 7)
  }

} 

plot_2D_contrast <- function(tmp, outdir, lab = NULL, name = NULL){
  

    p <- ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = lab, 
           y= 'Age at transmission source',size='Pairs\ncount') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) 
    
    if('LABEL_DIRECTION' %in% names(tmp)){
      p <- p + facet_grid(LABEL_DIRECTION~.)
    }
    
    ggsave(p, file = paste0(outdir, '-output-contrast_2D_', name, '.png'), w = 7, h = 7)
  
} 

plot_source_contrast <- function(tmp, outdir, lab = NULL, name = NULL){
  
  p <- ggplot(tmp, aes(x = AGE_TRANSMISSION.SOURCE)) + 
    geom_line(aes(y = M)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
    theme_bw() + 
    labs( y = lab, 
         x= 'Age at transmission source') +
    # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_fill_viridis_c() 
  
  if('LABEL_DIRECTION' %in% names(tmp)){
    p <- p + facet_grid(LABEL_DIRECTION~.)
  }
  
  ggsave(p, file = paste0(outdir, '-output-contrast_source_', name, '.png'), w = 7, h = 7)
  
} 

plot_recipient_contrast <- function(tmp, outdir, lab = NULL, name = NULL){
  
  p <- ggplot(tmp, aes(x = AGE_INFECTION.RECIPIENT)) + 
    geom_line(aes(y = M)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
    theme_bw() + 
    labs( y = lab, 
          x= 'Age at infection recipient') +
    # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_fill_viridis_c() 
  
  if('LABEL_DIRECTION' %in% names(tmp)){
    p <- p + facet_grid(LABEL_DIRECTION~.)
  }
  
  ggsave(p, file = paste0(outdir, '-output-contrast_recipient_', name, '.png'), w = 7, h = 7)
  
} 

plot_force_infection <- function(force_infection, outdir, lab = NULL){
  
  communities <- force_infection[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- force_infection[ COMM == communities[i]]

    p <- ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = paste0('Estimated median\n', lab, ' force infection'), 
           y= 'Age at transmission source',size='Pairs\ncount') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
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
    
    ggsave(p, file = paste0(outdir, '-output-', lab, 'force_infection_2D_',  communities[i], '.png'), w = 7, h = 7)
  }
  
} 

plot_force_infection_age_source <- function(force_infection_age_source, outdirL){
  
  tmp <- copy(force_infection_age_source)
  
  tmp[, Age := AGE_TRANSMISSION.SOURCE]
  
  tmp[, Direction :=LABEL_DIRECTION]

  
  p <- ggplot(tmp, aes(x = Age)) + 
    geom_bar(aes(y = M, fill = Direction), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = Direction), position = position_dodge()) + 
    labs(x = 'Age', y = 'Force of infection exerted', fill = '') + 
    theme_bw() +
    facet_grid(PERIOD~LABEL_COMMUNITY, scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age.png'), w = 12, h = 9)
  
}

plot_force_infection_sex_source <- function(force_infection_sex_source, outdir){
  
  tmp <- copy(force_infection_sex_source)
  
  tmp[, Direction :=LABEL_DIRECTION]
  
  p <- ggplot(tmp, aes(x = PERIOD)) + 
    geom_bar(aes(y = M, fill = Direction), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = Direction), position = position_dodge()) + 
    labs(x = '', y = 'Force of infection exerted', fill = '') + 
    theme_bw() +
    facet_grid(LABEL_COMMUNITY~., scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_sex.png'), w = 6, h = 7)
  
}


plot_force_infection_age_group <- function(force_infection_aggregated_age_group, outdir){
  
  communities <- force_infection_aggregated_age_group[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    tmp <- force_infection_aggregated_age_group[ COMM == communities[i]]
    tmp[, `Direction` := LABEL_DIRECTION]
    tmp[, `Age Recipient` := AGE_GROUP_INFECTION.RECIPIENT]
    
    p <- ggplot(tmp, aes(x = AGE_GROUP_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, fill = PERIOD), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = PERIOD), position = position_dodge()) + 
      labs(x = 'Age source', y = 'Force of infection exerted', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg() 
    ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_group_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}

plot_force_infection_age_classification <- function(force_infection_aggregated_age_classification, outdir){
  
  communities <- force_infection_aggregated_age_classification[, unique(COMM)]
  force_infection_aggregated_age_classification[, AGE_CLASSIFICATION.SOURCE := factor(AGE_CLASSIFICATION.SOURCE, levels  = c('Younger', 'Same age', 'Older'))]
  
  for(i in seq_along(communities)){
    
    tmp <- force_infection_aggregated_age_classification[ COMM == communities[i]]
    tmp[, `Direction` := LABEL_DIRECTION]
    tmp[, `Age Recipient` := AGE_GROUP_INFECTION.RECIPIENT]
    
    p <- ggplot(tmp, aes(x = AGE_CLASSIFICATION.SOURCE)) + 
      geom_bar(aes(y = M, fill = PERIOD), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = PERIOD), position = position_dodge()) + 
      labs(x = 'Age source', y = 'Force of infection exerted', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg() 
    ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_classification_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}

plot_contribution_sex_source <- function(contribution_sex_source, eligible_prop, outdir, lab = NULL){
  
  tmp <- copy(contribution_sex_source)
  tmp[, type  := 'Contribution to HIV infection']
  tmp <- rbind(tmp, eligible_prop, fill=TRUE)
  tmp[, type := factor(type , levels = c('Contribution to HIV infection','Share in the census eligible individuals'))]
  
  if(is.null(lab)) lab =  'Contribution to infection'
  
 p <- ggplot(tmp, aes(x = LABEL_DIRECTION)) + 
    geom_bar(aes(y = M, alpha = type, col = type), stat = 'identity', position = "identity", fill = 'cornflowerblue') + 
    geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.2, col = 'grey40') + 
    labs(x = '', y = 'Percent of census eligible individuals', 
         fill = lab, col = '', alpha = '') + 
    theme_bw() +
    facet_grid(LABEL_COMMUNITY~PERIOD, scale = 'free')+
    scale_color_manual(values = c('cornflowerblue', 'black')) + 
    scale_alpha_manual(values = c(1,0)) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') 
  ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_sex.png'), w = 6, h = 7)
  
}

plot_contribution_age_source <- function(contribution_age_source, eligible_prop, outdir, lab = NULL){
  
  tmp <- copy(contribution_age_source)
  tmp[, type  := 'Contribution to HIV infection']
  
  tmp1 <- copy(eligible_prop)
  setnames(tmp1, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
  tmp <- rbind(tmp, tmp1, fill=TRUE)
  tmp[, type := factor(type , levels = c('Contribution to HIV infection','Share in the census eligible individuals'))]
  
  communities <- tmp[, unique(COMM)]
  
  if(is.null(lab)) lab =  'Contribution to infection'
  
  for(i in seq_along(communities)){
    
    tmp1 <- tmp[COMM == communities[i]]
    
    p <- ggplot(tmp1, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, alpha = type, col = type), stat = 'identity', position = "identity", fill = 'cornflowerblue') + 
      geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.5, col = 'grey40') + 
      labs(x = 'Age', y = 'Percent of census eligible individuals', 
           fill = lab, col = '', alpha = '') + 
      theme_bw() +
      facet_grid(LABEL_DIRECTION~PERIOD, scale = 'free')+
      scale_color_manual(values = c('cornflowerblue', 'black')) + 
      scale_alpha_manual(values = c(1,0)) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') 
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_', communities[i], '.png'), w = 9, h = 7)
    
  }
}

plot_contribution_age_group <- function(contribution_age_group_source, outdir, lab = NULL){
  
  communities <- contribution_age_group_source[, unique(COMM)]
  
  if(is.null(lab)) lab =  'Contribution to infection'
  
  for(i in seq_along(communities)){
    
    tmp <- contribution_age_group_source[ COMM == communities[i]]
    tmp[, `Direction` := LABEL_DIRECTION]
    tmp[, `Age Recipient` := AGE_GROUP_INFECTION.RECIPIENT]
    
    p <- ggplot(tmp, aes(x = AGE_GROUP_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, fill = PERIOD), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = PERIOD), position = position_dodge()) + 
      labs(x = 'Age source', y = lab, fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg()  +
      ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_group_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}

plot_contribution_age_classification <- function(contribution_age_classification_source, outdir, lab = NULL){
  
  communities <- contribution_age_classification_source[, unique(COMM)]
  contribution_age_classification_source[, AGE_CLASSIFICATION.SOURCE := factor(AGE_CLASSIFICATION.SOURCE, levels  = c('Younger', 'Same age', 'Older'))]
  
  if(is.null(lab)) lab =  'Contribution to infection'
  
  for(i in seq_along(communities)){
    
    tmp <- contribution_age_classification_source[ COMM == communities[i]]
    tmp[, `Direction` := LABEL_DIRECTION]
    tmp[, `Age Recipient` := AGE_GROUP_INFECTION.RECIPIENT]
    
    p <- ggplot(tmp, aes(x = AGE_CLASSIFICATION.SOURCE)) + 
      geom_bar(aes(y = M, fill = PERIOD), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = PERIOD), position = position_dodge()) + 
      labs(x = 'Age source', y = lab, fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg()  +
      ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_classification_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}

plot_contribution_sex_source_by_round <- function(expected_contribution_sex_source_round, unsuppressed_prop_sex, outdir, lab = NULL){
  
  tmp <- copy(expected_contribution_sex_source_round)
  tmp[, type  := 'Contribution to HIV infection']
  tmp[, SEX := substr(LABEL_DIRECTION, 1, 1) ]
  tmp <- rbind(tmp, unsuppressed_prop_sex, fill=TRUE)
  tmp[, type := factor(type , levels = c('Contribution to HIV infection','Share in the unsuppressed HIV+ census eligible individuals'))]

  if(is.null(lab)) lab =  'Contribution to infection'
  
  p <- ggplot(tmp, aes(x = ROUND)) + 
    geom_bar(aes(y = M, alpha = type, col = type), stat = 'identity', position = "identity", fill = 'cornflowerblue') + 
    geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.2, col = 'grey40') + 
    labs(x = '', y = 'Percent', 
         fill = lab, col = '', alpha = '') + 
    theme_bw() +
    facet_grid(LABEL_COMMUNITY~SEX, scale = 'free')+
    scale_color_manual(values = c('cornflowerblue', 'black')) + 
    scale_alpha_manual(values = c(1,0)) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') 
  ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_sex_byround.png'), w = 6, h = 7)
  
}

plot_contribution_age_source_by_round <- function(contribution_age_source, unsuppressed_prop_age, outdir, lab = NULL){
  
  tmp <- copy(contribution_age_source)
  tmp[, type  := 'Contribution to HIV infection']
  tmp[, SEX := substr(LABEL_DIRECTION, 1, 1) ]
  
  tmp1 <- copy(unsuppressed_prop_age)
  setnames(tmp1, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
  tmp <- rbind(tmp, tmp1, fill=TRUE)
  tmp[, type := factor(type , levels = c('Contribution to HIV infection','Share in the unsuppressed HIV+ census eligible individuals'))]
  
  communities <- tmp[, unique(COMM)]
  
  if(is.null(lab)) lab =  'Contribution to infection'
  
  for(i in seq_along(communities)){
    
    tmp1 <- tmp[COMM == communities[i]]
    
    p <- ggplot(tmp1, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, alpha = type, col = type), stat = 'identity', position = "identity", fill = 'cornflowerblue') + 
      geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.5, col = 'grey40') + 
      labs(x = 'Age', y = 'Percent', 
           fill = lab, col = '', alpha = '') + 
      theme_bw() +
      facet_grid(ROUND~SEX, scale = 'free')+
      scale_color_manual(values = c('cornflowerblue', 'black')) + 
      scale_alpha_manual(values = c(1,0)) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') 
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_by_round_', communities[i], '.png'), w = 8, h = 9)
    
  }
}

plot_contribution_age_group_by_round <- function(contribution_age_group_source, outdir, lab = NULL){
  
  communities <- contribution_age_group_source[, unique(COMM)]
  
  if(is.null(lab)) lab =  'Contribution to HIV infection'
  
  for(i in seq_along(communities)){
    
    tmp <- contribution_age_group_source[ COMM == communities[i]]
    tmp[, `Direction` := LABEL_DIRECTION]
    tmp[, `Age Recipient` := AGE_GROUP_INFECTION.RECIPIENT]
    
    p <- ggplot(tmp, aes(x = AGE_GROUP_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, fill = ROUND), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = ROUND), position = position_dodge()) + 
      labs(x = 'Age source', y = lab, fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg() +
      ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_group_by_round_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}

plot_contribution_age_classification_by_round <- function(contribution_age_classification_source, outdir, lab = NULL){
  
  communities <- contribution_age_classification_source[, unique(COMM)]
  contribution_age_classification_source[, AGE_CLASSIFICATION.SOURCE := factor(AGE_CLASSIFICATION.SOURCE, levels  = c('Younger', 'Same age', 'Older'))]
  
  if(is.null(lab)) lab =  'Contribution to infection'
  
  for(i in seq_along(communities)){
    
    tmp <- contribution_age_classification_source[ COMM == communities[i]]
    tmp[, `Direction` := LABEL_DIRECTION]
    tmp[, `Age Recipient` := AGE_GROUP_INFECTION.RECIPIENT]
    
    p <- ggplot(tmp, aes(x = AGE_CLASSIFICATION.SOURCE)) + 
      geom_bar(aes(y = M, fill = ROUND), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = ROUND), position = position_dodge()) + 
      labs(x = 'Age source', y = lab, fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg()  +
      ggtitle(contribution_age_classification_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_classification_by_round_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}

plot_transmission_risk_sex_source_by_round <- function(transmission_risk_sex_source_round, outdir){
  
  tmp <- copy(transmission_risk_sex_source_round)
  tmp[, SEX := substr(LABEL_DIRECTION, 1, 1) ]
  
  p <- ggplot(tmp, aes(x = ROUND)) + 
    geom_bar(aes(y = M, fill = SEX), stat = 'identity', position = "identity") + 
    geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.2, col = 'grey40') + 
    labs(x = '', y = 'Exerted transmission risk per year', 
         col = '', alpha = '', fill ='') + 
    theme_bw() +
    facet_grid(LABEL_COMMUNITY~SEX, scale = 'free')+
    scale_fill_manual(values = c('pink', 'blue')) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    guides(alpha = 'none')
  ggsave(p, file = paste0(outdir, '-output-', 'transmisison_risk', '_sex_byround.png'), w = 6, h = 7)
  
}

plot_transmission_risk_age_source_by_round <- function(transmission_risk_round, outdir){
  
  tmp <- copy(transmission_risk_round)
  tmp[, SEX := substr(LABEL_DIRECTION, 1, 1) ]
  
  communities <- transmission_risk_round[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    tmp1 <- tmp[ COMM == communities[i] & AGE_INFECTION.RECIPIENT %in% c(25, 35, 45)]
    tmp1[, `Age recipient` := AGE_INFECTION.RECIPIENT]

    p <- ggplot(tmp1, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_line(aes(y = M, col = SEX)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX), alpha = 0.5 ) + 
      labs(x = 'Age source', y = 'Exerted transmission risk per year', 
           col = '',  fill ='') + 
      theme_bw() +
      facet_grid(`Age recipient`~ROUND, scale = 'free', label = 'label_both')+
      scale_fill_manual(values = c('pink', 'blue')) + 
      scale_color_manual(values = c('pink', 'blue')) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[ COMM == communities[i], unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-', 'transmisison_risk', '_age_byround_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}


plot_median_age_source <- function(median_age_source, outdir){
  
  cat("\nPlot median age at transmission of the source by age at infection of recipient\n")
  
  p <- ggplot(median_age_source) + 
    geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = PERIOD)) + 
    geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, fill = PERIOD), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient', y = 'Median age at transmission source') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +
    facet_grid(LABEL_DIRECTION~LABEL_COMMUNITY) 
  
  ggsave(p, file = paste0(outdir, '-MedianAgeSource_ByAgeRecipient.png'), w = 7, h = 6)
  

  
}

plot_PPC_augmented_recipient <- function(predict_z, incidence_cases_recipient, susceptible_recipient, outdir){
  
  predict_z <- merge(predict_z, incidence_cases_recipient[, .(INDEX_DIRECTION, INDEX_COMMUNITY, INDEX_TIME, AGE_INFECTION.RECIPIENT, INCIDENT_CASES, INCIDENT_CASES_UB, INCIDENT_CASES_LB)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'))
  predict_z <- merge(predict_z, susceptible_recipient[, .(INDEX_DIRECTION, INDEX_COMMUNITY, INDEX_TIME, AGE_INFECTION.RECIPIENT, SUSCEPTIBLE)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_INFECTION.RECIPIENT'))
  
  communities <- predict_z[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_z[ COMM == communities[i]]

    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
      geom_point(aes(y = INCIDENT_CASES), col = 'darkred') +
      geom_errorbar(aes(ymax = INCIDENT_CASES_UB, ymin = INCIDENT_CASES_LB), col = 'darkred', width = 0.2) +
      theme_bw() + 
      labs(x = 'Age at infection recipient', y = 'Augmented transmission events (Z)') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-PPC_augmented_recipient_', communities[i], '.png'), w = 7, h = 7)
    
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M/(PERIOD_SPAN*SUSCEPTIBLE))) +
      geom_ribbon(aes(ymin = CL/(PERIOD_SPAN*SUSCEPTIBLE), ymax = CU/(PERIOD_SPAN*SUSCEPTIBLE)), alpha = 0.5) +
      geom_point(aes(y = INCIDENT_CASES/(PERIOD_SPAN*SUSCEPTIBLE)), col = 'darkred') +
      geom_errorbar(aes(ymax = INCIDENT_CASES_UB/(PERIOD_SPAN*SUSCEPTIBLE), ymin = INCIDENT_CASES_LB/(PERIOD_SPAN*SUSCEPTIBLE)), col = 'darkred', width = 0.2) +
      theme_bw() + 
      labs(x = 'Age at infection recipient', y = 'Augmented transmission events (Z) per PY') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-PPC_augmented_perPY_recipient_', communities[i], '.png'), w = 7, h = 7)
  }
  
}

plot_PPC_observed_recipient <- function(predict_y, count_data, outdir){
  
  data <- count_data[, list(count = sum(count)), by = c('LABEL_DIRECTION', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_INFECTION.RECIPIENT', 'PERIOD_SPAN')]
  
  communities <- predict_y[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_y[ COMM == communities[i]]
    tmp1 <- data[ LABEL_COMMUNITY == tmp[, unique(LABEL_COMMUNITY)]]
    
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
      geom_point(data = tmp1, aes(y = count), col = 'darkred') +
      theme_bw() + 
      labs(x = 'Age at infection', y = 'Observed transmission events (Y)') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-PPC_observed_recipient_', communities[i], '.png'), w = 7, h = 7)
    
  }
  
}

plot_PPC_observed_source <- function(predict_y, count_data, outdir){
  
  data <- count_data[, list(count = sum(count)), by = c('LABEL_DIRECTION', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_TRANSMISSION.SOURCE', 'PERIOD_SPAN')]
  
  communities <- predict_y[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_y[ COMM == communities[i]]
    tmp1 <- data[ LABEL_COMMUNITY == tmp[, unique(LABEL_COMMUNITY)]]
    
    p <- ggplot(tmp, aes( x = AGE_TRANSMISSION.SOURCE)) + 
      geom_line(aes(y = M)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
      geom_point(data = tmp1, aes(y = count), col = 'darkred') +
      theme_bw() + 
      labs(x = 'Age at transmission', y = 'Observed transmission events (Y)') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    
    ggsave(p, file = paste0(outdir, '-output-PPC_observed_source_', communities[i], '.png'), w = 7, h = 7)
    
   }
  
}

plot_observed_to_augmented <- function(predict_y, predict_z, unsuppressed_count, outdir){
  
  predict_z[, type := 'Augmented (Z)']
  predict_y[, type := 'Observed (Y)']
  
  predict_df <- rbind(predict_z, predict_y)
  
  predict_df <- merge(predict_df, unsuppressed_count[, .(INDEX_DIRECTION, INDEX_COMMUNITY, INDEX_TIME, AGE_TRANSMISSION.SOURCE, count)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'AGE_TRANSMISSION.SOURCE'))
  
  communities <- predict_df[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_df[ COMM == communities[i]]

    p <- ggplot(tmp, aes( x = AGE_TRANSMISSION.SOURCE)) + 
      geom_line(aes(y = M, col = type)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.5) + 
      theme_bw() + 
      labs(x = 'Age at transmission', y = 'Transmission events') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-observed_vs_augmented_source_', communities[i], '.png'), w = 7, h = 7)
    
    p <- ggplot(tmp, aes( x = AGE_TRANSMISSION.SOURCE)) + 
      geom_line(aes(y = M/(count*PERIOD_SPAN), col = type)) + 
      geom_ribbon(aes(ymin = CL/(count*PERIOD_SPAN), ymax = CU/(count*PERIOD_SPAN), fill = type), alpha = 0.5) + 
      theme_bw() + 
      labs(x = 'Age at transmission', y = 'Transmission events per PY') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-observed_vs_augmented_perPY_source_', communities[i], '.png'), w = 7, h = 7)
  }
  
}
