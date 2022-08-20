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

plot_intensity_PP_by_round <- function(intensity_PP, outdir){
  
  communities <- intensity_PP[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- intensity_PP[ COMM == communities[i]]
    
    p <- ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = 'Estimated median transmission rate\nper year', 
           y= 'Age at transmission source',size='Pairs\ncount') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~LABEL_ROUND) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) + 
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    
    ggsave(p, file = paste0(outdir, '-output-intensity_transmission_by_round_',  communities[i], '.png'), w = 12, h = 7)
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
    
    if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
      p <- p + facet_grid(LABEL_DIRECTION~ROUND)
    }else if('LABEL_DIRECTION' %in% names(tmp)){
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
  
  if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~ROUND)
  }else if('LABEL_DIRECTION' %in% names(tmp)){
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
  
  if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~ROUND)
  }else if('LABEL_DIRECTION' %in% names(tmp)){
    p <- p + facet_grid(LABEL_DIRECTION~.)
  }
  
  ggsave(p, file = paste0(outdir, '-output-contrast_recipient_', name, '.png'), w = 7, h = 7)
  
} 

plot_force_infection <- function(force_infection, outdir, lab = NULL){
  
  communities <- force_infection[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- force_infection[ COMM == communities[i]]
    tmp <- tmp[ROUND != 'R014']
    p <- ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = paste0('Estimated median\n', lab, ' force infection'), 
           y= 'Age at transmission source',size='Pairs\ncount') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~ROUND) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) + 
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    
    ggsave(p, file = paste0(outdir, '-output-', lab, 'force_infection_2D_',  communities[i], '.png'), w = 10, h = 7)
  }
  
} 

plot_force_infection_age_source <- function(force_infection_age_source, outdir){
  
  tmp <- copy(force_infection_age_source)
  
  tmp[, Age := AGE_TRANSMISSION.SOURCE]
  tmp[, Direction :=LABEL_DIRECTION]
  
  communities <- tmp[, unique(COMM)]
  
  tmp <- tmp[ROUND != 'R014']
  
  p <- ggplot(tmp, aes(x = Age)) + 
    geom_bar(aes(y = M, fill = Direction), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = Direction), position = position_dodge()) + 
    labs(x = 'Age', y = 'Force of infection exerted', fill = '') + 
    theme_bw() +
    facet_grid(ROUND~LABEL_COMMUNITY, scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_source.png'), w = 12, h = 9)
  
  p <- ggplot(tmp, aes(x = Age)) + 
    geom_bar(aes(y = M, fill = ROUND), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = ROUND), position = position_dodge()) + 
    labs(x = 'Age', y = 'Force of infection exerted', fill = '') + 
    theme_bw() +
    facet_grid(Direction~LABEL_COMMUNITY, scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_source2.png'), w = 12, h = 9)
  
}

plot_force_infection_age_recipient <- function(force_infection_age_recipient, outdir){
  
  tmp <- copy(force_infection_age_recipient)
  
  tmp[, Age := AGE_INFECTION.RECIPIENT]
  tmp[, Direction :=LABEL_DIRECTION]
  
  communities <- tmp[, unique(COMM)]
  
  tmp <- tmp[ROUND != 'R014']
  
  p <- ggplot(tmp, aes(x = Age)) + 
    geom_bar(aes(y = M, fill = Direction), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = Direction), position = position_dodge()) + 
    labs(x = 'Age', y = 'Force of infection received', fill = '') + 
    theme_bw() +
    facet_grid(ROUND~LABEL_COMMUNITY, scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_recipient.png'), w = 12, h = 9)
  
  p <- ggplot(tmp, aes(x = Age)) + 
    geom_bar(aes(y = M, fill = ROUND), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = ROUND), position = position_dodge()) + 
    labs(x = 'Age', y = 'Force of infection received', fill = '') + 
    theme_bw() +
    facet_grid(Direction~LABEL_COMMUNITY, scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_recipient2.png'), w = 12, h = 9)
  
}


plot_force_infection_sex_source <- function(force_infection_sex_source, outdir){
  
  tmp <- copy(force_infection_sex_source)
  
  tmp[, Direction :=LABEL_DIRECTION]
  
  tmp <- tmp[ROUND != 'R014']
  
  p <- ggplot(tmp, aes(x = ROUND)) + 
    geom_bar(aes(y = M, fill = Direction), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = Direction), position = position_dodge()) + 
    labs(x = '', y = 'Force of infection exerted', fill = '') + 
    theme_bw() +
    facet_grid(LABEL_COMMUNITY~., scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_sex.png'), w = 10, h = 7)
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

plot_contribution_sex_source <- function(contribution_sex_source, unsuppressed_prop_sex, prevalence_prop_sex,outdir, lab = NULL){
  
  type_cont <- 'Contribution from male sources\nto HIV-1 infection'
  
  tmp <- copy(contribution_sex_source)
  tmp[, type  := type_cont]
  # tmp <- rbind(tmp, unsuppressed_prop_sex, fill=TRUE)
  # tmp[, type := factor(type , levels = c(type_cont,unique(unsuppressed_prop_sex$type), unique(prepare_prevalence_sex$type)))]
  # 
  tmp <- tmp[LABEL_DIRECTION == 'Male -> Female']
  # tmp <- tmp[ROUND != 'R014' & LABEL_DIRECTION == 'Male -> Female']
  

  prevalence_prop_sex[, type  := 'Share of males among HIV+\nindividuals']
  unsuppressed_prop_sex[, type  := 'Share of males among HIV+\nunsuppressed individuals']
  tmp1 <- rbind(prevalence_prop_sex, unsuppressed_prop_sex, fill=TRUE)
  tmp1[, type := factor(type , levels = c(unique(prevalence_prop_sex$type), unique(unsuppressed_prop_sex$type)))]
  tmp1 <- tmp1[LABEL_DIRECTION == 'Male -> Female']
  tmp1[, round :=as.numeric(gsub('R0(.+)', '\\1', ROUND)) ]
  
  p <-  ggplot(tmp, aes(x = round)) + 
     geom_line(aes(y = M, col = type)) + 
     geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5, fill = 'lightblue3') + 
    geom_errorbar(data = tmp1, aes(ymin = CL, ymax = CU), col = 'grey50', width = 0.1, size = 0.5)  +
    geom_line(data = tmp1, aes(y = M, linetype = type), col = 'black') +
    geom_point(data = tmp1, aes(y = M, shape = type), col = 'black', size= 2) + 
     labs(x = '', y = 'Percent', col = '', fill = '') + 
     theme_bw() +
     facet_grid(.~LABEL_COMMUNITY) + 
     scale_color_manual(values = 'lightblue3') + 
     scale_linetype_manual(values = c('dotted','solid')) + 
     scale_shape_manual(values = c(17, 16)) + 
     theme(strip.background = element_rect(colour="white", fill="white"),
           # axis.text.x = element_text(angle= 70, hjust = 1),
          strip.text = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 70, hjust = 1),
          axis.title.x = element_blank(), 
          # legend.justification = 'bottom',
          # legend.position='right',
          # legend.direction='vertical',
          legend.position = c(0.84,0.87),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          legend.margin = margin(),
          legend.spacing = unit(0.03, "cm"),
          legend.title = element_blank())  + 
     scale_x_continuous(labels = tmp[order(round), unique(LABEL_ROUND)], breaks = tmp[order(round), unique(round)]) + 
     scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0)), limits = c(0, 1)) + 
     guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 2))
   
   if(is.null(lab)) lab =  'Contribution'
  ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_sex.png'), w = 8, h = 5.3)
  
  # just inland for andrea's poster
  p <-  ggplot(tmp[COMM == 'inland'], aes(x = round)) + 
    geom_errorbar(data = tmp1[COMM == 'inland'], aes(ymin = CL, ymax = CU), col = 'grey50', width = 0.1, size = 0.5)  +
    geom_line(data = tmp1[COMM == 'inland'], aes(y = M, linetype = type), col = 'black') +
    geom_point(data = tmp1[COMM == 'inland'], aes(y = M, shape = type), col = 'black', size= 2) + 
    geom_line(aes(y = M, col = type)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU), alpha = 0.5, fill = 'lightblue3') + 
    labs(x = '', y = 'Percent', col = '', fill = '') + 
    theme_bw() +
    scale_color_manual(values = 'lightblue3') + 
    scale_linetype_manual(values = c('dotted','solid')) + 
    scale_shape_manual(values = c(17, 16)) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          # axis.text.x = element_text(angle= 70, hjust = 1),
          strip.text = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 70, hjust = 1),
          # legend.justification = 'bottom',
          # legend.position='right',
          # legend.direction='vertical',
          legend.position = c(0.80,0.87),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          legend.margin = margin(),
          legend.spacing = unit(0.03, "cm"),
          legend.title = element_blank())  + 
    scale_x_continuous(labels = tmp[order(round), unique(LABEL_ROUND)], breaks = tmp[order(round), unique(round)]) + 
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0)), limits = c(0, 1)) + 
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), shape = guide_legend(order = 2))
  
  if(is.null(lab)) lab =  'Contribution'
  ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_sex_inland.png'), w = 6, h = 5.3)
}

plot_contribution_age_source <- function(contribution_age_source, unsuppressed_prop_age, outdir, lab = NULL){
  
  type_cont <- 'Contribution to HIV-1 infection'
  
  tmp <- copy(contribution_age_source)
  tmp[, type  := type_cont]
  tmp[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' sources')]
  # tmp[, AGE_TRANSMISSION.SOURCE:= AGE_TRANSMISSION.SOURCE - 0.5] # account for the step function
  # tmpp <- copy(tmp[AGE_TRANSMISSION.SOURCE == max(AGE_TRANSMISSION.SOURCE)])
  # tmpp[, AGE_TRANSMISSION.SOURCE := AGE_TRANSMISSION.SOURCE + 1]
  # tmp <- rbind(tmp, tmpp)
  
  tmp1 <- copy(unsuppressed_prop_age)
  setnames(tmp1, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
  tmp1[, type := 'Share among HIV+\nunsuppressed individuals']
  tmp1[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' sources')]

  # tmp1[, AGE_TRANSMISSION.SOURCE:= AGE_TRANSMISSION.SOURCE - 0.5]
  # tmp11 <- copy(tmp1[AGE_TRANSMISSION.SOURCE == max(AGE_TRANSMISSION.SOURCE)])
  # tmp11[, AGE_TRANSMISSION.SOURCE := AGE_TRANSMISSION.SOURCE + 1]
  # tmp1 <- rbind(tmp1, tmp11)
  tmp1 <- merge(tmp1, df_round, by.x = 'ROUND', by.y = 'round')
  
  tmp2 <- tmp[round == min(round)]
  tmp2[, type := paste0(type, '\nin ', gsub('\n', ', ', LABEL_ROUND))]
  tmp2 <- select(tmp2, -LABEL_ROUND)
  # tmp2[, AGE_TRANSMISSION.SOURCE:= AGE_TRANSMISSION.SOURCE - 0.5] # account for the step function
  # tmp22 <- copy(tmp2[AGE_TRANSMISSION.SOURCE == max(AGE_TRANSMISSION.SOURCE)])
  # tmp22[, AGE_TRANSMISSION.SOURCE := AGE_TRANSMISSION.SOURCE + 1]
  # tmp2 <- rbind(tmp2, tmp22)
  
  communities <- tmp[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    tmp.p <- tmp[COMM == communities[i]]
    tmp1.p <- tmp1[COMM == communities[i]]
    tmp2.p <- tmp2[COMM == communities[i]]
    
    p <- ggplot(tmp.p, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX, size = type), alpha = 0.5) + 
      geom_line(aes(y = M, col = SEX, size = type), stat = 'identity', position = "identity") + 
      scale_color_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) +
      new_scale_color() +
      geom_line(data = tmp1.p, aes(y = M,alpha = type)) + 
      geom_line(data = tmp2.p, aes(y = M, col = SEX, linetype=type)) + 
      scale_color_manual(values = c('Male sources'='royalblue3','Female sources'='deeppink')) + 
      labs(x = 'Age of the source', y = 'Percent') + 
      theme_bw() +
      facet_grid(LABEL_ROUND~SEX)+
      # scale_color_manual(values = c('#FF4949')) + 
      scale_alpha_manual(values = 1) + 
      scale_size_manual(values = 0.5) + 
      scale_fill_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) + 
      # scale_fill_manual(values = c('Male source'='#C6DCE4','Female source'='#F2D1D1')) + 
      scale_linetype_manual(values  = 'dashed') +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'none', 
            legend.title = element_blank(), 
            panel.grid.minor = element_blank()) + 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05)), limits = c(0,NA))+ 
      scale_x_continuous(breaks = c(seq(min(tmp.p[, unique(AGE_TRANSMISSION.SOURCE)]), max(tmp.p[, unique(AGE_TRANSMISSION.SOURCE)]), 5), 
                         max(tmp.p[, unique(AGE_TRANSMISSION.SOURCE)]))) 
      
    p_legend <- ggplot(tmp.p[SEX == 'Female sources'], aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.5) + 
      geom_line(aes(y = M, col = type), stat = 'identity', position = "identity") + 
      geom_line(data = tmp1.p, aes(y = M, alpha = type), col = 'black') + 
      geom_line(data = tmp2.p, aes(y = M, linetype=type), col = 'deeppink') + 
      theme_bw() +
      scale_alpha_manual(values = 1) + 
      scale_color_manual(values = 'lightpink1') +
      scale_fill_manual(values = 'lightpink1') + 
      scale_linetype_manual(values  = 'dashed') +
      theme(legend.position = 'bottom', 
            legend.title = element_blank()) + 
      guides(fill = guide_legend(order = 1), color=guide_legend(order = 1),
             alpha = guide_legend(order = 3), 
             linetype = guide_legend(order = 2))
    
    pp <- ggarrange(p, legend.grob = get_legend(p_legend), legend = 'bottom') + 
      theme(panel.background = element_rect(fill='white'))
    if(is.null(lab)) lab =  'Contribution'
    ggsave(pp, file = paste0(outdir, '-output-', lab, '_age_', communities[i], '.png'), w = 8, h = 9)
    
    ### with ci on the unsuppressed 
    
    p <- p +  geom_errorbar(data = tmp1.p, aes(ymin = CL, ymax = CU,alpha = type), width = 0.5, col='grey50')
    pp <- ggarrange(p, legend.grob = get_legend(p_legend), legend = 'bottom')+ 
      theme(panel.background = element_rect(fill='white'))
    if(is.null(lab)) lab =  'Contribution'
    ggsave(pp, file = paste0(outdir, '-output-', lab, '_age_', communities[i], '2.png'), w = 8, h = 9)

    
  }
}

plot_contribution_age_group <- function(contribution_age_group_source, outdir, lab = NULL){
  
  communities <- contribution_age_group_source[, unique(COMM)]
  
  tmp <- copy(contribution_age_group_source)
  tmp[, AGE_LABEL := paste0('Age recipient: ', AGE_GROUP_INFECTION.RECIPIENT)]
  tmp[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' sources')]
  
  
  for(i in seq_along(communities)){
    
    tmp1 <- tmp[ COMM == communities[i]]

    p <- ggplot(tmp1, aes(x = AGE_GROUP_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, fill = LABEL_ROUND), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = ROUND), position = position_dodge()) + 
      labs(x = 'Age source', y = 'Share in HIV-1 transmissions', fill = '') + 
      theme_bw() +
      facet_grid(AGE_LABEL~SEX)+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) +
      ggsci::scale_fill_lancet()  +
      # ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])+ 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) + 
      guides(fill = guide_legend(byrow= T, nrow = 2))
    
    if(is.null(lab)) lab =  'Contribution'
    ggsave(p, file = paste0(outdir, '-output-', lab, '_age_group_',  communities[i], '.png'), w = 7,h =8)
  }
  
}

plot_contribution_age_classification <- function(contribution_age_classification_source, outdir, lab = NULL){
  
  communities <- contribution_age_classification_source[, unique(COMM)]
  
  tmp <- copy(contribution_age_classification_source)
  tmp[, AGE_CLASSIFICATION.SOURCE := factor(AGE_CLASSIFICATION.SOURCE, levels  = c('Younger', 'Same age', 'Older'))]
  tmp[, AGE_LABEL := paste0('Age recipient: ', AGE_GROUP_INFECTION.RECIPIENT)]
  tmp[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' sources')]
  
  for(i in seq_along(communities)){
    
    tmp1 <- tmp[ COMM == communities[i]]
    
    p <- ggplot(tmp1, aes(x = AGE_CLASSIFICATION.SOURCE)) + 
      geom_bar(aes(y = M, fill = LABEL_ROUND), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = LABEL_ROUND), position = position_dodge()) + 
      labs(x = 'Age classification of the source', y = 'HIV-1 transmission flows', fill = '') + 
      theme_bw() +
      facet_grid(AGE_LABEL~SEX)+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) +
      ggsci::scale_fill_lancet()  +
      # ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])+ 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) + 
      guides(fill = guide_legend(byrow= T, nrow = 2))
    
    if(is.null(lab)) lab =  'Contribution'
    ggsave(p, file = paste0(outdir, '-output-', lab, '_age_classification_',  communities[i], '.png'), w = 7, h = 8)
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
      labs(x = 'Age source', y = 'Contribution to infection\nadjusted by the number of HIV+ unsuppressed', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg() +
      ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])+ 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) 
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_group_by_round_',  communities[i], '.png'), w = 8, h = 7)
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
      labs(x = 'Age source', y = 'Contribution to infection\nadjusted by the number of HIV+ unsuppressed', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg()  +
      ggtitle(contribution_age_classification_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])+ 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) 
    ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_age_classification_by_round_',  communities[i], '.png'), w = 9, h = 7)
  }
  
}

plot_transmission_risk_age_source <- function(transmission_risk_age_source, outdir){
  
  tmp <- copy(transmission_risk_age_source)
  tmp[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' source')]
  
  tmp2 <- tmp[round == min(round)]
  tmp2[, type := paste0('In ', gsub('\n', '', LABEL_ROUND))]
  tmp2 <- select(tmp2, -LABEL_ROUND)
  tmp2[, AGE_TRANSMISSION.SOURCE:= AGE_TRANSMISSION.SOURCE - 0.5] # account for the step function
  tmp22 <- copy(tmp2[AGE_TRANSMISSION.SOURCE == max(AGE_TRANSMISSION.SOURCE)])
  tmp22[, AGE_TRANSMISSION.SOURCE := AGE_TRANSMISSION.SOURCE + 1]
  tmp2 <- rbind(tmp2, tmp22)
  
  communities <- tmp[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    tmp.p <- tmp[COMM == communities[i]]
    tmp2.p <- tmp2[COMM == communities[i]]
    
    p <- ggplot(tmp.p, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M,fill = SEX), stat = 'identity', position = "identity") + 
      geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.5, col = 'grey40') + 
      geom_step(data = tmp2.p, aes(y = M, color = type), linetype = 'solid') + 
      labs(x = 'Age of the source', y = 'Transmission risk exerted per unsuppressed per year') + 
      theme_bw() +
      facet_grid(LABEL_ROUND~SEX)+
      scale_color_manual(values = c('#FF4949')) + 
      scale_fill_manual(values = c('Male source'='lightblue3','Female source'='lightpink1')) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            legend.title = element_blank()) + 
      scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
      guides(fill = 'none')
    
    ggsave(p, file = paste0(outdir, '-output-', 'transmission_risk_age_', communities[i], '.png'), w = 8, h = 9)
    
  }
}


plot_transmission_risk_sex_source <- function(transmission_risk_sex_source_round, outdir){
  
  tmp <- copy(transmission_risk_sex_source_round)
  tmp[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' source')  ]
  
  p <- ggplot(tmp, aes(x = LABEL_ROUND)) + 
    geom_bar(aes(y = M, fill = SEX), stat = 'identity', position = "identity") + 
    geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.2, col = 'grey40') + 
    labs(x = '', y = 'Transmission risk per unsuppressed per year', 
         col = '', alpha = '', fill ='') + 
    theme_bw() +
    facet_grid(LABEL_COMMUNITY~SEX)+
    scale_fill_manual(values = c('Male source'='lightblue3','Female source'='lightpink1')) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'none', 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    guides(alpha = 'none') + 
    scale_y_continuous(expand = expansion(mult = c(0, .05)))
  ggsave(p, file = paste0(outdir, '-output-', 'transmission_risk', '_sex.png'), w = 8, h = 6)
  
}



plot_median_age_source <- function(median_age_source, outdir){
  
  p <- ggplot(median_age_source) + 
    geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = ROUND)) + 
    geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, fill = ROUND), alpha = 0.15) + 
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


plot_PPC_augmented_recipient_round <- function(predict_z_recipient_round, incidence_cases_recipient_round, susceptible_recipient_count, outdir){
  
  predict_z <- merge(predict_z_recipient_round, incidence_cases_recipient_round[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, AGE_INFECTION.RECIPIENT, INCIDENT_CASES, INCIDENT_CASES_UB, INCIDENT_CASES_LB)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'))
  predict_z <- merge(predict_z, susceptible_recipient_count[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, AGE_INFECTION.RECIPIENT, SUSCEPTIBLE)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'))
  
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
      facet_grid(LABEL_DIRECTION~ROUND) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-PPC_augmented_recipient_byround_', communities[i], '.png'), w = 10, h = 7)
    
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M/(ROUND_SPANYRS*SUSCEPTIBLE))) +
      geom_ribbon(aes(ymin = CL/(ROUND_SPANYRS*SUSCEPTIBLE), ymax = CU/(ROUND_SPANYRS*SUSCEPTIBLE)), alpha = 0.5) +
      geom_point(aes(y = INCIDENT_CASES/(ROUND_SPANYRS*SUSCEPTIBLE)), col = 'darkred') +
      geom_errorbar(aes(ymax = INCIDENT_CASES_UB/(ROUND_SPANYRS*SUSCEPTIBLE), ymin = INCIDENT_CASES_LB/(ROUND_SPANYRS*SUSCEPTIBLE)), col = 'darkred', width = 0.2) +
      theme_bw() + 
      labs(x = 'Age at infection recipient', y = 'Augmented transmission events (Z) per PY') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~ROUND) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-PPC_augmented_perPY_recipient_byround_', communities[i], '.png'), w = 10, h = 7)
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
      labs(x = 'Age at transmission', y = 'Transmission events per unsuppressed per year') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-observed_vs_augmented_perPY_source_', communities[i], '.png'), w = 7, h = 7)
  }
  
}


