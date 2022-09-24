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
    
    if(all(c('ROUND', 'LABEL_DIRECTION', 'COMM') %in% names(tmp))){
      p <- p + facet_grid(LABEL_DIRECTION+COMM~ROUND)
    }else if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
      p <- p + facet_grid(LABEL_DIRECTION~ROUND)
    }else if(all(c('COMM', 'LABEL_DIRECTION') %in% names(tmp))){
      p <- p + facet_grid(LABEL_DIRECTION~COMM)
    }else if('LABEL_DIRECTION' %in% names(tmp)){
      p <- p + facet_grid(LABEL_DIRECTION~.)
    }
    
    
    ggsave(p, file = paste0(outdir, '-output-contrast_2D_', name, '.png'), w = 9, h = 9)
  
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

  if(all(c('ROUND', 'LABEL_DIRECTION', 'COMM') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION+COMM~ROUND)
  }else if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~ROUND)
  }else if(all(c('COMM', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~COMM)
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
  
  if(all(c('ROUND', 'LABEL_DIRECTION', 'COMM') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION+COMM~ROUND)
  }else if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~ROUND)
  }else if(all(c('COMM', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~COMM)
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
  
  type_cont <- 'Contribution from male sources\nto HIV infection'
  
  tmp <- copy(contribution_sex_source)
  tmp[, type  := type_cont]

  tmp <- tmp[LABEL_DIRECTION == 'Male -> Female']
  tmp[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 7, 0)]
  tmp[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]

  prevalence_prop_sex[, type  := 'Share of males among HIV-positive\nindividuals']
  unsuppressed_prop_sex[, type  := 'Share of males among HIV-positive\nunsuppressed individuals']
  tmp1 <- rbind(prevalence_prop_sex, unsuppressed_prop_sex, fill=TRUE)
  tmp1[, type := factor(type , levels = c(unique(prevalence_prop_sex$type), unique(unsuppressed_prop_sex$type)))]
  tmp1 <- tmp1[LABEL_DIRECTION == 'Male -> Female']
  tmp1[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 7, 0)]
  
  tmp2 <- rbind(tmp, tmp1, fill=TRUE)
  tmp2[, type := factor(type , levels = c(unique(tmp$type), levels(tmp1$type)))]
  
  p <-  ggplot(tmp2, aes(x = INDEX_ROUND2, group = type)) +
      geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0, size = 0.5, position = position_dodge(width = 0.2))  +
    geom_line(aes(y = M, linetype = type), col = 'black', position = position_dodge(width = 0.2)) +
    geom_point(aes(y = M, shape = type, col = type), size= 1.7, position = position_dodge(width = 0.2)) + 
     labs(x = '', y = 'Percent', col = '', shape = '', linetype = '') + 
     theme_bw() +
     facet_grid(.~LABEL_COMMUNITY, scale = 'free_x') + 
     scale_linetype_manual(values = c('solid', 'dashed', 'dotted')) + 
     scale_shape_manual(values = c(16, 17, 15)) + 
      scale_color_manual(values = c("black", "#E69F00",  "#009E73")) + 
     theme(strip.background = element_rect(colour="white", fill="white"),
           # axis.text.x = element_text(angle= 70, hjust = 1),
          strip.text = element_text(size = rel(1)),
          axis.text.x = element_text(angle = 40, hjust = 1),
          axis.title.x = element_blank(), 
          # legend.justification = 'bottom',
          # legend.position='right',
          # legend.direction='vertical',
          legend.position = c(0.81,0.14),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          legend.margin = margin(),
          legend.title = element_blank())  + 
     scale_x_continuous(labels = tmp[order(INDEX_ROUND2), (LABEL_ROUND2)], breaks = tmp[order(INDEX_ROUND2), unique(INDEX_ROUND2)]) + 
     scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0)), limits = c(0, 1)) + 
     guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1))
   
   if(is.null(lab)) lab =  'Contribution'
  ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_sex.png'), w = 7, h = 4.5)
  
}

plot_contribution_age_source_unsuppressed <- function(contribution_age_source, unsuppressed_prop_age, outdir, lab = NULL){
  
  # restricted rounds
  Rounds <- paste0('R0', c(15,18))
  
  # prepare dataset
  type_cont <- 'Contribution to HIV infection'
  
  tmp <- copy(contribution_age_source)
  tmp[, type  := type_cont]
  tmp[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' sources')]
  
  tmp1 <- copy(unsuppressed_prop_age)
  setnames(tmp1, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
  tmp1[, type := 'Share among HIV-positive\nunsuppressed individuals']
  tmp1[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' sources')]
  
  tmp2 <- tmp[ROUND %in% Rounds]
  tmp2[, min_INDEX_ROUND :=  min(INDEX_ROUND), by = 'COMM']
  tmp2 <- tmp2[INDEX_ROUND == min_INDEX_ROUND]
  tmp2[, type := paste0(type, '\nin ', gsub('\n', ', ', LABEL_ROUND))]
  tmp2 <- select(tmp2, -LABEL_ROUND)
  
  communities <- tmp[, unique(COMM)]
  
  # prepare function
  plot.p <- function(tmp.p, tmp1.p, tmp2.p){
    ggplot(tmp.p, aes(x = AGE_TRANSMISSION.SOURCE)) +
      geom_ribbon(data = tmp1.p, aes(ymin = CL, ymax = CU, linetype=type), alpha = 0.3, fill='grey50')+
      geom_line(data = tmp1.p, aes(y = M,linetype = type)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX, size = type), alpha = 0.6) + 
      geom_line(aes(y = M, col = SEX, size = type), stat = 'identity', position = "identity") + 
      scale_color_manual(values = c('Male sources'='royalblue3','Female sources'='deeppink')) + 
      new_scale_color() +
      geom_line(data = tmp2.p, aes(y = M, col = SEX, alpha=type), linetype ='solid') +  
      scale_color_manual(values = c('Male sources'='paleturquoise4','Female sources'='pink4')) +
      labs(x = 'Age of the source', y = 'Percent') + 
      theme_bw() +
      facet_grid(LABEL_ROUND~SEX) +
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
  }
  
  
  # make plots
  for(i in seq_along(communities)){
    
    tmp.p <- tmp[COMM == communities[i]]
    tmp1.p <- tmp1[COMM == communities[i]]
    tmp2.p <- tmp2[COMM == communities[i]]
    
    
    ## all rounds
    p.all <- plot.p(tmp.p, tmp1.p, tmp2.p)
      
    ## round 15 and 18
    p <- plot.p(tmp.p[ROUND %in% Rounds], tmp1.p[ROUND %in% Rounds], tmp2.p)
    
    ## legend
    p_legend <- ggplot(tmp.p[SEX == 'Female sources'], aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.6) + 
      geom_line(aes(y = M, col = type), stat = 'identity', position = "identity") + 
      geom_line(data = tmp1.p, aes(y = M, size = type), col = 'black', linetype = 'dashed') + 
      geom_ribbon(data = tmp1.p, aes(ymin = CL, ymax = CU, size = type), alpha = 0.3, fill='grey50') + 
      geom_line(data = tmp2.p, aes(y = M, linetype=type), col = 'pink4') + 
      theme_bw() +
      scale_alpha_manual(values = 1) + 
      scale_size_manual(values = 0.5) + 
      scale_color_manual(values = 'deeppink') +
      scale_fill_manual(values = 'lightpink1') + 
      scale_linetype_manual(values  = 'solid') +
      theme(legend.position = 'bottom', 
            legend.title = element_blank()) + 
      guides(fill = guide_legend(order = 1), color=guide_legend(order = 1),
             alpha = guide_legend(order = 3), 
             linetype = guide_legend(order = 2))

    pp.all <- ggarrange(p.all, legend.grob = get_legend(p_legend), legend = 'bottom') + 
      theme(panel.background = element_rect(fill='white'))
    
    pp <- ggarrange(p, legend.grob = get_legend(p_legend), legend = 'bottom') + 
      theme(panel.background = element_rect(fill='white'))
    
    if(is.null(lab)) lab =  'Contribution'
    
    if(communities[i] == 'inland'){
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_entended_', communities[i], '.png'), w = 9.5, h = 12.5)
    }else{
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_entended_', communities[i], '.png'), w = 9.5, h = 9)
      
    }
    
    ggsave(pp, file = paste0(outdir, '-output-', lab, '_age_', communities[i], '.png'), w = 8, h = 6)
    
  }
}

plot_contribution_age_source <- function(contribution_age_source, outdir, lab = NULL){
  
  # restricted rounds
  Rounds <- paste0('R0', c(15,18))
  
  # prepare dataset
  tmp <- copy(contribution_age_source)
  tmp[, SEX := paste0(gsub('(.+) ->.*', '\\1', LABEL_DIRECTION), ' sources')]
  
  communities <- tmp[, unique(COMM)]
  
  # prepare function
  plot.p <- function(tmp.p){
    ggplot(tmp.p, aes(x = AGE_TRANSMISSION.SOURCE)) +
      geom_bar(aes(y = M, fill = SEX), stat = 'identity', position = position_dodge(0.9)) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = SEX), alpha = 0.6, position = position_dodge(0.9), width = 0.5) +
      # scale_fill_manual(values = c('Male sources'='royalblue3','Female sources'='deeppink')) + 
      scale_fill_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) + 
      labs(x = 'Age of the source', y = 'Contribution to HIV infection') + 
      theme_bw() +
      facet_grid(LABEL_ROUND~.) +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.title = element_blank(), 
            panel.grid.minor = element_blank()) + 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05)), limits = c(0,NA))+ 
      scale_x_continuous(breaks = c(seq(min(tmp.p[, unique(AGE_TRANSMISSION.SOURCE)]), max(tmp.p[, unique(AGE_TRANSMISSION.SOURCE)]), 5), 
                                    max(tmp.p[, unique(AGE_TRANSMISSION.SOURCE)]))) 
  }
  
  
  # make plots
  for(i in seq_along(communities)){
    
    tmp.p <- tmp[COMM == communities[i]]
    
    ## all rounds
    pp.all <- plot.p(tmp.p) + theme(legend.position = 'bottom')
    
    ## round 15 and 18
    pp <- plot.p(tmp.p[ROUND %in% Rounds])+ theme(legend.position =  c(0.85,0.93))

    if(is.null(lab)) lab =  'Contribution_Sex'
    
    if(communities[i] == 'inland'){
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_entended_', communities[i], '.png'), w = 7, h = 12.5)
    }else{
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_entended_', communities[i], '.png'), w = 7, h = 9)
      
    }
    
    ggsave(pp, file = paste0(outdir, '-output-', lab, '_age_', communities[i], '.png'), w = 6, h = 6)
    
  }
}


plot_contribution_age_group <- function(contribution_age_group_source, outdir, lab = NULL){
  
  communities <- contribution_age_group_source[, unique(COMM)]
  
  tmp <- copy(contribution_age_group_source)
  tmp[, AGE_LABEL := paste0('Age recipient: ', AGE_GROUP_INFECTION.RECIPIENT)]
  tmp[, SEX := paste0(gsub('.* -> (.+)', '\\1', LABEL_DIRECTION), ' recipients')]
  
  
  for(i in seq_along(communities)){
    
    tmp1 <- tmp[ COMM == communities[i]]
    
    if(communities[i] == 'inland'){
      colors <- palette_round_inland
    }else{
      colors <- palette_round_fishing
    }

    p <- ggplot(tmp1, aes(x = AGE_GROUP_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, fill = LABEL_ROUND), stat = 'identity', position =position_dodge(width = 0.9)) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = LABEL_ROUND), position =position_dodge(width = 0.9), width = 0.3) + 
      labs(x = 'Age source', y = 'Share in HIV transmissions', fill = '') + 
      theme_bw() +
      facet_grid(AGE_LABEL~SEX)+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) +
      scale_fill_manual(values = colors) +
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
  tmp[, SEX := paste0(gsub('.* -> (.+)', '\\1', LABEL_DIRECTION), ' recipients')]
  
  for(i in seq_along(communities)){
    
    if(communities[i] == 'inland'){
      colors <- palette_round_inland
    }else{
      colors <- palette_round_fishing
    }
    
    tmp1 <- tmp[ COMM == communities[i]]
    
    p <- ggplot(tmp1, aes(x = AGE_CLASSIFICATION.SOURCE, group = LABEL_ROUND)) + 
      geom_bar(aes(y = M, fill = LABEL_ROUND), stat = 'identity', position =position_dodge(width = 0.9)) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = LABEL_ROUND), position =position_dodge(width = 0.9), width = 0.3) + 
      labs(x = 'Age classification of the source', y = 'HIV transmission flows', fill = '') + 
      theme_bw() +
      facet_grid(AGE_LABEL~SEX)+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) +
      scale_fill_manual(values = colors) +
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
      labs(x = 'Age source', y = 'Contribution to infection\nadjusted by the number of HIV-positive unsuppressed', fill = '') + 
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
      labs(x = 'Age source', y = 'Contribution to infection\nadjusted by the number of HIV-positive unsuppressed', fill = '') + 
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
    facet_grid(SEX~LABEL_COMMUNITY, scale = 'free_x')+
    scale_fill_manual(values = c('Male source'='lightblue3','Female source'='lightpink1')) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'none', 
          axis.text.x = element_text(angle = 40, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    guides(alpha = 'none') + 
    scale_y_continuous(expand = expansion(mult = c(0, .05)))
  ggsave(p, file = paste0(outdir, '-output-', 'transmission_risk', '_sex.png'), w = 8, h = 6.5)
  
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

plot_median_age_source_group <- function(median_age_source_group, outdir){
  
  communities <- median_age_source_group[, unique(COMM)]
  median_age_source_group[, SEX_LABEL := paste0(gsub('.* -> (.+)', '\\1', LABEL_DIRECTION), ' recipients')]
  median_age_source_group <- median_age_source_group[ROUND %in% paste0('R0', c(15, 18))]
  median_age_source_group[, mean_age_group := mean(c(as.numeric(gsub('(.+)-.*', '\\1', AGE_GROUP_INFECTION.RECIPIENT)), 
                                                    as.numeric(gsub('.*-(.+)', '\\1', AGE_GROUP_INFECTION.RECIPIENT)))), by = 'AGE_GROUP_INFECTION.RECIPIENT']
  cols <- palette_round_inland[c(4, 7)]
  
  for(i in seq_along(communities)){
    tmp <- median_age_source_group[COMM == communities[i]]
    
    p <- ggplot(tmp) + 
      geom_point(aes(x = mean_age_group, y = M, col = LABEL_ROUND), position = position_dodge(2)) + 
      geom_errorbar(aes(x = mean_age_group, ymin= CL, ymax = CU, col = LABEL_ROUND), alpha = 0.5, position = position_dodge(2), width = 1) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
      theme_bw() + 
      labs(x = 'Age recipient', y = 'Median age source') +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.minor.x = element_blank(), 
            legend.title = element_blank()) +
      facet_grid(.~SEX_LABEL) + 
      scale_color_manual(values = cols) + 
      scale_y_continuous(expand = c(0,0), 
                         breaks = c(seq(min(range_age_non_extended), max(range_age_non_extended), 5), 
                           max(range_age_non_extended))) + 
      coord_cartesian(ylim = range_age_non_extended, xlim= range_age_non_extended) +
      scale_x_continuous(breaks = median_age_source_group[order(median_age_source_group), unique(mean_age_group)], 
                         labels = median_age_source_group[order(median_age_source_group), unique(AGE_GROUP_INFECTION.RECIPIENT)], 
                         expand = c(0,0))
    
    
    ggsave(p, file = paste0(outdir, '-MedianAgeSource_ByAgeGroupRecipient_', communities[i], '.png'), w = 7, h = 4.5)
    
  }
}

plot_PPC_augmented_recipient_round <- function(predict_z_recipient_round, incidence_cases_recipient_round, eligible_count_recipient, outdir){
  
  predict_z <- merge(predict_z_recipient_round, incidence_cases_recipient_round[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, AGE_INFECTION.RECIPIENT, INCIDENT_CASES, INCIDENT_CASES_UB, INCIDENT_CASES_LB)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'))
  predict_z <- merge(predict_z, eligible_count_recipient[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, AGE_INFECTION.RECIPIENT, ELIGIBLE)], 
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
      geom_line(aes(y = M/(ROUND_SPANYRS*ELIGIBLE))) +
      geom_ribbon(aes(ymin = CL/(ROUND_SPANYRS*ELIGIBLE), ymax = CU/(ROUND_SPANYRS*ELIGIBLE)), alpha = 0.5) +
      geom_point(aes(y = INCIDENT_CASES/(ROUND_SPANYRS*ELIGIBLE)), col = 'darkred') +
      geom_errorbar(aes(ymax = INCIDENT_CASES_UB/(ROUND_SPANYRS*ELIGIBLE), ymin = INCIDENT_CASES_LB/(ROUND_SPANYRS*ELIGIBLE)), col = 'darkred', width = 0.2) +
      theme_bw() + 
      labs(x = 'Age at infection recipient', y = 'Augmented transmission events (Z) per PY') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~ROUND) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-PPC_augmented_perPY_recipient_byround_', communities[i], '.png'), w = 10, h = 7)
  
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M/(ROUND_SPANYRS*ELIGIBLE), col = LABEL_ROUND)) +
      # geom_ribbon(aes(ymin = CL/(ROUND_SPANYRS*ELIGIBLE), ymax = CU/(ROUND_SPANYRS*ELIGIBLE), fill = LABEL_ROUND), alpha = 0.5) +
      theme_bw() + 
      labs(x = 'Age at infection recipient', y = 'Augmented transmission events (Z) per PY') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~.) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-augmented_perPY_recipient_byround_', communities[i], '.png'), w = 7, h = 10)
    
    
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

plot_counterfactual <- function(counterfactuals_p_f, counterfactuals_p_a, eligible_count_round, lab, outdir){
  
  # target label
  df_target = data.table(variable = c('TREATED.SPREADERS', 'TREATED.NONCOMPLIER', 'TREATED.RANDOM'), 
                       counterfactual_index = 1:3, 
                       LABEL_TARGET = c('Main spreaders', 'Greatest difference with\nfemale ART uptake', 'Random'))
  df_target[, LABEL_TARGET := factor(LABEL_TARGET, levels = df_target[, LABEL_TARGET])]
  
  # unlist as much as female
  label.f = 'ART same as female'
  budget.counterfactual <- counterfactuals_p_f$budget 
  relative_incidence_counterfactual <- counterfactuals_p_f$relative_incidence_counterfactual 
  incidence_counterfactual <- counterfactuals_p_f$incidence_counterfactual 
  budget.counterfactual[, label := label.f]
  relative_incidence_counterfactual[, label := label.f]
  incidence_counterfactual[, label := label.f]
  
  # unlist all
  label.a = 'ART all'
  budget.counterfactual.a <- counterfactuals_p_a$budget 
  relative_incidence_counterfactual.a <- counterfactuals_p_a$relative_incidence_counterfactual 
  incidence_counterfactual.a <- counterfactuals_p_a$incidence_counterfactual 
  budget.counterfactual.a[, label := label.a]
  relative_incidence_counterfactual.a[, label := label.a]
  incidence_counterfactual.a[, label := label.a]
  
  # combine
  budget.counterfactual <- rbind(budget.counterfactual, budget.counterfactual.a)
  relative_incidence_counterfactual <- rbind(relative_incidence_counterfactual, relative_incidence_counterfactual.a)
  incidence_counterfactual <- rbind(incidence_counterfactual, incidence_counterfactual.a)
  
  # restrict to one round
  Round <- 'R018'
  budget.counterfactual <- budget.counterfactual[ROUND == Round & SEX == 'M']
  relative_incidence_counterfactual <- relative_incidence_counterfactual[ROUND == Round & IS_MF == T]
  incidence_counterfactual <- incidence_counterfactual[ROUND == Round & IS_MF == T]
  
  # format budget
  bc <- melt.data.table(budget.counterfactual, id.vars = c('ROUND', 'SEX', 'COMM', 'label'))
  bc <- merge(bc, df_target, by = 'variable')
  bc[, lab := lab]
  
  # find incidence rate
  tmp <- copy(eligible_count_round[, .(ROUND, SEX, COMM, AGEYRS, SUSCEPTIBLE)])
  tmp[, IS_MF := ifelse(SEX == 'F', 1, 0)]
  ic <- merge(incidence_counterfactual, tmp, 
                                    by.x = c('ROUND', 'COMM', 'AGE_INFECTION.RECIPIENT', 'IS_MF'), 
                                    by.y = c('ROUND', 'COMM', 'AGEYRS', 'IS_MF'))
  ic[, M := M / SUSCEPTIBLE]
  ic[, CL := CL / SUSCEPTIBLE]
  ic[, CU := CU / SUSCEPTIBLE]
  ic <- merge(ic, df_target, by = 'counterfactual_index')
  
  # reduction in incidence
  ric <- merge(relative_incidence_counterfactual, df_target, by = 'counterfactual_index')
    
  # sex label
  # icf[, SEX := 'Female']
  # icf[LABEL_DIRECTION == 'Female -> Male', SEX := 'Male']
  # 
  
  communities <- ric[, unique(COMM)]
  cols <- c('#CC3636', '#F57328', '#367E18')
  
  for(i in seq_along(communities)){

      tmp1 <- bc[COMM == communities[i]]
      tmp2 <- ic[COMM == communities[i]]
      tmp3 <- ric[COMM == communities[i]]
      
      # budget
      p1 <- ggplot(tmp1, aes(x = label)) + 
        geom_bar(aes(fill = LABEL_TARGET, y = value), stat = 'identity', position = position_dodge(0.7), width = 0.6) + 
        labs(x = '', y = paste0('\nNumber of male treated'), fill = '') + 
        theme_bw() +
        theme(strip.background = element_rect(colour="white", fill="white"),
              strip.text = element_text(size = rel(1)),
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(), 
              # legend.position = c(0.85, 0.83), 
              legend.key.size = unit(0.4, 'cm'),
              legend.position = 'none', 
              axis.title.x = element_blank(),
              legend.title = element_blank()
        ) +
        scale_fill_manual(values = cols) +
        scale_y_continuous(expand = expansion(mult = c(0, .05))) +
        facet_grid(.~lab)
      
      # incidence rate
      p2 <- ggplot(tmp2) + 
        geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M*100, col = LABEL_TARGET)) + 
        geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL*100, ymax = CU*100, group = interaction(LABEL_TARGET), fill = LABEL_TARGET), alpha = 0.25) + 
        labs(x = 'Age female recipient', y = 'Incidence rate per 100 person-year\nin female recipient') + 
        theme_bw() +
        theme(strip.background = element_rect(colour="white", fill="white"),
              strip.text = element_text(size = rel(1)),
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(), 
              legend.position = 'none',
              legend.title = element_blank()) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        scale_y_continuous() + 
        scale_x_continuous(expand = c(0,0)) + 
        facet_grid(.~label)

      # reduction infection
      p3 <- ggplot(tmp3, aes(x = AGEYRS)) + 
        geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = LABEL_TARGET)) + 
        geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, group = interaction(LABEL_TARGET), fill = LABEL_TARGET), alpha = 0.25) + 
        labs(x = 'Age female recipient', y = '% reduction in incidence infections\namong female recipients', 
             fill = '', col = '') + 
        theme_bw() +
        theme(strip.background = element_rect(colour="white", fill="white"),
              strip.text = element_text(size = rel(1)),
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(), 
              legend.title = element_blank(), 
              # legend.direction = 'vertical',
              # axis.title.x = element_blank(), 
              # axis.text.x = element_blank(), 
              legend.position = 'bottom') +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
        scale_x_continuous(expand = c(0,0)) + 
        facet_grid(.~label)

      p <- grid.arrange(p1, p2, p3, layout_matrix = rbind(c(NA, 1), c(NA, 2), c(3,3)), widths = c(0.02, 0.95), heights = c(0.3, 0.3, 0.42))
      
      file = paste0(outdir, '-output-counterfactual_incidence_panel_', gsub(' ' , '', lab), '_', communities[i], '.png')
      ggsave(p, file = file, w = 7, h = 9)
      

  }
}


plot_counterfactual_relative_incidence_old2 <- function(eligible_count_round.counterfactual, relative_incidence_counterfactual, 
                                                   incidence_factual, incidence_counterfactual, outdir, only_participant= F){
  
  # restrict to one round
  Round <- 'R018'
  icf <- incidence_factual[ROUND == Round]
  icc <- incidence_counterfactual[ROUND == Round & IS_MF == T]
  rei <- relative_incidence_counterfactual[ROUND == Round & IS_MF == T]
  ecr <- eligible_count_round.counterfactual[ROUND == Round & SEX == 'M']

  # sex label
  icf[, SEX := 'Female']
  icf[LABEL_DIRECTION == 'Female -> Male', SEX := 'Male']
  
  communities <- rei[, unique(COMM)]
  cols <- c('#CC3636', '#F57328', '#367E18')
  counterfactual_compared <- list(c(1,4,3), c(2,5,3))
  
  # label male
  male_label = ''
  if(only_participant) male_label = 'participants'
  
  for(i in seq_along(communities)){
    for(j in 1:length(counterfactual_compared)){
      Counterfactual <- counterfactual_compared[[j]]
      
      tmp <- rei[COMM == communities[i] & counterfactual_index %in% Counterfactual]
      tmp1 <- ecr[COMM == communities[i] & counterfactual_index %in% Counterfactual]
      tmp2 <- icf[COMM == communities[i]]
      tmp3 <- icc[COMM == communities[i] & counterfactual_index %in% Counterfactual]
      
      # counterfactual label
      counterfactual_label <- paste0('Counterfactual with higher ART\ncoverage among male ', male_label)
      contributions <- spreaders[spreader_category %in% 1:2, sort(unique(label))]
      artdiff <- noncomplier[ROUND == Round & COMM == communities[i] & spreader_category %in% 4:5, round(sort(unique(label), decreasing = T), 2)]
      counterfactual_labels <- c(paste0('Males contributing to ', contributions, '% of transmissions'),
                                 paste0('All male ', male_label), 
                                 paste0('Males with ART coverage difference compared to\nfemale of at least ', artdiff))
      counterfactual_labels_index <- c(1, 2, 4, 5, 3)
      
      tmp[, COUNTERFACTUAL_LABEL := factor(counterfactual_labels[counterfactual_index], levels = counterfactual_labels[counterfactual_labels_index])]
      tmp1[, COUNTERFACTUAL_LABEL := factor(counterfactual_labels[counterfactual_index], levels = (counterfactual_labels[counterfactual_labels_index]))]
      tmp3[, COUNTERFACTUAL_LABEL := factor(counterfactual_labels[counterfactual_index], levels = (counterfactual_labels[counterfactual_labels_index]))]
      
      # percentage point increase
      tmp1[, index_plot := which(Counterfactual == counterfactual_index), by = 'counterfactual_index']
      p1 <- ggplot(tmp1, aes(x = AGEYRS)) + 
        geom_step(aes(col = COUNTERFACTUAL_LABEL, y = INCREASE_ART_COVERAGE + 0.005*index_plot)) + 
        labs(x = 'Age source', y = paste0('Percentage point increase in ART\ncoverage among male ', male_label), col = counterfactual_label) + 
        theme_bw() +
        theme(strip.background = element_rect(colour="white", fill="white"),
              strip.text = element_text(size = rel(1)),
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(), 
              legend.position = c(0.78, 0.76), 
              legend.key.size = unit(0.4, 'cm'),
              legend.direction = 'vertical', 
              # legend.title = element_blank()
        ) +
        scale_color_manual(values = cols) +
        scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
        scale_x_continuous(expand = c(0,0)) 
      
      p2 <- ggplot() + 
        geom_line(data = tmp2, aes(x = AGE_INFECTION.RECIPIENT, y = M, linetype = SEX), col = 'black') + 
        geom_ribbon(data = tmp2, aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, group = SEX), alpha = 0.15) + 
        geom_line(data = tmp3, aes(x = AGE_INFECTION.RECIPIENT, y = M, col = COUNTERFACTUAL_LABEL)) + 
        geom_ribbon(data = tmp3, aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, fill = COUNTERFACTUAL_LABEL), alpha = 0.15) + 
        labs(x = 'Age', y = '\nHIV incidence infections', 
             fill = counterfactual_label, col = counterfactual_label) + 
        theme_bw() +
        theme(strip.background = element_rect(colour="white", fill="white"),
              strip.text = element_text(size = rel(1)),
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(), 
              axis.title.x = element_blank(), 
              # axis.text.x = element_blank(), 
              legend.position = c(0.9, 0.86),
              legend.title = element_blank()) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        scale_y_continuous() + 
        scale_x_continuous(expand = c(0,0)) +
        guides(color = 'none', fill = 'none')

      p3 <- ggplot(tmp, aes(x = AGEYRS)) + 
        geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = COUNTERFACTUAL_LABEL)) + 
        geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, fill = COUNTERFACTUAL_LABEL), alpha = 0.5) + 
        labs(x = 'Age', y = '% reduction in HIV incidence infections\namong female recipients', 
             fill = counterfactual_label, col = counterfactual_label) + 
        theme_bw() +
        theme(strip.background = element_rect(colour="white", fill="white"),
              strip.text = element_text(size = rel(1)),
              panel.grid.major.x = element_blank(), 
              panel.grid.minor.x = element_blank(), 
              # axis.title.x = element_blank(), 
              # axis.text.x = element_blank(), 
              legend.position = 'none') +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
        scale_x_continuous(expand = c(0,0)) 

      p <- grid.arrange(p1, p2, p3, layout_matrix = rbind(c(NA, 1), c(NA, 2), c(3,3)), widths = c(0.02, 0.95), heights = c(0.33, 0.33, 0.33))
      
      file = paste0(outdir, '-output-counterfactual_incidence_panel_', j, '_', communities[i], '.png')
      if(only_participant){
        file = paste0(outdir, '-output-counterfactual_incidence_panel_', j, '_', 'only_participant', '_', communities[i], '.png')
      }
      ggsave(p, file = file, w = 8, h = 9)
      
      
    }
    
    
  }
}

plot_counterfactual_relative_incidence_old <- function(eligible_count_round.counterfactual, relative_incidence_counterfactual, 
                                                   incidence_factual, incidence_counterfactual, outdir){
  
  # restrict to one round
  Round <- 'R018'
  icf <- incidence_factual[ROUND == Round]
  icc <- incidence_counterfactual[ROUND == Round & IS_MF == T]
  rei <- relative_incidence_counterfactual[ROUND == Round & IS_MF == T]
  ecr <- eligible_count_round.counterfactual[ROUND == Round & SEX == 'M']
  
  # sex label
  icf[, SEX := 'Female']
  icf[LABEL_DIRECTION == 'Female -> Male', SEX := 'Male']
  
  communities <- rei[, unique(COMM)]
  cols <- c('#F08A5D', '#B83B5E', '#6A2C70')
  
  for(i in seq_along(communities)){
    
    tmp <- rei[COMM == communities[i]]
    tmp1 <- ecr[COMM == communities[i]]
    tmp2 <- icf[COMM == communities[i]]
    tmp3 <- icc[COMM == communities[i]]
    
    # counterfactual label
    counterfactual_label <- 'Counterfactual with higher ART\ncoverage among male sources'
    
    contributions <- spreaders[spreader_category %in% 1:2, sort(unique(label))]
    artdiff <- noncomplier[ROUND == Round & COMM == communities[i] & spreader_category %in% 4:5, round(sort(unique(label))* 100)]
    counterfactual_labels <- c(paste0('Male sources contributing to ', contributions, '% of transmissions'),
                               'All male sources', 
                               paste0('Males with ART coverage difference compared to female at least ', artdiff))
    
    tmp[, COUNTERFACTUAL_LABEL := factor(counterfactual_labels[counterfactual_index], levels = counterfactual_labels)]
    tmp1[, COUNTERFACTUAL_LABEL := factor(counterfactual_labels[counterfactual_index], levels = (counterfactual_labels))]
    tmp3[, COUNTERFACTUAL_LABEL := factor(counterfactual_labels[counterfactual_index], levels = (counterfactual_labels))]
    
    # percentage point increase
    tmp1[is.na(spreader_category), spreader_category := 4]
    tmp1 <- tmp1[!(spreader_category == 1 & counterfactual_index %in% 2:3) ]
    tmp1 <- tmp1[!(spreader_category == 2 & counterfactual_index %in% 3) ]
    
    p1 <- ggplot(tmp1, aes(x = AGEYRS, y = INCREASE_ART_COVERAGE)) + 
      geom_bar(aes(fill = COUNTERFACTUAL_LABEL), stat = 'identity', position = 'identity') + 
      labs(x = 'Age source', y = 'Percentage point increase in ART\ncoverage among male sources', fill = counterfactual_label) + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.position = c(0.82, 0.78), 
            legend.key.size = unit(0.4, 'cm'),
            legend.direction = 'vertical', 
            # legend.title = element_blank()
      ) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
      scale_x_continuous(expand = c(0,0)) 
    
    p2 <- ggplot() + 
      geom_line(data = tmp2, aes(x = AGE_INFECTION.RECIPIENT, y = M, linetype = SEX), col = 'black') + 
      geom_ribbon(data = tmp2, aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, group = SEX), alpha = 0.5) + 
      geom_line(data = tmp3, aes(x = AGE_INFECTION.RECIPIENT, y = M, col = COUNTERFACTUAL_LABEL)) + 
      geom_ribbon(data = tmp3, aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, fill = COUNTERFACTUAL_LABEL), alpha = 0.15) + 
      labs(x = 'Age', y = '\nHIV incidence infections', 
           fill = counterfactual_label, col = counterfactual_label) + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            axis.title.x = element_blank(), 
            # axis.text.x = element_blank(), 
            legend.position = c(0.9, 0.86),
            legend.title = element_blank()) +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous() + 
      scale_x_continuous(expand = c(0,0)) +
      guides(color = 'none', fill = 'none')
    
    p3 <- ggplot(tmp, aes(x = AGEYRS)) + 
      geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = COUNTERFACTUAL_LABEL)) + 
      geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, fill = COUNTERFACTUAL_LABEL), alpha = 0.5) + 
      labs(x = 'Age', y = '% reduction in HIV incidence infections\namong female recipients', 
           fill = counterfactual_label, col = counterfactual_label) + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            # axis.title.x = element_blank(), 
            # axis.text.x = element_blank(), 
            legend.position = 'none') +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
      scale_x_continuous(expand = c(0,0)) 
    
    p <- grid.arrange(p1, p2, p3, layout_matrix = rbind(c(NA, 1), c(NA, 2), c(3,3)), widths = c(0.02, 0.95), heights = c(0.33, 0.33, 0.33))
    ggsave(p, file = paste0(outdir, '-output-counterfactual_incidence_', communities[i], '.png'), w = 8, h = 9)
    
  }
}

plot_incidence_transmission <- function(incidence_tranmission, outdir){
  tmp <- copy(incidence_tranmission)
  
  tmp[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 7, 0)]
  tmp[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  
  
  p <- vector(mode = 'list', length = length(tmp[, unique(INDEX_DIRECTION)]))
  for(i in seq_along(tmp[, unique(INDEX_DIRECTION)])){
    if(i == 1){
      cols <- grDevices::colorRampPalette(c("#4C0033", '#790252', '#AF0171', '#E80F88', '#EE6983', "#FFC4C4"))(tmp[, length(unique(AGE_GROUP_TRANSMISSION.SOURCE))])
      # cols <- gradient_color(c("#4C0033", "#E80F88"))(tmp[, length(unique(AGE_GROUP_TRANSMISSION.SOURCE))])
      lab <- 'Female sources'
    }else{
      cols <- grDevices::colorRampPalette(c("#002B5B", '#0080bf', '#00acdf', '#55d0ff', '#7ce8ff'))(tmp[, length(unique(AGE_GROUP_TRANSMISSION.SOURCE))])
      lab <- 'Male sources'
    }
    cols <- rev(cols)
    
    tmp1 <- unique(tmp[, .(INDEX_ROUND2, LABEL_ROUND2)])
    
    p[[i]] <- ggplot(tmp[INDEX_DIRECTION == i], aes(x = INDEX_ROUND2, group = AGE_GROUP_TRANSMISSION.SOURCE)) +
      geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0, size = 0.5, position = position_dodge(width = 0.2))  +
      geom_line(aes(y = M, col = AGE_GROUP_TRANSMISSION.SOURCE), position = position_dodge(width = 0.2)) +
      geom_point(aes(y = M, col = AGE_GROUP_TRANSMISSION.SOURCE), size= 1.7, position = position_dodge(width = 0.2)) + 
      labs(x = '', y = 'HIV incident transmissions per person-year\nrelative to first round', col = lab) + 
      theme_bw() +
      facet_grid(.~LABEL_COMMUNITY, scale = 'free_x') + 
      scale_color_manual(values = cols) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            axis.text.x = element_text(angle = 40, hjust = 1),
            axis.title.x = element_blank(), 
            legend.position = 'bottom',
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.margin = margin())  + 
      scale_x_continuous(labels = tmp1[order(INDEX_ROUND2), (LABEL_ROUND2)], breaks = tmp1[order(INDEX_ROUND2), unique(INDEX_ROUND2)]) + 
      scale_y_continuous(labels = scales::percent, limits =  tmp[, range(CL, CU)])  
      # coord_cartesian(ylim = c(NA, 3.50))
    
  }
  
  pp <- ggarrange(plotlist = p, ncol = 1, legend = 'bottom')
  ggsave(pp, file = paste0(outdir, '-output-', 'Incidencetransmission', '_sex.png'), w = 7, h = 8)
  
}

plot_incidence_infection <- function(incidence_infection, outdir){
  tmp <- copy(incidence_infection)
  
  tmp[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 7, 0)]
  tmp[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  
  p <- vector(mode = 'list', length = length(tmp[, unique(INDEX_DIRECTION)]))
  for(i in seq_along(tmp[, unique(INDEX_DIRECTION)])){
    if(i == 2){
      cols <- grDevices::colorRampPalette(c("#4C0033", '#790252', '#AF0171', '#E80F88', '#EE6983', "#FFC4C4"))(tmp[, length(unique(AGE_GROUP_INFECTION.RECIPIENT))])
      lab <- 'Female recipients'
    }else{
      cols <- grDevices::colorRampPalette(c("#002B5B", '#0080bf', '#00acdf', '#55d0ff', '#7ce8ff'))(tmp[, length(unique(AGE_GROUP_INFECTION.RECIPIENT))])
      lab <- 'Male recipients'
    }
    cols <- rev(cols)
    
    tmp1 <- unique(tmp[, .(INDEX_ROUND2, LABEL_ROUND2)])
    
    p[[i]] <- ggplot(tmp[INDEX_DIRECTION == i], aes(x = INDEX_ROUND2, group = AGE_GROUP_INFECTION.RECIPIENT)) +
      geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0, size = 0.5, position = position_dodge(width = 0.2))  +
      geom_line(aes(y = M, col = AGE_GROUP_INFECTION.RECIPIENT), position = position_dodge(width = 0.2)) +
      geom_point(aes(y = M, col = AGE_GROUP_INFECTION.RECIPIENT), size= 1.7, position = position_dodge(width = 0.2)) + 
      labs(x = '', y = 'HIV incident infections per person-year\nrelative to first round', col = lab) + 
      theme_bw() +
      facet_grid(.~LABEL_COMMUNITY, scale = 'free_x') + 
      scale_color_manual(values = cols) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            axis.text.x = element_text(angle = 40, hjust = 1),
            axis.title.x = element_blank(), 
            legend.position = 'bottom',
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.margin = margin())  + 
      scale_x_continuous(labels = tmp1[order(INDEX_ROUND2), (LABEL_ROUND2)], breaks = tmp1[order(INDEX_ROUND2), unique(INDEX_ROUND2)]) + 
      scale_y_continuous(labels = scales::percent, limits =  tmp[, range(CL, CU)])  
      # coord_cartesian(ylim = c(NA, 2.00))
    
   
  }
  
  pp <- ggarrange(plotlist = p, ncol = 1, legend = 'bottom')
  ggsave(pp, file = paste0(outdir, '-output-', 'Incidenceinfection', '_sex.png'), w = 7, h = 8)
  
}

plot_NNT <- function(NNT, outdir){
  
  # select round 
  Round <- 'R018'
  nnt <- NNT[ROUND == Round & IS_MF == T]
  
  # select age
  Ageyrs <- seq(20, 49, 5)
  nnt <- nnt[AGE_TRANSMISSION.SOURCE %in% Ageyrs]
  nnt[, AGE_TRANSMISSION.SOURCE := as.character(AGE_TRANSMISSION.SOURCE)]
  
  # label
  counterfactual_label <- 'Age of treated males'
  
  # plot
  p <- ggplot(nnt, aes(x = AGE_INFECTION.RECIPIENT)) + 
    geom_line(aes(y = M, col = AGE_TRANSMISSION.SOURCE)) + 
    geom_ribbon(aes(ymin= CL, ymax = CU, fill = AGE_TRANSMISSION.SOURCE), alpha = 0.1) + 
    labs(x = 'Age recipient', y = 'Number of HIV-positive male needed to treat\nto prevent one infection in female', 
         fill = counterfactual_label, col = counterfactual_label) + 
    theme_bw() +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          # axis.title.x = element_blank(), 
          # axis.text.x = element_blank(), 
          legend.position = 'bottom') +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    scale_x_continuous(expand = c(0,0)) + 
    facet_grid(LABEL_COMMUNITY~.)
  ggsave(p, file = paste0(outdir, '-output-counterfactual_NNT.png'), w = 6, h = 7)
  
  p <- p + scale_y_log10() 
  ggsave(p, file = paste0(outdir, '-output-counterfactual_NNT_log.png'), w = 6, h = 7)
}

plot_NNT_group <- function(NNT_grouped, outdir){
  
  # select round 
  Round <- 'R018'
  nnt <- NNT_grouped[ROUND == Round & IS_MF == T]

  # select age
  Ageyrs <- seq(20, 49, 5)
  nnt <- nnt[AGE_TRANSMISSION.SOURCE %in% Ageyrs]
  nnt[, AGE_TRANSMISSION.SOURCE := as.character(AGE_TRANSMISSION.SOURCE)]
  
  # label
  counterfactual_label <- 'Age of treated males'
  
  # plot
  p <- ggplot(nnt, aes(x = AGE_GROUP_INFECTION.RECIPIENT)) + 
    geom_errorbar(aes(ymin= CL, ymax = CU, group = AGE_TRANSMISSION.SOURCE), alpha = 0.5, position = position_dodge(0.5), width = 0.2) + 
    geom_point(aes(y = M, col = AGE_TRANSMISSION.SOURCE), position = position_dodge(0.5)) + 
    labs(x = 'Age recipient', y = 'Number of HIV-positive male needed to treat\nto prevent one infection in female', 
         col = counterfactual_label) + 
    theme_bw() +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          # axis.title.x = element_blank(), 
          # axis.text.x = element_blank(), 
          legend.position = 'bottom') +
    scale_color_viridis_d() +
    facet_grid(LABEL_COMMUNITY~.)
  ggsave(p, file = paste0(outdir, '-output-counterfactual_NNT_group.png'), w = 6, h = 7)
  
  p <- p + scale_y_log10() 
  ggsave(p, file = paste0(outdir, '-output-counterfactual_NNT_group_log.png'), w = 6, h = 7)
}

