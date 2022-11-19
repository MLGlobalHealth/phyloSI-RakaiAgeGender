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
    w = 20
  }else if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~ROUND)
    w = 20
  }else if(all(c('COMM', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~COMM)
    w = 9
  }else if('LABEL_DIRECTION' %in% names(tmp)){
    p <- p + facet_grid(LABEL_DIRECTION~.)
    w = 9
  }
  
  
  ggsave(p, file = paste0(outdir, '-output-contrast_2D_', name, '.png'), w = w, h = 9)
  
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
    w = 20
  }else if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~ROUND)
    w = 20
  }else if(all(c('COMM', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~COMM)
    w = 9
  }else if('LABEL_DIRECTION' %in% names(tmp)){
    p <- p + facet_grid(LABEL_DIRECTION~.)
    w = 9
  }
  
  ggsave(p, file = paste0(outdir, '-output-contrast_source_', name, '.png'), w = w, h = 7)
  
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
    w = 20
  }else if(all(c('ROUND', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~ROUND)
    w = 20
  }else if(all(c('COMM', 'LABEL_DIRECTION') %in% names(tmp))){
    p <- p + facet_grid(LABEL_DIRECTION~COMM)
    w = 9
  }else if('LABEL_DIRECTION' %in% names(tmp)){
    p <- p + facet_grid(LABEL_DIRECTION~.)
    w = 9
  }
  
  ggsave(p, file = paste0(outdir, '-output-contrast_recipient_', name, '.png'), w = w, h = 7)
  
} 

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
    
    p <- ggplot(tmp, aes(x = AGE_TRANSMISSION.SOURCE, y = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      geom_point(data = tmp1[count > 0], aes(size = count), col = 'grey50') +
      theme_bw() + 
      labs(y = 'Age at infection recipient', fill = 'Estimated median\nobserved transmission\nevent rate', 
           x= 'Age at transmission source',size='Pairs\ncount') +
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
    
    p <- ggplot(tmp, aes(x = AGE_TRANSMISSION.SOURCE, y = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(y = 'Age at infection recipient', fill = 'Estimated median transmission rate\nper year', 
           x= 'Age at transmission source',size='Pairs\ncount') +
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
    
    if(communities[i] == 'inland'){
      # inland has more round and need a wider figure
      ggsave(p, file = paste0(outdir, '-output-intensity_transmission_by_round_',  communities[i], '.png'), w = 18, h = 7)
      
    }else{
      ggsave(p, file = paste0(outdir, '-output-intensity_transmission_by_round_',  communities[i], '.png'), w = 12, h = 7)
      
    }
  }
  
} 

plot_PPC_observed_source <- function(predict_y, count_data, outdir){
  
  # sum count of observed transmission events across recipient 
  # to find observed transmission events by age of the source
  data <- count_data[, list(count = sum(count)), by = c('LABEL_SOURCE', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_TRANSMISSION.SOURCE', 'PERIOD_SPAN', 'LABEL_DIRECTION')]
  
  communities <- predict_y[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_y[ COMM == communities[i]]
    tmp1 <- data[ LABEL_COMMUNITY == tmp[, unique(LABEL_COMMUNITY)]]
    
    # plot1
    p <- ggplot(tmp, aes( x = AGE_TRANSMISSION.SOURCE)) + 
      geom_line(aes(y = M, linetype = 'Fit')) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, linetype = 'Fit'), alpha = 0.5) + 
      geom_point(data = tmp1, aes(y = count, size = 'Data'), col = 'darkred') +
      theme_bw() + 
      labs(x = 'Age', y = 'Observed transmission events\nfrom phylogenetic analysis', 
           size = '', linetype = '') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(PERIOD~LABEL_SOURCE) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      guides(size = guide_legend(order = 1))
    
    ggsave(p, file = paste0(outdir, '-output-PPC_observed_source_', communities[i], '.png'), w = 7, h = 5)
    ggsave(p, file = paste0(outdir, '-output-PPC_observed_source_', communities[i], '.pdf'), w = 7, h = 5)
    
    # plot2
    df <- merge(tmp, tmp1, by = c('LABEL_DIRECTION', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_TRANSMISSION.SOURCE', 'PERIOD_SPAN'))
    df <- merge(df, unique(df_age_aggregated[, .(AGE_TRANSMISSION.SOURCE, AGE_GROUP_TRANSMISSION.SOURCE)]), by = c('AGE_TRANSMISSION.SOURCE'))
    df[, AGE_GROUP_TRANSMISSION.SOURCE := paste0('Age: ', AGE_GROUP_TRANSMISSION.SOURCE)]
    set.seed(12)
    df[, jitter := runif(length(count), 0, 1), by= c('AGE_TRANSMISSION.SOURCE', 'LABEL_DIRECTION', 'PERIOD', 'count')]
    df[, count_jitter := count + jitter]
    df[, CL_jitter := CL + jitter]
    df[, CU_jitter := CU + jitter]
    df[, M_jitter := M + jitter]
    df[, LABEL_DIRECTION2 := 'Men -> Women']
    df[LABEL_DIRECTION == 'Female -> Male', LABEL_DIRECTION2 := 'Women -> Men']
    
    p <- ggplot(df) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
      geom_errorbar(aes(x=count_jitter, ymin=CL_jitter, ymax=CU_jitter),  color = 'grey50', width = 0, size = 0.5)+
      geom_point(aes(y=M_jitter, x=count_jitter, color=PERIOD), size = 1) + 
      theme_bw() + 
      facet_grid(.~LABEL_DIRECTION2) +
      labs(x = 'Observed transmission events\nin RCCS participants',
           y = 'Predicted transmission events\nin RCCS participants', 
           col ='', fill = '') +
      scale_color_viridis_d(option = 'A', end = 0.9, begin = 0.1) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(0.9)),
            legend.position = 'bottom') +
      guides(color = guide_legend(byrow = T, nrow = 4))
    ggsave(p, file = paste0(outdir, '-output-PPC_observed_point_source_', communities[i], '.pdf'), w = 6, h = 4.8)
    
  }
  
}

plot_PPC_observed_recipient <- function(predict_y, count_data, outdir){
  
  # sum count of observed transmission events across source 
  # to find observed transmission events by age of the recipient
  data <- count_data[, list(count = sum(count)), by = c('LABEL_RECIPIENT', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_INFECTION.RECIPIENT', 'PERIOD_SPAN')]
  
  communities <- predict_y[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_y[ COMM == communities[i]]
    tmp1 <- data[ LABEL_COMMUNITY == tmp[, unique(LABEL_COMMUNITY)]]
    
    # plot 1
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M, linetype = 'Fit')) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, linetype = 'Fit'), alpha = 0.5) + 
      geom_point(data = tmp1, aes(y = count, size = 'Data'), col = 'darkred') +
      theme_bw() + 
      labs(x = 'Age', y = 'Observed transmission events\nfrom phylogenetic analysis', 
           size ='', linetype = '') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_RECIPIENT~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      guides(size = guide_legend(order = 1))
    ggsave(p, file = paste0(outdir, '-output-PPC_observed_recipient_', communities[i], '.png'), w = 7, h = 5)

    # plot 2
    df <- merge(tmp, tmp1, by = c('LABEL_RECIPIENT', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_INFECTION.RECIPIENT', 'PERIOD_SPAN'))
    df <- merge(df, unique(df_age_aggregated[, .(AGE_INFECTION.RECIPIENT, AGE_GROUP_INFECTION.RECIPIENT)]), by = c('AGE_INFECTION.RECIPIENT'))
    df[, AGE_GROUP_INFECTION.RECIPIENT := paste0('Age: ', AGE_GROUP_INFECTION.RECIPIENT)]
    set.seed(12)
    df[, jitter := runif(length(count), 0, 0.5), by= c('AGE_INFECTION.RECIPIENT', 'LABEL_RECIPIENT', 'PERIOD', 'count')]
    df[, count_jitter := count + jitter]
    df[, CL_jitter := CL + jitter]
    df[, CU_jitter := CU + jitter]
    df[, M_jitter := M + jitter]
    
    p <- ggplot(df) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
      geom_errorbar(aes(x=count_jitter, ymin=CL_jitter, ymax=CU_jitter),  color = 'grey50', width = 0, size = 0.5)+
      geom_point(aes(y=M_jitter, x=count_jitter, color=PERIOD), size = 1) + 
      theme_bw() + 
      labs(x = 'Observed transmission events\nin RCCS participants',
           y = 'Predicted transmission events\nin RCCS participants', 
           col ='', fill = '') +
        scale_color_viridis_d(option = 'A', end = 0.9, begin = 0.1) + 
      facet_wrap(AGE_GROUP_INFECTION.RECIPIENT~LABEL_RECIPIENT, scale = 'free', ncol = 3) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(0.9)),
            legend.position = 'bottom') +
      guides(size = guide_legend(order = 1))
    ggsave(p, file = paste0(outdir, '-output-PPC_observed_point_recipient_', communities[i], '.pdf'), w = 6, h = 5.2)
    
  }
  
}

plot_PPC_augmented_recipient_round <- function(predict_z_recipient_round, incidence_cases_recipient_round, outdir){
  
  predict_z <- merge(predict_z_recipient_round, incidence_cases_recipient_round[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, AGE_INFECTION.RECIPIENT, INCIDENT_CASES, INCIDENT_CASES_UB, INCIDENT_CASES_LB)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'))
  
  communities <- predict_z[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_z[ COMM == communities[i]]
    
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M, linetype = 'Fit')) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, linetype = 'Fit'), alpha = 0.5) + 
      geom_point(aes(y = INCIDENT_CASES, size = 'Data'), col = 'darkred') +
      geom_errorbar(aes(ymax = INCIDENT_CASES_UB, ymin = INCIDENT_CASES_LB), col = 'darkred', width = 0.2) +
      theme_bw() + 
      labs(x = 'Age', y = 'Transmission events', size = '', linetype = '') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_ROUND~LABEL_RECIPIENT) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      guides(size = guide_legend(order = 1))
    
    if(communities[i] == 'inland'){
      # inland has more round than fishing
      ggsave(p, file = paste0(outdir, '-output-PPC_augmented_recipient_byround_', communities[i], '.png'), w = 7, h = 14)
      ggsave(p, file = paste0(outdir, '-output-PPC_augmented_recipient_byround_', communities[i], '.pdf'), w = 7, h = 14)
      
    }else{
      ggsave(p, file = paste0(outdir, '-output-PPC_augmented_recipient_byround_', communities[i], '.png'), w = 6, h = 10)
      ggsave(p, file = paste0(outdir, '-output-PPC_augmented_recipient_byround_', communities[i], '.pdf'), w = 6, h = 10)
    }
    
  }
  
}

plot_PPC_incidence_rate_round <- function(predict_incidence_rate_round, incidence_cases_recipient_round, outdir){
  
  predict_z <- merge(predict_incidence_rate_round, incidence_cases_recipient_round[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, 
                                                                                       AGE_INFECTION.RECIPIENT, INCIDENCE, LB, UB)], 
                     by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'))
  
  communities <- predict_z[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_z[ COMM == communities[i]]
    
    # plot 1
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M*100, linetype = 'Fit')) +
      geom_ribbon(aes(ymin = CL*100, ymax = CU*100, linetype = 'Fit'), alpha = 0.5) +
      geom_point(aes(y = INCIDENCE*100, size= 'Data'), col = 'darkred') +
      # geom_errorbar(aes(ymax = UB*100, ymin = LB*100), col = 'darkred', width = 0.2) +
      theme_bw() + 
      labs(x = 'Age', y = 'Incidence rate per 100 PY', linetype = '', size = '') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_ROUND~LABEL_RECIPIENT) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      guides(size=guide_legend(order=1))
    
    if(communities[i] == 'inland'){
      # inland has more round than fishing
      ggsave(p, file = paste0(outdir, '-output-PPC_incidencerate_perPY_recipient_byround_', communities[i], '.png'), w = 7, h = 14)

    }else{
      ggsave(p, file = paste0(outdir, '-output-PPC_incidencerate_perPY_recipient_byround_', communities[i], '.png'), w = 6, h = 10)
    }
    
    # plot 2
    df <- merge(tmp, unique(df_age_aggregated[, .(AGE_INFECTION.RECIPIENT, AGE_GROUP_INFECTION.RECIPIENT)]), by = c('AGE_INFECTION.RECIPIENT'))
    df[, AGE_GROUP_INFECTION.RECIPIENT := paste0('Age: ', AGE_GROUP_INFECTION.RECIPIENT)]
    
    set.seed(12)
    df[, jitter := runif(length(INCIDENCE), 0, 0.005), by= c('AGE_GROUP_INFECTION.RECIPIENT', 'LABEL_RECIPIENT', 'ROUND', 'INCIDENCE')]
    df[, INCIDENCE_jitter := (INCIDENCE + jitter)*100]
    df[, CL_jitter := (CL + jitter)*100]
    df[, CU_jitter := (CU + jitter)*100]
    df[, M_jitter := (M + jitter)*100]
    df[, SEX_LABEL := gsub('(.+) recipients', '\\1', LABEL_RECIPIENT)]
    df[, SEX_LABEL := ifelse(SEX_LABEL=='Male', 'Men', 'Women')]
    
    p <- ggplot(df) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
      geom_errorbar(aes(x=INCIDENCE_jitter, ymin=CL_jitter, ymax=CU_jitter),  color = 'grey70', width = 0, size = 0.2)+
      geom_point(aes(y=M_jitter, x=INCIDENCE_jitter, 
                     color=LABEL_ROUND), size = 1) + 
      theme_bw() + 
      labs(x = 'Prior mean incidence rates\nper 100 person-years',
           y = 'Estimated incidence rates\nper 100 person-years', 
           col ='', fill = '') +
      scale_color_manual(values = palette_round_inland) + 
      facet_grid(.~SEX_LABEL) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(0.9)),
            legend.position = 'bottom')  +
      guides(color = guide_legend(byrow = T, nrow = 3))
    ggsave(p, file = paste0(outdir, '-output-PPC_incidencerate_perPY_point_recipient_byround_', communities[i], '.pdf'), w = 5.5, h =4.5)
    
  }
  
}

plot_observed_to_augmented <- function(predict_y, predict_z, outdir){
  
  predict_z[, type := 'All']
  predict_y[, type := 'Observed']
  
  predict_df <- rbind(predict_z, predict_y)
  
  communities <- predict_df[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_df[ COMM == communities[i]]
    
    p <- ggplot(tmp, aes( x = AGE_TRANSMISSION.SOURCE)) + 
      geom_line(aes(y = M, col = type)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.5) + 
      theme_bw() + 
      labs(x = 'Age', y = 'Transmission events') +
      facet_grid(PERIOD~LABEL_SOURCE) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    ggsave(p, file = paste0(outdir, '-output-observed_vs_augmented_source_', communities[i], '.png'), w = 7, h = 7)
  }
  
}

plot_force_infection <- function(force_infection, outdir, lab = NULL){
  
  communities <- force_infection[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- force_infection[ COMM == communities[i]]
    p <- ggplot(tmp, aes(x = AGE_TRANSMISSION.SOURCE, y = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(y = 'Age at infection recipient', fill = paste0('Estimated median\n', lab, ' force infection'), 
           x= 'Age at transmission source',size='Pairs\ncount') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~LABEL_ROUND) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) 
    
    if(communities[i] == 'inland'){
      # inland has more round than fishing
      ggsave(p, file = paste0(outdir, '-output-', lab, 'FOI_2D_',  communities[i], '.png'), w = 18, h = 7)
      
    }else{
      ggsave(p, file = paste0(outdir, '-output-', lab, 'FOI_2D_',  communities[i], '.png'), w = 10, h = 7)
      
    }
  }
  
} 

plot_force_infection_sex_source <- function(force_infection_sex_source, outdir){
  
  communities <- force_infection_sex_source[, unique(COMM)]
  
  for(i in seq_along(communities)){
    tmp <- force_infection_sex_source[COMM == communities[i]]
    
    p <- ggplot(tmp, aes(x = LABEL_ROUND)) + 
      geom_bar(aes(y = M, fill = LABEL_SOURCE), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = LABEL_SOURCE), position = position_dodge()) + 
      labs(x = '', y = 'Force of infection exerted', fill = '') + 
      theme_bw() +
      facet_grid(LABEL_COMMUNITY~., scale = 'free')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.text.x = element_text(angle= 50, hjust = 1)) +
      ggsci::scale_fill_npg() 
    
    if(communities[i] == 'inland'){
      # inalnd has more round than fishing
      ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_sex_', communities[i], '.png'), w = 7, h = 5)
      
    }else{
      ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_sex_', communities[i], '.png'), w = 5, h = 5)
    }

  }
}

plot_force_infection_sex_age_source <- function(force_infection_age_source, outdir){

  communities <- force_infection_age_source[, unique(COMM)]
  
  for(i in seq_along(communities)){
    tmp <- force_infection_age_source[COMM == communities[i]]
    
    p <- ggplot(tmp, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M, fill = LABEL_SOURCE), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = LABEL_DIRECTION), position = position_dodge()) + 
      labs(x = 'Age', y = 'Force of infection exerted', fill = '') + 
      theme_bw() +
      facet_grid(ROUND~., scale = 'free')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg() 
    
    if(communities[i] == 'inland'){
      ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_source_', communities[i], '.png'), w = 7, h = 14)
    }else{
      ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_source_', communities[i], '.png'), w = 7, h = 9)
      
    }

  }

}

plot_force_infection_sex_age_recipient <- function(force_infection_age_recipient, outdir){
  
  communities <- force_infection_age_recipient[, unique(COMM)]
  
  for(i in seq_along(communities)){
    tmp <- force_infection_age_recipient[COMM == communities[i]]
    
    p <- ggplot(tmp, aes(x = AGE_INFECTION.RECIPIENT)) + 
      geom_bar(aes(y = M, fill = LABEL_RECIPIENT), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = LABEL_RECIPIENT), position = position_dodge()) + 
      labs(x = 'Age', y = 'Force of infection received', fill = '') + 
      theme_bw() +
      facet_grid(ROUND~., scale = 'free')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggsci::scale_fill_npg() 
    
    if(communities[i] == 'inland'){
      ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_recipient_', communities[i], '.png'), w = 7, h = 14)
    }else{
      ggsave(p, file = paste0(outdir, '-output-', 'FOI', '_age_recipient_', communities[i], '.png'), w = 7, h = 9)
      
    }
    
  }

}

plot_contribution_sex_source <- function(contribution_sex_source, unsuppressed_prop_sex, prevalence_prop_sex,outdir, lab = NULL){
  
  # y axis label
  type_cont <- 'Contribution from male sources\nto HIV infection'
  
  # contribution from output
  tmp <- copy(contribution_sex_source)
  tmp[, type  := type_cont]

  # keep only male contributions
  tmp <- tmp[LABEL_DIRECTION == 'Male -> Female']
  
  # add share of male among infected and infected unsuppressed
  prevalence_prop_sex[, type  := 'Share of males among HIV-positive\nindividuals']
  unsuppressed_prop_sex[, type  := 'Share of males among HIV-positive\nunsuppressed individuals']
  tmp1 <- rbind(prevalence_prop_sex, unsuppressed_prop_sex, fill=TRUE)
  tmp1[, type := factor(type , levels = c(unique(prevalence_prop_sex$type), unique(unsuppressed_prop_sex$type)))]
  tmp1 <- tmp1[LABEL_DIRECTION == 'Male -> Female']
  
  # index round such that x axis of inland and fishing are not the same
  tmp[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 9, 0)]
  tmp[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  tmp1[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 9, 0)]
  
  # combine
  tmp2 <- rbind(tmp, tmp1, fill=TRUE)
  tmp2[, type := factor(type , levels = c(unique(tmp$type), levels(tmp1$type)))]
  
  # plot
  p <-   ggplot(tmp2, aes(x = INDEX_ROUND2, group = type)) +
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
     scale_x_continuous(labels = tmp[order(INDEX_ROUND2), (LABEL_ROUND2)], 
                        breaks = tmp[order(INDEX_ROUND2), unique(INDEX_ROUND2)]) + 
     scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0)), limits = c(0, 1)) + 
     guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1), shape = guide_legend(order = 1))
   
   if(is.null(lab)) lab =  'Contribution'
  ggsave(p, file = paste0(outdir, '-output-', gsub(' ', '_', lab), '_sex.png'), w = 7, h = 4.5)
  
}

plot_contribution_age_source_unsuppressed <- function(contribution_age_source, unsuppressed_prop_age, outdir, lab = NULL){
  
  # restricted rounds for main 
  Rounds.all <- list('inland' = paste0('R0', c(15,18)), 'fishing' = paste0('R0', c(15,18)))
  
  # y label
  type_cont <- 'Contribution to HIV incidence'
  
  # contribution output
  cont <- copy(contribution_age_source)
  cont[, type  := type_cont]

  # find median age contribution for triangles
  median_age <- copy(cont)
  median_age[, WEIGHT_CONTRIBUTION := M / sum(M), by = c('LABEL_SOURCE', 'LABEL_ROUND', 'COMM', 'INDEX_ROUND', 'ROUND')]
  median_age <- median_age[, list(AGE_MEDIAN_CONTRIBUTION = matrixStats::weightedMedian(AGE_TRANSMISSION.SOURCE, WEIGHT_CONTRIBUTION )), 
                       by = c('LABEL_SOURCE', 'LABEL_ROUND', 'COMM', 'INDEX_ROUND', 'ROUND')]
    
  # find unsuppressed share by sex and age
  uns <- copy(unsuppressed_share_sex_age)
  setnames(uns, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
  uns[, type := 'Share among HIV-positive\nunsuppressed individuals']

  # prepare plotting function
  plot.p <- function(cont.p, uns.p, cont_b.p, median_age.p, level_y){
    ggplot(cont.p, aes(x = AGE_TRANSMISSION.SOURCE)) +
      geom_ribbon(data = uns.p, aes(ymin = CL, ymax = CU, linetype=type), alpha = 0.3, fill='grey50')+
      geom_line(data = uns.p, aes(y = M,linetype = type)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = LABEL_SOURCE, size = type), alpha = 0.6) + 
      geom_line(aes(y = M, col = LABEL_SOURCE, size = type), stat = 'identity', position = "identity") + 
      scale_color_manual(values = c('Male sources'='royalblue3','Female sources'='deeppink')) + 
      scale_fill_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) + 
      new_scale_fill() +
      new_scale_color() +
      geom_point(data = select(median_age.p[INDEX_ROUND == min(INDEX_ROUND)], -'LABEL_ROUND'), aes(x = AGE_MEDIAN_CONTRIBUTION, y = level_y, col = LABEL_SOURCE, fill = LABEL_SOURCE), size = 3,  shape = 25) + 
      scale_color_manual(values = c('Male sources'='paleturquoise4','Female sources'='pink4')) +
      scale_fill_manual(values = c('Male sources'='paleturquoise4','Female sources'='pink4')) + 
      new_scale_color() +
      new_scale_fill() +
      geom_point(data = median_age.p, aes(x = AGE_MEDIAN_CONTRIBUTION, y = level_y, col = LABEL_SOURCE, fill = LABEL_SOURCE), size = 3,  shape = 25) + 
      scale_color_manual(values = c('Male sources'='royalblue3','Female sources'='deeppink')) + 
      scale_fill_manual(values = c('Male sources'='royalblue3','Female sources'='deeppink')) + 
      new_scale_color() +
      new_scale_fill() +
      geom_point(data = median_age.p[INDEX_ROUND == min(INDEX_ROUND)], aes(x = AGE_MEDIAN_CONTRIBUTION, y = level_y, col = LABEL_SOURCE, fill = LABEL_SOURCE), size = 3,  shape = 25) + 
      scale_color_manual(values = c('Male sources'='paleturquoise4','Female sources'='pink4')) +
      scale_fill_manual(values = c('Male sources'='paleturquoise4','Female sources'='pink4')) + 
      new_scale_color() +
      geom_line(data = cont_b.p, aes(y = M, col = LABEL_SOURCE, alpha=type), linetype ='solid') +  
      scale_color_manual(values = c('Male sources'='paleturquoise4','Female sources'='pink4')) +
      labs(x = 'Age of the source', y = 'Percent') + 
      theme_bw() +
      facet_grid(LABEL_ROUND~LABEL_SOURCE) +
      # scale_color_manual(values = c('#FF4949')) + 
      scale_alpha_manual(values = 1) + 
      scale_size_manual(values = 0.5) + 
      # scale_fill_manual(values = c('Male source'='#C6DCE4','Female source'='#F2D1D1')) + 
      scale_linetype_manual(values  = 'dashed') +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'none', 
            legend.title = element_blank(), 
            panel.grid.minor = element_blank()) + 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05)), limits = c(0,NA))+ 
      scale_x_continuous(expand = c(0,0),
                         breaks = c(seq(min(cont.p[, unique(AGE_TRANSMISSION.SOURCE)]), max(cont.p[, unique(AGE_TRANSMISSION.SOURCE)]), 5))) 
  }
  
  
  communities <- cont[, unique(COMM)]
  
  # make plots
  for(i in seq_along(communities)){
    
    cont.p <- cont[COMM == communities[i]]
    uns.p <- uns[COMM == communities[i]]
    median_age.p <- median_age[COMM == communities[i]]
    Rounds <- Rounds.all[[communities[i]]]
    
    # find contribution in the first round 
    cont_b.p <- cont.p[ROUND %in% Rounds]
    cont_b.p[, min_INDEX_ROUND :=  min(INDEX_ROUND), by = 'COMM']
    cont_b.p <- cont_b.p[INDEX_ROUND == min_INDEX_ROUND]
    cont_b.p[, type := paste0(type, '\nin ', gsub('\n', ', ', LABEL_ROUND))]
    cont_b.p <- select(cont_b.p, -LABEL_ROUND)
    
    # plot all rounds
    p.all <- plot.p(cont.p, uns.p, cont_b.p, median_age.p, 0.003) # last argument tune the level of the triangle
      
    # plot selected rounds
    p <- plot.p(cont.p[ROUND %in% Rounds], uns.p[ROUND %in% Rounds], cont_b.p, median_age.p[ROUND %in% Rounds], 0.003)
    
    # plot legend
    p_legend <- ggplot(cont.p[LABEL_SOURCE == 'Female sources'], aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = type), alpha = 0.6) + 
      geom_line(aes(y = M, col = type), stat = 'identity', position = "identity") + 
      geom_line(data = uns.p, aes(y = M, size = type), col = 'black', linetype = 'dashed') + 
      geom_ribbon(data = uns.p, aes(ymin = CL, ymax = CU, size = type), alpha = 0.3, fill='grey50') + 
      geom_line(data = cont_b.p, aes(y = M, linetype=type), col = 'pink4') + 
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
             linetype = 'none')

    # add legends
    pp.all <- ggarrange(p.all, legend.grob = get_legend(p_legend), legend = 'bottom') + 
      theme(panel.background = element_rect(fill='white'))
    
    pp <- ggarrange(p, legend.grob = get_legend(p_legend), legend = 'bottom') + 
      theme(panel.background = element_rect(fill='white'))
    
    # save
    if(is.null(lab)) lab =  'Contribution'
    
    if(communities[i] == 'inland'){
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_extended_', communities[i], '.pdf'), w = 9.5, h = 14.5)
    }else{
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_extended_', communities[i], '.pdf'), w = 9.5, h = 9)
      
    }
    
    ggsave(pp, file = paste0(outdir, '-output-', lab, '_age_', communities[i], '.pdf'), w = 5.5, h = 4.2)
    
  }
}

plot_contribution_age_source <- function(contribution_age_source, outdir, lab = NULL){
  
  # selected rounds for main
  Rounds.all <- list('inland' = paste0('R0', c(10, 12, 14, 16,18)), 'fishing' = paste0('R0', c(15,18)))
  
  # contribution
  cont <- copy(contribution_age_source)

  # median age of contributors
  median_age <- copy(cont)
  median_age[, WEIGHT_CONTRIBUTION := M / sum(M), by = c('LABEL_SOURCE', 'LABEL_ROUND', 'COMM', 'INDEX_ROUND', 'ROUND')]
  median_age <- median_age[, list(AGE_MEDIAN_CONTRIBUTION = matrixStats::weightedMedian(AGE_TRANSMISSION.SOURCE, WEIGHT_CONTRIBUTION )), 
                       by = c('LABEL_SOURCE', 'LABEL_ROUND', 'COMM', 'INDEX_ROUND', 'ROUND')]
  
  # prepare plotting function
  plot.p <- function(cont.p, median_age.p, level_y){
    ggplot(cont.p, aes(x = AGE_TRANSMISSION.SOURCE)) +
      geom_line(aes(y = M, col = LABEL_SOURCE)) + 
      geom_point(data = median_age.p, aes(x = AGE_MEDIAN_CONTRIBUTION, y = level_y, col = LABEL_SOURCE, fill = LABEL_SOURCE), size = 3, shape = 25) + 
      scale_fill_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) +
      new_scale_fill() + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = LABEL_SOURCE), alpha = 0.5) +
      scale_color_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) +
      scale_fill_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) + 
      labs(x = 'Age', y = '\nContribution to HIV incidence') + 
      theme_bw() +
      facet_grid(LABEL_ROUND~.) +
      theme(strip.background = element_rect(colour="white", fill="white"),
            legend.title = element_blank(), 
            panel.grid.minor = element_blank(), 
            strip.text = element_text(size = 9.3), 
            axis.title = element_text(size = 12)) + 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05)), limits = c(0,NA))+ 
      scale_x_continuous(expand = c(0,0), 
                         breaks = c(seq(min(cont.p[, unique(AGE_TRANSMISSION.SOURCE)]), max(cont.p[, unique(AGE_TRANSMISSION.SOURCE)]), 5))) 
  }
  
  
  # make plots
  communities <- cont[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    cont.p <- cont[COMM == communities[i]]
    median_age.p <- median_age[COMM == communities[i]]
    Rounds <- Rounds.all[[communities[i]]]
    
    # all rounds
    pp.all <- plot.p(cont.p, median_age.p, 0.003) + #last argument tune the level of the triangles
      theme(legend.position =  'bottom')
    
    # selected rounds
    pp <- plot.p(cont.p[ROUND %in% Rounds], median_age.p[ROUND %in% Rounds], 0.0035)+ 
      theme(legend.position =  c(0.78,0.95))
    
    # selected rounds horizontal
    pp.hori <-  plot.p(cont.p[ROUND %in% Rounds], median_age.p[ROUND %in% Rounds], 0.003) + 
      facet_grid(.~LABEL_ROUND) +
      theme(legend.position =  'none', 
            strip.text = element_blank())

    if(is.null(lab)) lab =  'Contribution_sex'
    
    if(communities[i] == 'inland'){
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_extended_', communities[i], '.pdf'), w = 7, h = 14.5)
      ggsave(pp, file = paste0(outdir, '-output-', lab, '_age_', communities[i], '.pdf'), w = 4.5, h = 7.8)
      ggsave(pp.hori, file = paste0(outdir, '-output-', lab, '_age_horizontal_', communities[i], '.pdf'), w = 8, h = 2.6)
    }else{
      ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_extended_', communities[i], '.pdf'), w = 7, h = 9)
      
    }

  }
}

plot_contribution_age_source_sex_ratio <- function(expected_contribution_age_source_sex_ratio, outdir, lab = NULL){
  
  # selected rounds
  Rounds <- list('fishing' = paste0('R0', c(15:18)), 'inland'= paste0('R0', c(10, 12, 14, 16, 18)))

  # plot function
  plot.p <- function(cont.p){
    ggplot(cont.p, aes(x = AGE_TRANSMISSION.SOURCE)) +
      geom_hline(yintercept = 1, linetype = 'dashed', alpha = 0.5) +
      geom_line(aes(y = M, col = LABEL_ROUND)) + 
      # geom_ribbon(aes(ymin = CL, ymax = CU, fill = LABEL_ROUND), alpha = 0.5) +
      labs(x = 'Age of the source', y = 'Male-female ratio contribution to HIV incidence') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.title = element_blank(), 
            panel.grid.minor = element_blank()) + 
      scale_y_log10(expand = expansion(mult = c(0, .05)))+ 
      scale_x_continuous(expand = c(0,0), 
                         breaks = c(seq(min(cont.p[, unique(AGE_TRANSMISSION.SOURCE)]), max(cont.p[, unique(AGE_TRANSMISSION.SOURCE)]), 5), 
                                    max(cont.p[, unique(AGE_TRANSMISSION.SOURCE)]))) 
  }
  
  
  # make plots
  
  communities <- expected_contribution_age_source_sex_ratio[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    cont.p <- expected_contribution_age_source_sex_ratio[COMM == communities[i]]
    Rounds.c <- Rounds[[communities[i]]]
    
    # color palette
    if(communities[i] == 'inland'){
      colors <- palette_round_inland[df_round[COMM == 'inland' ,INDEX_ROUND]]
      colors_reduced <- palette_round_inland[df_round[COMM == 'inland' & ROUND %in% Rounds.c,INDEX_ROUND]]
    }else{
      colors <- palette_round_fishing[df_round[COMM == 'fishing' ,INDEX_ROUND]]
      colors_reduced <- palette_round_fishing[df_round[COMM == 'fishing' & ROUND %in% Rounds.c,INDEX_ROUND]]
    }
    
    # all rounds
    pp.all <- plot.p(cont.p) + 
      theme(legend.position = 'bottom') + 
      scale_color_manual(values =colors) +
      scale_fill_manual(values = colors) +
      guides(color = guide_legend(byrow = T, nrow = 3), fill = guide_legend(byrow = T, nrow = 3))
    
    # all round with CI
    pp.allCI <- pp.all + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = LABEL_ROUND), alpha = 0.5) +
      facet_grid(LABEL_ROUND~.) + 
      theme(legend.position = 'bottom')
    
    # selected rounds
    pp <- plot.p(cont.p[ROUND %in% Rounds.c])+ 
      theme(legend.position =  c(0.8,0.2)) + 
      scale_color_manual(values =colors_reduced) +
      scale_fill_manual(values = colors_reduced) 
    
    #save
    if(is.null(lab)) lab =  'Contribution_sex_ratio'

    if(communities[i] == 'inland'){
      ggsave(pp.allCI, file = paste0(outdir, '-output-', lab, '_age_extended_CI_', communities[i], '.pdf'), w = 6, h = 15)
    }else{
      ggsave(pp.allCI, file = paste0(outdir, '-output-', lab, '_age_extended_CI_', communities[i], '.pdf'), w = 6, h = 9)
      
    }

    ggsave(pp.all, file = paste0(outdir, '-output-', lab, '_age_extended_', communities[i], '.pdf'), w = 6, h = 6)
    ggsave(pp, file = paste0(outdir, '-output-', lab, '_age_', communities[i], '.pdf'), w = 5, h = 4.5)
    
  }
}

plot_contribution_age_ungroup <- function(expected_contribution_age_group_source, outdir, lab = NULL){
  
  # select rounds
  Rounds <- c('R015', 'R018')
  
  # plot
  communities <- expected_contribution_age_group_source[, unique(COMM)]
  
  tmp <- copy(expected_contribution_age_group_source)
  tmp[, AGE_LABEL := paste0('Age recipient:\n', AGE_GROUP_INFECTION.RECIPIENT)]
  
  for(i in seq_along(communities)){
    
    tmp1 <- tmp[ COMM == communities[i] & ROUND %in% Rounds]
    
    if(communities[i] == 'inland'){
      colors <- palette_round_inland[c(6, 9)]
    }else{
      colors <- palette_round_fishing[c(1, 5)]
    }
    
    p <- ggplot(tmp1, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_line(aes(y = M, col = LABEL_ROUND)) + 
      geom_ribbon(aes(ymin = CL, ymax = CU, fill = LABEL_ROUND), alpha = 0.5) + 
      labs(x = 'Age source', y = 'Share in HIV incidence', fill = '', col ='') + 
      theme_bw() +
      facet_grid(AGE_LABEL~LABEL_RECIPIENT)+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      # ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])+ 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) 
    
    if(is.null(lab)) lab =  'Contribution'
    ggsave(p, file = paste0(outdir, '-output-', lab, '_age_ungroup_',  communities[i], '.png'), w = 7,h =9)
  }
  
}


plot_contribution_age_group <- function(contribution_age_group_source, outdir, lab = NULL){
  
  communities <- contribution_age_group_source[, unique(COMM)]
  
  tmp <- copy(contribution_age_group_source)
  tmp[, AGE_LABEL := paste0('Age recipient: ', AGE_GROUP_INFECTION.RECIPIENT)]
  
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
      labs(x = 'Age source', y = 'Share in HIV incidence', fill = '') + 
      theme_bw() +
      facet_grid(AGE_LABEL~LABEL_RECIPIENT)+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) +
      scale_fill_manual(values = colors) +
      # ggtitle(contribution_age_group_source[ COMM == communities[i], unique(LABEL_COMMUNITY)])+ 
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) + 
      guides(fill = guide_legend(byrow= T, nrow = 3))
    
    if(is.null(lab)) lab =  'Contribution'
    ggsave(p, file = paste0(outdir, '-output-', lab, '_age_group_',  communities[i], '.png'), w = 7,h =8)
  }
  
}

plot_contribution_age_classification <- function(contribution_age_classification_source, outdir, lab = NULL){
  
  communities <- contribution_age_classification_source[, unique(COMM)]
  
  tmp <- copy(contribution_age_classification_source)
  tmp[, AGE_CLASSIFICATION.SOURCE := factor(AGE_CLASSIFICATION.SOURCE, levels  = c('Younger', 'Same age', 'Older'))]
  tmp[, AGE_LABEL := paste0('Age recipient: ', AGE_GROUP_INFECTION.RECIPIENT)]

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
      labs(x = 'Age classification of the source', y = 'Share in HIV incidence', fill = '') + 
      theme_bw() +
      facet_grid(AGE_LABEL~LABEL_RECIPIENT)+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank()) +
      scale_fill_manual(values = colors) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05))) + 
      guides(fill = guide_legend(byrow= T, nrow = 3))
    
    if(is.null(lab)) lab =  'Contribution'
    ggsave(p, file = paste0(outdir, '-output-', lab, '_age_classification_',  communities[i], '.png'), w = 7, h = 8)
  }
  
}

plot_transmission_risk_sex_source <- function(transmission_risk_sex_source_round, outdir){
  
  tmp <- copy(transmission_risk_sex_source_round)
  
  p <- ggplot(tmp, aes(x = LABEL_ROUND)) + 
    geom_bar(aes(y = M, fill = LABEL_SOURCE), stat = 'identity', position = "identity") + 
    geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.2, col = 'grey40') + 
    labs(x = '', y = 'Transmission risk per unsuppressed per year', 
         col = '', alpha = '', fill ='') + 
    theme_bw() +
    facet_grid(LABEL_SOURCE~LABEL_COMMUNITY, scale = 'free_x')+
    scale_fill_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) + 
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

plot_transmission_risk_sex_age_source <- function(transmission_risk_age_source, outdir){
  
  tran <- copy(transmission_risk_age_source)

  tranf <- tran[round == min(round)]
  tranf[, type := paste0('In ', gsub('\n', ' ', LABEL_ROUND))]
  tranf <- select(tranf, -LABEL_ROUND)
  tranf[, AGE_TRANSMISSION.SOURCE:= AGE_TRANSMISSION.SOURCE - 0.5] # account for the step function
  tmp <- copy(tranf[AGE_TRANSMISSION.SOURCE == max(AGE_TRANSMISSION.SOURCE)])
  tmp[, AGE_TRANSMISSION.SOURCE := AGE_TRANSMISSION.SOURCE + 1]
  tranf <- rbind(tranf, tmp)
  
  communities <- tran[, unique(COMM)]
  
  for(i in seq_along(communities)){
    
    tran.p <- tran[COMM == communities[i]]
    tranf.p <- tranf[COMM == communities[i]]
    
    p <- ggplot(tran.p, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_bar(aes(y = M,fill = LABEL_SOURCE), stat = 'identity', position = "identity") + 
      geom_errorbar(aes(ymin = CL, ymax = CU), width = 0.5, col = 'grey40') + 
      geom_step(data = tranf.p, aes(y = M, color = type), linetype = 'solid') + 
      labs(x = 'Age', y = 'Transmission risk exerted per unsuppressed per year') + 
      theme_bw() +
      facet_grid(LABEL_ROUND~LABEL_SOURCE)+
      scale_color_manual(values = c('#FF4949')) + 
      scale_fill_manual(values = c('Male sources'='lightblue3','Female sources'='lightpink1')) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            legend.title = element_blank()) + 
      scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
      guides(fill = 'none')
    
    if( communities[i] == 'inland'){
      ggsave(p, file = paste0(outdir, '-output-', 'transmission_risk_sex_age_', communities[i], '.png'), w = 8, h = 14)
    }else{
      ggsave(p, file = paste0(outdir, '-output-', 'transmission_risk_sex_age_', communities[i], '.png'), w = 8, h = 9)
      
    }

  }
}

plot_incidence_transmission <- function(incidence_tranmission, outdir){
  
  inc <- copy(incidence_tranmission)
  
  # start index round fishing after inland so they can be plotted side by side
  inc[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 9, 0)]
  inc[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  
  # plot
  directions <- inc[, unique(INDEX_DIRECTION)]
  p <- list()
  
  for(i in seq_along(directions)){
    
    # find color palette
    if(i == 1){
      cols <- grDevices::colorRampPalette(c("#4C0033", '#790252', '#AF0171', '#E80F88', '#EE6983', "#FFC4C4"))(inc[, length(unique(AGE_GROUP_TRANSMISSION.SOURCE))])
      lab <- 'Female sources'
    }else{
      cols <- grDevices::colorRampPalette(c("#002B5B", '#0080bf', '#00acdf', '#55d0ff', '#7ce8ff'))(inc[, length(unique(AGE_GROUP_TRANSMISSION.SOURCE))])
      lab <- 'Male sources'
    }
    cols <- rev(cols)
    
    tmp1 <- unique(inc[, .(INDEX_ROUND2, LABEL_ROUND2)])
    
    p[[i]] <- ggplot(inc[INDEX_DIRECTION == i], aes(x = INDEX_ROUND2, group = AGE_GROUP_TRANSMISSION.SOURCE)) +
      geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0, size = 0.5, position = position_dodge(width = 0.2))  +
      geom_line(aes(y = M, col = AGE_GROUP_TRANSMISSION.SOURCE), position = position_dodge(width = 0.2)) +
      geom_point(aes(y = M, col = AGE_GROUP_TRANSMISSION.SOURCE), size= 1.7, position = position_dodge(width = 0.2)) + 
      labs(x = '', y = 'HIV incidence rate per person-year\nrelative to first round', col = lab) + 
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
      scale_y_continuous(labels = scales::percent, limits =  inc[, range(CL, CU)])  
    # coord_cartesian(ylim = c(NA, 3.50))
    
  }
  
  pp <- ggarrange(plotlist = p, ncol = 1, legend = 'bottom')
  ggsave(pp, file = paste0(outdir, '-output-', 'Incidence_transmission', '_sex.png'), w = 7, h = 8)
  
}

plot_incidence_infection <- function(incidence_infection, outdir){
  
  inc <- copy(incidence_infection)
  
  # start index round fishing after inland so they can be plotted side by side
  inc[, INDEX_ROUND2 := INDEX_ROUND + ifelse(COMM == 'fishing', 9, 0)]
  inc[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
  
  # plot
  p <- list()
  
  for(i in seq_along(inc[, unique(INDEX_DIRECTION)])){
    
    if(i == 2){
      cols <- grDevices::colorRampPalette(c("#4C0033", '#790252', '#AF0171', '#E80F88', '#EE6983', "#FFC4C4"))(inc[, length(unique(AGE_GROUP_INFECTION.RECIPIENT))])
      lab <- 'Female recipients'
    }else{
      cols <- grDevices::colorRampPalette(c("#002B5B", '#0080bf', '#00acdf', '#55d0ff', '#7ce8ff'))(inc[, length(unique(AGE_GROUP_INFECTION.RECIPIENT))])
      lab <- 'Male recipients'
    }
    
    cols <- rev(cols)
    
    tmp1 <- unique(inc[, .(INDEX_ROUND2, LABEL_ROUND2)])
    
    p[[i]] <- ggplot(inc[INDEX_DIRECTION == i], aes(x = INDEX_ROUND2, group = AGE_GROUP_INFECTION.RECIPIENT)) +
      geom_errorbar(aes(ymin = CL, ymax = CU), col = 'grey50', width = 0, size = 0.5, position = position_dodge(width = 0.2))  +
      geom_line(aes(y = M, col = AGE_GROUP_INFECTION.RECIPIENT), position = position_dodge(width = 0.2)) +
      geom_point(aes(y = M, col = AGE_GROUP_INFECTION.RECIPIENT), size= 1.7, position = position_dodge(width = 0.2)) + 
      labs(x = '', y = 'HIV incidence rate per person-year\nrelative to first round', col = lab) + 
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
      scale_y_continuous(labels = scales::percent, limits =  inc[, range(CL, CU)])  
    # coord_cartesian(ylim = c(NA, 2.00))
    
    
  }
  
  pp <- ggarrange(plotlist = p, ncol = 1, legend = 'bottom')
  ggsave(pp, file = paste0(outdir, '-output-', 'Incidence_infection', '_sex.png'), w = 7, h = 8)
  
}

plot_median_age_source <- function(median_age_source, outdir){
  
  # select median age source
  tmp <- median_age_source[quantile == 'C50']
  
  # plot
  p <- ggplot(tmp) + 
    geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = ROUND)) + 
    geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, fill = ROUND), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    labs(x = 'Age', y = 'Median age at transmission source') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +
    facet_grid(LABEL_COMMUNITY~LABEL_RECIPIENT) 
  
  ggsave(p, file = paste0(outdir, '-output-MedianAgeSource_ByAgeRecipient.png'), w = 7, h = 6)
  
}

plot_median_age_source_group <- function(median_age_source_group, expected_contribution_age_group_source2, 
                                         reported_contact, outdir){
  
  mas <- copy(median_age_source_group)
  
  # select rounds to be plotted
  mas <- mas[(COMM == 'inland' & ROUND %in% paste0('R0', c(15, 18)))|(COMM == 'fishing' & ROUND %in% paste0('R0', c(15, 18)))]
  
  # dcast
  mag <- dcast.data.table(mas, LABEL_ROUND + COMM + AGE_GROUP_INFECTION.RECIPIENT + LABEL_RECIPIENT ~ quantile, value.var = 'M')
  
  # find expected transmission flows towards each age group of the recipient for the width of the boxplot
  eca <- expected_contribution_age_group_source2[, .(LABEL_ROUND, COMM, AGE_GROUP_INFECTION.RECIPIENT, LABEL_RECIPIENT, M)]
  setnames(eca, 'M', "M_CONTRIBUTION")
  
  # merge
  mac <- merge(mag, eca, by =  c('LABEL_ROUND', 'COMM', 'AGE_GROUP_INFECTION.RECIPIENT', 'LABEL_RECIPIENT'))

  # find boxplot quantiles for reported contact
  rec <- copy(reported_contact)
  rec <- rec[, list(C10= quantile(part.age, 0.1), 
                    C25 = quantile(part.age, 0.25), 
                    C50 = quantile(part.age, 0.5), 
                    C75 = quantile(part.age, 0.75), 
                    C90 = quantile(part.age, 0.9), 
                    CONTRIBUTION = .N), by = c('COMM', 'ROUND', 'LABEL_RECIPIENT', 'AGE_GROUP')]
  rec[, M_CONTRIBUTION := CONTRIBUTION / sum(CONTRIBUTION), by = c('COMM', 'ROUND')]
  set(rec, NULL, 'CONTRIBUTION', NULL)
  setnames(rec, 'AGE_GROUP', 'AGE_GROUP_INFECTION.RECIPIENT')
  
  # make labels
  mac[, LABEL := paste0('Estimated transmission\nsources, ', LABEL_ROUND)]
  rec[, LABEL := paste0('Estimated sexual contacts\nin last year, ', mac[grepl('Round 15', LABEL_ROUND), as.character(unique(LABEL_ROUND))])]
  dt <- rbind(mac, rec, fill=TRUE)
  
  # plot
  communities <- mac[, unique(COMM)]
  cols <- palette_round_inland[c(4, 7)] # color are the same for inland and fishing ( because same round selected)
  
  for(i in seq_along(communities)){
    
    tmp <- dt[COMM == communities[i]]
    tmp <- tmp[order(LABEL, AGE_GROUP_INFECTION.RECIPIENT, LABEL_RECIPIENT)]
    
    # width of the boxplot need to be proportional to a max (just for the figure to look pretty)
    widths <- tmp[order(LABEL, AGE_GROUP_INFECTION.RECIPIENT, LABEL_RECIPIENT), M_CONTRIBUTION]
    widths <- widths * 1.2 / max(widths)
    
    p <- ggplot(tmp) + 
      geom_boxplot(stat = "identity",
                   aes(x = AGE_GROUP_INFECTION.RECIPIENT, fill = LABEL,
                       lower  = C25,upper = C75, middle = C50, ymin = C10, ymax = C90, 
                       group= interaction(LABEL, AGE_GROUP_INFECTION.RECIPIENT)),
                   width = widths, varwidth = T, col = 'black', outlier.shape = NA, 
                   size = 0.2, position = position_dodge2(width = widths)) +
      facet_grid(.~LABEL_RECIPIENT) + 
      theme_bw() +
      labs(x = 'Age recipient', y = 'Age source') +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.text.x = element_text(angle= 30, hjust =1),
            panel.grid.minor.x = element_blank(), 
            legend.spacing.x = unit(0.1, 'cm'),
            legend.title = element_blank()) +
      # scale_color_manual(values = cols)  +
      scale_fill_manual(values = c('grey50', cols))  +
      scale_y_continuous(expand = c(0,0), limits = range_age_non_extended,
                         breaks = c(seq(min(range_age_non_extended), max(range_age_non_extended), 5), 
                           max(range_age_non_extended))) + 
      guides(color = guide_legend(order = 1))
    ggsave(p, file = paste0(outdir, '-output-MedianAgeSource_ByAgeGroupRecipient_', communities[i], '.pdf'), w = 5.5, h = 3.7)
    
  }
}

plot_counterfactual <- function(counterfactuals_p_f, counterfactuals_p_f05, 
                                counterfactuals_p_959595, counterfactuals_p_909090, 
                                incidence_factual, lab, outdir){
  
  #
  # labels
  #
  
  # label in scenario where treated male take up art as much as female
  label.f = 'ART coverage in men as in women\nSuppression rate in men as in women (95%)'
  if(!grepl('Diagnosed', lab)){
    label.f <- paste0('Diagnosed rate in men as in women (95%)\n', label.f)
  }
  
  # label in scenario where treated male take up art half as much as female
  label.f05 = paste0('Half-way to\n', label.f)
  
  # label in scenario where treated male are treated at 95%
  label.959595 = '95% receiving ART\n95% suppression rate'
  if(!grepl('Diagnosed', lab)){
    label.959595 <- paste0('95% diagnosed\n', label.959595)
  }
  
  # label in scenario where treated male are treated at half 90/90/90
  label.909090 = '90% receiving ART\n90% suppression rate'
  if(!grepl('Diagnosed', lab)){
    label.909090 <- paste0('90% diagnosed\n', label.909090)
  }
  
  
  #
  # unlist objects in scenario where male are diagnosed/treated/suppressed as much as female
  #

  # number of male treated regardless of age
  budget.counterfactual <- counterfactuals_p_f$budget 
  budget.counterfactual[, label := label.f]
  
  # number of male treated by age
  eligible_count_round.counterfactual <- counterfactuals_p_f$eligible_count_round.counterfactual
  eligible_count_round.counterfactual[, label := label.f]
  
  # relative incident cases counterfactual compared to factual by age
  relative_incidence_counterfactual <- counterfactuals_p_f$relative_incidence_counterfactual 
  relative_incidence_counterfactual[, label := label.f]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all <- counterfactuals_p_f$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all[, label := label.f]
  
  # incident cases counterfactual 
  incidence_counterfactual <- counterfactuals_p_f$incidence_counterfactual 
  incidence_counterfactual[, label := label.f]
  

  #
  # unlist objects in scenario where male are diagnosed/treated/suppressed half as much as female
  #
  
  # number of male treated regardless of age
  budget.counterfactual.f05 <- counterfactuals_p_f05$budget 
  budget.counterfactual.f05[, label := label.f05]
  
  # number of male treated by age
  eligible_count_round.counterfactual.f05 <- counterfactuals_p_f05$eligible_count_round.counterfactual
  eligible_count_round.counterfactual.f05[, label := label.f05]
  
  # relative incident cases counterfactual compared to factual by age
  relative_incidence_counterfactual.f05 <- counterfactuals_p_f05$relative_incidence_counterfactual 
  relative_incidence_counterfactual.f05[, label := label.f05]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all.f05 <- counterfactuals_p_f05$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all.f05[, label := label.f05]
  
  # incident cases counterfactual 
  incidence_counterfactual.f05 <- counterfactuals_p_f05$incidence_counterfactual 
  incidence_counterfactual.f05[, label := label.f05]
  
  
  #
  # unlist objects in scenario where treated male are treated at 95-95-95
  #
  
  # number of male treated regardless of age
  budget.counterfactual.959595 <- counterfactuals_p_959595$budget 
  budget.counterfactual.959595[, label := label.959595]
  
  # number of male treated regardless by age
  eligible_count_round.counterfactual.959595 <- counterfactuals_p_959595$eligible_count_round.counterfactual
  eligible_count_round.counterfactual.959595[, label := label.959595]
  
  # relative incident cases counterfactual compared to factual by age
  relative_incidence_counterfactual.959595 <- counterfactuals_p_959595$relative_incidence_counterfactual 
  relative_incidence_counterfactual.959595[, label := label.959595]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all.959595 <- counterfactuals_p_959595$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all.959595[, label := label.959595]
  
  # incident cases counterfactual 
  incidence_counterfactual.959595 <- counterfactuals_p_959595$incidence_counterfactual 
  incidence_counterfactual.959595[, label := label.959595]


  #
  # unlist objects in scenario where treated male are treated 90 90 90
  #
  
  # number of male treated regardless of age
  budget.counterfactual.909090 <- counterfactuals_p_909090$budget 
  budget.counterfactual.909090[, label := label.909090]

  # number of male treated regardless by age
  eligible_count_round.counterfactual.909090 <- counterfactuals_p_909090$eligible_count_round.counterfactual
  eligible_count_round.counterfactual.909090[, label := label.909090]
  
  # relative incident cases counterfactual compared to factual by age
  relative_incidence_counterfactual.909090 <- counterfactuals_p_909090$relative_incidence_counterfactual 
  relative_incidence_counterfactual.909090[, label := label.909090]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all.909090 <- counterfactuals_p_909090$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all.909090[, label := label.909090]
  
  # incident cases counterfactual 
  incidence_counterfactual.909090 <- counterfactuals_p_909090$incidence_counterfactual 
  incidence_counterfactual.909090[, label := label.909090]
  
  
  #
  # combine both scenarios
  #
  
  budget.counterfactual <- do.call('rbind', list(budget.counterfactual, budget.counterfactual.f05, 
                                                 budget.counterfactual.959595, budget.counterfactual.909090))
  budget.counterfactual[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  eligible_count_round.counterfactual <-  do.call('rbind', list(eligible_count_round.counterfactual, eligible_count_round.counterfactual.f05,
                                                                eligible_count_round.counterfactual.959595, eligible_count_round.counterfactual.909090))
  eligible_count_round.counterfactual[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  relative_incidence_counterfactual <- do.call('rbind', list(relative_incidence_counterfactual, relative_incidence_counterfactual.f05,
                                                             relative_incidence_counterfactual.959595, relative_incidence_counterfactual.909090))
  relative_incidence_counterfactual[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  relative_incidence_counterfactual_all <- do.call('rbind', list(relative_incidence_counterfactual_all, relative_incidence_counterfactual_all.f05,
                                                                 relative_incidence_counterfactual_all.959595, relative_incidence_counterfactual_all.909090))
  relative_incidence_counterfactual_all[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  incidence_counterfactual <- do.call('rbind', list(incidence_counterfactual, incidence_counterfactual.f05,
                                                    incidence_counterfactual.959595, incidence_counterfactual.909090))
  incidence_counterfactual[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  

  #
  # restrict to one round and to male to female direction
  #
  
  Round <- 'R018'
  budget.counterfactual <- budget.counterfactual[ROUND == Round & SEX == 'M']
  relative_incidence_counterfactual <- relative_incidence_counterfactual[ROUND == Round & IS_MF == T]
  incidence_counterfactual <- incidence_counterfactual[ROUND == Round& IS_MF == T]
  icf <- incidence_factual[ROUND == Round ]
  ecf <- eligible_count_round.counterfactual[ROUND == Round & SEX == 'M']
  relative_incidence_counterfactual_all <- relative_incidence_counterfactual_all[ROUND == Round & IS_MF == T]
  
  
  #
  # Clean and merge to target labels
  #
  
  # merge udget by by age to target labels
  ecf[, INFECTED_SUPPRESSED := INFECTED - INFECTED_NON_SUPPRESSED]
  ecf[, INFECTED_ALREADY_SUPPRESSED := INFECTED_SUPPRESSED - TREATED]
  ecf <- ecf[, .(ROUND, SEX, AGEYRS, COMM, label, INFECTED_NON_SUPPRESSED, INFECTED_ALREADY_SUPPRESSED, TREATED)]
  ecf <- melt.data.table(ecf, id.vars = c('ROUND', 'SEX', 'AGEYRS', 'COMM', 'label'))
  
  # make labels
  label.suppressed = 'Virally suppressed in round 18'; label.unsuppressed = 'Virally unsuppressed'; 
  label.new.suppressed = 'Additionally suppressed\nin intervention'
  ecf[, VARIABLE_LABEL := label.suppressed]
  ecf[variable == 'TREATED', VARIABLE_LABEL := label.new.suppressed]
  ecf[variable == 'INFECTED_NON_SUPPRESSED', VARIABLE_LABEL := label.unsuppressed]
  ecf[, VARIABLE_LABEL := factor(VARIABLE_LABEL, levels = c(label.unsuppressed, label.new.suppressed, label.suppressed))]
  ecf[, REDUCTION_UNSUPPRESSED := value[variable == 'INFECTED_NON_SUPPRESSED'] / (value[variable == 'INFECTED_NON_SUPPRESSED'] + value[variable == 'TREATED']), by= c('ROUND', 'SEX', 'AGEYRS', 'COMM', 'label')]
  ecf[, VARIABLE_LABEL2 := VARIABLE_LABEL]
  ecf[VARIABLE_LABEL == label.new.suppressed, VARIABLE_LABEL2 := label]
  ecf[, VARIABLE_LABEL2 := factor(VARIABLE_LABEL2, levels = c(label.unsuppressed, 
                                                              ecf[!VARIABLE_LABEL2 %in% c(label.suppressed, label.unsuppressed), as.character(unique(VARIABLE_LABEL2))], 
                                                              label.suppressed))]
  
  # format budget regardless of age, merge to target labels and add total number of unsuppressed in factual
  bc <- melt.data.table(budget.counterfactual, id.vars = c('ROUND', 'SEX', 'COMM', 'label'))
  # bc <- ecf[, .(value = sum(value)), by = c('ROUND', 'SEX', 'COMM', 'label', 'counterfactual_index', 'VARIABLE_LABEL', 'variable')]
  # bc <- merge(bc, tmp, by = c('ROUND', 'SEX', 'COMM', 'label', 'counterfactual_index'), allow.cartesian = T)
  # bc <- merge(bc, df_target, by = 'counterfactual_index')
  bc[, lab := lab]
  
  # add sex label
  icf[, SEX_LABEL := 'In women']
  icf[IS_MF == F, SEX_LABEL := 'In men']
  icf[, SEX_LABEL := factor(SEX_LABEL, levels = c('In women', 'In men'))]
  
  # merge incidence cases counterfactual by ageto target labels
  # ic <- merge(incidence_counterfactual, df_target, by = 'counterfactual_index')
  ic <- copy(incidence_counterfactual)
  
  # merge incidence cases counterfactual relative to factual to target labels
  # ric <- merge(relative_incidence_counterfactual, df_target, by = 'counterfactual_index')
  ric <- copy(relative_incidence_counterfactual)
  # ric.all <- merge(relative_incidence_counterfactual_all, df_target, by = 'counterfactual_index')
  ric.all <- copy(relative_incidence_counterfactual_all)
  
  communities <- ric[, unique(COMM)]
  cols <- c('#F1A661', '#C55300', '#425F57')
  
  for(i in seq_along(communities)){
    
    bc.c <- bc[COMM == communities[i]]
    ic.c <- ic[COMM == communities[i]]
    ric.c <- ric[COMM == communities[i]]
    icf.c <- icf[COMM == communities[i]]
    ecf.c <- ecf[COMM == communities[i]]
    ric.all.c <- ric.all[COMM == communities[i]]
    
    # remove 90-90-90
    bc.c <- bc.c[label != label.909090]
    ic.c <- ic.c[label != label.909090]
    ric.c <- ric.c[label != label.909090]
    ecf.c <- ecf.c[label != label.909090]
    ric.all.c <- ric.all.c[label != label.909090]
    
    # budget by age group
    p <- ggplot(ecf.c, aes(x = AGEYRS)) +
      geom_bar(aes(y = value, fill = VARIABLE_LABEL2), stat = 'identity') +  
      facet_wrap(~label, ncol = 1) + 
      theme_bw() + 
      labs(y = 'Number of infected men', x = 'Age', col = 'Male treated') + 
      theme(legend.title = element_blank(),
            legend.position = 'non',
            strip.background = element_rect(colour="white", fill="white"), 
            strip.text = element_blank()) + 
      scale_fill_manual(values = c('grey50', cols[c(2,1,3)], 'grey80')) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
      guides(fill = guide_legend(byrow = T, nrow = 6))
    file = paste0(outdir, '-output-counterfactual_strategy_budget_age_', gsub(' ' , '', lab), '_', communities[i], '.pdf')
    ggsave(p, file = file, w = 3.9, h = 3.7)
    
    # for the legend
    p <- ggplot(ecf.c, aes(x = AGEYRS)) +
      geom_bar(aes(y = value, fill = VARIABLE_LABEL), stat = 'identity') +  
      facet_wrap(~label, ncol = 1) + 
      theme_bw() + 
      labs(y = 'Number of infected men', x = 'Age', col = 'Male treated') + 
      theme(legend.title = element_blank(),
            legend.position = 'bottom',
            strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1))) + 
      scale_fill_manual(values = c('grey50', cols[c(3,1,4,2)], 'grey80')) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
      guides(fill = guide_legend(byrow = T, nrow = 4))
    file = paste0(outdir, '-output-counterfactual_strategy_budget_legend_age_', gsub(' ' , '', lab), '_', communities[i], '.pdf')
    ggsave(p, file = file, w = 7, h = 10)
    
    # reduction in unsuppressed
    p <- ggplot(ecf.c, aes(x = AGEYRS)) +
      geom_line(aes(y = 1 - REDUCTION_UNSUPPRESSED, col = label), stat = 'identity') +  
      theme_bw() + 
      labs(y = '% reduction in number\nof unsuppressed', x = 'Age', col = 'Male treated') + 
      theme(legend.title = element_blank(),
            legend.position = 'bottom',
            strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1))) +
      scale_y_continuous(labels = scales::percent, limits = c(0, NA)) + 
      scale_color_manual(values = cols) + 
      guides(color = guide_legend(byrow = T, nrow = 4))
    file = paste0(outdir, '-output-counterfactual_strategy_budget_reduction_age_', gsub(' ' , '', lab), '_', communities[i], '.png')
    ggsave(p, file = file, w = 5, h = 6)
    
    
    # budget regardless of age
    p1 <- ggplot(bc.c) + 
      geom_bar(aes(y = value, x = label, fill = label), stat = 'identity') +  
      labs(x = '', y = paste0('Additional number\nof men suppressed'), fill = '') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x = element_blank(), 
            # axis.text.x = element_text(angle = 90, hjust = 1),
            # legend.position = c(0.85, 0.83), 
            legend.key.size = unit(0.4, 'cm'),
            legend.position = c(0.9,0.9), 
            axis.title.x = element_blank(),
            legend.title = element_blank(), 
            panel.grid.minor.y = element_blank(), 
      ) +
      scale_fill_manual(values = cols) +
      guides(fill = 'none') + 
      scale_y_continuous(expand = expansion(mult = c(0, .05))) 

    # incidence cases by age
    p2 <- ggplot(ic.c) + 
      geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, group = interaction(label), fill = label), alpha = 0.2) + 
      # geom_ribbon(data = icf.c, aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, linetype = SEX_LABEL), alpha = 0.2) +
      geom_line(data = icf.c, aes(x = AGE_INFECTION.RECIPIENT, y = M, linetype = SEX_LABEL), col = 'black') +
      geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = label)) + 
      labs(x = 'Age', y = '\nIncidence cases', linetype = 'No intervention', fill = 'Intervention', color = 'Intervention') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            axis.title.x = element_blank(), 
            legend.direction = 'vertical',
            legend.spacing.y = unit(0.25, 'cm'),
            legend.position = 'bottom') +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous() + 
      scale_x_continuous(expand = c(0,0)) + 
      guides(fill = guide_legend(order = 2, byrow = T), col = guide_legend(order = 2, byrow = T), linetype = guide_legend(order = 1, byrow = T))
    
    # reduction infection cases by age
    p3 <- ggplot(ric.c, aes(x = AGEYRS)) + 
      geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = label)) + 
      geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, group = interaction(label), fill = label), alpha = 0.25) + 
      labs(x = 'Age', y = '% reduction in\nincidence in women', 
           fill = '', col = '') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.title = element_blank(), 
            # legend.direction = 'vertical',
            # axis.title.x = element_blank(), 
            # axis.text.x = element_blank(), 
            legend.position = 'none') +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(labels = scales::percent, limits = c(0,0.8), expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0)) 

    # reduction infection cases regardless of age
    p4 <- ggplot(ric.all.c) + 
      geom_bar(aes(x = label, y = M, fill = label), stat = 'identity') + 
      geom_errorbar(aes(x = label, ymin= CL, ymax = CU, group = interaction(label)), alpha = 0.5, width = 0.15) + 
      labs(x = '', y = '% reduction in\nincidence in women', 
           fill = '', col = '') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.title = element_blank(), 
            axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.position = 'none') +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(labels = scales::percent, limits = c(0,0.75), expand = c(0,0)) 

    # arrange
    
    p3 <- ggarrange(p3, legend.grob = get_legend(p2), legend = 'bottom')
    p2 <- p2 + theme(legend.position= 'none')
    
    p <- grid.arrange(p1, p4, p2, p3, layout_matrix = rbind(c(NA, 1, 1), c(2, 2, 2), c(NA,NA, 3), c(4, 4, 4)),
                      widths = c(0.007, 0.022, 0.95), heights = c(0.15, 0.15, 0.2, 0.39))

    # save
    file = paste0(outdir, '-output-counterfactual_strategy_incidence_panel_', gsub(' ' , '', lab), '_', communities[i], '.pdf')
    ggsave(p, file = file, w = 5.5, h = 9.5)
    
    file = paste0(outdir, '-output-counterfactual_strategy_incidence_panel_plot1_', gsub(' ' , '', lab), '_', communities[i], '.pdf')
    ggsave(p1, file = file, w = 4, h = 1.85)

    file = paste0(outdir, '-output-counterfactual_strategy_incidence_panel_plot2_', gsub(' ' , '', lab), '_', communities[i], '.pdf')
    ggsave(p4, file = file, w = 4, h = 1.85)
    
    p2 <- p2 + theme(axis.title.x = element_text())
    file = paste0(outdir, '-output-counterfactual_strategy_incidence_panel_plot3_', gsub(' ' , '', lab), '_', communities[i], '.pdf')
    ggsave(p2, file = file, w = 4, h = 2.9)
    
  }
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

plot_counterfactual_one <- function(counterfactuals_p_a, incidence_factual, lab, outdir){
  
  #
  # target label
  #
  
  df_target = data.table(variable = c('TREATED.SPREADERS', 'TREATED.NONCOMPLIER', 'TREATED.RANDOM'), 
                         counterfactual_index = 1:3, 
                         type = c('main spreaders', 'non compliers', 'random'),
                         LABEL_TARGET = c('Top 1/3 men sources', 
                                          'Men age groups with greatest\ndifference in ART coverage\ncompared to women', 
                                          'Untargeted intervention'))
  df_target[, LABEL_TARGET := factor(LABEL_TARGET, levels = df_target[, LABEL_TARGET])]
  
  #
  # unlist objects from counterfactual
  #
  
  # number of male treated by age
  eligible_count_round.counterfactual <- counterfactuals_p_a$eligible_count_round.counterfactual
  
  # number of male treated regardless of age
  budget.counterfactual <- counterfactuals_p_a$budget 
  
  # relative incidence cases counterfactual compared to factual by age
  relative_incidence_counterfactual <- counterfactuals_p_a$relative_incidence_counterfactual 
  
  # relative incidence cases counterfactual compared to factual regardless age
  relative_incidence_counterfactual_all <- counterfactuals_p_a$relative_incidence_counterfactual_all 
  
  # incidence cases counterfactual by age
  incidence_counterfactual <- counterfactuals_p_a$incidence_counterfactual 
  
  
  #
  # restrict to one round and to male to female direction
  #
  
  Round <- 'R018'
  budget.counterfactual <- budget.counterfactual[ROUND == Round & SEX == 'M']
  relative_incidence_counterfactual <- relative_incidence_counterfactual[ROUND == Round & IS_MF == T]
  incidence_counterfactual <- incidence_counterfactual[ROUND == Round & IS_MF == T]
  icf <- incidence_factual[ROUND == Round & IS_MF == T]
  relative_incidence_counterfactual_all <- relative_incidence_counterfactual_all[ROUND == Round & IS_MF == T]
  ecr <- eligible_count_round.counterfactual[ROUND == Round & SEX == 'M']
  
  
  #
  # Clean objects
  #
  
  # format budget and merge to target label
  bc <- melt.data.table(budget.counterfactual, id.vars = c('ROUND', 'SEX', 'COMM'))
  bc <- merge(bc, df_target, by = 'variable')
  bc[, lab := lab]
  
  # merge incidence cases counterfactual to target label 
  ic <- merge(incidence_counterfactual, df_target, by = 'counterfactual_index')
  ecr <- merge(ecr, df_target, by = 'counterfactual_index')
  
  # merge relative incidence cases counterfactual to target label 
  ric <- merge(relative_incidence_counterfactual, df_target, by = 'counterfactual_index')
  ric.all <- merge(relative_incidence_counterfactual_all, df_target, by = 'counterfactual_index')
  
  
  #
  # Plot
  #
  
  communities <- ric[, unique(COMM)]
  cols <- c('#CC3636', '#F57328', '#367E18')
  
  for(i in seq_along(communities)){
    
    bc.c <- bc[COMM == communities[i]]
    ic.c <- ic[COMM == communities[i]]
    ric.c <- ric[COMM == communities[i]]
    icf.c <- icf[COMM == communities[i]]
    ric.all.c <- ric.all[COMM == communities[i]]
    ecr.c <- ecr[COMM == communities[i]]
    
    # budget total
    p <- ggplot(ecr.c, aes(x = AGEYRS)) + 
      geom_bar(aes(fill = LABEL_TARGET, y = TREATED), stat = 'identity') + 
      labs(x = 'Age', y = paste0('Additional number\nof men treated'), fill = '') + 
      theme_bw() +
      facet_wrap(~LABEL_TARGET, ncol = 1) +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.title.x = element_blank(),
            legend.title = element_blank(), 
      ) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) 
    file = paste0(outdir, '-output-counterfactual_target_incidence_budget_', gsub(' ' , '', lab), '_all_', communities[i], '.pdf')
    ggsave(p, file = file, w = 5.5, h = 7.6)
    
    # budget total
    p1 <- ggplot(bc.c, aes(x = 1)) + 
      geom_bar(aes(fill = LABEL_TARGET, y = value), stat = 'identity', position = position_dodge(0.7), width = 0.6) + 
      labs(x = '', y = paste0('Additional number\nof men treated'), fill = '') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x = element_blank(), 
            # legend.position = c(0.85, 0.83), 
            legend.key.size = unit(0.4, 'cm'),
            legend.position = 'none', 
            axis.title.x = element_blank(),
            legend.title = element_blank(), 
            panel.grid.minor.y = element_blank(), 
      ) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(expand = expansion(mult = c(0, .05))) 
    
    # incidence cases by age
    fit_data_label <- 'No intervention\nRound 18'
    p2 <- ggplot(ic.c) + 
      geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, group = interaction(LABEL_TARGET), fill = LABEL_TARGET), alpha = 0.2) + 
      geom_ribbon(data = icf.c, aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, linetype = fit_data_label), alpha = 0.2) + 
      geom_line(data = icf.c, aes(x = AGE_INFECTION.RECIPIENT, y = M, linetype = fit_data_label), col = 'black') +
      geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = LABEL_TARGET)) + 
      labs(x = 'Age', y = 'Incidence cases') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            axis.title.x = element_blank(), 
            legend.position = 'bottom',
            legend.box= 'vertical', 
            legend.title = element_blank()) +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous() + 
      scale_x_continuous(expand = c(0,0)) 
    
    # reduction infection by age
    p3 <- ggplot(ric.c, aes(x = AGEYRS)) + 
      geom_line(aes(x = AGE_INFECTION.RECIPIENT, y = M, col = LABEL_TARGET)) + 
      geom_ribbon(aes(x = AGE_INFECTION.RECIPIENT, ymin= CL, ymax = CU, group = interaction(LABEL_TARGET), fill = LABEL_TARGET), alpha = 0.25) + 
      labs(x = 'Age', y = '% reduction in\nincidence in women', 
           fill = '', col = '') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.title = element_blank(), 
            # legend.direction = 'vertical',
            # axis.title.x = element_blank(), 
            # axis.text.x = element_blank(), 
            legend.position = 'none') +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(labels = scales::percent, limits = c(0,0.25), expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0)) 
    
    # ireduction infection regardless of age
    p4 <- ggplot(ric.all.c, aes(x = AGEYRS)) + 
      geom_bar(aes(x = counterfactual_index, y = M, fill = LABEL_TARGET), stat = 'identity') + 
      geom_errorbar(aes(x = counterfactual_index, ymin= CL, ymax = CU, group = interaction(LABEL_TARGET)), alpha = 0.5, width = 0.15) + 
      labs(x = '', y = '% reduction in\nincidence in women', 
           fill = '', col = '') + 
      theme_bw() +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            legend.title = element_blank(), 
            axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.position = 'none') +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      scale_y_continuous(labels = scales::percent, limits = c(0,0.65), expand = c(0,0)) 
    
    # arrange
    p3 <- ggarrange(p3, legend.grob = get_legend(p2), legend = 'bottom')
    p2 <- p2 + theme(legend.position= 'none')
    
    p <- grid.arrange(p1, p4, p2, p3, layout_matrix = rbind(c(NA, 1, 1), c(2, 2, 2), c(NA,NA, 3), c(4, 4, 4)), 
                      widths = c(0.02, 0.012, 0.95), heights = c(0.16, 0.16, 0.2, 0.34))
    
    # save
    file = paste0(outdir, '-output-counterfactual_target_incidence_panel_', gsub(' ' , '', lab), '_all_', communities[i], '.pdf')
    ggsave(p, file = file, w = 5.5, h = 7.6)
    
  }
}

