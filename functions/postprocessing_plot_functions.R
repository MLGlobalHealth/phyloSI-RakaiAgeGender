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
    
    ggsave(p, file = paste0(outdir, '-intensity_transmission_',  communities[i], '.png'), w = 7, h = 7)
  }

} 

plot_force_infection <- function(force_infection, outdir){
  
  communities <- force_infection[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- force_infection[ COMM == communities[i]]

    p <- ggplot(tmp, aes(y = AGE_TRANSMISSION.SOURCE, x = AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = 'Estimated median\nforce infection', 
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
    
    ggsave(p, file = paste0(outdir, '-force_infection_',  communities[i], '.png'), w = 7, h = 7)
  }
  
} 

plot_force_infection_age_source <- function(force_infection_age_source, outdir){
  
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
  ggsave(p, file = paste0(outdir, '-force_infection_exerted_age.png'), w = 6, h = 7)
  
}

plot_force_infection_sex_source <- function(force_infection_sex_source, outdir){
  
  tmp <- copy(force_infection_sex_source)
  
  tmp[, Direction :=LABEL_DIRECTION]
  
  p <- ggplot(tmp, aes(x = Direction)) + 
    geom_bar(aes(y = M, fill = PERIOD), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = PERIOD), position = position_dodge()) + 
    labs(x = '', y = 'Force of infection exerted', fill = '') + 
    theme_bw() +
    facet_grid(LABEL_COMMUNITY~., scale = 'free')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    ggsci::scale_fill_npg() 
  ggsave(p, file = paste0(outdir, '-force_infection_exerted_sex.png'), w = 6, h = 7)
  
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

plot_PPC <- function(predict_y, count_data, outdir){
  data <- count_data[, list(count = sum(count)), by = c('LABEL_DIRECTION', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_INFECTION.RECIPIENT')]
  
  communities <- predict_y[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_y[ COMM == communities[i]]
    tmp1 <- data[ LABEL_COMMUNITY == tmp[, unique(LABEL_COMMUNITY)]]
    
    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M)) + 
      geom_ribbon(aes(ymin = CL, ymax = CL)) + 
      geom_point(data = tmp1, aes(y = count), col = 'grey50') +
      theme_bw() + 
      labs(x = 'Age at infection recipient', y = 'Total transmission events') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    
    ggsave(p, file = paste0(outdir, '-PPC_', communities[i], '.png'), w = 7, h = 7)
  }
  
}

plot_observed_to_augmented <- function(predict_y, predict_z, outdir){
  
  predict_z[, type := 'Augmented (Z)']
  predict_y[, type := 'Observed (Y)']
  
  predict_df <- rbind(predict_z, predict_y)
  
  communities <- predict_y[, unique(COMM)]
  for(i in seq_along(communities)){
    
    tmp <- predict_df[ COMM == communities[i]]

    p <- ggplot(tmp, aes( x = AGE_INFECTION.RECIPIENT)) + 
      geom_line(aes(y = M, col = type)) + 
      geom_ribbon(aes(ymin = CL, ymax = CL, fill = type), alpha = 0.5) + 
      theme_bw() + 
      labs(x = 'Age at infection', y = 'Total infected') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(LABEL_DIRECTION~PERIOD) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      ggtitle(tmp[,unique(LABEL_COMMUNITY)])
    
    ggsave(p, file = paste0(outdir, '-predict_total_infected_', communities[i], '.png'), w = 7, h = 7)
  }
  
}
