plot_intensity_PP <- function(intensity_PP, count_data, direction, outdir){
  
  intensity_PP <- intensity_PP[label_direction == direction]
  count_data <- count_data[label_direction == direction]
  
  count_data_reduced <- count_data[count > 0]
  
  dates = unique(count_data$date_infection_reduced.RECIPIENT)
  mid_date = floor(seq(1, length(dates), ceiling( length(dates)/3)))
  
  idx = list(1:(mid_date[2]-1), 
             mid_date[2]:(mid_date[3]-1),
             mid_date[3]:length(dates))
  
  p = list()
  for(i in 1:length(idx)){
    Date = dates[idx[[i]]]
    tmp <- intensity_PP[date_infection_reduced.RECIPIENT %in% Date]
    tmp1 <- subset(count_data_reduced, date_infection_reduced.RECIPIENT %in% Date)
    p[[i]] <- ggplot(tmp, aes(y = age_infection_reduced.SOURCE, x = age_infection_reduced.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      geom_point(data = tmp1, aes(size = count), col = 'grey50') +
      theme_bw() + 
      coord_fixed() +
      labs(x = 'Age at infection recipient', fill = 'transmission rate') +
      geom_contour(aes(z = M), col = 'red', alpha = 0.8) + 
      facet_grid(.~date_infection_reduced_name.RECIPIENT) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.title.y = element_blank()) +
      scale_fill_viridis_c(limits = range(intensity_PP$M)) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_size_continuous(range = c(1, 3))
    
    if(i != length(idx)){
      p[[i]] <-  p[[i]] + theme(axis.title.x = element_blank())
    }
  }

  
  p <- ggarrange(plotlist = p, nrow = 3, common.legend = T, legend= 'bottom')
  p <- grid.arrange(p, top = text_grob(direction, size = 16), 
                    left = text_grob('Age at infection source', rot = 90, size = 13))
  ggsave(p, file = paste0(outdir, '-intensity_transmission_', gsub(' -> ', '_', direction), '.png'), w = 8, h = 10)
}


plot_mean_age_source_old <- function(age_source, direction, outdir){
  
  age_source <- age_source[label_direction == direction]
  
  dates = unique(age_source$date_infection_reduced.RECIPIENT)
  mid_date = floor(seq(1, length(dates), ceiling( length(dates)/3)))
  
  idx = list(1:(mid_date[2]-1), 
             mid_date[2]:(mid_date[3]-1),
             mid_date[3]:length(dates))
  
  p = list()
  for(i in 1:length(idx)){
    Date = dates[idx[[i]]]
    tmp <- age_source[date_infection_reduced.RECIPIENT %in% Date]
    p[[i]] <- ggplot(tmp, aes(x = age_infection_reduced.RECIPIENT)) + 
      geom_point(aes(y = M)) + 
      geom_errorbar(aes(ymin= CL, ymax = CU)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
      theme_bw() + 
      coord_fixed() +
      labs(x = 'Age at infection recipient') +
      facet_grid(.~date_infection_reduced_name.RECIPIENT) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.title.y = element_blank()) +
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) 
    
    if(i != length(idx)){
      p[[i]] <-  p[[i]] + theme(axis.title.x = element_blank())
    }
  }
  
  
  p <- ggarrange(plotlist = p, nrow = 3, common.legend = T, legend= 'bottom')
  p <- grid.arrange(p, top = text_grob(direction, size = 16), 
                    left = text_grob('Mean age at infection source', rot = 90, size = 13))
  ggsave(p, file = paste0(outdir, '-MeanAgeSource_ByAgeRecipient_', gsub(' -> ', '_', direction), '.png'), w = 8, h = 10)
}

plot_mean_age_source <- function(age_source, direction, outdir){
  
  age_source <- age_source[label_direction == direction]

  p <- ggplot(age_source, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_line(aes(y = M, col = date_infection_reduced_name.RECIPIENT)) + 
    geom_ribbon(aes(ymin= CL, ymax = CU, fill = date_infection_reduced_name.RECIPIENT), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    labs(x = 'Age at infection recipient', y = 'Mean age at infection source',
         col = 'Date infection recipient', fill = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0), limits = range(age_source$age_infection_reduced.RECIPIENT)) + 
    scale_color_viridis_d(option = 'A', end = 0.9)+ 
    scale_fill_viridis_d(option = 'A', end = 0.9) +
    guides(fill=guide_legend(nrow=4,byrow=TRUE), col=guide_legend(nrow=4,byrow=TRUE))+
    ggtitle(direction)

  ggsave(p, file = paste0(outdir, '-MeanAgeSource_ByAgeRecipient_', gsub(' -> ', '_', direction), '.png'), w = 6, h = 6)
}

plot_mean_age_source_overall <- function(age_source_overall, direction, outdir){
  
  age_source_overall <- age_source_overall[label_direction == direction]
  
  p <- ggplot(age_source_overall, aes(x = date_infection_reduced.RECIPIENT)) + 
    geom_line(aes(y = M)) + 
    geom_ribbon(aes(ymin= CL, ymax = CU), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    labs(y = 'Mean age at infection source', x = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_date(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    ggtitle(direction)
  
  ggsave(p, file = paste0(outdir, '-MeanAgeSource_', gsub(' -> ', '_', direction), '.png'), w = 6, h = 4)
}


