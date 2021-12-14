plot_intensity_reduced_PP <- function(intensity_PP, count_data, outdir){
  
  count_data_reduced <- count_data[count > 0]
  
  fct <- function(intensity_PP, count_data_reduced){
    ggplot(intensity_PP, aes(y = age_infection_reduced.SOURCE, x = age_infection_reduced.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      geom_point(data = count_data_reduced, aes(size = count), col = 'grey50') +
      theme_bw() + 
      coord_fixed() +
      labs(x = 'Age at infection recipient', fill = 'transmission rate', y= 'Age at infection source') +
      geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(label_direction~date_infection_reduced_name.RECIPIENT) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c(limits = range(intensity_PP$M)) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_size_continuous(range = c(1, 3))
  }
  
  p <- fct(intensity_PP, count_data_reduced)
  ggsave(p, file = paste0(outdir, '-intensity_transmission', '.png'), w = 10, h = 8)
  
  p = list(); p1 = list()
  Dates = unique(intensity_PP$date_infection_reduced_name.RECIPIENT)
  Directions = unique(intensity_PP$label_direction)
  for(j in 1:length(Directions)){
    for(i in 1:length(Dates)){
      Date = Dates[i]; Direction = Directions[j]
      tmp <- subset(intensity_PP, date_infection_reduced_name.RECIPIENT == Date & label_direction == Direction)
      tmp1 <- subset(count_data_reduced, date_infection_reduced_name.RECIPIENT == Date & label_direction == Direction)
      
      p[[i]] <- fct(tmp, tmp1) + theme(legend.position='none')
      
      if(i != length(Dates)){
        p[[i]] <- p[[i]] + 
          theme(strip.text.y = element_blank())
      }
      if(i != 1){
        p[[i]] <- p[[i]] + 
          theme(axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank())
      }
      
      if(j != 1){
        p[[i]] <- p[[i]] + 
          theme(strip.text.x = element_blank())
      }
      if(j != length(Directions)){
        p[[i]] <- p[[i]] + 
          theme(axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank())
      }
      
    }
    p1[[j]] <- ggarrange(plotlist = p, nrow = 1, widths = c(1, 0.87, 0.87, 1))
    
  }
  
  p <- ggarrange(plotlist = p1, nrow = 2, heights = c(1, 1))
  ggsave(p, file = paste0(outdir, '-intensity_transmission_idt', '.png'), w = 10, h = 5)
  
}

plot_intensity_PP <- function(intensity_PP, outdir){
  
  Directions <- unique(intensity_PP$label_direction)
  
  for(i in 1:length(Directions)){
    Direction <- Directions[i]
    tmp <- subset(intensity_PP, label_direction == Direction)
    
    p <- ggplot(intensity_PP, aes(y = age_infection_evaluated.SOURCE, x = age_infection_evaluated.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      theme_bw() + 
      coord_fixed() +
      labs(x = 'Age at infection recipient', fill = 'transmission rate', y= 'Age at infection source') +
      geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_wrap(.~date_infection_evaluated_name.RECIPIENT) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c(limits = range(intensity_PP$M)) + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_size_continuous(range = c(1, 3)) + 
      ggtitle(Direction)
    
    ggsave(p, file = paste0(outdir, '-intensity_transmission_all_', gsub(' -> ', '_', Direction), '.png'), w = 10, h = 8)
  }

  
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

plot_mean_age_source <- function(age_source, outdir){

  p <- ggplot(age_source, aes(x = age_infection_reduced.RECIPIENT)) + 
    geom_line(aes(y = M, col = date_infection_reduced_name.RECIPIENT)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    labs(x = 'Age at infection recipient', y = 'Mean age at infection source',
         col = 'Date infection recipient',) +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0), limits = range(age_source$age_infection_reduced.RECIPIENT)) + 
    scale_color_viridis_d(option = 'A', end = 0.9)+ 
    guides(col=guide_legend(nrow=2,byrow=TRUE))+
    facet_grid(.~label_direction)
  
  ggsave(p, file = paste0(outdir, '-MeanAgeSource_ByAgeRecipient.png'), w = 8, h = 5)
  
  Dates <- age_source$date_infection_reduced.RECIPIENT
  Dates <- Dates[seq(1, length(Dates), length.out = 3)]
  tmp <- subset(age_source, date_infection_reduced.RECIPIENT %in% Dates)
  p <- ggplot(tmp, aes(x = age_infection_reduced.RECIPIENT)) + 
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
    guides(fill=guide_legend(nrow=2,byrow=TRUE), col=guide_legend(nrow=2,byrow=TRUE))+
    facet_grid(.~label_direction)

  ggsave(p, file = paste0(outdir, '-MeanAgeSource_ByAgeRecipient_withCI.png'), w = 8, h = 5)
}

plot_mean_age_source_overall <- function(age_source_overall, outdir){
  
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
    facet_grid(label_direction~.)
  
  ggsave(p, file = paste0(outdir, '-MeanAgeSource.png'), w = 7, h = 5)
}


