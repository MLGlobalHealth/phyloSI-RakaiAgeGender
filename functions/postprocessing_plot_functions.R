plot_intensity_PP <- function(intensity_PP, count_data, outdir){
  
  count_data_reduced <- count_data[count > 0]
  
  dates = unique(count_data$date_infection_reduced.RECIPIENT)
  mid_date = floor(length(dates) / 2)
  
  tmp <- subset(intensity_PP, date_infection_reduced.RECIPIENT %in% dates[1:mid_date])
  tmp1 <- subset(count_data_reduced, date_infection_reduced.RECIPIENT %in% dates[1:mid_date])
  p1 <- ggplot(tmp, aes(x = age_infection_reduced.SOURCE, y = age_infection_reduced.RECIPIENT)) + 
    geom_raster(aes(fill = M)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
    geom_point(data = tmp1, aes(size = count), col = 'grey50') +
    theme_bw() + 
    coord_fixed() +
    labs(x = 'Age at infection source', y = 'Age at infection recipient', fill = 'transmission rate') +
    geom_contour(aes(z = M), col = 'red') + 
    facet_grid(label_direction~date_infection_reduced_name.RECIPIENT) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_fill_viridis_c(limits = range(intensity_PP$M)) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) 
  
  tmp <- subset(intensity_PP, date_infection_reduced.RECIPIENT %in% dates[(mid_date+1):length(dates)])
  tmp1 <- subset(count_data_reduced, date_infection_reduced.RECIPIENT %in% dates[(mid_date+1):length(dates)])
  p2 <- ggplot(tmp, aes(x = age_infection_reduced.SOURCE, y = age_infection_reduced.RECIPIENT)) + 
    geom_raster(aes(fill = M)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
    geom_point(data = tmp1, aes(size = count), col = 'grey50') +
    theme_bw() + 
    coord_fixed() +
    labs(x = 'Age at infection source', y = 'Age at infection recipient', fill = 'transmission rate') +
    geom_contour(aes(z = M), col = 'red') + 
    facet_grid(label_direction~date_infection_reduced_name.RECIPIENT) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_fill_viridis_c(limits = range(intensity_PP$M)) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) 
  
  p <- ggarrange(p1, p2, nrow = 2, common.legend = T, legend= 'bottom')
  ggsave(p, file = paste0(outdir, '-intensity_transmission.png'), w = 10, h = 8)
}
