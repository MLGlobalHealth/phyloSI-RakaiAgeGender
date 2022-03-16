# plot_intensity_reduced_PP <- function(intensity_PP, count_data, outdir){
#   
#   count_data_reduced <- count_data[count > 0]
#   
#   fct <- function(intensity_PP, count_data_reduced){
#     ggplot(intensity_PP, aes(y = age_infection_reduced.SOURCE, x = age_infection_reduced.RECIPIENT)) + 
#       geom_raster(aes(fill = M)) + 
#       geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
#       geom_point(data = count_data_reduced, aes(size = count), col = 'grey50') +
#       theme_bw() + 
#       coord_fixed() +
#       labs(x = 'Age at infection recipient', fill = 'transmission rate', y= 'Age at transmission source') +
#       geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
#       facet_grid(label_direction~date_infection_reduced_name.RECIPIENT) + 
#       theme(strip.background = element_rect(colour="white", fill="white"),
#             strip.text = element_text(size = rel(1)),
#             legend.position = 'bottom') +
#       scale_fill_viridis_c(limits = range(intensity_PP$M)) + 
#       scale_x_continuous(expand = c(0,0)) + 
#       scale_y_continuous(expand = c(0,0)) + 
#       scale_size_continuous(range = c(1, 3))
#   }
#   
#   p <- fct(intensity_PP, count_data_reduced)
#   ggsave(p, file = paste0(outdir, '-intensity_transmission', '.png'), w = 6, h = 6)
#   
#   p = list(); p1 = list()
#   Dates = unique(intensity_PP$date_infection_reduced_name.RECIPIENT)
#   Directions = unique(intensity_PP$label_direction)
#   for(j in 1:length(Directions)){
#     for(i in 1:length(Dates)){
#       Date = Dates[i]; Direction = Directions[j]
#       tmp <- subset(intensity_PP, date_infection_reduced_name.RECIPIENT == Date & label_direction == Direction)
#       tmp1 <- subset(count_data_reduced, date_infection_reduced_name.RECIPIENT == Date & label_direction == Direction)
#       
#       p[[i]] <- fct(tmp, tmp1) + theme(legend.position='none')
#       
#       if(i != length(Dates)){
#         p[[i]] <- p[[i]] + 
#           theme(strip.text.y = element_blank())
#       }
#       if(i != 1){
#         p[[i]] <- p[[i]] + 
#           theme(axis.title.y = element_blank(),
#                 axis.ticks.y = element_blank(),
#                 axis.text.y = element_blank())
#       }
#       
#       if(j != 1){
#         p[[i]] <- p[[i]] + 
#           theme(strip.text.x = element_blank())
#       }
#       if(j != length(Directions)){
#         p[[i]] <- p[[i]] + 
#           theme(axis.title.x = element_blank(),
#                 axis.ticks.x = element_blank(),
#                 axis.text.x = element_blank())
#       }
#       
#     }
#     p1[[j]] <- ggarrange(plotlist = p, nrow = 1, widths = c(1, 0.87, 0.87, 1))
#     
#   }
#   
#   p <- ggarrange(plotlist = p1, nrow = 2, heights = c(1, 1))
#   ggsave(p, file = paste0(outdir, '-intensity_transmission_idt', '.png'), w = 8, h = 5)
#   
# }

plot_intensity_PP <- function(intensity_PP, count_data, outdir){
  
  p <- ggplot(intensity_PP, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_raster(aes(fill = M)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
    geom_point(data = count_data[count > 0], aes(size = count), col = 'grey50') +
    theme_bw() + 
    labs(x = 'Age at infection recipient', fill = 'Estimated median\ntransmission rate', 
         y= 'Age at transmission source',size='Pairs\ncount') +
    geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
    facet_grid(label_direction~label_time) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_fill_viridis_c() + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    scale_size_continuous(range = c(1, 3), breaks = sort(unique(count_data[count > 0]$count))) +
    coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) + 
    guides(fill = guide_colorbar(order = 1), 
           shape = guide_legend(order = 2)) 
  ggsave(p, file = paste0(outdir, '-intensity_transmission.png'), w = 7, h = 7)
  
  Dates <- unique(intensity_PP$label_time)
  p = list(); i = 1
  for(Date in Dates){
    
    tmp <- subset(intensity_PP, label_time == Date)
    tmp1 <- subset(count_data, label_time == Date)
    
    p[[i]] <- ggplot(tmp, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      geom_point(data = tmp1[count > 0], aes(size = count), col = 'grey50') +
      theme_bw() + 
      labs(fill = 'Estimated\nmedian\ntransmission\nrate', y= 'Age at transmission source',size='Pairs\ncount') +
      geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(label_direction~label_time) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'right',
            axis.title.x = element_blank()) +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_size_continuous(range = c(1, 3), breaks = sort(unique(tmp1[count > 0]$count))) +
      coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) 
    
    if(i != 1){
      p[[i]] <- p[[i]] + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
    }
    
    # if(i != length(Dates)){
    #   p[[i]] <- p[[i]] + theme(strip.text.y = element_blank())
    # }
    
    i = i + 1
  }
  
  p <- grid.arrange(grobs = p, ncol =2, bottom = text_grob('Age at infection recipient'))
  ggsave(p, file = paste0(outdir, '-intensity_transmission_scaled.png'), w = 10, h = 7)
} 


plot_incident_cases <- function(incident_cases, outdir){
  
  p <- ggplot(incident_cases, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_raster(aes(fill = M)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient', fill = 'Estimated incident cases\nper year', 
         y= 'Age at transmission source',size='Pairs\ncount') +
    geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
    facet_grid(label_direction~label_time) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_fill_viridis_c() + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    guides(fill = guide_colorbar(order = 1), 
           shape = guide_legend(order = 2)) 
  ggsave(p, file = paste0(outdir, '-incident_cases.png'), w = 7, h = 7)
  
} 

plot_relative_intensity_PP <- function(relative_intensity_PP, outdir){
  
  p <- ggplot(relative_intensity_PP, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
    geom_raster(aes(fill = M)) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient', fill = 'Transmission flow', 
         y= 'Age at transmission source',size='Pairs\ncount') +
    geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
    facet_grid(label_direction~label_time) + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_fill_viridis_c() + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    guides(fill = guide_colorbar(order = 1), 
           shape = guide_legend(order = 2)) 
  ggsave(p, file = paste0(outdir, '-relative_intensity_PP.png'), w = 7, h = 7)
  
} 

plot_relative_intensity_PP_standardised <- function(relative_intensity_PP_aggregated, relative_incident_cases_aggregated, df_age_aggregated, outdir){
  
  relative_intensity_PP_aggregated[, type := 'Unstandardised']
  relative_incident_cases_aggregated[, type := 'Standardised']
  
  tmp <- rbind(relative_intensity_PP_aggregated, relative_incident_cases_aggregated)
  tmp <- merge(tmp, unique(df_age_aggregated[, .(age_group_transmission.SOURCE, age_group_infection.RECIPIENT, 
                                                 age_from.RECIPIENT, age_from.SOURCE)]), by = c('age_group_transmission.SOURCE', 'age_group_infection.RECIPIENT'))
  tmp  <- tmp[age_from.RECIPIENT != 5 & age_from.SOURCE != 5]
  tmp[, age_age_cat := paste0(age_group_transmission.SOURCE, '->', age_group_infection.RECIPIENT)]
  tmp[, `Age Recipient` := age_group_infection.RECIPIENT]
  tmp[, `Age Source` := age_group_transmission.SOURCE]
  
  tmp[age_from.RECIPIENT == age_from.SOURCE, age_source := 'Same age']
  tmp[age_from.RECIPIENT > age_from.SOURCE, age_source := 'Younger']
  tmp[age_from.RECIPIENT < age_from.SOURCE, age_source := 'Older']
  tmp[, age_source := factor(age_source, levels = c('Younger', 'Same age', 'Older'))]
  
  p <- ggplot(tmp[is_before_cutoff_date == 1 & is_mf == 1], aes(x = `Age Source`)) + 
    geom_bar(aes(y = M, fill = type), stat = 'identity', position = "dodge") + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = type), position = "dodge") + 
    labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
    theme_bw() +
    facet_grid(.~`Age Recipient`, scale = 'free', label = 'label_both')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom', 
          axis.text.x = element_text(angle = 70,hjust =1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) + 
    ggtitle('Male --> Female')
  ggsave(p, file = paste0(outdir, '-relative_intensity_PP_comparison_mf.png'), w = 7, h = 5.5)
  
  p <- ggplot(tmp[is_before_cutoff_date == 1 & is_mf == 0], aes(x = `Age Source`)) + 
    geom_bar(aes(y = M, fill = type), stat = 'identity', position = "dodge") + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = type), position = "dodge") + 
    labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
    theme_bw() +
    facet_grid(.~`Age Recipient`, scale = 'free', label = 'label_both')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom', 
          axis.text.x = element_text(angle = 70,hjust =1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) + 
    ggtitle('Female --> Male')
  ggsave(p, file = paste0(outdir, '-relative_intensity_PP_comparison_fm.png'), w = 7, h = 5.5)
  
  p <- ggplot(tmp[type == 'Standardised' & is_mf == 1], aes(x = `Age Source`)) + 
    geom_bar(aes(y = M, fill = label_time), stat = 'identity', position = "dodge") + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = label_time), position = "dodge") + 
    labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
    theme_bw() +
    facet_grid(.~`Age Recipient`, scale = 'free', label = 'label_both')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom', 
          axis.text.x = element_text(angle = 70,hjust =1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) + 
    ggtitle('Male --> Female') + 
    scale_fill_viridis_d()
  ggsave(p, file = paste0(outdir, '-relative_intensity_PP_standardised_mf.png'), w = 7, h = 5.5)
  
  p <- ggplot(tmp[type == 'Standardised' & is_mf == 0], aes(x = `Age Source`)) + 
    geom_bar(aes(y = M, fill = label_time), stat = 'identity', position = "dodge") + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = label_time), position = "dodge") + 
    labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
    theme_bw() +
    facet_grid(.~`Age Recipient`, scale = 'free', label = 'label_both')+
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom', 
          axis.text.x = element_text(angle = 70,hjust =1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) + 
    ggtitle('Female --> Male') + 
    scale_fill_viridis_d()
  ggsave(p, file = paste0(outdir, '-relative_intensity_PP_standardised_fm.png'), w = 7, h = 5.5)
}

plot_sex_source_standardised <- function(sex_source, sex_source_standardised, outdir){
  sex_source[, type := 'Unstandardised']
  sex_source_standardised[, type := 'Standardised']
  
  tmp <- rbind(sex_source, sex_source_standardised)
  tmp[is_mf == 1, sex_source := 'Male']
  tmp[is_mf == 0, sex_source := 'Female']
  
  p <- ggplot(tmp[is_before_cutoff_date == 1], aes(x = sex_source)) + 
    geom_bar(aes(y = M, fill = type), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = type), position = position_dodge()) + 
    theme_bw() + 
    labs(x = 'Sex source', y = 'Transmission flows', fill = '') + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom')
  ggsave(p, file = paste0(outdir, '-sex_source_comparison.png'), w = 6, h = 4)
  
  p <- ggplot(tmp[type == 'Standardised'], aes(x = sex_source)) + 
    geom_bar(aes(y = M, fill = label_time), stat = 'identity', position = position_dodge()) + 
    geom_errorbar(aes(ymin = CL, ymax = CU, group = label_time), position = position_dodge()) + 
    theme_bw() + 
    labs(x = 'Sex source', y = 'Transmission flows', fill = '') + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  + 
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') + 
    scale_fill_viridis_d()
  ggsave(p, file = paste0(outdir, '-sex_source_standardised.png'), w = 6, h = 4)
  
}

plot_median_age_source <- function(age_source, outdir){
  
  cat("\nPlot mean age at transmission of the source by age at infection of recipient\n")
  
  p <- ggplot(age_source) + 
    geom_line(aes(x = age_infection.RECIPIENT, y = M, col = label_time)) + 
    geom_ribbon(aes(x = age_infection.RECIPIENT, ymin= CL, ymax = CU, fill = label_time), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    geom_rect(data = range_age_observed, aes(xmin=max_age_infection.RECIPIENT, xmax=Inf, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    geom_rect(data = range_age_observed, aes(xmin=-Inf, xmax=min_age_infection.RECIPIENT, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    theme_bw() + 
    labs(x = 'Age at infection recipient', y = 'Median age at transmission source',
         col = 'Date infection recipient', fill = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +
    facet_grid(.~label_direction)
  
  ggsave(p, file = paste0(outdir, '-MedianAgeSource_ByAgeRecipient.png'), w = 7, h = 5)
  
}

plot_incident_cases_age_source <- function(incident_cases_age_source, outdir){
  
  p <- ggplot(incident_cases_age_source) + 
    geom_line(aes(x = age_transmission.SOURCE, y = M, col = label_time)) + 
    geom_ribbon(aes(x = age_transmission.SOURCE, ymin= CL, ymax = CU, fill = label_time), alpha = 0.15) + 
    geom_rect(data = range_age_observed, aes(xmin=max_age_infection.RECIPIENT, xmax=Inf, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    geom_rect(data = range_age_observed, aes(xmin=-Inf, xmax=min_age_infection.RECIPIENT, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    theme_bw() + 
    labs(x = 'Age at transmission source', y = 'Median number of incident cases per year',
         col = 'Date trasmission', fill = 'Date trasmission') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    coord_cartesian(xlim = range_age_non_extended) +
    facet_grid(.~label_direction)
  
  ggsave(p, file = paste0(outdir, '-incident_cases_age_source.png'), w = 7, h = 5)
  
}

plot_median_age_source_difference <- function(age_source_difference, outdir){

  cat("\nPlot difference age at transmission of the source to age at infection of recipient\n")
  
  p <- ggplot(age_source_difference) + 
    geom_line(aes(x = age_infection.RECIPIENT, y = M, col = label_time)) + 
    geom_ribbon(aes(x = age_infection.RECIPIENT, ymin= CL, ymax = CU, fill = label_time), alpha = 0.15) + 
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey50') + 
    geom_rect(data = range_age_observed, aes(xmin=max_age_infection.RECIPIENT, xmax=Inf, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    geom_rect(data = range_age_observed, aes(xmin=-Inf, xmax=min_age_infection.RECIPIENT, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    theme_bw() + 
    labs(x = 'Age at infection recipient', y = 'Median difference age at transmission source\nto age at infection recipient',
         col = 'Date infection recipient', fill = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    coord_cartesian(xlim = range_age_non_extended) +
    facet_grid(.~label_direction)
  
  ggsave(p, file = paste0(outdir, '-MedianAgeSourceDifference_ByAgeRecipient.png'), w = 7, h = 5)
  
}

plot_median_age_source_overall <- function(age_source_overall, outdir){
  
  cat("\nPlot mean age at transmission of the source overall\n")
  
  p <- ggplot(age_source_overall, aes(x = label_time)) + 
    geom_point(aes(y = M)) + 
    geom_errorbar(aes(ymin= CL, ymax = CU), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    labs(y = "Median age at transmission source\nadjusted by the recipient's incidence rate", x = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    facet_grid(.~label_direction) 
  
  ggsave(p, file = paste0(outdir, '-MedianAgeSource.png'), w = 7, h = 5)
}

plot_median_age_recipient <- function(age_recipient, outdir){
  
  cat("\nPlot mean age at infection of the recipient by age at transmission of source\n")
  
  p <- ggplot(age_recipient) + 
    geom_line(aes(x = age_transmission.SOURCE, y = M, col = label_time)) + 
    geom_ribbon(aes(x = age_transmission.SOURCE, ymin= CL, ymax = CU, fill = label_time), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    geom_rect(data = range_age_observed, aes(xmin=max_age_transmission.SOURCE, xmax=Inf, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    geom_rect(data = range_age_observed, aes(xmin=-Inf, xmax=min_age_transmission.SOURCE, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    theme_bw() + 
    labs(x = 'Age at transmission source', y = 'Median age at infection recipient',
         col = 'Date infection recipient', fill = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +    
    facet_grid(.~label_direction)
  
  ggsave(p, file = paste0(outdir, '-MedianAgeRecipient_ByAgeSource.png'), w = 7, h = 5)
  
}

plot_median_age_recipient_difference <- function(age_recipient_difference, outdir){
  
  cat("\nPlot difference age at transmission of the source to age at infection of recipient\n")
  
  p <- ggplot(age_recipient_difference) + 
    geom_line(aes(x = age_transmission.SOURCE, y = M, col = label_time)) + 
    geom_ribbon(aes(x = age_transmission.SOURCE, ymin= CL, ymax = CU, fill = label_time), alpha = 0.15) + 
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey50') + 
    geom_rect(data = range_age_observed, aes(xmin=max_age_transmission.SOURCE, xmax=Inf, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    geom_rect(data = range_age_observed, aes(xmin=-Inf, xmax=min_age_transmission.SOURCE, 
                                             ymin=-Inf, ymax=Inf), alpha=.2) +
    theme_bw() + 
    labs(x = 'Age at transmission source', y = 'Median difference age at infection recipient\nto age at transmission source',
         col = 'Date infection recipient', fill = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    facet_grid(.~label_direction) +
    coord_cartesian(xlim = range_age_non_extended) 
  
  ggsave(p, file = paste0(outdir, '-MedianAgeRecipientDifference_ByAgeSource.png'), w = 7, h = 5)
  
}

plot_median_age_recipient_overall <- function(age_recipient_overall, outdir){
  
  cat("\nPlot mean age at transmission of the source overall\n")
  
  p <- ggplot(age_recipient_overall, aes(x = label_time)) + 
    geom_point(aes(y = M)) + 
    geom_errorbar(aes(ymin= CL, ymax = CU), alpha = 0.15) + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    labs(y = "Median age at infection recipient\nadjusted by the source's incidence rate", x = 'Date infection recipient') +
    theme(strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)),
          legend.position = 'bottom') +
    facet_grid(.~label_direction)
  
  ggsave(p, file = paste0(outdir, '-MedianAgeRecipient.png'), w = 7, h = 5)
}
