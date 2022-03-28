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

plot_intensity_PP <- function(intensity_PP, count_data, range_age_observed, outdir){
  
  intensity_PP <- intensity_PP[age_transmission.SOURCE >= range_age_observed[, min_age] & age_transmission.SOURCE <= range_age_observed[, max_age]]
  intensity_PP <- intensity_PP[age_infection.RECIPIENT >= range_age_observed[, min_age] & age_infection.RECIPIENT <= range_age_observed[, max_age]]
  
  communities <- intensity_PP[, unique(comm)]
  for(i in seq_along(communities)){
    tmp <- intensity_PP[ comm == communities[i]]
    tmp1 <- count_data[ comm == communities[i]]
    p <- ggplot(tmp, aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
      geom_raster(aes(fill = M)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
      geom_point(data = tmp1[count > 0], aes(size = count), col = 'grey50') +
      theme_bw() + 
      labs(x = 'Age at infection recipient', fill = 'Estimated median\ntransmission rate', 
           y= 'Age at transmission source',size='Pairs\ncount') +
      # geom_contour(aes(z = M), col = 'red', alpha = 0.8, bins = 5) + 
      facet_grid(label_direction~label_time) + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_fill_viridis_c() + 
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_size_continuous(range = c(1, 3), breaks = sort(unique(tmp1[count > 0]$count))) +
      coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) + 
      guides(fill = guide_colorbar(order = 1), 
             shape = guide_legend(order = 2)) + 
      ggtitle(tmp[,unique(label_community)])
    
    ggsave(p, file = paste0(outdir, '-intensity_transmission_', communities[i], '.png'), w = 7, h = 7)
  }

} 


plot_transmission_flows <- function(transmission_flows, range_age_observed, outdir, lab=NULL, count_data = NULL, with_contour = F){
  
  communities <- transmission_flows[, unique(comm)]
  
  for(i in seq_along(communities)){
    
    tmp <- transmission_flows[ comm == communities[i]]
    
    index_groups <- tmp[, sort(unique(index_group))]
    
    levels <- tmp[, {
      prob = c(0.5, 0.8, 0.9)
      level = getLevel(age_infection.RECIPIENT, age_transmission.SOURCE, M, prob)
      list(prob = prob, level = level)
    }, by = c('index_group')]
    
    plots <- list()
    
    for(j in seq_along(index_groups)){
      
      p <- ggplot(tmp[index_group == index_groups[j]], aes(y = age_transmission.SOURCE, x = age_infection.RECIPIENT)) + 
        geom_raster(aes(fill = M)) + 
        geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'white') + 
        theme_bw() + 
        labs(fill = paste0(lab, '\nTransmission\nflows'), 
             size='Pairs\ncount') +
        facet_grid(label_direction~label_time) + 
        theme(strip.background = element_rect(colour="white", fill="white"),
              strip.text = element_text(size = rel(1)),
              legend.position = 'right', 
              axis.title.y= element_blank(), 
              axis.title.x= element_blank()) +
        scale_fill_viridis_c() + 
        scale_x_continuous(expand = c(0,0)) + 
        scale_y_continuous(expand = c(0,0)) + 
        guides(fill = guide_colorbar(order = 1), 
               shape = guide_legend(order = 2)) 
      
      
      if(tmp[index_group == index_groups[j], unique(is_mf)]){
        p <- p + 
          theme(strip.text.x = element_blank())
      }
      
      if(tmp[index_group == index_groups[j], unique(is_before_cutoff_date)]){
        p <- p + 
          theme(strip.text.y = element_blank())
      }
      
      if(!tmp[index_group == index_groups[j], unique(is_mf)]){
        p <- p + 
          theme(axis.text.x  = element_blank(), 
                axis.title.x = element_blank(), 
                axis.ticks.x = element_blank())
      }
      
      if(!tmp[index_group == index_groups[j], unique(is_before_cutoff_date)]){
        p <- p + 
          theme(axis.text.y = element_blank(), 
                axis.title.y = element_blank(), 
                axis.ticks.y = element_blank())
      }
      
      if(with_contour){

        tmp1 <- merge(levels[index_group == index_groups[j]], tmp, by = 'index_group', allow.cartesian=TRUE)
        tmp1[, diff := abs(M - level)]
        mindiff <- tmp1[, list(diff = min(diff)), by = 'prob']
        tmp1 <- merge(tmp1, mindiff, by = c('prob', 'diff'))
        tmp1[, prob_label := paste0(prob * 100, '%')]

        p <- p +
          geom_contour(aes(z=M, col = ..level..), breaks = tmp1[, level]) +
          geom_text(data = tmp1, aes(label = prob_label, col = level), size = 4) + 
            scale_color_gradient(low = 'darkred', high = 'lightpink2') + 
            guides(col="none")
      }
      
      if(!is.null(count_data)){
        tmp1 <- count_data[ comm == communities[i] & index_group == index_groups[j]]
        p <- p + 
          geom_point(data = tmp1[count > 0], aes(size = count), col = 'grey50', alpha = 0.5) + 
          scale_size_continuous(limits = count_data[  count > 0, range(count)], 
                                breaks = count_data[  count > 0, sort(unique(count))])
      }
      
      plots[[j]] <- p
      
    }

     p <- ggarrange(plotlist = plots[c(1,3,2,4)], common.legend = T, legend = 'right')
     pf <- grid.arrange(p, bottom = 'Age at infection recipient                     ', 
                  left = 'Age at transmission source',
                  top = tmp[,unique(label_community)])

    ggsave(pf, file = paste0(outdir, '-transmission_flows', lab, '_', communities[i], '.png'), w = 8, h = 6)
  }
  
} 

getLevel <- function(x,y,z,prob=0.95) {
  dx <- diff(unique(x)[1:2])
  dy <- diff(unique(y)[1:2])
  sz <- sort(z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

plot_transmission_flows_aggregated <- function(transmission_flows_aggregated, standardised_transmission_flows_aggregated, df_age_aggregated, outdir)
{
  
  transmission_flows_aggregated[, type := 'Unstandardised']
  standardised_transmission_flows_aggregated[, type := 'Standardised']
  
  tmp <- rbind(transmission_flows_aggregated, standardised_transmission_flows_aggregated)
  
  tmp[, `Age Recipient` := age_group_infection.RECIPIENT]
  tmp[, `Age Source` := age_group_transmission.SOURCE]
  
  tmp[, Direction :=label_direction]
  
  communities <- tmp[, unique(comm)]
  for(i in seq_along(communities)){
    tmp1 <- tmp[ comm == communities[i]]
    
    p <- ggplot(tmp1[is_before_cutoff_date == 1], aes(x = `Age Source`)) + 
      geom_bar(aes(y = M, fill = type), stat = 'identity', position = "dodge") + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = type), position = "dodge") + 
      labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, scale = 'free', label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.text.x = element_text(angle = 70,hjust =1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) +
      ggtitle(tmp1[,unique(label_community)])
    
    ggsave(p, file = paste0(outdir, '-transmission_flows_aggregated_comparison_',  communities[i], '.png'), w = 7, h = 5.5)
    
    p <- ggplot(tmp1[type == 'Standardised'], aes(x = `Age Source`)) + 
      geom_bar(aes(y = M, fill = label_time), stat = 'identity', position = "dodge") + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = label_time), position = "dodge") + 
      labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, scale = 'free', label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.text.x = element_text(angle = 70,hjust =1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) +
      scale_fill_viridis_d() +
      ggtitle(tmp1[,unique(label_community)])
    ggsave(p, file = paste0(outdir, '-transmission_flows_aggregated_',  communities[i], '.png'), w = 7, h = 5.5)
  }
}

plot_transmission_flows_aggregated2 <- function(transmission_flows_aggregated, standardised_transmission_flows_aggregated, df_age_aggregated, outdir)
{
  
  transmission_flows_aggregated[, type := 'Unstandardised']
  standardised_transmission_flows_aggregated[, type := 'Standardised']
  
  tmp <- rbind(transmission_flows_aggregated, standardised_transmission_flows_aggregated)
  
  tmp[, `Age Recipient` := age_group_infection.RECIPIENT]
  tmp[, `Age Source` := factor(age_classification.SOURCE, levels = c('Younger', 'Same age', 'Older'))]
  
  tmp[, Direction := label_direction]
  
  communities <- tmp[, unique(comm)]
  for(i in seq_along(communities)){
    tmp1 <- tmp[ comm == communities[i]]
    
    p <- ggplot(tmp1[is_before_cutoff_date == 1], aes(x = `Age Source`)) + 
      geom_bar(aes(y = M, fill = type), stat = 'identity', position = "dodge") + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = type), position = "dodge") + 
      labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, scale = 'free', label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.text.x = element_text(angle = 70,hjust =1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) +
      ggtitle(tmp1[,unique(label_community)])
    ggsave(p, file = paste0(outdir, '-transmission_flows_aggregated2_comparison_', communities[i], '.png'), w = 7, h = 5.75)
    
    p <- ggplot(tmp1[type == 'Standardised'], aes(x = `Age Source`)) + 
      geom_bar(aes(y = M, fill = label_time), stat = 'identity', position = "dodge") + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = label_time), position = "dodge") + 
      labs(x = 'Age source', y = 'Transmission flows', fill = '') + 
      theme_bw() +
      facet_grid(Direction~`Age Recipient`, scale = 'free', label = 'label_both')+
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom', 
            axis.text.x = element_text(angle = 70,hjust =1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = tmp[, range(c(CU, 0))]) +
      scale_fill_viridis_d()+
      ggtitle(tmp1[,unique(label_community)])
    ggsave(p, file = paste0(outdir, '-transmission_flows_aggregated2_', communities[i], '.png'), w = 7, h = 5.75)
  }
}

plot_sex_source_standardised <- function(sex_source, sex_source_standardised, outdir){
  
  sex_source[, type := 'Unstandardised']
  sex_source_standardised[, type := 'Standardised']
  
  tmp <- rbind(sex_source, sex_source_standardised)
  tmp[is_mf == 1, sex_source := 'Male']
  tmp[is_mf == 0, sex_source := 'Female']
  
  communities <- tmp[, unique(comm)]
  for(i in seq_along(communities)){
    tmp1 <- tmp[ comm == communities[i]]
    
    p <- ggplot(tmp1[is_before_cutoff_date == 1], aes(x = sex_source)) + 
      geom_bar(aes(y = M, fill = type), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = type), position = position_dodge()) + 
      theme_bw() + 
      labs(x = 'Sex source', y = 'Transmission flows', fill = '') + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') + 
      ggtitle(tmp1[, unique(label_community)])
    ggsave(p, file = paste0(outdir, '-sex_source_comparison_', communities[i], '.png'), w = 6, h = 4)
    
    p <- ggplot(tmp1[type == 'Standardised'], aes(x = sex_source)) + 
      geom_bar(aes(y = M, fill = label_time), stat = 'identity', position = position_dodge()) + 
      geom_errorbar(aes(ymin = CL, ymax = CU, group = label_time), position = position_dodge()) + 
      theme_bw() + 
      labs(x = 'Sex source', y = 'Transmission flows', fill = '') + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))  + 
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') + 
      scale_fill_viridis_d() + 
      ggtitle(tmp1[, unique(label_community)])
    ggsave(p, file = paste0(outdir, '-sex_source_standardised_', communities[i], '.png'), w = 6, h = 4)
  }
  
  
}

plot_median_age_source <- function(age_source, outdir){
  
  cat("\nPlot mean age at transmission of the source by age at infection of recipient\n")
  
  communities <- age_source[, unique(comm)]
  
  for(i in seq_along(communities)){
    tmp <- age_source[comm == communities[i]]
    
    p <- ggplot(tmp) + 
      geom_line(aes(x = age_infection.RECIPIENT, y = M)) + 
      geom_ribbon(aes(x = age_infection.RECIPIENT, ymin= CL, ymax = CU), alpha = 0.15) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
      theme_bw() + 
      labs(x = 'Age at infection recipient', y = 'Median age at transmission source',
           col = 'Date infection recipient', fill = 'Date infection recipient') +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +
      facet_grid(label_time~label_direction) + 
      ggtitle(tmp[, unique(label_community)])
    
    ggsave(p, file = paste0(outdir, '-MedianAgeSource_ByAgeRecipient_', communities[i], '.png'), w = 7, h = 6)
  }

  
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

plot_median_age_recipient <- function(age_recipient, age_recipient_standardised, outdir){
  
  cat("\nPlot mean age at infection of the recipient by age at transmission of source\n")
  
  age_recipient[, type := 'Unstandardised']
  age_recipient_standardised[, type := 'Standardised']
  
  tmp <- rbind(age_recipient, age_recipient_standardised)
  
  communities <- tmp[, unique(comm)]
  
  for(i in seq_along(communities)){
    tmp1 <- tmp[comm == communities[i]]
    
    p <- ggplot(tmp1) + 
      geom_line(aes(x = age_transmission.SOURCE, y = M, col = type)) + 
      geom_ribbon(aes(x = age_transmission.SOURCE, ymin= CL, ymax = CU, fill = type), alpha = 0.15) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
      theme_bw() + 
      labs(x = 'Age at transmission source', y = 'Median age at infection recipient',
           col = 'Date infection recipient', fill = 'Date infection recipient') +
      theme(strip.background = element_rect(colour="white", fill="white"),
            strip.text = element_text(size = rel(1)),
            legend.position = 'bottom') +
      scale_x_continuous(expand = c(0,0)) + 
      scale_y_continuous(expand = c(0,0)) + 
      coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +    
      facet_grid(label_time~label_direction) + 
      ggtitle(tmp1[, unique(label_community)])
    
    ggsave(p, file = paste0(outdir, '-MedianAgeRecipient_ByAgeSource_', communities[i], '.png'), w = 7, h = 5)
    
  }

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
