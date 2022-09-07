

save_statistics_transmission_events <- function(pairs, outdir){
  
  stat <- list()
  
  stat[[1]] <- nrow(pairs)
  stat[[2]] <- format(pairs[, range(DATE_INFECTION.RECIPIENT)], '%B %d, %Y')
  
  print(stat)
  
  file = paste0(outdir, '-data-pairs_stat.rds')
  saveRDS(stat, file)
}
