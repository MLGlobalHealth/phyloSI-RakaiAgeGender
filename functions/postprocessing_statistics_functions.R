save_statistics_expected_contribution <- function(expected_contribution_sex_source, expected_contribution_age_source2, outdir){
  
  # contribution of men in first and last round 
  ecm <- expected_contribution_sex_source[LABEL_SOURCE == 'Male sources']
  ecm[, MAX_INDEX_ROUND := max(INDEX_ROUND), by = c('COMM')]
  ecm <- ecm[INDEX_ROUND == 1 | INDEX_ROUND == MAX_INDEX_ROUND, .(M, CL, CU), by = c('COMM', 'ROUND', 'MAX_INDEX_ROUND', 'INDEX_ROUND')]
  ecm <- ecm[, .(M = round(M*100, 2), CL = round(CL * 100, 2), CU = round(CU * 100, 2)), by = c('COMM', 'ROUND')]
  ecm <- ecm[order(COMM, ROUND)]
  
  # age median of contribution 
  median_age <- copy(expected_contribution_age_source2)
  median_age[, WEIGHT_CONTRIBUTION := M / sum(M), by = c('LABEL_SOURCE', 'LABEL_ROUND', 'COMM', 'INDEX_ROUND', 'ROUND')]
  median_age <- median_age[, list(AGE_MEDIAN_CONTRIBUTION = matrixStats::weightedMedian(AGE_TRANSMISSION.SOURCE, WEIGHT_CONTRIBUTION )), 
                           by = c('LABEL_SOURCE', 'LABEL_ROUND', 'COMM', 'INDEX_ROUND', 'ROUND')]
  median_age[, MAX_INDEX_ROUND := max(INDEX_ROUND), by = c('COMM')]
  median_age <- median_age[INDEX_ROUND == 1 | INDEX_ROUND == MAX_INDEX_ROUND, .(AGE_MEDIAN_CONTRIBUTION), 
                           by = c('COMM', 'ROUND', 'MAX_INDEX_ROUND', 'INDEX_ROUND', 'LABEL_SOURCE')]
  median_age <- median_age[, .(AGE_MEDIAN_CONTRIBUTION = round(AGE_MEDIAN_CONTRIBUTION, 2)), by = c('COMM', 'ROUND', 'LABEL_SOURCE')]
  median_age <- median_age[order(COMM, ROUND, LABEL_SOURCE)]
  
  # save
  stats <- list()
  stats$contribution_male <- ecm
  stats$median_age <- median_age
  
  saveRDS(stats, paste0(outdir, '-output-contribution.rds'))
  
}

save_statistics_PPC <- function(predict_y_source_recipient, count_data, predict_incidence_rate_round, incidence_cases_recipient_round, outdir){
  
  stats <- list()
  
  # pairs
  data <- count_data[, list(count = sum(count)), by = c('LABEL_SOURCE', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_TRANSMISSION.SOURCE', 'PERIOD_SPAN', 'AGE_INFECTION.RECIPIENT')]
  
  tmp <- merge(predict_y_source_recipient, data, by = c('LABEL_SOURCE', 'LABEL_COMMUNITY', 'PERIOD', 'AGE_TRANSMISSION.SOURCE', 'PERIOD_SPAN', 'AGE_INFECTION.RECIPIENT'))
  tmp[, within.CI := count >= CL & count <= CU]
  tmp[, MAE := abs(count - M)]
  
  stats[['pairs_WCI']] <- tmp[, round(mean(within.CI) *100, 2)]
  stats[['pairs_MAE']] <- tmp[, round(mean(MAE) , 2)]
  
  # incidence rate
  tmp <- merge(predict_incidence_rate_round, incidence_cases_recipient_round[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, 
                                                                                 AGE_INFECTION.RECIPIENT, INCIDENCE, LB, UB)], 
               by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'))
  tmp[, within.CI := INCIDENCE >= CL & INCIDENCE <= CU]
  tmp[, MAE := abs(INCIDENCE - M)]
  
  stats[['incidence_WCI']] <- tmp[, round(mean(within.CI) *100, 2)]
  stats[['incidence_MAE']] <- tmp[, round(mean(MAE) , 2)]
  
  saveRDS(stats, file = paste0(outdir, '-statistics_prediction.RDS'))
}
