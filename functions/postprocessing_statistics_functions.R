save_statistics_expected_contribution <- function(expected_contribution_sex_source, median_age_source, outdir){
  
  n_digits <- 1
  
  # contribution of men in first and last round 
  ecm <- expected_contribution_sex_source
  ecm[, MAX_INDEX_ROUND := max(INDEX_ROUND), by = c('COMM')]
  ecm <- ecm[INDEX_ROUND == 1 | INDEX_ROUND == MAX_INDEX_ROUND, .(M, CL, CU), by = c('COMM', 'LABEL_SOURCE', 'ROUND', 'MAX_INDEX_ROUND', 'INDEX_ROUND')]
  ecm <- ecm[, .(M = round(M*100, n_digits), CL = round(CL * 100, n_digits), CU = round(CU * 100, n_digits)), by = c('COMM', 'LABEL_SOURCE', 'ROUND')]
  ecm <- ecm[order(COMM, LABEL_SOURCE, ROUND)]
  
  # age median of contribution 
  median_age <- copy(median_age_source)
  median_age[, MAX_INDEX_ROUND := max(INDEX_ROUND), by = c('COMM')]
  median_age <- median_age[quantile == 'C50']
  median_age <- median_age[INDEX_ROUND == 1 | INDEX_ROUND == MAX_INDEX_ROUND, .(M, CL, CU), by = c('COMM', 'LABEL_SOURCE', 'ROUND', 'MAX_INDEX_ROUND', 'INDEX_ROUND')]
  median_age <- median_age[order(COMM, LABEL_SOURCE, ROUND)]
  
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
  stats[['pairs_MAE']] <- tmp[, mean(MAE) ]
  
  # incidence rate
  tmp <- merge(predict_incidence_rate_round, incidence_cases_recipient_round[, .(INDEX_DIRECTION, INDEX_COMMUNITY, ROUND, 
                                                                                 AGE_INFECTION.RECIPIENT, INCIDENCE, LB, UB)], 
               by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'ROUND', 'AGE_INFECTION.RECIPIENT'))
  tmp[, within.CI := INCIDENCE >= CL & INCIDENCE <= CU]
  tmp[, MAE := abs(INCIDENCE - M)]
  
  stats[['incidence_WCI']] <- tmp[, round(mean(within.CI) *100, 2)]
  stats[['incidence_MAE']] <- tmp[, mean(MAE) ]
  
  saveRDS(stats, file = paste0(outdir, '-statistics_prediction.RDS'))
}

make_transmission_flows_table <- function(expected_contribution_age_classification_male, expected_contribution_sex_age_group_recipient, 
                                          expected_contribution_sex_age_classification_male, expected_contribution_age_group_recipient, 
                                          expected_contribution_age_classification_male_total, expected_contribution_sex_source, outdir){
  
  
  table <- list()
  
  Rounds <- c('R010', 'R015', 'R018')
  
  n_digits <- 1
  
  #
  # transmission flows by gender-direction age group recipient and age classification male
  #
  
  saa <- copy(expected_contribution_age_classification_male)
  saa <- saa[, .(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT, AGE_CLASSIFICATION.MALE, M, CL, CU)]
  saa[, M := format(round(M * 100, n_digits),  nsmall = n_digits) ]
  saa[, CI := paste0('[', format(round(CL * 100, n_digits),  nsmall = n_digits), '-', format(round(CU * 100, n_digits),  nsmall = n_digits), ']') ]
  saa[, CI := gsub(' ', '', CI)]
  set(saa, NULL, "CL", NULL)
  set(saa, NULL, "CU", NULL)
  
  saa <- saa[order(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT, AGE_CLASSIFICATION.MALE)]
  
  tmp <- dcast.data.table(saa, COMM + ROUND + LABEL_DIRECTION + AGE_CLASSIFICATION.MALE~ AGE_GROUP_INFECTION.RECIPIENT, value.var = 'M')
  names(tmp)[grepl('-', names(tmp))] <- paste0('M_', names(tmp)[grepl('-', names(tmp))])
  saa <- dcast.data.table(saa, COMM + ROUND + LABEL_DIRECTION + AGE_CLASSIFICATION.MALE~ AGE_GROUP_INFECTION.RECIPIENT, value.var = 'CI')
  names(saa)[grepl('-', names(tmp))] <- paste0('CI_', names(tmp)[grepl('-', names(tmp))])
  saa <- merge(tmp, saa, by = c('COMM', 'ROUND', 'LABEL_DIRECTION', 'AGE_CLASSIFICATION.MALE'))
  
  saa <- saa[ROUND %in% Rounds]
  
  table[['flows_gender_ager_agem']] <- saa
  
  
  #
  # transmission flows by gender-direction age group recipient
  #
  
  saa <- copy(expected_contribution_sex_age_group_recipient)
  saa <- saa[, .(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT, M, CL, CU)]
  saa[, M := format(round(M * 100, n_digits),  nsmall = n_digits) ]
  saa[, CI := paste0('[', format(round(CL * 100, n_digits),  nsmall = n_digits), '-', format(round(CU * 100, n_digits),  nsmall = n_digits), ']') ]
  saa[, CI := gsub(' ', '', CI)]
  set(saa, NULL, "CL", NULL)
  set(saa, NULL, "CU", NULL)
  
  saa <- saa[order(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT)]
  
  tmp <- dcast.data.table(saa, COMM + ROUND + LABEL_DIRECTION ~ AGE_GROUP_INFECTION.RECIPIENT, value.var = 'M')
  names(tmp)[grepl('-', names(tmp))] <- paste0('M_', names(tmp)[grepl('-', names(tmp))])
  saa <- dcast.data.table(saa, COMM + ROUND + LABEL_DIRECTION ~ AGE_GROUP_INFECTION.RECIPIENT, value.var = 'CI')
  names(saa)[grepl('-', names(tmp))] <- paste0('CI_', names(tmp)[grepl('-', names(tmp))])
  saa <- merge(tmp, saa, by = c('COMM', 'ROUND', 'LABEL_DIRECTION'))
  
  saa <- saa[ROUND %in% Rounds]
  
  table[['flows_gender_ager']] <- saa
  
  
  #
  # transmission flows by gender-direction age classification male 
  #
  
  saa <- copy(expected_contribution_sex_age_classification_male)
  saa <- saa[, .(COMM, ROUND, LABEL_DIRECTION, AGE_CLASSIFICATION.MALE, M, CL, CU)]
  saa[, M := format(round(M * 100, n_digits),  nsmall = n_digits) ]
  saa[, CI := paste0('[', format(round(CL * 100, n_digits),  nsmall = n_digits), '-', format(round(CU * 100, n_digits),  nsmall = n_digits), ']') ]
  saa[, CI := gsub(' ', '', CI)]
  set(saa, NULL, "CL", NULL)
  set(saa, NULL, "CU", NULL)
  
  saa <- saa[order(COMM, ROUND, LABEL_DIRECTION, AGE_CLASSIFICATION.MALE)]
  
  saa <- saa[ROUND %in% Rounds]
  
  table[['flows_gender_agem']] <- saa
  
  
  #
  # transmission flows by age group recipient 
  #
  
  saa <- copy(expected_contribution_age_group_recipient)
  saa <- saa[, .(COMM, ROUND, AGE_GROUP_INFECTION.RECIPIENT, M, CL, CU)]
  saa[, M := format(round(M * 100, n_digits),  nsmall = n_digits) ]
  saa[, CI := paste0('[', format(round(CL * 100, n_digits),  nsmall = n_digits), '-', format(round(CU * 100, n_digits),  nsmall = n_digits), ']') ]
  saa[, CI := gsub(' ', '', CI)]
  set(saa, NULL, "CL", NULL)
  set(saa, NULL, "CU", NULL)
  
  saa <- saa[order(COMM, ROUND, AGE_GROUP_INFECTION.RECIPIENT)]
  
  tmp <- dcast.data.table(saa, COMM + ROUND ~ AGE_GROUP_INFECTION.RECIPIENT, value.var = 'M')
  names(tmp)[grepl('-', names(tmp))] <- paste0('M_', names(tmp)[grepl('-', names(tmp))])
  saa <- dcast.data.table(saa, COMM + ROUND ~ AGE_GROUP_INFECTION.RECIPIENT, value.var = 'CI')
  names(saa)[grepl('-', names(tmp))] <- paste0('CI_', names(tmp)[grepl('-', names(tmp))])
  saa <- merge(tmp, saa, by = c('COMM', 'ROUND'))
  
  saa <- saa[ROUND %in% Rounds]
  
  table[['flows_ager']] <- saa
  
  
  #
  # transmission flows by age classification male
  #
  
  saa <- copy(expected_contribution_age_classification_male_total)
  saa <- saa[, .(COMM, ROUND, AGE_CLASSIFICATION.MALE, M, CL, CU)]
  saa[, M := format(round(M * 100, n_digits),  nsmall = n_digits) ]
  saa[, CI := paste0('[', format(round(CL * 100, n_digits),  nsmall = n_digits), '-', format(round(CU * 100, n_digits),  nsmall = n_digits), ']') ]
  saa[, CI := gsub(' ', '', CI)]
  set(saa, NULL, "CL", NULL)
  set(saa, NULL, "CU", NULL)
  
  saa <- saa[order(COMM, ROUND, AGE_CLASSIFICATION.MALE)]
  
  saa <- saa[ROUND %in% Rounds]
  
  table[['flows_agem']] <- saa
  
  
  #
  # transmission flows by gender-direction 
  #
  
  saa <- copy(expected_contribution_sex_source)
  saa <- saa[, .(COMM, ROUND, LABEL_DIRECTION, M, CL, CU)]
  saa[, M := format(round(M * 100, n_digits),  nsmall = n_digits) ]
  saa[, CI := paste0('[', format(round(CL * 100, n_digits),  nsmall = n_digits), '-', format(round(CU * 100, n_digits),  nsmall = n_digits), ']') ]
  saa[, CI := gsub(' ', '', CI)]
  set(saa, NULL, "CL", NULL)
  set(saa, NULL, "CU", NULL)
  
  saa <- saa[order(COMM, ROUND, LABEL_DIRECTION)]
  
  saa <- saa[ROUND %in% Rounds]
  
  table[['flows_gender']] <- saa
  
  saveRDS(table, file = paste0(outdir, '-transmission_flows_table.rds'))
  
  return(table)
}

save_median_age_diff <- function(median_age_diff_group, outdir){
  
  table <- list()
  
  #
  # median age difference male to female
  #
  
  mad <- copy(median_age_diff_group)
  mad <- mad[ROUND == "R018" & quantile == 'C50']
  mad <- mad[, .(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT, quantile, M, CL, CU)]
  mad <- mad[order(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT)]
  table[['median_age_difference']] <- mad
  
  #
  #Save
  #
  
  saveRDS(table, paste0(outdir, '-output-median_age_difference.rds'))
}


