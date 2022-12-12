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

save_median_age_source <- function(median_age_source_group, median_age_source_group2, 
                                   median_age_source, outdir){
  table <- list()
  
  #
  # Median age source
  #
  
  mas <- copy(median_age_source_group)
  mas <- mas[ROUND %in% c('R010', 'R018') & quantile == 'C50']
  mas <- mas[, .(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT, M, CL, CU)]
  mas <- mas[order(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT)]
  
  table[['median_age_source_agg']] <- mas
  
  #
  # Median age source other aggregate
  #
  
  mas <- copy(median_age_source_group2)
  mas <- mas[ROUND %in% c('R010', 'R018') & quantile == 'C50']
  mas <- mas[, .(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT, M, CL, CU)]
  mas <- mas[order(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT)]
  
  table[['median_age_source_agg2']] <- mas
  
  #
  # Median age source total
  #
  
  mas <- copy(median_age_source)
  mas <- mas[ROUND %in% c('R010', 'R018') & quantile == 'C50']
  mas <- mas[, .(COMM, ROUND, LABEL_SOURCE, M, CL, CU)]
  mas <- mas[order(COMM, ROUND, LABEL_SOURCE)]
  
  table[['median_age_source']] <- mas
  
  #
  # Save
  #
  
  saveRDS(table, paste0(outdir, '-output-median_age_source.rds'))
  
}

save_counterfactual_results <- function(counterfactuals_p_f, counterfactuals_p_f05, 
                                        counterfactuals_p_959595, counterfactuals_p_909090, 
                                        lab, outdir){
  #
  # labels
  #
  
  # label in scenario where treated male take up art as much as female
  label.f = 'ART coverage in men as in women\nSuppression rate in men as in women'
  if(!grepl('Diagnosed', lab)){
    label.f <- paste0('Diagnosed rate in men as in women\n', label.f)
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
  
  # relative incident cases counterfactual compared to factual by age group 
  relative_incidence_counterfactual_group  <- counterfactuals_p_f$relative_incidence_counterfactual_group2
  relative_incidence_counterfactual_group[, label := label.f]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all <- counterfactuals_p_f$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all[, label := label.f]
  
  # sex ratio incident cases counterfactual  by age group
  sex_ratio_incidence_counterfactual_group <- counterfactuals_p_f$sex_ratio_incidence_counterfactual_group2 
  sex_ratio_incidence_counterfactual_group[, label := label.f]
  
  # sex ratio incident cases counterfactual  by age group regardless of age
  sex_ratio_incidence_counterfactual_all <- counterfactuals_p_f$sex_ratio_incidence_counterfactual_all 
  sex_ratio_incidence_counterfactual_all[, label := label.f]
  
  
  #
  # unlist objects in scenario where male are diagnosed/treated/suppressed half as much as female
  #
  
  # number of male treated regardless of age
  budget.counterfactual.f05 <- counterfactuals_p_f05$budget 
  budget.counterfactual.f05[, label := label.f05]
  
  # relative incident cases counterfactual compared to factual by age group 
  relative_incidence_counterfactual_group.f05  <- counterfactuals_p_f05$relative_incidence_counterfactual_group2
  relative_incidence_counterfactual_group.f05[, label := label.f05]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all.f05 <- counterfactuals_p_f05$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all.f05[, label := label.f05]
  
  # sex ratio incident cases counterfactual  by age group
  sex_ratio_incidence_counterfactual_group.f05 <- counterfactuals_p_f05$sex_ratio_incidence_counterfactual_group2 
  sex_ratio_incidence_counterfactual_group.f05[, label := label.f05]
  
  # sex ratio incident cases counterfactual  by age group regardless of age
  sex_ratio_incidence_counterfactual_all.f05 <- counterfactuals_p_f05$sex_ratio_incidence_counterfactual_all 
  sex_ratio_incidence_counterfactual_all.f05[, label := label.f05]
  
  
  #
  # unlist objects in scenario where treated male are treated at 95-95-95
  #
  
  # number of male treated regardless of age
  budget.counterfactual.959595 <- counterfactuals_p_959595$budget 
  budget.counterfactual.959595[, label := label.959595]
  
  # relative incident cases counterfactual compared to factual by age group 
  relative_incidence_counterfactual_group.959595  <- counterfactuals_p_959595$relative_incidence_counterfactual_group2
  relative_incidence_counterfactual_group.959595[, label := label.959595]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all.959595 <- counterfactuals_p_959595$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all.959595[, label := label.959595]
  
  # sex ratio incident cases counterfactual  by age group
  sex_ratio_incidence_counterfactual_group.959595 <- counterfactuals_p_959595$sex_ratio_incidence_counterfactual_group2 
  sex_ratio_incidence_counterfactual_group.959595[, label := label.959595]
  
  # sex ratio incident cases counterfactual  by age group regardless of age
  sex_ratio_incidence_counterfactual_all.959595 <- counterfactuals_p_959595$sex_ratio_incidence_counterfactual_all 
  sex_ratio_incidence_counterfactual_all.959595[, label := label.959595]
  
  
  #
  # unlist objects in scenario where treated male are treated 90 90 90
  #
  
  # number of male treated regardless of age
  budget.counterfactual.909090 <- counterfactuals_p_909090$budget 
  budget.counterfactual.909090[, label := label.909090]
  
  # relative incident cases counterfactual compared to factual by age group 
  relative_incidence_counterfactual_group.909090  <- counterfactuals_p_909090$relative_incidence_counterfactual_group2
  relative_incidence_counterfactual_group.909090[, label := label.909090]
  
  # relative incident cases counterfactual compared to factual regardless of age
  relative_incidence_counterfactual_all.909090 <- counterfactuals_p_909090$relative_incidence_counterfactual_all 
  relative_incidence_counterfactual_all.909090[, label := label.909090]
  
  # sex ratio incident cases counterfactual  by age group
  sex_ratio_incidence_counterfactual_group.909090 <- counterfactuals_p_909090$sex_ratio_incidence_counterfactual_group2 
  sex_ratio_incidence_counterfactual_group.909090[, label := label.909090]
  
  # sex ratio incident cases counterfactual  by age group regardless of age
  sex_ratio_incidence_counterfactual_all.909090 <- counterfactuals_p_909090$sex_ratio_incidence_counterfactual_all 
  sex_ratio_incidence_counterfactual_all.909090[, label := label.909090]
  
  
  #
  # combine both scenarios
  #
  
  budget.counterfactual <- do.call('rbind', list(budget.counterfactual, budget.counterfactual.f05, 
                                                 budget.counterfactual.959595, budget.counterfactual.909090))
  budget.counterfactual[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  relative_incidence_counterfactual_group <- do.call('rbind', list(relative_incidence_counterfactual_group, relative_incidence_counterfactual_group.f05, 
                                                                   relative_incidence_counterfactual_group.959595, relative_incidence_counterfactual_group.909090))
  relative_incidence_counterfactual_group[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  relative_incidence_counterfactual_all <- do.call('rbind', list(relative_incidence_counterfactual_all, relative_incidence_counterfactual_all.f05, 
                                                                 relative_incidence_counterfactual_all.959595, relative_incidence_counterfactual_all.909090))
  relative_incidence_counterfactual_all[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  sex_ratio_incidence_counterfactual_group <- do.call('rbind', list(sex_ratio_incidence_counterfactual_group, sex_ratio_incidence_counterfactual_group.f05, 
                                                                    sex_ratio_incidence_counterfactual_group.959595, sex_ratio_incidence_counterfactual_group.909090))
  sex_ratio_incidence_counterfactual_group[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  sex_ratio_incidence_counterfactual_all <- do.call('rbind', list(sex_ratio_incidence_counterfactual_all, sex_ratio_incidence_counterfactual_all.f05, 
                                                                  sex_ratio_incidence_counterfactual_all.959595, sex_ratio_incidence_counterfactual_all.909090))
  sex_ratio_incidence_counterfactual_all[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
  
  #
  # restrict to one round and to male to female direction
  #
  
  Round <- 'R018'
  budget.counterfactual <- budget.counterfactual[ROUND == Round & SEX == 'M']
  relative_incidence_counterfactual_group <- relative_incidence_counterfactual_group[ROUND == Round & IS_MF == T]
  relative_incidence_counterfactual_all <- relative_incidence_counterfactual_all[ROUND == Round & IS_MF == T]
  sex_ratio_incidence_counterfactual_group <- sex_ratio_incidence_counterfactual_group[ROUND == Round]
  sex_ratio_incidence_counterfactual_all <- sex_ratio_incidence_counterfactual_all[ROUND == Round]
  
  
  #
  # budget
  #
  
  bc <- copy(budget.counterfactual[order(COMM, ROUND, SEX, label)])
  bc[, `:=` (TREATED = round(TREATED), TREATED_CL = round(TREATED_CL), TREATED_CU = round(TREATED_CU))]
  
  #
  # relative incidence rate
  #
  
  ric <- copy(relative_incidence_counterfactual_group)
  ric <-ric[order(COMM, ROUND, LABEL_RECIPIENT, AGE_GROUP_INFECTION.RECIPIENT, label), .(COMM, ROUND, LABEL_RECIPIENT, AGE_GROUP_INFECTION.RECIPIENT, label, M, CL, CU)]
  ric.all <- copy(relative_incidence_counterfactual_all)
  ric.all <-ric.all[order(COMM, ROUND, LABEL_RECIPIENT, label), .(COMM, ROUND, LABEL_RECIPIENT, label, M, CL, CU)]
  ric.all[, AGE_GROUP_INFECTION.RECIPIENT := 'All']
  ric <- rbind(ric, ric.all)
  ric[, `:=` (M = format(round(M*100, 1), nsmall = 1), CL = format(round(CL*100, 1), nsmall = 1), 
              CU = format(round(CU*100, 1), nsmall = 1))]
  
  #
  # sex incidence rate ratio
  #
  
  ser <- copy(sex_ratio_incidence_counterfactual_group)
  ser <-ser[order(COMM, ROUND, AGE_GROUP_INFECTION.RECIPIENT, label), .(COMM, ROUND, AGE_GROUP_INFECTION.RECIPIENT, label, M, CL, CU)]
  ser.all <- copy(sex_ratio_incidence_counterfactual_all)
  ser.all <-ser.all[order(COMM, ROUND, label), .(COMM, ROUND, label, M, CL, CU)]
  ser.all[, AGE_GROUP_INFECTION.RECIPIENT := 'All']
  ser <- rbind(ser, ser.all)
  ser[, `:=` (M = format(round(M, 1), nsmall = 1), CL = format(round(CL, 1), nsmall = 1), 
              CU = format(round(CU, 1), nsmall = 1))]
  
  
  #
  # save
  #
  table <- list('budget' = bc, 
                'relative_incidence' = ric, 
                'indidence_rate_ratio' = ser)
  file = paste0(outdir, '-output-statistics_budget_counterfactual_', gsub(' ' , '', lab), '.rds')
  saveRDS(table, file)
  return(table)
}


