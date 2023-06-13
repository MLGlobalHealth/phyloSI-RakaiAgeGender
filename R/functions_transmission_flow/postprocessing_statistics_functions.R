save_statistics_expected_contribution <- function(expected_contribution_sex_source, expected_contribution_age_group_source_total, outdir){
  
  n_digits <- 1
  
  # contribution of men in first and last round 
  ecm <- copy(expected_contribution_sex_source)
  ecm[, MAX_INDEX_ROUND := max(INDEX_ROUND), by = c('COMM')]
  ecm <- ecm[INDEX_ROUND == 1 | INDEX_ROUND == MAX_INDEX_ROUND, .(M, CL, CU), by = c('COMM', 'LABEL_SOURCE', 'ROUND', 'MAX_INDEX_ROUND', 'INDEX_ROUND')]
  ecm <- ecm[, .(M = round(M*100, n_digits), CL = round(CL * 100, n_digits), CU = round(CU * 100, n_digits)), by = c('COMM', 'LABEL_SOURCE', 'ROUND')]
  ecm <- ecm[order(COMM, LABEL_SOURCE, ROUND)]
  
  # contribution of men by age in last round 
  ecma <- copy(expected_contribution_age_group_source_total)
  ecma[, MAX_INDEX_ROUND := max(INDEX_ROUND), by = c('COMM')]
  ecma <- ecma[ INDEX_ROUND == MAX_INDEX_ROUND, .(M, CL, CU), by = c('COMM', 'ROUND', 'LABEL_SOURCE', 'AGE_GROUP_TRANSMISSION.SOURCE')]
  ecma <- ecma[, .(M = round(M*100, n_digits), CL = round(CL * 100, n_digits), CU = round(CU * 100, n_digits)), by = c('COMM', 'ROUND', 'LABEL_SOURCE', 'AGE_GROUP_TRANSMISSION.SOURCE')]
  ecma <- ecma[order(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_TRANSMISSION.SOURCE)]
  
  # save
  stats <- list()
  stats$contribution_male <- ecm
  stats$contribution_male_age <- ecma
  
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
  
  saveRDS(table, file = paste0(outdir, '-output-transmission_flows_table.rds'))
  
  return(table)
}

save_median_age_diff <- function(median_age_diff_group1, median_age_diff_group2, outdir){
  
  table <- list()
  n_digits <- 1
  
  #
  # median age difference male to female age group 1
  #
  
  mad <- copy(median_age_diff_group1)
  mad <- mad[ROUND == "R018" & quantile == 'C50']
  mad <- mad[, .(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT, quantile, M, CL, CU)]
  mad <- mad[order(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT)]
  mad[, M := format(round(M, n_digits), nsmall=n_digits)]
  mad[, CL := format(round(CL, n_digits), nsmall=n_digits)]
  mad[, CU := format(round(CU, n_digits), nsmall=n_digits)]
  mad[, M := gsub(' ', '', M)]
  mad[, CL := gsub(' ', '', CL)]
  mad[, CU := gsub(' ', '', CU)]
  table[['median_age_difference_group1']] <- mad
  
  #
  # median age difference male to female age group 2
  #
  
  mad <- copy(median_age_diff_group2)
  mad <- mad[ROUND == "R018" & quantile == 'C50']
  mad <- mad[, .(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT, quantile, M, CL, CU)]
  mad <- mad[order(COMM, ROUND, LABEL_DIRECTION, AGE_GROUP_INFECTION.RECIPIENT)]
  mad[, M := format(round(M, n_digits), nsmall=n_digits)]
  mad[, CL := format(round(CL, n_digits), nsmall=n_digits)]
  mad[, CU := format(round(CU, n_digits), nsmall=n_digits)]
  mad[, M := gsub(' ', '', M)]
  mad[, CL := gsub(' ', '', CL)]
  mad[, CU := gsub(' ', '', CU)]
  table[['median_age_difference_group2']] <- mad
  
  #
  #Save
  #
  
  saveRDS(table, paste0(outdir, '-output-median_age_difference.rds'))
  

}

save_median_age_source <- function(median_age_source_group, median_age_source_group2, 
                                   median_age_source, outdir){
  table <- list()
  n_digits <- 1
  
  #
  # Median age source
  #
  
  mas <- copy(median_age_source_group)
  mas <- mas[ROUND %in% c('R010', 'R018') & quantile == 'C50']
  mas <- mas[, .(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT, M, CL, CU)]
  mas <- mas[order(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT)]
  mas[, M := format(round(M, n_digits), nsmall = n_digits)]
  mas[, CL := format(round(CL, n_digits), nsmall = n_digits)]
  mas[, CU := format(round(CU, n_digits), nsmall = n_digits)]
  
  table[['median_age_source_agg']] <- mas
  
  #
  # Median age source other aggregate
  #
  
  mas <- copy(median_age_source_group2)
  mas <- mas[ROUND %in% c('R010', 'R018') & quantile == 'C50']
  mas <- mas[, .(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT, M, CL, CU)]
  mas <- mas[order(COMM, ROUND, LABEL_SOURCE, AGE_GROUP_INFECTION.RECIPIENT)]
  mas[, M := format(round(M, n_digits), nsmall = n_digits)]
  mas[, CL := format(round(CL, n_digits), nsmall = n_digits)]
  mas[, CU := format(round(CU, n_digits), nsmall = n_digits)]
  
  table[['median_age_source_agg2']] <- mas
  
  #
  # Median age source total
  #
  
  mas <- copy(median_age_source)
  mas <- mas[ROUND %in% c('R010', 'R018') & quantile == 'C50']
  mas <- mas[, .(COMM, ROUND, LABEL_SOURCE, M, CL, CU)]
  mas <- mas[order(COMM, ROUND, LABEL_SOURCE)]
  mas[, M := format(round(M, n_digits), nsmall = n_digits)]
  mas[, CL := format(round(CL, n_digits), nsmall = n_digits)]
  mas[, CU := format(round(CU, n_digits), nsmall = n_digits)]
  
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
  
  # number of male treated by age
  budget_age.counterfactual <- counterfactuals_p_f$budget_age
  budget_age.counterfactual[, label := label.f]
  
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
  
  # number of male treated by age
  budget_age.counterfactual.f05 <- counterfactuals_p_f05$budget_age
  budget_age.counterfactual.f05[, label := label.f05]
  
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
  
  # number of male treated by age
  budget_age.counterfactual.959595 <- counterfactuals_p_959595$budget_age
  budget_age.counterfactual.959595[, label := label.959595]
  
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
  
  # number of male treated by age
  budget_age.counterfactual.909090 <- counterfactuals_p_909090$budget_age
  budget_age.counterfactual.909090[, label := label.909090]
  
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
  
  budget_age.counterfactual <- do.call('rbind', list(budget_age.counterfactual, budget_age.counterfactual.f05, 
                                                     budget_age.counterfactual.959595, budget_age.counterfactual.909090))
  budget_age.counterfactual[, label := factor(label, levels = c(label.f05, label.909090, label.f, label.959595))]
  
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
  budget_age.counterfactual <- budget_age.counterfactual[ROUND == Round & SEX == 'M']
  relative_incidence_counterfactual_group <- relative_incidence_counterfactual_group[ROUND == Round & IS_MF == T]
  relative_incidence_counterfactual_all <- relative_incidence_counterfactual_all[ROUND == Round & IS_MF == T]
  sex_ratio_incidence_counterfactual_group <- sex_ratio_incidence_counterfactual_group[ROUND == Round]
  sex_ratio_incidence_counterfactual_all <- sex_ratio_incidence_counterfactual_all[ROUND == Round]
  
  
  #
  # budget
  #
  n_digits <- 1
  bc <- copy(budget.counterfactual[order(COMM, ROUND, SEX, label)])
  bc[, `:=` (TREATED = format(round(TREATED, n_digits), nsmall = n_digits), 
             TREATED_CL = format(round(TREATED_CL, n_digits), nsmall = n_digits), 
             TREATED_CU = format(round(TREATED_CU), nsmall=n_digits))]
  bc[, TREATED := gsub(' ', '', TREATED)]
  bc[, TREATED_CL := gsub(' ', '', TREATED_CL)]
  bc[, TREATED_CU := gsub(' ', '', TREATED_CU)]
  
  
  #
  # relative incidence 
  #
  
  n_digits <- 1
  
  ric.all <- copy(relative_incidence_counterfactual_all)
  ric.all <-ric.all[order(COMM, ROUND, LABEL_RECIPIENT, label), .(COMM, ROUND, LABEL_RECIPIENT, label, M, CL, CU)]
  ric.all[, AGE_GROUP_INFECTION.RECIPIENT := 'All']
  ric.all[, `:=` (M = format(round(M*100, n_digits), nsmall = n_digits), CL = format(round(CL*100, n_digits), nsmall = n_digits), 
              CU = format(round(CU*100, n_digits), nsmall = n_digits))]
  
 
  #
  # sex incidence rate ratio
  #
  
  n_digits <- 2
  
  ser.all <- copy(sex_ratio_incidence_counterfactual_all)
  ser.all <-ser.all[order(COMM, ROUND, label), .(COMM, ROUND, label, M, CL, CU)]
  ser.all[, `:=` (M = format(round(M, n_digits), nsmall = n_digits), CL = format(round(CL, n_digits), nsmall = n_digits), 
              CU = format(round(CU, n_digits), nsmall = n_digits))]
  
  #
  # budget by age
  #
  
  n_digits <- 1
  
  bca <- copy(budget_age.counterfactual[order(label, COMM, ROUND, SEX, AGE_GROUP)])
  bca[PROP_TREATED_CL < 0, PROP_TREATED_CL := 0]
  bca[, `:=` (TREATED = format(round(TREATED, n_digits), nsmall = n_digits), 
             TREATED_CL = format(round(TREATED_CL, n_digits), nsmall = n_digits), 
             TREATED_CU = format(round(TREATED_CU), nsmall=n_digits), 
             PROP_TREATED = format(round(PROP_TREATED * 100, n_digits), n_digits),
             PROP_TREATED_CL = format(round(PROP_TREATED_CL * 100, n_digits), n_digits),
             PROP_TREATED_CU = format(round(PROP_TREATED_CU * 100, n_digits), n_digits))]
  
  #
  # relative incidence by age
  #
  
  n_digits <- 1
  
  ric <- copy(relative_incidence_counterfactual_group)
  ric <-ric[order(label, COMM, ROUND, LABEL_RECIPIENT, AGE_GROUP_INFECTION.RECIPIENT), .(COMM, ROUND, LABEL_RECIPIENT, AGE_GROUP_INFECTION.RECIPIENT, label, M, CL, CU)]
  ric[, `:=` (M = format(round(M*100, n_digits), nsmall = n_digits), CL = format(round(CL*100, n_digits), nsmall = n_digits), 
              CU = format(round(CU*100, n_digits), nsmall = n_digits))]
  
  
  #
  # sex incidence rate ratio by age
  #
  
  n_digits <- 2
  
  ser <- copy(sex_ratio_incidence_counterfactual_group)
  ser <-ser[order(label, COMM, ROUND, AGE_GROUP_INFECTION.RECIPIENT), .(COMM, ROUND, AGE_GROUP_INFECTION.RECIPIENT, label, M, CL, CU)]
  ser[, `:=` (M = format(round(M, n_digits), nsmall = n_digits), CL = format(round(CL, n_digits), nsmall = n_digits), 
              CU = format(round(CU, n_digits), nsmall = n_digits))]
  
  
  #
  # save
  #
  
  table <- list('budget' = bc, 
                'relative_incidence' = ric.all, 
                'indidence_rate_ratio' = ser.all, 
                'budget_age' = bca, 
                'relative_incidence_age' = ric, 
                'indidence_rate_ratio_age' = ser)
  file = paste0(outdir, '-output-statistics_budget_counterfactual_', gsub(' ' , '', lab), '.rds')
  saveRDS(table, file)
  return(table)
}

save_counterfactual_results_for_UNAIDS <- function(  counterfactuals_a_f,
                                                     counterfactuals_a_f05,
                                                     counterfactuals_a_959595,
                                                     counterfactuals_a_909090,
                                                     incidence_factual,
                                                     lab, 
                                                     outdir.table)
{
  Date <- Sys.Date()
  
  # labels
  lab_table <- data.table(label = c('Diagnosed rate in men as in women\nART coverage in men as in women\nSuppression rate in men as in women', 
                                    'Half-way to\nDiagnosed rate in men as in women\nART coverage in men as in women\nSuppression rate in men as in women',
                                    '95% diagnosed\n95% receiving ART\n95% suppression rate'), 
                          intervention = c('Closing the suppression gap in men relative to women', 
                                           'Closing half the suppression gap in men relative to women', 
                                           '95-95-95 in men'))
  
  var_table <- data.table(variable = c('INFECTED_NON_SUPPRESSED', 'INFECTED_ALREADY_SUPPRESSED', 'TREATED'), 
                          category = c('Remaining virally unsuppressed', 'Already virally suppressed', 'Additional virally suppressed'))
  
  # combine tables
  cct <- combine_counterfactual_tables(counterfactuals_a_f,
                                       counterfactuals_a_f05,
                                       counterfactuals_a_959595,
                                       counterfactuals_a_909090,
                                       incidence_factual,
                                       lab = lab, 
                                       Round = 'R018',
                                       include_909090 = F)
  
  # subfigure a
  ecf = copy(cct$eligible_count_round.counterfactual)
  ecf <- merge(ecf, var_table, by = 'variable')
  ecf <- merge(ecf, lab_table, by = 'label')
  set(ecf, NULL, 'variable', NULL)
  set(ecf, NULL, 'label', NULL)
  set(ecf, NULL, 'ROUND', NULL)
  set(ecf, NULL, 'SEX', NULL)
  set(ecf, NULL, 'COMM', NULL)
  ecf[, label := 'Number_Men_With_HIV']
  write.csv(ecf, row.names = F, 
            file.path(dirname(outdir.table),
                      paste0('UNAIDSGlobalReport-Number_Men_With_HIV-Figure4a-Monodetal_', Date,'.csv')), 
  )
  
  # subfigure b
  bc = copy(cct$budget.counterfactual)
  bc <- merge(bc, lab_table, by = 'label')
  set(bc, NULL, 'label', NULL)
  set(bc, NULL, 'ROUND', NULL)
  set(bc, NULL, 'SEX', NULL)
  set(bc, NULL, 'COMM', NULL)
  setnames(bc, c('TREATED', 'TREATED_CL' , 'TREATED_CU'), c('M', 'CL', 'CU'))
  bc[, label := 'Additional_Number_Men_With_Suppressed_Virus']
  write.csv(bc, row.names = F, 
            file.path(dirname(outdir.table),
                      paste0('UNAIDSGlobalReport-Additional_Number_Men_With_Suppressed_Virus-Figure4b-Monodetal_', Date,'.csv')), 
  )
  
  # subfigure c 
  ric.all = copy(cct$relative_incidence_counterfactual_all)
  ric.all <- merge(ric.all, lab_table, by = 'label')
  ric.all <- ric.all[, .(M, CL, CU, intervention)]
  ric.all[, label := '%_Reduction_In_Incidence_In_Women']
  write.csv(ric.all, row.names = F, 
            file.path(dirname(outdir.table),
                      paste0('UNAIDSGlobalReport-Percent_Reduction_In_Incidence_In_Women-Figure4c-Monodetal_', Date,'.csv')), 
  )
  
  # subfigure d - no intervention
  icf = copy(cct$incidence_factual)
  icf <- icf[, .(LABEL_GENDER_RECIPIENT, AGE_INFECTION.RECIPIENT, M, CL, CU)]
  icf[, M := M * 100]
  icf[, CL := CL * 100]
  icf[, CU := CU * 100]
  icf[, intervention := 'No intervention']
  icf[, label := 'Incidence_Rate_per_100_person_year']
  write.csv(icf, row.names = F, 
            file.path(dirname(outdir.table),
                      paste0('UNAIDSGlobalReport-Incidence_Rate_per_100_person_year-NoIntervention-Figure4d-Monodetal_', Date,'.csv')), 
  )
  
  # subfigure d - intervention
  ic = copy(cct$incidence_counterfactual)
  ic <- merge(ic, lab_table, by = 'label')
  ic <- ic[, .(intervention, AGE_INFECTION.RECIPIENT, M, CL, CU)]
  ic[, M := M * 100]
  ic[, CL := CL * 100]
  ic[, CU := CU * 100]
  ic[, label := 'Incidence_Rate_per_100_person_year']
  write.csv(ic, row.names = F, 
            file.path(dirname(outdir.table),
                      paste0('UNAIDSGlobalReport-Incidence_Rate_per_100_person_year-WithIntervention-Figure4d-Monodetal_', Date,'.csv')), 
  )
}
