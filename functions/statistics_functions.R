

save_statistics_transmission_events <- function(pairs, outdir){
  
  stat <- list()
  
  # number of pairs
  stat[[1]] <- nrow(pairs)
  
  # range of the infection between source and recipient
  stat[[2]] <- format(pairs[, range(DATE_INFECTION.RECIPIENT)], '%B %d, %Y')
  
  print(stat)
  
  file = paste0(outdir, '-data-pairs_stat.rds')
  saveRDS(stat, file)
}

save_statistics_incidence_rate_ratio_trends <- function(ic, outdir.table){
  
  # select relative incident rate ratio and unsuppressed rate ratio for inland across all the rounds 
  tmp <- ic[COMM== 'inland' & INDEX_ROUND != min(INDEX_ROUND)]
  
  # fit a one way anova
  one.way <- aov(INCIDENT_RATE_RATIO_REF ~ UNSUPPRESSION_RATE_RATIO_RATIO_M, data = tmp)
  summary(one.way)
  p_value <- round(summary(one.way)[[1]][["Pr(>F)"]][1], 4) # pvalue of the anova
  f_value <- round(summary(one.way)[[1]][["F value"]][1], 4) # f value of the anova
  df <- summary(one.way)[[1]][["Df"]][1] # df of the anova
  
  # fit a linear regression
  fit.lm <- lm(INCIDENT_RATE_RATIO_REF ~ UNSUPPRESSION_RATE_RATIO_RATIO_M, data = tmp)
  summary(fit.lm)
  
  # save
  stats = list()
  stats$p_value <- p_value
  stats$f_value <- f_value
  stats$df <- nrow(tmp) -1
  
  saveRDS(stats, paste0(outdir.table, '-data-incidence_rate_ratio_trends.rds'))
}

save_statistics_incidence_rate_trends <- function(icrr, icr, median_age, icrrs){
  
  ps <- c(0.5, 0.2, 0.8)
  p_labs <- c('M','CL','CU')
  
  # relative incidence by gender and age group
  df_incidence_rel_age = icrr[, list(q= quantile(INCIDENCE_REL, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'age_group', 'SEX')]	
  df_incidence_rel_age = dcast(df_incidence_rel_age, ... ~ q_label, value.var = "q")
  df_incidence_age = icrr[, list(q= quantile(INCIDENT_RATE_SUSCEPTIBLE, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'SEX', 'age_group')]	
  df_incidence_age = dcast(df_incidence_age, ... ~ q_label, value.var = "q")
  df_incidence_rel_age <- df_incidence_rel_age[, .(M= round((1-M)*100, 1), CL = round((1-CU)*100, 1), CU = round((1-CL) * 100, 1)), by=c('COMM', 'ROUND', 'age_group', 'SEX')]
  df_incidence_age <- df_incidence_age[, .(M= round((M)*100, 2), CL = round((CL)*100, 2), CU = round((CU) * 100, 2)), by=c('COMM', 'ROUND', 'age_group','SEX')]
  df_incidence_rel_age <- df_incidence_rel_age[order(COMM, ROUND, SEX, age_group)]
  df_incidence_age <- df_incidence_age[order(COMM, ROUND, SEX, age_group)]
  
  # relative incidence by gender regardless of age
  tmp <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
  tmp[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
  tmp <- tmp[, list(INCIDENT_RATE_SUSCEPTIBLE = sum(INCIDENT_CASES) / sum(SUSCEPTIBLE)), by = c('SEX', 'ROUND', 'COMM', 'REF.ROUND', 'iterations')] 
  tmp[, INCIDENCE_REL := INCIDENT_RATE_SUSCEPTIBLE / INCIDENT_RATE_SUSCEPTIBLE[ROUND == REF.ROUND], by = c('COMM', 'SEX', 'iterations')]
  df_incidence_rel = tmp[, list(q= quantile(INCIDENCE_REL, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'SEX')]	
  df_incidence_rel = dcast(df_incidence_rel, ... ~ q_label, value.var = "q")
  df_incidence = tmp[, list(q= quantile(INCIDENT_RATE_SUSCEPTIBLE, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'SEX')]	
  df_incidence = dcast(df_incidence, ... ~ q_label, value.var = "q")
  df_incidence_rel <- df_incidence_rel[, .(M= round((1-M)*100, 1), CL = round((1-CU)*100, 1), CU = round((1-CL) * 100, 1)), by=c('COMM', 'ROUND', 'SEX')]
  df_incidence <- df_incidence[, .(M= round((M)*100, 2), CL = round((CL)*100, 2), CU = round((CU) * 100, 2)), by=c('COMM', 'ROUND', 'SEX')]
  df_incidence_rel <- df_incidence_rel[order(COMM, ROUND, SEX)]
  df_incidence <- df_incidence[order(COMM, ROUND, SEX)]
  
  # relative incidence ratio by age
  inc_rel_ratio_age <- icrrs[,.(M= round((M), 2), CL = round((CL), 2), CU = round((CU), 2)), by=c('COMM', 'ROUND', 'age_group')]
  inc_rel_ratio_age <- inc_rel_ratio_age[order(COMM, ROUND, age_group)]
  
  # median age at infection
  median_age[, MEDIAN_AGEYRS := round(MEDIAN_AGEYRS, 2)]
  median_age <- median_age[order(COMM, ROUND, SEX_LABEL), .(COMM, ROUND, SEX_LABEL, MEDIAN_AGEYRS)]
  
  #save
  stats <- list()
  
  stats$inc_rel_age <- df_incidence_rel_age[(COMM == 'inland' & ROUND %in% c('R018')) | (COMM == 'fishing' & ROUND %in% c('R018'))]
  stats$inc_rel <- df_incidence_rel[(COMM == 'inland' & ROUND %in% c('R018')) | (COMM == 'fishing' & ROUND %in% c('R018'))]
  stats$inc_age <- df_incidence_age[(COMM == 'inland' & ROUND %in% c('R010', 'R018')) | (COMM == 'fishing' & ROUND %in% c('R015', 'R018'))]
  stats$inc <- df_incidence[(COMM == 'inland' & ROUND %in% c('R010','R018')) | (COMM == 'fishing' & ROUND %in% c('R015','R018'))]
  stats$median_age <- median_age
  stats$inc_rel_ratio_age <- inc_rel_ratio_age
  
  saveRDS(stats, paste0(outdir.table, '-data-incidence_rate_trends.rds'))
  
}
