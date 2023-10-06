

save_statistics_transmission_events <- function(pairs, pairs.all, outdir){
  
  stat <- list()
  
  # number of pairs
  stat[[1]] <- nrow(pairs)
  
  # range of the infection between source and recipient
  stat[[2]] <- format(pairs[, range(DATE_INFECTION.RECIPIENT)], '%B %d, %Y')
  
  # number of pairs MM and FF
  tmp <- pairs.all[COMM.RECIPIENT == "inland" & COMM.SOURCE == 'inland']
  tmp <- tmp[COMM.SOURCE != 'neuro' & COMM.RECIPIENT != "neuro"]
  tmp[, DIRECTION := paste0(SEX.SOURCE, SEX.RECIPIENT)]
  stat[[3]] <- tmp[, list(N = .N), by = 'DIRECTION']
  
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

save_statistics_incidence_rate_trends <- function(icrr, icr, icrrs, icrrt, medage, fmr){
  
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
  
  # relative incidence ratio total
  inc_rel_ratio <- icrrt[,.(M= round((M), 2), CL = round((CL), 2), CU = round((CU), 2)), by=c('COMM', 'ROUND')]
  inc_rel_ratio <- inc_rel_ratio[order(COMM, ROUND)]
  
  # median age at infection
  n_digit <- 1
  median_age <- copy(medage[order(COMM,SEX, ROUND)])
  median_age[, M_roundn := format(round(M, n_digit), nsmall = n_digit)]
  median_age[, CL_roundn  := format(round(CL, n_digit), nsmall = n_digit)]
  median_age[, CU_roundn  := format(round(CU, n_digit), nsmall = n_digit)]
  median_age[, M_round0 := round(M)]
  median_age[, CL_round0  := round(CL)]
  median_age[, CU_round0  := round(CU)]

  # incidence ratio total
  fm_incidence_rate_ratio <- fmr[,.(M= round((M), 2), CL = round((CL), 2), CU = round((CU), 2)), by=c('COMM', 'ROUND')]
  fm_incidence_rate_ratio <- fm_incidence_rate_ratio[order(COMM, ROUND)]
  
  #save
  stats <- list()
  
  stats$inc_rel_age <- df_incidence_rel_age[(COMM == 'inland' & ROUND %in% c('R018')) | (COMM == 'fishing' & ROUND %in% c('R018'))]
  stats$inc_rel <- df_incidence_rel[(COMM == 'inland' & ROUND %in% c('R018')) | (COMM == 'fishing' & ROUND %in% c('R018'))]
  stats$inc_age <- df_incidence_age[(COMM == 'inland' & ROUND %in% c('R010', 'R018')) | (COMM == 'fishing' & ROUND %in% c('R015', 'R018'))]
  stats$inc <- df_incidence[(COMM == 'inland' & ROUND %in% c('R010','R018')) | (COMM == 'fishing' & ROUND %in% c('R015','R018'))]
  stats$median_age <- median_age
  stats$inc_rel_ratio_age <- inc_rel_ratio_age
  stats$fm_incidence_rate_ratio <- fm_incidence_rate_ratio
  
  saveRDS(stats, paste0(outdir.table, '-data-incidence_rate_trends.rds'))
  
}

save_number_pairs_round16_18 <- function(pars, df_round, outdir){
  
  # find number of pairs observed
  dp <- copy(pairs)
  
  # we add pairs to the corresponding rounds (move their date of infection so that it falls within the observational period)
  tmp <- df_round[, .(ROUND, MIN_SAMPLE_DATE, MAX_SAMPLE_DATE, COMM)]
  tmp[, MIN_SAMPLE_DATE_NEXT_ROUND := c(MIN_SAMPLE_DATE[2:nrow(tmp)], MAX_SAMPLE_DATE[nrow(tmp)])]
  stopifnot(tmp[, all(MAX_SAMPLE_DATE <= MIN_SAMPLE_DATE_NEXT_ROUND)])
  
  dp <- merge(dp, tmp,  by.x = 'COMM.RECIPIENT', by.y = 'COMM', allow.cartesian=TRUE)
  dp <- dp[DATE_INFECTION.RECIPIENT >= MIN_SAMPLE_DATE & DATE_INFECTION.RECIPIENT <= MIN_SAMPLE_DATE_NEXT_ROUND]
  stopifnot(nrow(dp) == nrow(pairs))
  set(dp, NULL, 'MIN_SAMPLE_DATE', NULL)
  set(dp, NULL, 'MAX_SAMPLE_DATE', NULL)
  set(dp, NULL, 'MIN_SAMPLE_DATE_NEXT_ROUND', NULL)
  
  # number pairs inround 16 to 18
  n <- nrow(dp[ROUND %in% c('R016', 'R017', 'R018')])
  
  # save
  saveRDS(n, paste0(outdir, '-data-number_pairs_rounds_16-18.rds'))
}

