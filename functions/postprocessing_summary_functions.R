make_convergence_diagnostics_stats = function(fit, re, outdir)
{
  
  stopifnot(!is.null(fit))
  
  summary = rstan::summary(fit)$summary
  eff_sample_size_cum = summary[,9][!is.na(summary[,9])]
  Rhat_cum = summary[,10][!is.na(summary[,10])]
  cat("the minimum and maximum effective sample size are ", range(eff_sample_size_cum), "\n")
  cat("the minimum and maximum Rhat are ", range(Rhat_cum), "\n")
  if(min(eff_sample_size_cum) < 500) cat('\nEffective sample size smaller than 500 \n')
  
  tryCatch({
    sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
    sampler_diagnostics <- data.table()
    for (i in colnames(sampler_params[[1]])) {
      tmp <- data.table(t(sapply(sampler_params, function(x) quantile(x[, i],probs = c(0.025,0.5,0.975)))))
      tmp[, diagnostics:=i ]
      tmp[, chain:= seq_len(length(sampler_params))]
      sampler_diagnostics <- rbind(sampler_diagnostics, tmp)
    }
    print(sampler_diagnostics)
  }, error = function(e) e)
  
  
  
  check_all_diagnostics(fit, outdir)
  
  # compute WAIC and LOO
  tryCatch({
    
    if('log_lik' %in% names(re)){
      log_lik <- loo::extract_log_lik(fit)
      log_lik = log_lik[!is.na(log_lik[,1]),]
      .WAIC = loo::waic(log_lik)
      .LOO = loo::loo(log_lik)
      print(.WAIC); print(.LOO)
      WAIC = .WAIC$pointwise
      LOO = .LOO$pointwise
    }} , error = function(e) e)
  
  # time of execution
  time = sum(rstan::get_elapsed_time(fit))
  
  # save
  saveRDS(eff_sample_size_cum, file = paste0(outdir, "-eff_sample_size_cum.rds"))
  saveRDS(Rhat_cum, file = paste0(outdir,  "-Rhat_cum.rds"))
  saveRDS(.WAIC, file = paste0(outdir, "-WAIC.rds"))
  saveRDS(.LOO, file = paste0(outdir, "-LOO.rds"))
  saveRDS(sampler_diagnostics, file = paste0(outdir, "-sampler_diagnostics.rds"))
  saveRDS(time, file = paste0(outdir, "-time_elapsed.rds"))
  
 
  # return(summary)
}

check_all_diagnostics <- function(fit, outdir) {
  check_n_eff(fit)
  check_rhat(fit)
  n_div <- check_div(fit)
  n_treedepth <- check_treedepth(fit,15)
  check_energy(fit)
  
  cat('\nn_div',n_div, '\n' )
  cat('\nn_treedepth',n_treedepth, '\n' )
  saveRDS(list(n_div,n_treedepth), file=paste0(outdir,'-diagnostics.rds'))
}

check_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                n, N, max_depth, 100 * n / N))
  if (n > 0)
    print('  Run again with max_depth set to a larger value to avoid saturation')
  return(n)
}

check_energy <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)*2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('E-BFMI indicated no pathological behavior')
  else
    print('  E-BFMI below 0.2 indicates you may need to reparameterize your model')
}

check_n_eff <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  # remove na neff
  fit_summary <- fit_summary[!is.na(fit_summary[,5]),]
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(rstan::extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      print(sprintf('n_eff / iter for parameter %s is %s!',
                    rownames(fit_summary)[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('n_eff / iter looks reasonable for all parameters')
  else
    print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
}

check_div <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                n, N, 100 * n / N))
  if (n > 0)
    print('  Try running with larger adapt_delta to remove the divergences')
  # return iterations ended with a divergence
  return(n)
}

check_rhat <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  # remove na neff
  fit_summary <- fit_summary[!is.na(fit_summary[,6]),]
  N <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      print(sprintf('Rhat for parameter %s is %s!',
                    rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('Rhat looks reasonable for all parameters')
  else
    print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
}

find_summary_output <- function(samples, output, vars, transform = NULL, standardised.vars = NULL, names = NULL, operation = NULL){
  
  # summarise outputs by period
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[[output]]) )
  if(!is.null(names)){
    setnames(tmp1, 2:(length(names) + 1), names)
    }else if(tmp1[, max(Var2)] == df_age[, max(INDEX_AGE)]){
    setnames(tmp1, 2:5, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'))
    }else{
    setnames(tmp1, 2:5, c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'))
  }

  if('INDEX_AGE' %in% names(tmp1)){
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
    tmp1 <- merge(tmp1, df_age_aggregated, by = c('AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'))
  }

  if(!is.null(transform)){
    tmp1[, value := sapply(value, transform)]
  }
  
  #  sum force of infection
  if(is.null(operation)){
    tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', vars)]
  } else{
    tmp1 <- tmp1[, list(value = sapply(value, operation)), by = c('iterations', vars)]
  }

  # standardised
  if(!is.null(standardised.vars)){
    tmp1[, total_value := sum(value), by = c('iterations', standardised.vars)]
    tmp1[, value := value / total_value]
  }

  #summarise
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=vars]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")

  
  if('INDEX_DIRECTION' %in% vars)
    tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  if('INDEX_COMMUNITY' %in% vars)
    tmp1 <- merge(tmp1, df_community, by = 'INDEX_COMMUNITY')
  if('INDEX_TIME' %in% vars)
    tmp1 <- merge(tmp1, df_period, by = c('INDEX_TIME', 'COMM'))
  if('INDEX_AGE' %in% vars)
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  
  file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')))
  if(!is.null(standardised.vars)){
    file = paste0(file, 'standardisedby_', tolower(paste0(gsub('INDEX_', '', standardised.vars), collapse = '_')))
  }
  file = paste0(file, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
}

find_summary_output_by_round <- function(samples, output, vars, 
                                         transform = NULL, standardised.vars = NULL, names = NULL, operation = NULL, log_offset_round = NULL, 
                                         log_offset_formula = 'LOG_OFFSET', per_unsuppressed = F, per_susceptible = F, posterior_samples = F, relative_baseline = F, 
                                         invert = F, median_age_source = F, quantile_age_source = F, sex_ratio = F){
  
  # summarise outputs by round
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[[output]]) )
  if(!is.null(names)){
    setnames(tmp1, 2:(length(names) + 1), names)
  }else if(tmp1[, max(Var2)] == df_age[, max(INDEX_AGE)]){
    setnames(tmp1, 2:5, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'))
  }else{
    setnames(tmp1, 2:5, c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'INDEX_AGE'))
  }
  
  if('INDEX_AGE' %in% names(tmp1)){
    # merge to map of age and age agregated
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
    tmp1 <- merge(tmp1, df_age_aggregated, by = c('AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'))
  }
  
  if('INDEX_AGE_TRANSMISSION.SOURCE' %in% names(tmp1)){
    # find age source if output is stratified by the index of age of source
    tmp1[, AGE_TRANSMISSION.SOURCE:=unique(df_age[order(AGE_TRANSMISSION.SOURCE), .(AGE_TRANSMISSION.SOURCE)])$AGE_TRANSMISSION.SOURCE[INDEX_AGE_TRANSMISSION.SOURCE]]
  }
  
  if('INDEX_AGE_INFECTION.RECIPIENT' %in% names(tmp1)){
    # find age recipient if output is stratified by the index of age of recipient
    tmp1[, AGE_INFECTION.RECIPIENT:=unique(df_age[order(AGE_INFECTION.RECIPIENT), .(AGE_INFECTION.RECIPIENT)])$AGE_INFECTION.RECIPIENT[INDEX_AGE_INFECTION.RECIPIENT]]
  }
  
  if('INDEX_ROUND' %in% names(tmp1)){
    #  merge to map rounds and community
    tmp <- merge(df_round, df_community, by = c('COMM'))
    tmp1 <- merge(tmp1, tmp, by = c('INDEX_ROUND', 'INDEX_COMMUNITY'))
  }
  
  if(!is.null(log_offset_round)){
    # add a log offset specified by the formula
    tmp1 <- merge(tmp1, log_offset_round, by = c('ROUND', 'AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE', 'INDEX_DIRECTION', 'INDEX_COMMUNITY'))
    tmp1[, value := value + eval(rlang::parse_expr(log_offset_formula))]
  }
  
  if(!is.null(transform)){
    # transform the value
    tmp1[, value := sapply(value, transform)]
  }
  
  # sum value across pre-specified vars
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', vars)]
  if(!is.null(operation)){
    # transform the value
    tmp1 <- tmp1[, list(value = sapply(value, operation)), by = c('iterations', vars)]
  }
  
  if(!is.null(standardised.vars)){
    # standardised the value by pre-specified vars
    tmp1[, total_value := sum(value), by = c('iterations', standardised.vars)]
    tmp1[, value := value / total_value]
  }
  
  if(median_age_source){
    # take median age source
    setnames(tmp1, 'value', 'delta')
    vars <- standardised.vars
    tmp1 <- tmp1[, list(value = matrixStats::weightedMedian(AGE_TRANSMISSION.SOURCE, delta ), 
                        quantile = c('C50')), by = c('iterations', vars)]
    vars = c(vars, 'quantile')
  }
  
  if(quantile_age_source){
    # take quantile of the age of the source
    setnames(tmp1, 'value', 'delta')
    vars <- standardised.vars
    tmp1 <- tmp1[, list(value = Hmisc::wtd.quantile(x = AGE_TRANSMISSION.SOURCE, weight = delta,
                                                        probs = c(0.1, 0.25, 0.5, 0.75, 0.9), normwt = TRUE),
                        quantile = c('C10', 'C25', 'C50', 'C75', 'C90')), by = c('iterations', vars)]
    vars = c(vars, 'quantile')
  }

  if(per_unsuppressed){
    # divide by the number of unsuppressed
    tmp <- copy(eligible_count_round)
    if('AGE_TRANSMISSION.SOURCE' %in% vars)
      setnames(tmp, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
    if('INDEX_DIRECTION' %in% vars)
      tmp[, INDEX_DIRECTION := ifelse(SEX == 'M', df_direction[IS_MF == 1, INDEX_DIRECTION], df_direction[IS_MF == 0, INDEX_DIRECTION])]
    if('INDEX_COMMUNITY' %in% vars)
      tmp <- merge(tmp, df_community, by = 'COMM')
    if('INDEX_ROUND' %in% vars)
      tmp <- merge(tmp, df_round, by = c('COMM', 'ROUND'))
    
    tmp <- tmp[,list(TOTAL_INFECTED_NON_SUPPRESSED = sum(INFECTED_NON_SUPPRESSED)), by = vars]
    
    tmp1 <- merge(tmp, tmp1, by = vars)
    tmp1[, value := value / TOTAL_INFECTED_NON_SUPPRESSED]
  }
  
  if(per_susceptible){
    # divide by the number of susceptible
    tmp <- copy(eligible_count_round)
    if('AGE_GROUP_INFECTION.RECIPIENT' %in% vars){
      setnames(tmp, 'AGEYRS', 'AGE_INFECTION.RECIPIENT')
      tmp <- merge(tmp, unique(df_age_aggregated[, .(AGE_INFECTION.RECIPIENT, AGE_GROUP_INFECTION.RECIPIENT)]), by = c('AGE_INFECTION.RECIPIENT'))
    }
    if('AGE_GROUP_TRANSMISSION.SOURCE' %in% vars){
      setnames(tmp, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
      tmp <- merge(tmp, unique(df_age_aggregated[, .(AGE_TRANSMISSION.SOURCE, AGE_GROUP_TRANSMISSION.SOURCE)]), by = c('AGE_TRANSMISSION.SOURCE'))
    }
    if('INDEX_DIRECTION' %in% vars){
      if('AGE_GROUP_TRANSMISSION.SOURCE' %in% vars)
        tmp[, INDEX_DIRECTION := ifelse(SEX == 'M', df_direction[IS_MF == 1, INDEX_DIRECTION], df_direction[IS_MF == 0, INDEX_DIRECTION])]
      if('AGE_GROUP_INFECTION.RECIPIENT' %in% vars)
        tmp[, INDEX_DIRECTION := ifelse(SEX == 'M', df_direction[IS_MF == 0, INDEX_DIRECTION], df_direction[IS_MF == 1, INDEX_DIRECTION])]
    }
    if('INDEX_COMMUNITY' %in% vars)
      tmp <- merge(tmp, df_community, by = 'COMM')
    if('INDEX_ROUND' %in% vars)
      tmp <- merge(tmp, df_round, by = c('COMM', 'ROUND'))
    
    tmp <- tmp[,list(TOTAL_SUSCEPTIBLE = sum(SUSCEPTIBLE)), by = vars]
    
    tmp1 <- merge(tmp, tmp1, by = vars)
    tmp1[, value := value / TOTAL_SUSCEPTIBLE]
  }
  
  if(relative_baseline){
    # take value relative to first round
    vars.without.index.round <- vars[which(vars != 'INDEX_ROUND')]
    tmp1[, value_baseline := value[INDEX_ROUND == min(INDEX_ROUND)], by = c('iterations', vars.without.index.round)]
    tmp1[, value := value / value_baseline]
  }
  
  if(invert){
    # invert
    tmp1[, value := 1 / value ]
  }
  
  if(sex_ratio){
    # take ratio of the value by sex
    tmp1 <- select(tmp1, - 'total_value')
    tmp1 <- dcast(tmp1, ... ~ INDEX_DIRECTION, value.var = 'value')
    setnames(tmp1, c('1', '2'), c('value_FM', 'value_MF'))
    tmp1[, value := value_MF / value_FM]
    vars <- vars[-which(vars == 'INDEX_DIRECTION')]
  }
  
  if(posterior_samples == T){
    # return the posterior samples
    return(tmp1)
  }
  
  #summarise
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=vars]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")
  
  # merge by all the maps
  if('INDEX_DIRECTION' %in% vars)
    tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  if('INDEX_COMMUNITY' %in% vars)
    tmp1 <- merge(tmp1, df_community, by = 'INDEX_COMMUNITY')
  if('INDEX_ROUND' %in% vars)
    tmp1 <- merge(tmp1, df_round, by = c('INDEX_ROUND', 'COMM'))
  if('INDEX_AGE' %in% vars)
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  if('INDEX_TIME' %in% vars)
    tmp1 <- merge(tmp1, df_period, by = c('INDEX_TIME', 'COMM'))
  
  file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')))
  if(!is.null(standardised.vars)){
    file = paste0(file, 'standardisedby_', tolower(paste0(gsub('INDEX_', '', standardised.vars), collapse = '_')))
  }
  file = paste0(file, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
}

find_eligible_count_round_95suppression_given_ART <- function(eligible_count_smooth, proportion_prevalence, proportion_unsuppressed){
  
  eligible_count_round_95suppression_given_ART <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence)
  eligible_count_round_95suppression_given_ART[, ROUND := paste0('R0', ROUND)]
  
  tmp <- copy(proportion_unsuppressed)
  # prop_unsuppressed = 1 - prop_suppressed
  # prop_unsuppressed = 1 - prop_art * 0.95
  # prop_unsuppressed = 1 - (1- prop_not_art) * 0.95
  # note that in our central def PROP_UNSUPPRESSED_M = prop_not_art
  tmp[, PROP_UNSUPPRESSED_M := 1 - (1-PROP_UNSUPPRESSED_M)*0.95 ]  
  
  tmp <- merge(eligible_count_round_95suppression_given_ART, tmp, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  if(only.participant.treated){
    # assuming that only participants are treated and non-participants are all unsuppressed
    tmp[, INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_M + INFECTED * (1-PARTICIPATION) * 1]
  }else{
    # assuming that participant and non-participant are treated at the same proportion
    tmp[, INFECTED_NON_SUPPRESSED := INFECTED * PROP_UNSUPPRESSED_M]
  }
  return(tmp)
}

find_spreaders <- function(expected_contribution_age_source, outdir){
  

  find.index.mass <- function(x, p){
    # iterate over the index of x until the sum of x at the indices sum to p
    # start with the index with the greatest x value and enlarge on the left or on the right
    # by the following indices with the greatest x value

    stopifnot(p <= sum(x))
    
    n <- length(x)
    
    index.max <- which.max(x)
    index.groups <- index.max
    enter <- T
    
    while(enter){
      if(max(index.groups) == n){
        new.index = min(index.groups) - 1
      }else if(min(index.groups) == 1){
        new.index = max(index.groups) + 1
      }else{
        proposed.indices <- c(min(index.groups) - 1, max(index.groups) + 1)
        new.index <- proposed.indices[which.max(x[proposed.indices])]
      }
      
      new.index.groups <- c(index.groups, new.index)
      if(sum(x[new.index.groups]) > p){
        enter = F
      }else if(sum(x[new.index.groups]) == p){
        index.groups <- c(index.groups, new.index)
        enter = F
      } else{
        index.groups <- c(index.groups, new.index)
      }
    }
    
    # check
    stopifnot(x[index.groups] <= p)
    
    
    return(sort(index.groups))
  }
  
  # age groups that contribute to 33%
  spreaders <- expected_contribution_age_source[, list(AGEYRS = AGE_TRANSMISSION.SOURCE[find.index.mass(M, 1/3)]), by = c('COMM', 'ROUND', 'LABEL_DIRECTION')]
  spreaders[, first_main_spreader := T]
  
  # age groups that contribute to 66%
  tmp <- expected_contribution_age_source[, list(AGEYRS = AGE_TRANSMISSION.SOURCE[find.index.mass(M, 2/3)]), by = c('COMM', 'ROUND', 'LABEL_DIRECTION')]
  tmp[, second_main_spreader := T]
  spreaders <- merge(spreaders, tmp, by = c('COMM', 'ROUND', 'LABEL_DIRECTION', 'AGEYRS'), all.x = T, all.y = T)
  spreaders[is.na(first_main_spreader), first_main_spreader := F]
  stopifnot(nrow( spreaders[first_main_spreader == T & second_main_spreader == F]) == 0)
  
  # age groups that contribute to 100%
  tmp <- expected_contribution_age_source[, list(AGEYRS = AGE_TRANSMISSION.SOURCE[find.index.mass(M, sum(M))]), by = c('COMM', 'ROUND', 'LABEL_DIRECTION')]
  tmp[, third_main_spreader := T]
  spreaders <- merge(spreaders, tmp, by = c('COMM', 'ROUND', 'LABEL_DIRECTION', 'AGEYRS'), all.x = T, all.y = T)
  spreaders[is.na(first_main_spreader), first_main_spreader := F]
  spreaders[is.na(second_main_spreader), second_main_spreader := F]
  
  # make sex label
  spreaders[, SEX := 'F']
  spreaders[LABEL_DIRECTION == 'Male -> Female', SEX := 'M']
  set(spreaders, NULL, 'LABEL_DIRECTION', NULL)
  
  # check
  stopifnot(nrow( spreaders[first_main_spreader == T & (second_main_spreader == F | third_main_spreader == F)]) == 0)
  stopifnot(nrow( spreaders[second_main_spreader == T & third_main_spreader == F]) == 0)
  
  # add category of spreader in one column
  spreaders[third_main_spreader == T, spreader_category := 3]
  spreaders[second_main_spreader == T, spreader_category := 2]
  spreaders[first_main_spreader == T, spreader_category := 1]
  
  # enlarge
  tmp <- spreaders[first_main_spreader == T]
  tmp[, spreader_category := 1]
  df_spreaders <- tmp
  tmp <- spreaders[second_main_spreader == T]
  tmp[, spreader_category := 2]
  df_spreaders <- rbind(df_spreaders, tmp)
  tmp <- spreaders[third_main_spreader == T]
  tmp[, spreader_category := 3]
  df_spreaders <- rbind(df_spreaders, tmp)
  
  set(df_spreaders, NULL, c('first_main_spreader', 'second_main_spreader', 'third_main_spreader'), NULL)
  
  # label
  df_spreaders[, type := 'main spreaders']
  df_spreaders[, label := 100]
  df_spreaders[spreader_category == 2, label := 66]
  df_spreaders[spreader_category == 1, label := 33]
  
  
  # save
  file = paste0(outdir, '-output-spreaders.rds')
  saveRDS(df_spreaders, file)
  
  
  return(df_spreaders)
}

make_counterfactual_target <- function(samples, spreaders, log_offset_round, stan_data, 
                                       eligible_count_smooth, proportion_unsuppressed, proportion_prevalence, 
                                       only_participant = F, art_up_to_female = F, outdir){
  
  # find treated male if they were all treated
  all.males <- spreaders[spreader_category == 3]
  all.males[, PROPORTION_TO_CONSIDER := 1]
  eligible_count_round.all <- find_counterfactual_unsuppressed_count_target(all.males, copy(eligible_count_smooth), copy(proportion_unsuppressed), 
                                                                            copy(proportion_prevalence), stan_data, only_participant, art_up_to_female)
  
  # select spreader that contributes to 33% 
  selected.spreaders <- spreaders[spreader_category == 1, .(COMM, ROUND, AGEYRS, SEX, type)]
  selected.spreaders[, PROPORTION_TO_CONSIDER := 1]
  # find budget 
  eligible_count_round.spreaders <- find_counterfactual_unsuppressed_count_target(copy(selected.spreaders), copy(eligible_count_smooth), copy(proportion_unsuppressed), 
                                                                                  copy(proportion_prevalence), stan_data, only_participant, art_up_to_female)
  selected.spreaders <- merge(selected.spreaders,  eligible_count_round.spreaders[, .(ROUND, SEX, COMM, AGEYRS, TREATED)], by= c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  budget <- selected.spreaders[, list(TREATED.SPREADERS = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM')]
  stopifnot(nrow(selected.spreaders[TREATED < 0 & ROUND == 'R018']) == 0)
  
  # find age groups with the largest difference in art uptake compared to frmale
  noncomplier <- find_male_with_greatest_art_diff(proportion_unsuppressed, eligible_count_round.all, only_participant, only.participant.treated, budget)
  eligible_count_round.noncomplier <- find_counterfactual_unsuppressed_count_target(copy(noncomplier), copy(eligible_count_smooth), copy(proportion_unsuppressed), 
                                                                                    copy(proportion_prevalence), stan_data, only_participant, art_up_to_female)
  tmp <- eligible_count_round.noncomplier[, list(TREATED.NONCOMPLIER = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM')]
  budget <- merge(budget, tmp, by = c('ROUND', 'SEX', 'COMM'))
  stopifnot(nrow(eligible_count_round.noncomplier[TREATED < 0 & ROUND == 'R018' & SEX == 'M']) == 0)
  
  # find random selection of age groups
  set.seed(12)
  random_male <- find_random_male(eligible_count_round.all, only_participant, only.participant.treated, budget)
  eligible_count_round.random <- find_counterfactual_unsuppressed_count_target(copy(random_male), copy(eligible_count_smooth), copy(proportion_unsuppressed), 
                                                                               copy(proportion_prevalence), stan_data, only_participant, art_up_to_female)
  tmp <- eligible_count_round.random[, list(TREATED.RANDOM = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM')]
  budget <- merge(budget, tmp, by = c('ROUND', 'SEX', 'COMM'))
  stopifnot(nrow(eligible_count_round.random[TREATED < 0 & ROUND == 'R018'  & SEX == 'M']) == 0)
  
  # check budget is the same for all scenario of counterfactul
  budget[SEX == 'M']
  stopifnot(budget[SEX == 'M'& ROUND == 'R018',all(abs(TREATED.SPREADERS - TREATED.NONCOMPLIER) <  1e-10)])
  stopifnot(budget[SEX == 'M'& ROUND == 'R018',all(abs(TREATED.SPREADERS - TREATED.RANDOM) <  1e-10)])
  
  # combine spreader and noncompliers
  male_to_treat <- do.call('rbind', list(selected.spreaders, noncomplier, random_male))
  eligible_count_round.counterfactual.list <- list(eligible_count_round.spreaders, eligible_count_round.noncomplier, eligible_count_round.random)
  
  # find unsuppressed and relative incidence under counterfactual scenarios
  n_counterfactual <- male_to_treat[, length(unique(type))]
  eligible_count_round.counterfactual <- incidence_counterfactual <- vector(mode = 'list', length = n_counterfactual)
  relative_incidence_counterfactual <- vector(mode = 'list', length = n_counterfactual)
  relative_incidence_counterfactual_all <- vector(mode = 'list', length = n_counterfactual)
  for(i in 1:n_counterfactual){
    
    Type = c("main spreaders", "non compliers", "random" )[i]
    
    # target
    target <- male_to_treat[type == Type]
    
    # find unsuppressed under counterfactual
    eligible_count_round.counterfactual[[i]] <- eligible_count_round.counterfactual.list[[i]]
    
    # find offset under counterfactual
    log_offset_round.counterfactual <- find_log_offset_by_round(stan_data, copy(eligible_count_round.counterfactual[[i]]))
    
    # find incidence counterfactual by age of the recipient
    incidence_counterfactual[[i]] <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), 
                                                                  transform = 'exp', 
                                                                  log_offset_round = log_offset_round.counterfactual, 
                                                                  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')
    # find relative difference incidence  by age of the recipient
    relative_incidence_counterfactual[[i]] <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                                                     log_offset_round, log_offset_round.counterfactual,
                                                                                     transform = 'exp')
    
    # find relative difference incidence 
    relative_incidence_counterfactual_all[[i]] <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                                                         log_offset_round, log_offset_round.counterfactual,
                                                                                         transform = 'exp')
    
    # tag with the index of the counterfactual
    incidence_counterfactual[[i]][, counterfactual_index := i]
    relative_incidence_counterfactual[[i]][, counterfactual_index := i]
    eligible_count_round.counterfactual[[i]][, counterfactual_index := i]
    relative_incidence_counterfactual_all[[i]][, counterfactual_index := i]
  }
  incidence_counterfactual <- do.call('rbind', incidence_counterfactual)
  relative_incidence_counterfactual <- do.call('rbind', relative_incidence_counterfactual)
  eligible_count_round.counterfactual <- do.call('rbind', eligible_count_round.counterfactual)
  relative_incidence_counterfactual_all <- do.call('rbind', relative_incidence_counterfactual_all)
  
  # group
  counterfactuals <- list(incidence_counterfactual = incidence_counterfactual, 
                          relative_incidence_counterfactual = relative_incidence_counterfactual, 
                          eligible_count_round.counterfactual = eligible_count_round.counterfactual, 
                          relative_incidence_counterfactual_all = relative_incidence_counterfactual_all,
                          budget = budget)
  
  # save
  file = paste0(outdir, '-output-counterfactuals-target')
  if(only_participant){
    file <- paste0(file, '_only_participant')
  }
  if(art_up_to_female){
    file <- paste0(file, '_art_up_to_female')
  }
  file <- paste0(file, '.rds')
  
  saveRDS(counterfactuals, file)
  
  return(counterfactuals)
}

make_counterfactual <- function(samples, targeted.males, log_offset_round, stan_data, 
                                eligible_count_smooth, proportion_unsuppressed, proportion_prevalence, 
                                only_participant = F, art_up_to_female = F, outdir){
  
  # find unsuppressed and relative incidence under counterfactual scenarios
  n_counterfactual <- targeted.males[, length(unique(category))]
  eligible_count_round.counterfactual <- incidence_counterfactual <- vector(mode = 'list', length = n_counterfactual)
  relative_incidence_counterfactual <- vector(mode = 'list', length = n_counterfactual)
  relative_incidence_counterfactual_all <- budget <- vector(mode = 'list', length = n_counterfactual)
  for(i in 1:n_counterfactual){
    
    # target
    selected_males <- targeted.males[category == i, .(COMM, ROUND, AGEYRS, SEX)]
    
    # find unsuppressed under counterfactual
    eligible_count_round.counterfactual[[i]] <- find_counterfactual_unsuppressed_count(copy(selected_males), copy(eligible_count_smooth), copy(proportion_unsuppressed), 
                                                                                       copy(proportion_prevalence), stan_data, only_participant, art_up_to_female)
    
    # budget
    selected_males <- merge(selected_males,  eligible_count_round.counterfactual[[i]][, .(ROUND, SEX, COMM, AGEYRS, TREATED)], by= c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
    budget[[i]] <- selected_males[, list(TREATED = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM')]
    print(selected_males[TREATED < 0 & ROUND == 'R018'])
    
    # find offset under counterfactual
    log_offset_round.counterfactual <- find_log_offset_by_round(stan_data, copy(eligible_count_round.counterfactual[[i]]))
    
    # find incidence counterfactual by age of the recipient
    incidence_counterfactual[[i]] <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), 
                                                                  transform = 'exp', 
                                                                  log_offset_round = log_offset_round.counterfactual, 
                                                                  log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED')
    # find relative difference incidence  by age of the recipient
    relative_incidence_counterfactual[[i]] <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                                                     log_offset_round, log_offset_round.counterfactual,
                                                                                     transform = 'exp')
    
    # find relative difference incidence 
    relative_incidence_counterfactual_all[[i]] <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND'),
                                                                                         log_offset_round, log_offset_round.counterfactual,
                                                                                         transform = 'exp')
    
    # tag with the index of the counterfactual
    incidence_counterfactual[[i]][, counterfactual_index := i]
    relative_incidence_counterfactual[[i]][, counterfactual_index := i]
    eligible_count_round.counterfactual[[i]][, counterfactual_index := i]
    relative_incidence_counterfactual_all[[i]][, counterfactual_index := i]
    budget[[i]][, counterfactual_index := i]
  }
  incidence_counterfactual <- do.call('rbind', incidence_counterfactual)
  relative_incidence_counterfactual <- do.call('rbind', relative_incidence_counterfactual)
  eligible_count_round.counterfactual <- do.call('rbind', eligible_count_round.counterfactual)
  relative_incidence_counterfactual_all <- do.call('rbind', relative_incidence_counterfactual_all)
  budget <- do.call('rbind', budget)
  
  # group
  counterfactuals <- list(incidence_counterfactual = incidence_counterfactual, 
                          relative_incidence_counterfactual = relative_incidence_counterfactual, 
                          eligible_count_round.counterfactual = eligible_count_round.counterfactual, 
                          relative_incidence_counterfactual_all = relative_incidence_counterfactual_all,
                          budget = budget)
  
  # save
  file = paste0(outdir, '-output-counterfactuals')
  if(only_participant){
    file <- paste0(file, '_only_participant')
  }
  if(art_up_to_female){
    file <- paste0(file, '_art_up_to_female')
  }
  file <- paste0(file, '.rds')
  
  saveRDS(counterfactuals, file)
  
  return(counterfactuals)
}


find_relative_incidence_counterfactual <- function(samples, output, vars, log_offset_round, log_offset_round.counterfactual,
                                         transform = NULL, log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED'){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # incidence rate under original scenario
  tmp1 <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round, log_offset_formula = log_offset_formula, posterior_samples = T)
  
  # incidence rate under countefactual scenario
  tmp <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round.counterfactual, log_offset_formula = log_offset_formula, posterior_samples = T)
  setnames(tmp, 'value', 'value_counterfactual')
  
  # find relative incidence
  tmp1 <- merge(tmp, tmp1, by = c(vars, 'iterations'))
  tmp1[, value := (value - value_counterfactual) / value]
  
  #summarise
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=vars]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")
  
  
  if('INDEX_DIRECTION' %in% vars)
    tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  if('INDEX_COMMUNITY' %in% vars)
    tmp1 <- merge(tmp1, df_community, by = 'INDEX_COMMUNITY')
  if('INDEX_ROUND' %in% vars)
    tmp1 <- merge(tmp1, df_round, by = c('INDEX_ROUND', 'COMM'))
  if('INDEX_AGE' %in% vars)
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  if('INDEX_TIME' %in% vars)
    tmp1 <- merge(tmp1, df_period, by = c('INDEX_TIME', 'COMM'))
  
  file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')), '_counterfactual')

  file = paste0(file, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
}

find_difference_incidence_counterfactual <- function(samples, output, vars, log_offset_round, log_offset_round.counterfactual,
                                                   transform = NULL, log_offset_formula = 'log_PROP_SUSCEPTIBLE + log_INFECTED_NON_SUPPRESSED'){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # incidence rate under original scenario
  tmp1 <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round, log_offset_formula = log_offset_formula, posterior_samples = T)
  
  # incidence rate under countefactual scenario
  tmp <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round.counterfactual, log_offset_formula = log_offset_formula, posterior_samples = T)
  setnames(tmp, 'value', 'value_counterfactual')
  
  # find difference incidence
  tmp1 <- merge(tmp, tmp1, by = c(vars, 'iterations'))
  tmp1[, value := (value - value_counterfactual) ]
  
  #summarise
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=vars]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")
  
  
  if('INDEX_DIRECTION' %in% vars)
    tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  if('INDEX_COMMUNITY' %in% vars)
    tmp1 <- merge(tmp1, df_community, by = 'INDEX_COMMUNITY')
  if('INDEX_ROUND' %in% vars)
    tmp1 <- merge(tmp1, df_round, by = c('INDEX_ROUND', 'COMM'))
  if('INDEX_AGE' %in% vars)
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  if('INDEX_TIME' %in% vars)
    tmp1 <- merge(tmp1, df_period, by = c('INDEX_TIME', 'COMM'))
  
  file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')), '_diff_counterfactual')
  
  file = paste0(file, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
}




find_counterfactual_unsuppressed_count_target <- function(selected.spreaders, eligible_count_smooth, proportion_unsuppressed, proportion_prevalence, stan_data, 
                                                   only_participant = F, art_up_to_female = F){
  
  # find proportion of unsuppressed female
  proportion_unsuppressed.counterfactual = copy(proportion_unsuppressed)
  proportion_unsuppressed.counterfactual[, PROP_UNSUPPRESSED_M.FEMALE := PROP_UNSUPPRESSED_M[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND')]
  
  # find age group for which the unsuppressed level should change
  selected.spreaders[, spreader := T]
  proportion_unsuppressed.counterfactual <- merge(proportion_unsuppressed.counterfactual, selected.spreaders, all.x = T, by = c('AGEYRS', 'COMM', 'ROUND', 'SEX'))
  proportion_unsuppressed.counterfactual[is.na(spreader), spreader := F]
  proportion_unsuppressed.counterfactual[is.na(PROPORTION_TO_CONSIDER), PROPORTION_TO_CONSIDER := 0]
  # proportion_unsuppressed.counterfactual[is.na(TREATED), TREATED := 0]
  
  # find percentage of reduction ART coverage
  proportion_unsuppressed.counterfactual[, INCREASE_ART_COVERAGE := 0]
  if(art_up_to_female){
    proportion_unsuppressed.counterfactual[spreader == T, INCREASE_ART_COVERAGE := (1 - PROP_UNSUPPRESSED_M.FEMALE) - (1 - PROP_UNSUPPRESSED_M ) ]
  }else{
    proportion_unsuppressed.counterfactual[spreader == T, INCREASE_ART_COVERAGE := (1 - 0) - (1 - PROP_UNSUPPRESSED_M ) ]
  }
  
  # set proportion unsuppressed of male to be the same as female for specific age groups
  proportion_unsuppressed.counterfactual[, PROP_UNSUPPRESSED_M.COUNTERFACTUAL := 1 - (1-PROP_UNSUPPRESSED_M)*0.95 ]
  if(art_up_to_female){
    proportion_unsuppressed.counterfactual[spreader == T, PROP_UNSUPPRESSED_M.COUNTERFACTUAL := 1 - (1-PROP_UNSUPPRESSED_M.FEMALE)*0.95] 
  } else{
    proportion_unsuppressed.counterfactual[spreader == T, PROP_UNSUPPRESSED_M.COUNTERFACTUAL := 0.05]
  }
  
  # find infected 
  eligible_count_round.counterfactual <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence)
  eligible_count_round.counterfactual[, ROUND := paste0('R0', ROUND)]
  df <- merge(eligible_count_round.counterfactual, proportion_unsuppressed.counterfactual, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # find unsuppressed factual but add the proportion of being unsuppressed given ART
  tmp <- eligible_count_round_95suppression_given_ART[, .(ROUND, SEX, COMM, AGEYRS, INFECTED_NON_SUPPRESSED)]
  setnames(tmp, 'INFECTED_NON_SUPPRESSED', 'INFECTED_NON_SUPPRESSED.FACTUAL')
  df <- merge(df, tmp, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # find unsuppressed counterfactual
  par <- copy(participation)
  par[, ROUND := paste0('R0', ROUND)]
  df <- merge(df, par, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  if(only_participant){# counterfactual
    
    if(only.participant.treated){ # baseline assumption
      df[spreader == T & SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_M.COUNTERFACTUAL + 
           INFECTED * PARTICIPATION * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_M + 
           INFECTED * (1-PARTICIPATION) * 1]
    }else{
      df[spreader == T& SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROPORTION_TO_CONSIDER  * PROP_UNSUPPRESSED_M.COUNTERFACTUAL + 
           INFECTED * PARTICIPATION * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_M + 
           INFECTED * (1-PARTICIPATION) * PROP_UNSUPPRESSED_M]
    }
    
  }else{
    
    df[spreader == T& SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_M.COUNTERFACTUAL + 
         INFECTED * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_M ]
    
  }
  
  df[spreader == F|SEX == 'F', INFECTED_NON_SUPPRESSED := INFECTED_NON_SUPPRESSED.FACTUAL]
  
  
  # find difference in unsuppressed
  df[, TREATED := INFECTED_NON_SUPPRESSED.FACTUAL - INFECTED_NON_SUPPRESSED]

  return(df)
}

find_counterfactual_unsuppressed_count <- function(targeted_males, eligible_count_smooth, 
                                                   proportion_unsuppressed, proportion_prevalence, stan_data, 
                                                   only_participant = F, art_up_to_female = F){
  
  # find proportion of unsuppressed female
  proportion_unsuppressed.counterfactual = copy(proportion_unsuppressed)
  proportion_unsuppressed.counterfactual[, PROP_UNSUPPRESSED_M.FEMALE := PROP_UNSUPPRESSED_M[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND')]
  
  # find age group for which the unsuppressed level should change
  targeted_males[, target := T]
  proportion_unsuppressed.counterfactual <- merge(proportion_unsuppressed.counterfactual, targeted_males, all.x = T, by = c('AGEYRS', 'COMM', 'ROUND', 'SEX'))
  proportion_unsuppressed.counterfactual[is.na(target), target := F]

  # find percentage of reduction ART coverage
  proportion_unsuppressed.counterfactual[, INCREASE_ART_COVERAGE := 0]
  if(art_up_to_female){
    proportion_unsuppressed.counterfactual[target == T, INCREASE_ART_COVERAGE := (1 - PROP_UNSUPPRESSED_M.FEMALE) - (1 - PROP_UNSUPPRESSED_M ) ]
  }else{
    proportion_unsuppressed.counterfactual[target == T, INCREASE_ART_COVERAGE := (1 - 0) - (1 - PROP_UNSUPPRESSED_M ) ]
  }
  
  # set proportion unsuppressed of male to be the same as female for specific age groups
  proportion_unsuppressed.counterfactual[, PROP_UNSUPPRESSED_M := 1 - (1-PROP_UNSUPPRESSED_M)*0.95 ]
  proportion_unsuppressed.counterfactual[, PROP_UNSUPPRESSED_M.COUNTERFACTUAL := PROP_UNSUPPRESSED_M]
  if(art_up_to_female){
    proportion_unsuppressed.counterfactual[target == T, PROP_UNSUPPRESSED_M.COUNTERFACTUAL := 1 - (1-PROP_UNSUPPRESSED_M.FEMALE)*0.95] 
  } else{
    proportion_unsuppressed.counterfactual[target == T, PROP_UNSUPPRESSED_M.COUNTERFACTUAL := 0.05]
  }
  
  # find infected 
  eligible_count_round.counterfactual <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence)
  eligible_count_round.counterfactual[, ROUND := paste0('R0', ROUND)]
  df <- merge(eligible_count_round.counterfactual, proportion_unsuppressed.counterfactual, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # find unsuppressed factual but add the proportion of being unsuppressed given ART
  tmp <- eligible_count_round_95suppression_given_ART[, .(ROUND, SEX, COMM, AGEYRS, INFECTED_NON_SUPPRESSED)]
  setnames(tmp, 'INFECTED_NON_SUPPRESSED', 'INFECTED_NON_SUPPRESSED.FACTUAL')
  df <- merge(df, tmp, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # find unsuppressed counterfactual
  par <- copy(participation)
  par[, ROUND := paste0('R0', ROUND)]
  df <- merge(df, par, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  if(only_participant){# counterfactual
    
    if(only.participant.treated){ # baseline assumption
      df[target == T & SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * (1 - 0.95^2 * (1-PROP_UNSUPPRESSED_M.COUNTERFACTUAL)) + 
           INFECTED * (1-PARTICIPATION) * 1]
    }else{
      df[target == T& SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION  * (1 - 0.95^2 * (1-PROP_UNSUPPRESSED_M.COUNTERFACTUAL)) + 
           INFECTED * (1-PARTICIPATION) * PROP_UNSUPPRESSED_M]
    }
    
  }else{
    
    df[target == T& SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED  * (1 - 0.95^2 * (1-PROP_UNSUPPRESSED_M.COUNTERFACTUAL)) ]
    
  }
  
  df[target == F|SEX == 'F', INFECTED_NON_SUPPRESSED := INFECTED_NON_SUPPRESSED.FACTUAL]
  
  
  # find difference in unsuppressed
  df[, TREATED := INFECTED_NON_SUPPRESSED.FACTUAL - INFECTED_NON_SUPPRESSED]
  
  return(df)
}




find_male_with_greatest_art_diff <- function(proportion_unsuppressed, eligible_count_round.all, only_participant, only.participant.treated, budget){
  
  # find difference in art uptake between male and female
  ppu <- copy(proportion_unsuppressed)
  ppu[, PROP_UNSUPPRESSED_M.FEMALE := PROP_UNSUPPRESSED_M[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND')]
  ppu[, INCREASE_ART_COVERAGE := (1 - PROP_UNSUPPRESSED_M.FEMALE) - (1 - PROP_UNSUPPRESSED_M )]
  ppu <- ppu[SEX == 'M']
  
  # get the budget of treating male
  eli <- copy(eligible_count_round.all[, .(SEX,COMM,AGEYRS,ROUND,TREATED, INFECTED, PROP_UNSUPPRESSED_M.COUNTERFACTUAL, INFECTED_NON_SUPPRESSED.FACTUAL, PARTICIPATION)])

  # find age groups to consider
  ppc <- merge(ppu, eli, by = c('SEX','COMM','AGEYRS','ROUND'))
  ppc <- merge(ppc, budget, by = c('SEX','COMM','ROUND'), allow.cartesian=TRUE)
  ppc <- ppc[order(COMM,ROUND,1-INCREASE_ART_COVERAGE)]
  ppc[, CUMSUM_TREATED := cumsum(TREATED), by = c('SEX', 'COMM', 'ROUND')]
  ppc[, proportion := CUMSUM_TREATED / TREATED.SPREADERS]
  ppc[, to_consider := proportion < 1 | proportion == min(proportion[proportion >= 1]), by = c('SEX', 'COMM', 'ROUND')]
  ppc[SEX == 'M', to_consider := proportion < 1 | proportion == min(proportion[proportion >= 1]), by = c('SEX', 'COMM', 'ROUND')]
  ppc[SEX == 'F', to_consider := 0]
  ppc <- ppc[ to_consider == 1]
  
  # find proportion to consider
  ppc[proportion > 1, UNSUPPRESSED.COUNTERFACTUAL := INFECTED_NON_SUPPRESSED.FACTUAL - (TREATED.SPREADERS - (CUMSUM_TREATED - TREATED))  ]
  if(only_participant){
    if(only.participant.treated){
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_M + (1-PARTICIPATION))) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else{
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_M) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }
  }else{
    ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_M) / (INFECTED *(PROP_UNSUPPRESSED_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
  }

  # keep variable of interest
  ppc <- ppc[, .(SEX, COMM, ROUND, AGEYRS, TREATED, PROPORTION_TO_CONSIDER)]

  # add tupe
  ppc[, type := 'non compliers']
  
  return(ppc)
}

find_random_male <- function(eligible_count_round.all,only_participant, only.participant.treated, budget){
  
  # find random ordering of males
  ppu <- copy(eligible_count_round.all[, .(ROUND, SEX, COMM, AGEYRS)])
  ppu[, ORDER := sample(1:length(AGEYRS)), by = c('SEX', 'COMM', 'ROUND')]

  # get the budget of treating male
  eli <- copy(eligible_count_round.all[, .(SEX,COMM,AGEYRS,ROUND,TREATED, INFECTED, PROP_UNSUPPRESSED_M, PROP_UNSUPPRESSED_M.COUNTERFACTUAL, INFECTED_NON_SUPPRESSED.FACTUAL, PARTICIPATION)])
  
  # find age groups to consider
  ppc <- merge(ppu, eli, by = c('SEX','COMM','AGEYRS','ROUND'))
  ppc <- merge(ppc, budget, by = c('SEX','COMM','ROUND'), allow.cartesian=TRUE)
  ppc <- ppc[order(COMM,ROUND,ORDER)]
  ppc[, CUMSUM_TREATED := cumsum(TREATED), by = c('SEX', 'COMM', 'ROUND')]
  ppc[, proportion :=  CUMSUM_TREATED / TREATED.SPREADERS]
  ppc[SEX == 'M', to_consider := proportion < 1 | proportion == min(proportion[proportion >= 1]), by = c('SEX', 'COMM', 'ROUND')]
  ppc[SEX == 'F', to_consider := 0]

  ppc <- ppc[ to_consider == 1]
  
  # find proportion to consider
  ppc[proportion > 1, UNSUPPRESSED.COUNTERFACTUAL := INFECTED_NON_SUPPRESSED.FACTUAL - (TREATED.SPREADERS - (CUMSUM_TREATED - TREATED))  ]
  if(only_participant){
    if(only.participant.treated){
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_M + (1-PARTICIPATION))) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else{
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_M) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }
  }else{
    ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_M) / (INFECTED *(PROP_UNSUPPRESSED_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
  }
  
  # keep variable of interest
  ppc <- ppc[, .(SEX, COMM, ROUND, AGEYRS, TREATED, PROPORTION_TO_CONSIDER)]
  
  # add tupe
  ppc[, type := 'random']
  
  return(ppc)
}



find_male_with_greatest_art_difference_category <- function(eligible_count_round_95suppression_given_ART, outdir){
  
  find.index.mass <- function(x, y, p){
    
    stopifnot(p <= sum(y))
    
    n <- length(x)
    
    index.max <- which.max(x)
    index.groups <- index.max
    enter <- T
    
    while(enter){
      if(max(index.groups) == n){
        new.index = min(index.groups) - 1
      }else if(min(index.groups) == 1){
        new.index = max(index.groups) + 1
      }else{
        proposed.indices <- c(min(index.groups) - 1, max(index.groups) + 1)
        new.index <- proposed.indices[which.max(x[proposed.indices])]
      }
      
      new.index.groups <- c(index.groups, new.index)
      if(sum(y[new.index.groups]) > p){
        enter = F
      }else if(sum(y[new.index.groups]) == p){
        index.groups <- c(index.groups, new.index)
        enter = F
      } else{
        index.groups <- c(index.groups, new.index)
      }
    }
    
    # check
    stopifnot(y[index.groups] <= p)
    
    return(sort(index.groups))
  }
  
  DF <- copy(eligible_count_round_95suppression_given_ART)
  DF[, PROP_UNSUPPRESSED_M.FEMALE := PROP_UNSUPPRESSED_M[SEX == 'F'], by = c('ROUND', 'COMM', 'AGEYRS')]
  DF[, INCREASE_ART_COVERAGE := (1 - PROP_UNSUPPRESSED_M.FEMALE) - (1 - PROP_UNSUPPRESSED_M ) ]
  DF[, PROP_INFECTED_NON_SUPPRESSED := INFECTED_NON_SUPPRESSED / sum(INFECTED_NON_SUPPRESSED),  by = c('ROUND', 'COMM', 'SEX')]
  
  # age groups of targeted males representing 33% of infected unsuppressed
  targeted_males <- DF[, list(AGEYRS = AGEYRS[find.index.mass(INCREASE_ART_COVERAGE, PROP_INFECTED_NON_SUPPRESSED, 1/3)]), by = c('COMM', 'ROUND', 'SEX')]
  targeted_males[, first_targeted_males := T]
  
  # age groups of targeted males representing 66% of infected unsuppressed
  tmp <- DF[, list(AGEYRS = AGEYRS[find.index.mass(INCREASE_ART_COVERAGE, PROP_INFECTED_NON_SUPPRESSED, 2/3)]), by = c('COMM', 'ROUND', 'SEX')]
  tmp[, second_targeted_males := T]
  targeted_males <- merge(targeted_males, tmp, by = c('COMM', 'ROUND', 'SEX', 'AGEYRS'), all.x = T, all.y = T)
  targeted_males[is.na(first_targeted_males), first_targeted_males := F]
  stopifnot(nrow( targeted_males[first_targeted_males == T & second_targeted_males == F]) == 0)
  
  # all males
  tmp <- DF[, list(AGEYRS = AGEYRS[find.index.mass(INCREASE_ART_COVERAGE, PROP_INFECTED_NON_SUPPRESSED, sum(PROP_INFECTED_NON_SUPPRESSED))]), by = c('COMM', 'ROUND', 'SEX')]
  tmp[, third_targeted_males := T]
  targeted_males <- merge(targeted_males, tmp, by = c('COMM', 'ROUND', 'SEX', 'AGEYRS'), all.x = T, all.y = T)
  targeted_males[is.na(first_targeted_males), first_targeted_males := F]
  targeted_males[is.na(second_targeted_males), second_targeted_males := F]

  # check
  stopifnot(nrow( targeted_males[first_targeted_males == T & (second_targeted_males == F | third_targeted_males == F)]) == 0)
  stopifnot(nrow( targeted_males[second_targeted_males == T & third_targeted_males == F]) == 0)
  
  # add category of spreader in one column
  targeted_males[third_targeted_males == T, category := 3]
  targeted_males[second_targeted_males == T, category := 2]
  targeted_males[first_targeted_males == T, category := 1]
  
  # enlarge
  tmp <- targeted_males[first_targeted_males == T]
  tmp[, category := 1]
  df_targeted_males <- tmp
  tmp <- targeted_males[second_targeted_males == T]
  tmp[, category := 2]
  df_targeted_males <- rbind(df_targeted_males, tmp)
  tmp <- targeted_males[third_targeted_males == T]
  tmp[, category := 3]
  df_targeted_males <- rbind(df_targeted_males, tmp)
  
  set(df_targeted_males, NULL, c('first_targeted_males', 'second_targeted_males', 'third_targeted_males'), NULL)
  
  # label
  df_targeted_males[, label := 100]
  df_targeted_males[category == 2, label := 66]
  df_targeted_males[category == 1, label := 33]
  
  
  # save
  file = paste0(outdir, '-output-targeted_males.rds')
  saveRDS(df_targeted_males, file)
  
  
  return(df_targeted_males)
}
