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
    setnames(tmp1, 2:4, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_TIME'))
    }else{
    setnames(tmp1, 2:4, c('INDEX_DIRECTION', 'INDEX_TIME', 'INDEX_AGE'))
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
  if('INDEX_TIME' %in% vars){
    tmp1 <- merge(tmp1, df_period, by = c('INDEX_TIME'))
    tmp1 <- merge(tmp1, df_community, by = 'COMM')
  }
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
                                         invert = F, median_age_source = F, quantile_age_source = F, sex_ratio = F, save_output = T){
  
  # summarise outputs by round
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[[output]]) )
  if(!is.null(names)){
    setnames(tmp1, 2:(length(names) + 1), names)
  }else if(tmp1[, max(Var2)] == df_age[, max(INDEX_AGE)]){
    setnames(tmp1, 2:4, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_ROUND'))
  }else{
    setnames(tmp1, 2:4, c('INDEX_DIRECTION', 'INDEX_ROUND', 'INDEX_AGE'))
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
    #  merge to map rounds 
    tmp1 <- merge(tmp1, df_round, by = c('INDEX_ROUND'))
  }
  
  if(!is.null(log_offset_round)){
    # add a log offset specified by the formula
    tmp1 <- merge(tmp1, log_offset_round, by = c('ROUND', 'AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE', 'INDEX_DIRECTION'))
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
    if('INDEX_ROUND' %in% vars)
      tmp <- merge(tmp, df_round, by = c('ROUND'))
    
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
    if('INDEX_ROUND' %in% vars)
      tmp <- merge(tmp, df_round, by = c('ROUND'))
    
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
  if('INDEX_ROUND' %in% vars){
    tmp1 <- merge(tmp1, df_round, by = c('INDEX_ROUND'))
    tmp1 <- merge(tmp1, df_community, by = 'COMM')
  }
  if('INDEX_AGE' %in% vars)
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  if('INDEX_TIME' %in% vars){
    tmp1 <- merge(tmp1, df_period, by = c('INDEX_TIME'))
    tmp1 <- merge(tmp1, df_community, by = 'COMM')
  }
  
  if(save_output){
    file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')))
    if(!is.null(standardised.vars)){
      file = paste0(file, 'standardisedby_', tolower(paste0(gsub('INDEX_', '', standardised.vars), collapse = '_')))
    }
    file = paste0(file, '.rds')
    saveRDS(tmp1, file)
  }
  
  return(tmp1)
}

find_relative_incidence_counterfactual <- function(samples, output, vars, log_offset_round, log_offset_round.counterfactual, transform = NULL){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # incidence rate under original scenario
  tmp1 <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round, posterior_samples = T)
  
  # incidence rate under countefactual scenario
  tmp <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round.counterfactual, posterior_samples = T)
  setnames(tmp, 'value', 'value_counterfactual')
  
  # find relative incidence
  tmp1 <- merge(tmp, tmp1, by = c(vars, 'iterations'))
  tmp1[, value := (value - value_counterfactual) / value]
  
  #summarise
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=vars]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")
  
  
  if('INDEX_DIRECTION' %in% vars)
    tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  if('INDEX_ROUND' %in% vars){
    tmp1 <- merge(tmp1, df_round, by = c('INDEX_ROUND'))
    tmp1 <- merge(tmp1, df_community, by = 'COMM')
  }
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
                                                   transform = NULL){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  # incidence rate under original scenario
  tmp1 <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round,  posterior_samples = T)
  
  # incidence rate under countefactual scenario
  tmp <- find_summary_output_by_round(samples, output, vars, transform = transform, log_offset_round = log_offset_round.counterfactual, posterior_samples = T)
  setnames(tmp, 'value', 'value_counterfactual')
  
  # find difference incidence
  tmp1 <- merge(tmp, tmp1, by = c(vars, 'iterations'))
  tmp1[, value := (value - value_counterfactual) ]
  
  #summarise
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=vars]	
  tmp1 = dcast(tmp1, ... ~ q_label, value.var = "q")
  
  
  if('INDEX_DIRECTION' %in% vars)
    tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  if('INDEX_ROUND' %in% vars){
    tmp1 <- merge(tmp1, df_round, by = c('INDEX_ROUND'))
    tmp1 <- merge(tmp1, df_community, by = 'COMM')
  }
  if('INDEX_AGE' %in% vars)
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  if('INDEX_TIME' %in% vars)
    tmp1 <- merge(tmp1, df_period, by = c('INDEX_TIME', 'COMM'))
  
  file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')), '_diff_counterfactual')
  
  file = paste0(file, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
}

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

find_spreaders <- function(expected_contribution_age_source, outdir){
  
  
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

find_male_with_greatest_supp_diff <- function(treatment_cascade, eligible_count_round.all, participation, 
                                             only_participant, only.participant.treated, budget){
  
  # find suppressed given infected in participants
  ppu <- copy(treatment_cascade[, .(ROUND, SEX, COMM, AGEYRS, PROP_UNSUPPRESSED_PARTICIPANTS_M, PROP_UNSUPPRESSED_NONPARTICIPANTS_M)])
  ppu[,PROP_SUPPRESSED_PARTICIPANTS_M := 1- PROP_UNSUPPRESSED_PARTICIPANTS_M ]
  ppu[,PROP_SUPPRESSED_NONPARTICIPANTS_M := 1- PROP_UNSUPPRESSED_NONPARTICIPANTS_M ]
  set(ppu, NULL, 'PROP_UNSUPPRESSED_NONPARTICIPANTS_M', NULL)
  set(ppu, NULL, 'PROP_UNSUPPRESSED_PARTICIPANTS_M', NULL)
  pa <- copy(participation)
  pa[, ROUND := paste0('R0', ROUND)]
  ppu <- merge(ppu, pa, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  if(nonparticipants.treated.like.participants){# baseline assumption
    ppu[, PROP_SUPPRESSED_M := PROP_SUPPRESSED_PARTICIPANTS_M * PARTICIPATION + PROP_SUPPRESSED_NONPARTICIPANTS_M * (1-PARTICIPATION)]
  }else if(nonparticipants.not.treated){
    ppu[, PROP_SUPPRESSED_M := PROP_SUPPRESSED_PARTICIPANTS_M * PARTICIPATION + 0 * (1-PARTICIPATION)]
  }else{
    ppu[, PROP_SUPPRESSED_M := PROP_SUPPRESSED_PARTICIPANTS_M * PARTICIPATION + PROP_SUPPRESSED_NONPARTICIPANTS_M * (1-PARTICIPATION)]
  }
  
  # find difference in suppression level uptake between male and female
  ppu[, PROP_SUPPRESSED_M.FEMALE := PROP_SUPPRESSED_M[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND')]
  ppu[, INCREASE_PROP_UNSUPPRESSED := PROP_SUPPRESSED_M.FEMALE - PROP_SUPPRESSED_M]
  ppu <- ppu[SEX == 'M']
  
  # get the budget of treating male
  eli <- copy(eligible_count_round.all[, .(SEX,COMM,AGEYRS,ROUND,TREATED, INFECTED, 
                                           PROP_UNSUPPRESSED_PARTICIPANTS_M, PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL, 
                                           PROP_UNSUPPRESSED_NONPARTICIPANTS_M, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL,
                                           INFECTED_NON_SUPPRESSED.FACTUAL)])
  
  # find age groups to consider by ordering the male by the different in suppression levl compared to female
  ppc <- merge(ppu, eli, by = c('SEX','COMM','AGEYRS','ROUND'))
  ppc <- merge(ppc, budget, by = c('SEX','COMM','ROUND'), allow.cartesian=TRUE)
  ppc <- ppc[order(COMM,ROUND,1-INCREASE_PROP_UNSUPPRESSED)]
  ppc[, CUMSUM_TREATED := cumsum(TREATED), by = c('SEX', 'COMM', 'ROUND')]
  ppc[, proportion := CUMSUM_TREATED / TREATED.SPREADERS]
  ppc[, to_consider := proportion < 1 | proportion == min(proportion[proportion >= 1]), by = c('SEX', 'COMM', 'ROUND')]
  
  # we keep all the age groups for which their sum < budget and the first one that goes above the budget
  ppc[SEX == 'M', to_consider := proportion < 1 | proportion == min(proportion[proportion >= 1]), by = c('SEX', 'COMM', 'ROUND')]
  ppc[SEX == 'F', to_consider := 0]
  ppc <- ppc[ to_consider == 1]
  
  # find proportion to consider
  # for the  age groups for which their sum < budget, we consider 100% of them
  # for the first age group that goes above the budget, we consider a proportion such that their additional contribution matchs the budget
  ppc[proportion > 1, UNSUPPRESSED.COUNTERFACTUAL := INFECTED_NON_SUPPRESSED.FACTUAL - (TREATED.SPREADERS - (CUMSUM_TREATED - TREATED))  ]
  if(only_participant){
    # counterfactual assuming onlu participants are targeted
    if(nonparticipants.treated.like.participants){# baseline assumption
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_PARTICIPANTS_M ) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else if(nonparticipants.not.treated){
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*1)) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else{
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*PROP_UNSUPPRESSED_NONPARTICIPANTS_M)) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }
  }else{
    # counterfactual assuming all population are targeted
    if(nonparticipants.treated.like.participants){# baseline assumption
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_PARTICIPANTS_M ) / (INFECTED *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else if(nonparticipants.not.treated){
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*1)) / (INFECTED*( PARTICIPATION * (PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) + 
                                                                                                                                                                                                  (1-PARTICIPATION)*(PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL -1) ) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else{
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*PROP_UNSUPPRESSED_NONPARTICIPANTS_M)) / (INFECTED*( PARTICIPATION * (PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) + 
                                                                                                                                                                                                                                    (1-PARTICIPATION)*(PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_NONPARTICIPANTS_M) ) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }  
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
  eli <- copy(eligible_count_round.all[, .(SEX,COMM,AGEYRS,ROUND,TREATED, INFECTED, 
                                           PROP_UNSUPPRESSED_PARTICIPANTS_M, PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL, 
                                           PROP_UNSUPPRESSED_NONPARTICIPANTS_M, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL,
                                           INFECTED_NON_SUPPRESSED.FACTUAL, PARTICIPATION)])
  
  # find age groups to consider
  ppc <- merge(ppu, eli, by = c('SEX','COMM','AGEYRS','ROUND'))
  ppc <- merge(ppc, budget, by = c('SEX','COMM','ROUND'), allow.cartesian=TRUE)
  ppc <- ppc[order(COMM,ROUND,ORDER)]
  ppc[, CUMSUM_TREATED := cumsum(TREATED), by = c('SEX', 'COMM', 'ROUND')]
  ppc[, proportion :=  CUMSUM_TREATED / TREATED.SPREADERS]
  
  # we keep all the age groups for which their sum < budget and the first one that goes above the budget
  ppc[SEX == 'M', to_consider := proportion < 1 | proportion == min(proportion[proportion >= 1]), by = c('SEX', 'COMM', 'ROUND')]
  ppc[SEX == 'F', to_consider := 0]
  ppc <- ppc[ to_consider == 1]
  
  # find proportion to consider
  # for the  age groups for which their sum < budget, we consider 100% of them
  # for the first age group that goes above the budget, we consider a proportion such that their additional contribution matchs the budget
  ppc[proportion > 1, UNSUPPRESSED.COUNTERFACTUAL := INFECTED_NON_SUPPRESSED.FACTUAL - (TREATED.SPREADERS - (CUMSUM_TREATED - TREATED))  ]
  if(only_participant){
    # counterfactual assuming onlu participants are targeted
    if(nonparticipants.treated.like.participants){# baseline assumption
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_PARTICIPANTS_M ) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else if(nonparticipants.not.treated){
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*1)) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else{
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*PROP_UNSUPPRESSED_NONPARTICIPANTS_M)) / (INFECTED*PARTICIPATION *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }
  }else{
    # counterfactual assuming all population are targeted
    if(nonparticipants.treated.like.participants){# baseline assumption
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *PROP_UNSUPPRESSED_PARTICIPANTS_M ) / (INFECTED *(PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else if(nonparticipants.not.treated){
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*1)) / (INFECTED*( PARTICIPATION * (PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) + 
                                                                                                                                                                                                  (1-PARTICIPATION)*(PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL -1) ) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }else{
      ppc[, PROPORTION_TO_CONSIDER := ifelse(proportion <= 1, 1, ( UNSUPPRESSED.COUNTERFACTUAL- INFECTED *(PARTICIPATION*PROP_UNSUPPRESSED_PARTICIPANTS_M + (1-PARTICIPATION)*PROP_UNSUPPRESSED_NONPARTICIPANTS_M)) / (INFECTED*( PARTICIPATION * (PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_PARTICIPANTS_M) + 
                                                                                                                                                                                                                                    (1-PARTICIPATION)*(PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL -PROP_UNSUPPRESSED_NONPARTICIPANTS_M) ) )), by = c('SEX', 'COMM', 'ROUND', 'AGEYRS')]
    }  
  }
  
  # keep variable of interest
  ppc <- ppc[, .(SEX, COMM, ROUND, AGEYRS, TREATED, PROPORTION_TO_CONSIDER)]
  
  # add tupe
  ppc[, type := 'random']
  
  return(ppc)
}


make_counterfactual_target <- function(samples, spreaders, log_offset_round, stan_data, 
                                       eligible_count_smooth, eligible_count_round, treatment_cascade, 
                                       proportion_prevalence, participation,
                                       only_participant = F, art_up_to_female = F, outdir){
  
  #
  # find number of treated male if all age groups were treated
  # 
  
  #not a compared scenario but used to find the budget of the 2nd and 3rd scneario
  all.males <- spreaders[spreader_category == 3] #spreader_category == 3 = all age groups
  all.males[, PROPORTION_TO_CONSIDER := 1]
  eligible_count_round.all <- find_counterfactual_unsuppressed_count_target(all.males, copy(eligible_count_smooth), copy(eligible_count_round),
                                                                            copy(treatment_cascade), copy(proportion_prevalence), copy(participation),
                                                                            only_participant, art_up_to_female)
  #
  # 1st scenario: main spreaders are treated
  # 
  
  # select age groups that contributes to 33% 
  selected.spreaders <- spreaders[spreader_category == 1, .(COMM, ROUND, AGEYRS, SEX, type)]
  selected.spreaders[, PROPORTION_TO_CONSIDER := 1]
  
  # find number of treated male by age group in this scenario 
  eligible_count_round.spreaders <- find_counterfactual_unsuppressed_count_target(copy(selected.spreaders), copy(eligible_count_smooth), copy(eligible_count_round),
                                                                                  copy(treatment_cascade), copy(proportion_prevalence), copy(participation),
                                                                                  only_participant, art_up_to_female)
  selected.spreaders <- merge(selected.spreaders,  eligible_count_round.spreaders[, .(ROUND, SEX, COMM, AGEYRS, TREATED)], by= c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  budget <- selected.spreaders[, list(TREATED.SPREADERS = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM')]
  stopifnot(nrow(selected.spreaders[TREATED < 0 & ROUND == 'R018']) == 0)
  
  
  #
  # 2nd scenario: male with greatest difference in art uptake compared to female
  # 
  
  # find age groups such that the sum of the unsuppressed matches the budget of the 1st scenario
  noncomplier <- find_male_with_greatest_supp_diff(treatment_cascade, eligible_count_round.all, participation, 
                                                  only_participant, only.participant.treated, budget)
  
  # find number of treated male by age group in this scenario 
  eligible_count_round.noncomplier <- find_counterfactual_unsuppressed_count_target(copy(noncomplier), copy(eligible_count_smooth), copy(eligible_count_round),
                                                                                    copy(treatment_cascade), copy(proportion_prevalence), copy(participation),
                                                                                    only_participant, art_up_to_female)
  tmp <- eligible_count_round.noncomplier[, list(TREATED.NONCOMPLIER = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM')]
  budget <- merge(budget, tmp, by = c('ROUND', 'SEX', 'COMM'))
  stopifnot(nrow(eligible_count_round.noncomplier[TREATED < 0 & ROUND == 'R018' & SEX == 'M']) == 0)
  
  
  #
  # 3rd scenario: random male 
  # 
  
  # find random age groups such that the sum of the unsuppressed matches the budget of the 1st (and 2nd) scenario
  set.seed(12)
  random_male <- find_random_male(eligible_count_round.all, only_participant, only.participant.treated, budget)
  
  # find number of treated male by age group in this scenario 
  eligible_count_round.random <- find_counterfactual_unsuppressed_count_target(copy(random_male), copy(eligible_count_smooth), copy(eligible_count_round),
                                                                               copy(treatment_cascade), copy(proportion_prevalence), copy(participation),
                                                                               only_participant, art_up_to_female)
  tmp <- eligible_count_round.random[, list(TREATED.RANDOM = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM')]
  budget <- merge(budget, tmp, by = c('ROUND', 'SEX', 'COMM'))
  stopifnot(nrow(eligible_count_round.random[TREATED < 0 & ROUND == 'R018'  & SEX == 'M']) == 0)
  
  #
  # Combine the scenarios
  #
  
  # check budget is the same for all scenarios 
  budget[SEX == 'M']
  stopifnot(budget[SEX == 'M'& ROUND == 'R018',all(abs(TREATED.SPREADERS - TREATED.NONCOMPLIER) < 1e-10)])
  stopifnot(budget[SEX == 'M'& ROUND == 'R018',all(abs(TREATED.SPREADERS - TREATED.RANDOM) <  1e-10)])
  
  # combine treated age groups in the scenarios
  male_to_treat <- do.call('rbind', list(selected.spreaders, noncomplier, random_male))
  eligible_count_round.counterfactual.list <- list(eligible_count_round.spreaders, eligible_count_round.noncomplier, eligible_count_round.random)
  
  # find incidence under counterfactual scenarios
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
    log_offset_round.counterfactual <- find_log_offset_by_round(stan_data, copy(eligible_count_round.counterfactual[[i]]), df_estimated_contact_rates,
                                                                use_number_susceptible_offset, use_contact_rates_prior)
    
    # find incidence counterfactual by age of the recipient
    incidence_counterfactual[[i]] <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), 
                                                                  transform = 'exp', 
                                                                  log_offset_round = log_offset_round.counterfactual, 
                                                                  save_output = F)
    # find relative difference incidence  by age of the recipient
    relative_incidence_counterfactual[[i]] <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                                                     log_offset_round, 
                                                                                     log_offset_round.counterfactual,
                                                                                     transform = 'exp')
    
    # find relative difference incidence 
    relative_incidence_counterfactual_all[[i]] <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                                                                         log_offset_round, 
                                                                                         log_offset_round.counterfactual,
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

find_counterfactual_unsuppressed_count_target <- function(selected.spreaders, eligible_count_smooth, eligible_count_round,
                                                          treatment_cascade, proportion_prevalence, participation, 
                                                          only_participant = F, art_up_to_female = F){
  
  # we will make a counterfactual treatment cascade
  treatment_cascade.counterfactual = copy(treatment_cascade)
  
  # merge to participation
  par <- copy(participation)
  par[, ROUND := paste0('R0', ROUND)]
  treatment_cascade.counterfactual <- merge(treatment_cascade.counterfactual, par, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # find age group for which the proportion of suppressed should change
  selected.spreaders[, spreader := T]
  treatment_cascade.counterfactual <- merge(treatment_cascade.counterfactual, selected.spreaders, all.x = T, by = c('AGEYRS', 'COMM', 'ROUND', 'SEX'))
  treatment_cascade.counterfactual[is.na(spreader), spreader := F]
  treatment_cascade.counterfactual[is.na(PROPORTION_TO_CONSIDER), PROPORTION_TO_CONSIDER := 0]
  
  # find proportion of suppressed
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_PARTICIPANTS_M := 1-PROP_UNSUPPRESSED_PARTICIPANTS_M]
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_NONPARTICIPANTS_M := 1-PROP_UNSUPPRESSED_NONPARTICIPANTS_M]
  
  # find suppressed proprotion in female 
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_PARTICIPANTS_M.FEMALE := PROP_SUPPRESSED_PARTICIPANTS_M[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND')]
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_NONPARTICIPANTS_M.FEMALE := PROP_SUPPRESSED_NONPARTICIPANTS_M[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND')]
  
  
  #
  # find proportion suppressed under the counterfactual 
  #
  
  # baseline 
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS_M]
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := PROP_SUPPRESSED_NONPARTICIPANTS_M]
  
  # spreaders
  if(art_up_to_female){
    # in the case where proportion suppressed in male is set to be the same as in female
    
    # for participants
    treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS_M.FEMALE]
    
    # for non-participants
    if(!only_participant){# counterfactual where only participants are treated
      treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := PROP_SUPPRESSED_NONPARTICIPANTS_M.FEMALE]
    }
    
  }else{
    # in the case where proportion suppressed in male is set to 100%
    
    # for participants
    treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL := 1]
    
    # for non-participants
    if(!only_participant){# counterfactual
      treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := 1]
    }
  }
  
  
  
  #
  # find proportion unsuppressed under the counterfactual for participants
  #
  
  # baseline
  treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL := PROP_UNSUPPRESSED_PARTICIPANTS_M]
  
  # spreaders                      
  treatment_cascade.counterfactual[spreader == T, PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL := 1 - PROP_SUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL]
  
  
  #
  # find proportion unsuppressed under the counterfactual for non-participants
  #
  
  # baseline
  if(nonparticipants.treated.like.participants){
    treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := PROP_UNSUPPRESSED_PARTICIPANTS_M]
  }else if(nonparticipants.not.treated){
    treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := 1]
  }else{
    treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := PROP_UNSUPPRESSED_NONPARTICIPANTS_M]
  }
  
  # spreaders
  if(!only_participant){# counterfactual
    if(nonparticipants.treated.like.participants){# baseline assumption
      treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL]
      treatment_cascade.counterfactual[spreader == T, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL]
    }else{
      treatment_cascade.counterfactual[spreader == T, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL]
    }
  }
  
  # recall the number of infected 
  eligible_count_round.counterfactual <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence, 
                                                                  participation, 
                                                                  nonparticipants.male.relative.infection, 
                                                                  nonparticipants.female.relative.infection)
  eligible_count_round.counterfactual[, ROUND := paste0('R0', ROUND)]
  df <- merge(select(eligible_count_round.counterfactual, -'PARTICIPATION'), treatment_cascade.counterfactual, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # find unsuppressed factual 
  tmp <- eligible_count_round[, .(ROUND, SEX, COMM, AGEYRS, INFECTED_NON_SUPPRESSED)]
  setnames(tmp, 'INFECTED_NON_SUPPRESSED', 'INFECTED_NON_SUPPRESSED.FACTUAL')
  df <- merge(df, tmp, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  #
  # find unsuppressed counterfactual for spreaders
  #
  
  if(nonparticipants.treated.like.participants){
    df[spreader == T & SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL + 
         INFECTED * (1-PARTICIPATION) * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL + 
         INFECTED * PARTICIPATION * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_PARTICIPANTS_M + 
         INFECTED * (1-PARTICIPATION) * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_PARTICIPANTS_M]
  }else if(nonparticipants.not.treated){
    df[spreader == T & SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL + 
         INFECTED * (1-PARTICIPATION) * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL + 
         INFECTED * PARTICIPATION * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_PARTICIPANTS_M + 
         INFECTED * (1-PARTICIPATION) * (1-PROPORTION_TO_CONSIDER) * 1]
  }else{
    df[spreader == T & SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL + 
         INFECTED * (1-PARTICIPATION) * PROPORTION_TO_CONSIDER * PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL + 
         INFECTED * PARTICIPATION * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_PARTICIPANTS_M + 
         INFECTED * (1-PARTICIPATION) * (1-PROPORTION_TO_CONSIDER) * PROP_UNSUPPRESSED_NONPARTICIPANTS_M]
  }
  
  
  #
  # find unsuppressed counterfactual for non-spreaders
  #
  df[spreader == F|SEX == 'F', INFECTED_NON_SUPPRESSED := INFECTED_NON_SUPPRESSED.FACTUAL]
  
  # find difference in unsuppressed
  df[, TREATED := INFECTED_NON_SUPPRESSED.FACTUAL - INFECTED_NON_SUPPRESSED]
  
  return(df)
}


make_counterfactual <- function(samples, log_offset_round, stan_data, 
                                eligible_count_smooth, eligible_count_round, 
                                treatment_cascade_samples, proportion_prevalence, participation,
                                only_participant = F, art_up_to_female = 1, s959595 = NULL, s909090 = NULL, outdir){
  
  stopifnot(is.null(art_up_to_female) + is.null(s959595) + is.null(s909090) == 2)
  
  # consider only all males
  selected_males <- eligible_count_round[, .(COMM, ROUND, AGEYRS, SEX)]
  
  # find unsuppressed under counterfactual
  eligible_count_round.counterfactual <- find_counterfactual_unsuppressed_count(copy(selected_males), copy(eligible_count_smooth), copy(eligible_count_round),
                                                                                copy(treatment_cascade_samples), copy(proportion_prevalence), copy(participation),
                                                                                only_participant, art_up_to_female, s959595, s909090)
  
  # find number of male treated under counterfactual
  budget <- eligible_count_round.counterfactual[AGEYRS == 0, .(ROUND, SEX, COMM, TREATED, TREATED_CL, TREATED_CU)]

  # find offset under counterfactual
  log_offset_round.counterfactual <- find_log_offset_by_round(stan_data, copy(eligible_count_round.counterfactual),df_estimated_contact_rates,
                                                              use_number_susceptible_offset, use_contact_rates_prior)
  
  # find incidence counterfactual by age of the recipient
  incidence_counterfactual <- find_summary_output_by_round(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'), 
                                                                transform = 'exp', 
                                                                log_offset_round = log_offset_round.counterfactual, 
                                                                save_output = F)
  # find relative difference incidence  by age of the recipient
  relative_incidence_counterfactual <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
                                                                                   log_offset_round, log_offset_round.counterfactual,
                                                                                   transform = 'exp')
  
  # find relative difference incidence 
  relative_incidence_counterfactual_all <- find_relative_incidence_counterfactual(samples, 'log_beta', c('INDEX_DIRECTION', 'INDEX_ROUND'),
                                                                                       log_offset_round, log_offset_round.counterfactual,
                                                                                       transform = 'exp')
  
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
  if(!is.null(art_up_to_female)){
    file <- paste0(file, '_art_up_to_female', art_up_to_female)
  }
  if(!is.null(s959595)){
    file <- paste0(file, '_959595', s959595)
  }
  if(!is.null(s909090)){
    file <- paste0(file, '_909090', s909090)
  }
  file <- paste0(file, '.rds')
  
  saveRDS(counterfactuals, file)
  
  return(counterfactuals)
}


find_counterfactual_unsuppressed_count <- function(selected_males,  eligible_count_smooth, eligible_count_round,
                                                   treatment_cascade_samples, proportion_prevalence, participation, 
                                                   only_participant = F, art_up_to_female = 1, s959595 = NULL, s909090 = NULL){
  
  # we will make a counterfactual treatment cascade
  treatment_cascade.counterfactual = copy(treatment_cascade_samples)
  
  # merge to participation
  par <- copy(participation)
  par[, ROUND := paste0('R0', ROUND)]
  treatment_cascade.counterfactual <- merge(treatment_cascade.counterfactual, par, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # find age group for which the art coverage should change
  selected_males[, spreader := T]
  treatment_cascade.counterfactual <- merge(treatment_cascade.counterfactual, selected_males, all.x = T, by = c('AGEYRS', 'COMM', 'ROUND', 'SEX'))
  treatment_cascade.counterfactual[is.na(spreader), spreader := F]
  
  # find unsuppressed proprotion in female 
  treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_PARTICIPANTS := 1 - PROP_SUPPRESSED_PARTICIPANTS]
  treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_NONPARTICIPANTS := 1 - PROP_SUPPRESSED_NONPARTICIPANTS]
  
  # find suppressed proprotion in female 
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_PARTICIPANTS.FEMALE := PROP_SUPPRESSED_PARTICIPANTS[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND', 'iterations')]
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_NONPARTICIPANTS.FEMALE := PROP_SUPPRESSED_NONPARTICIPANTS[SEX == 'F'], by = c('AGEYRS', 'COMM', 'ROUND', 'iterations')]
  
  
  
  #
  # find proportion suppressed under the counterfactual 
  #
  
  # baseline
  
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_PARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS]
  treatment_cascade.counterfactual[, PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_NONPARTICIPANTS]
  
  # spreaders     
  if(!is.null(art_up_to_female)){
    # in the case where suppression prop in male is set to be the same as in female
    
    # for participants
    treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_PARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS + (PROP_SUPPRESSED_PARTICIPANTS.FEMALE-PROP_SUPPRESSED_PARTICIPANTS)*art_up_to_female]
    
    # for non-participants
    if(!only_participant){# counterfactual where only participants are treated
      treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_NONPARTICIPANTS + (PROP_SUPPRESSED_NONPARTICIPANTS.FEMALE-PROP_SUPPRESSED_NONPARTICIPANTS)*art_up_to_female]
    }
    
    
  }else if(!is.null(s959595)){
    # in the case where suppression prop in male is set to be 95*0.95*0.95%
    
    # for participants
    treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_PARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS + (0.95*0.95*0.95-PROP_SUPPRESSED_PARTICIPANTS)*s959595]
    
    # for non-participants
    if(!only_participant){# counterfactual where only participants are treated
      treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_NONPARTICIPANTS + (0.95*0.95*0.95-PROP_SUPPRESSED_NONPARTICIPANTS)*s959595]
    }
    
  }else{
    # in the case where diganosed in male is set to be 90*90*90%
    
    # for participants
    treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_PARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS + (0.9*0.9*0.9-PROP_SUPPRESSED_PARTICIPANTS)*s909090]
    
    # for non-participants
    if(!only_participant){# counterfactual where only participants are treated
      treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_NONPARTICIPANTS + (0.9*0.9*0.9-PROP_SUPPRESSED_NONPARTICIPANTS)*s909090]
    }
    
  }
  
  
  
  #
  # find proportion unsuppressed under the counterfactual for participants
  #
  
  # baseline
  treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_PARTICIPANTS.COUNTERFACTUAL := PROP_UNSUPPRESSED_PARTICIPANTS]
  
  # spreadrs
  treatment_cascade.counterfactual[spreader == T, PROP_UNSUPPRESSED_PARTICIPANTS.COUNTERFACTUAL := 1 - PROP_SUPPRESSED_PARTICIPANTS.COUNTERFACTUAL]
  
  
  #
  # find proportion unsuppressed under the counterfactual for non-participants
  #
  
  # baseline
  if(nonparticipants.treated.like.participants){
    treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := PROP_UNSUPPRESSED_PARTICIPANTS]
  }else if(nonparticipants.not.treated){
    treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := 1]
  }else{
    treatment_cascade.counterfactual[, PROP_UNSUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := PROP_UNSUPPRESSED_NONPARTICIPANTS]
  }
  
  # spreaders
  if(!only_participant){# counterfactual
    if(nonparticipants.treated.like.participants){# baseline assumption
      treatment_cascade.counterfactual[spreader == T, PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := PROP_SUPPRESSED_PARTICIPANTS.COUNTERFACTUAL]
      treatment_cascade.counterfactual[spreader == T, PROP_UNSUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := 1 - PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL]
    }else{
      treatment_cascade.counterfactual[spreader == T, PROP_UNSUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL := 1 - PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL]
      
    }
  }
  
  # recall the number of infected 
  eligible_count_round.counterfactual <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence,
                                                                  participation, 
                                                                  nonparticipants.male.relative.infection, 
                                                                  nonparticipants.female.relative.infection)
  eligible_count_round.counterfactual[, ROUND := paste0('R0', ROUND)]
  df <- merge(select(eligible_count_round.counterfactual, -'PARTICIPATION'),treatment_cascade.counterfactual, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  

  #
  # find unsuppressed counterfactual for spreaders
  #
  
  df[spreader == T & SEX == "M", INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS.COUNTERFACTUAL +
       INFECTED * (1-PARTICIPATION) * PROP_UNSUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL]
  
  #
  # find unsuppressed factual
  if(nonparticipants.not.treated){
    # assuming that only participants are treated and non-participants are all unsuppressed
    df[, INFECTED_NON_SUPPRESSED.FACTUAL := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS + INFECTED * (1-PARTICIPATION) * 1]
  }else if(nonparticipants.treated.like.participants){
    # assuming that participant and non-participant are diagnosed and treated at the same proportion
    df[, INFECTED_NON_SUPPRESSED.FACTUAL := INFECTED * PROP_UNSUPPRESSED_PARTICIPANTS ]
  }else{
    # assuming that participant and non-participant are treated not at the same proportion
    # and using newly registered participants to inform suppression of non-participnats
    df[, INFECTED_NON_SUPPRESSED.FACTUAL := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS + INFECTED * (1-PARTICIPATION) * PROP_UNSUPPRESSED_NONPARTICIPANTS ]
    
  }
  
  
  #
  # find unsuppressed counterfactual for non-spreaders
  #
  
  df[spreader == F|SEX == 'F', INFECTED_NON_SUPPRESSED := INFECTED_NON_SUPPRESSED.FACTUAL]
  
  #
  # find difference in unsuppressed
  #
  
  df[, TREATED := INFECTED_NON_SUPPRESSED.FACTUAL - INFECTED_NON_SUPPRESSED]
  
  #
  # select variables
  #
  
  df1 <- df[, .(ROUND, SEX, COMM, AGEYRS, iterations, 
                PROP_UNSUPPRESSED_PARTICIPANTS, PROP_UNSUPPRESSED_NONPARTICIPANTS,
                PROP_SUPPRESSED_PARTICIPANTS.COUNTERFACTUAL, PROP_SUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL, 
                PROP_UNSUPPRESSED_PARTICIPANTS.COUNTERFACTUAL, PROP_UNSUPPRESSED_NONPARTICIPANTS.COUNTERFACTUAL, 
                INFECTED_NON_SUPPRESSED, INFECTED_NON_SUPPRESSED.FACTUAL, TREATED)]
  
  #
  # find total treated across ages
  #
  
  df2 <- df1[, list(TREATED = sum(TREATED)), by = c('ROUND', 'SEX', 'COMM', 'iterations')]
  df2[, AGEYRS := 0]
  df1 <- rbind(df1, df2,fill=TRUE)
  
  #
  # summarise
  #
  
  ps <- c(0.025,0.5,0.975)
  qlab <- c('CL','M','CU')
  
  df1 <- melt.data.table(df1, id.vars = c('ROUND', 'SEX', 'COMM', 'AGEYRS', 'iterations'))
  ns = df1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=paste0(variable, '_', qlab)), by=c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'variable')]
  ns = as.data.table(reshape2::dcast(ns, AGEYRS + SEX + COMM + ROUND  ~ q_label, value.var = "q"))
  ns <- merge(ns, unique(df[, .(ROUND, SEX, COMM, AGEYRS, ELIGIBLE, EMPIRICAL_PREVALENCE, PREVALENCE_CL, PREVALENCE_CU, PREVALENCE_M, 
                                PARTICIPANT.x, INFECTED, SUSCEPTIBLE, PARTICIPANT.y, PARTICIPATION, spreader)]), 
              by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'), all.x = T)
  setnames(ns, 'INFECTED_NON_SUPPRESSED_M', 'INFECTED_NON_SUPPRESSED')
  setnames(ns, 'INFECTED_NON_SUPPRESSED.FACTUAL_M', 'INFECTED_NON_SUPPRESSED.FACTUAL')
  setnames(ns, 'TREATED_M', 'TREATED')
  
  # subset to community
  ns <- ns[COMM %in% df_community[, unique(COMM)]]
  
  # if treeated negtive then do not do the counterfactual
  ns[TREATED < 0, INFECTED_NON_SUPPRESSED := INFECTED_NON_SUPPRESSED.FACTUAL]
  ns[TREATED < 0, PROP_UNSUPPRESSED_PARTICIPANTS_M.COUNTERFACTUAL :=  PROP_UNSUPPRESSED_PARTICIPANTS_M]
  ns[TREATED < 0, PROP_UNSUPPRESSED_NONPARTICIPANTS_M.COUNTERFACTUAL :=  PROP_UNSUPPRESSED_NONPARTICIPANTS_M]
  ns[TREATED < 0, TREATED := 0]
  ns[TREATED < 0, TREATED_CL := 0]
  ns[TREATED < 0, TREATED_CU := 0]
  
  return(ns)
}

