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
                                         log_offset_formula = 'LOG_OFFSET', per_unsuppressed = F){
  
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
    tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
    tmp1 <- merge(tmp1, df_age_aggregated, by = c('AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE'))
  }
  
  if('INDEX_ROUND' %in% names(tmp1)){
    tmp <- merge(df_round, df_community, by = c('COMM'))
    tmp1 <- merge(tmp1, tmp, by = c('INDEX_ROUND', 'INDEX_COMMUNITY'))
  }
  
  if(!is.null(log_offset_round)){
    tmp1 <- merge(tmp1, log_offset_round, by = c('ROUND', 'AGE_INFECTION.RECIPIENT', 'AGE_TRANSMISSION.SOURCE', 'INDEX_DIRECTION', 'INDEX_COMMUNITY'))
    tmp1[, value := value + eval(rlang::parse_expr(log_offset_formula))]
  }
  
  if(!is.null(transform)){
    tmp1[, value := sapply(value, transform)]
  }
  
  #  sum force of infection
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', vars)]
  if(!is.null(operation)){
    tmp1 <- tmp1[, list(value = sapply(value, operation)), by = c('iterations', vars)]
  }
  
  # standardised
  if(!is.null(standardised.vars)){
    tmp1[, total_value := sum(value), by = c('iterations', standardised.vars)]
    tmp1[, value := value / total_value]
  }
  
  # divide by the number of unsuppressed
  if(per_unsuppressed){
    tmp <- copy(eligible_count_round)
    if('AGE_TRANSMISSION.SOURCE' %in% vars)
      setnames(tmp, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
    if('INDEX_DIRECTION' %in% vars)
      tmp[, INDEX_DIRECTION := ifelse(SEX == 'M', df_direction[IS_MF == 1, INDEX_DIRECTION], df_direction[IS_MF == 0, INDEX_DIRECTION])]
    if('INDEX_COMMUNITY' %in% vars)
      tmp <- merge(tmp, df_community, by = 'COMM')
    
    tmp <- tmp[,list(TOTAL_INFECTED_NON_SUPPRESSED = sum(INFECTED_NON_SUPPRESSED)), by = vars]
    
    tmp1 <- merge(tmp, tmp1, by = vars)
    tmp1[, value := value / TOTAL_INFECTED_NON_SUPPRESSED]
  }
  
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
    tmp1 <- merge(tmp1, df_period, by = 'INDEX_TIME')
  
  file = paste0(outdir.table, '-output-', output, 'by_', tolower(paste0(gsub('INDEX_', '', vars), collapse = '_')))
  if(!is.null(standardised.vars)){
    file = paste0(file, 'standardisedby_', tolower(paste0(gsub('INDEX_', '', standardised.vars), collapse = '_')))
  }
  file = paste0(file, '.rds')
  saveRDS(tmp1, file)
  
  return(tmp1)
}


find_median_age_source <- function(samples, var){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[[var]]) )
  setnames(tmp1, 2:5, c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'INDEX_AGE'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'))
  tmp1[, delta := value / total_value]
  tmp1 <- tmp1[, list(value = matrixStats::weightedMedian(AGE_TRANSMISSION.SOURCE, delta )), by = c('iterations', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT')]
  
  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT')]	
  tmp1 = dcast(tmp1, INDEX_DIRECTION + INDEX_COMMUNITY + INDEX_ROUND + AGE_INFECTION.RECIPIENT ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  tmp1 <- merge(tmp1, df_community, by = 'INDEX_COMMUNITY')
  tmp1 <- merge(tmp1, df_round, by = 'INDEX_ROUND')
  
  return(tmp1)
}

prepare_count_data <- function(stan_data){
  
  tmp <- as.data.table(reshape2::melt(stan_data[['y']]))
  setnames(tmp, 1:4, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME'))
  tmp <- merge(tmp, df_direction, by = 'INDEX_DIRECTION')
  tmp <- merge(tmp, df_community, by = 'INDEX_COMMUNITY')
  tmp <- merge(tmp, df_period, by = c('INDEX_TIME', 'COMM'))
  tmp <- merge(tmp, df_age, by = 'INDEX_AGE')
  
  setnames(tmp, 'value', 'count')
  
  return(as.data.table(tmp))
}

prepare_incidence_cases <- function(incidence_cases){
  tmp <- copy(incidence_cases)
  setnames(tmp, 'AGEYRS', 'AGE_INFECTION.RECIPIENT')
  tmp[, IS_MF := as.numeric(SEX == 'F')]
  tmp <- merge(tmp, df_direction, by = 'IS_MF')
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp
}

prepare_susceptible_count <- function(eligible_count){

  tmp <- copy(eligible_count)
  if(!'SUSCEPTIBLE' %in%names(eligible_count)){
    tmp <- tmp[variable =='SUSCEPTIBLE']
    setnames(tmp, c('count'), c('SUSCEPTIBLE'))
  }

  setnames(tmp, c( 'AGEYRS'), c( 'AGE_INFECTION.RECIPIENT'))
  tmp[, IS_MF := as.numeric(SEX == 'F')]
  tmp <- merge(tmp, df_direction, by = 'IS_MF')
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp
}

prepare_eligible_proportion <- function(eligible_count, vars, standardised.vars){
  tmp1 <- eligible_count[variable == 'ELIGIBLE', list(count = sum(count)), by = vars]
  tmp1[, M := count / sum(count), by = standardised.vars]
  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  tmp1 <- merge(tmp1, df_period, by = 'INDEX_TIME')
  tmp1[, type := 'Share in the census eligible individuals']
}

prepare_unsuppressed_proportion <- function(eligible_count, vars, standardised.vars){
  tmp1 <- eligible_count[variable == 'INFECTED_NON_SUPPRESSED', list(count = sum(count)), by = vars]
  tmp1[, M := count / sum(count), by = standardised.vars]
  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  tmp1 <- merge(tmp1, df_period, by = 'INDEX_TIME')
  tmp1[, type := 'Share in the unsuppressed HIV + census eligible individuals']
}

prepare_unsuppressed_proportion_by_round <- function(file.unsuppressed.share, vars){
  
  tmp1 <- as.data.table(read.csv(file.unsuppressed.share))
  
  if(all(c('SEX', 'AGEYRS') %in%vars)){
    tmp1 <- tmp1[, .(ROUND, COMM, SEX, AGEYRS, UNSUPPRESSED_SHARE_AGE_AND_SEX_CL, UNSUPPRESSED_SHARE_AGE_AND_SEX_CU, UNSUPPRESSED_SHARE_AGE_AND_SEX_M)]
    setnames(tmp1, c('UNSUPPRESSED_SHARE_AGE_AND_SEX_CL', 'UNSUPPRESSED_SHARE_AGE_AND_SEX_CU', 'UNSUPPRESSED_SHARE_AGE_AND_SEX_M'), c('CL', "CU", 'M'))
  }else if(vars == 'SEX'){
    tmp1 <- unique(tmp1[, .(ROUND, COMM, SEX, UNSUPPRESSED_SHARE_SEX_CL, UNSUPPRESSED_SHARE_SEX_CU, UNSUPPRESSED_SHARE_SEX_M)])
    setnames(tmp1, c('UNSUPPRESSED_SHARE_SEX_CL', 'UNSUPPRESSED_SHARE_SEX_CU', 'UNSUPPRESSED_SHARE_SEX_M'), c('CL', "CU", 'M'))
  }else{
    stop()
  }

  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')

  # set round 14 to be the same as round 15 in inland
  tmp <- tmp1[ROUND == 15 & COMM == 'inland']
  tmp[, ROUND := 14]
  tmp1 <- rbind(tmp, tmp1)
  
  # set round 15S to be the same as round 15 in fishing
  tmp1[, ROUND := as.numeric(ROUND)]
  tmp <- tmp1[ROUND == 15 & COMM == 'fishing']
  tmp[, ROUND := 15.1]
  tmp1 <- rbind(tmp, tmp1)
  
  # merge to index round
  tmp1 <- merge(tmp1, df_round[, .(COMM,round, INDEX_ROUND, LABEL_ROUND)], by.x = c('COMM', 'ROUND'), by.y = c('COMM', 'round'))
  
  # type
  tmp1[, type := 'Share in the HIV+ unsuppressed census eligible individuals']
  
  return(tmp1)
}

prepare_prevalence_proportion_by_round <- function(file.prevalence.share, vars){

  tmp1 <- as.data.table(read.csv(file.prevalence.share))
  
  if(all(c('SEX', 'AGEYRS') %in%vars)){
    tmp1 <- tmp1[, .(ROUND, COMM, SEX, AGEYRS, PREVALENCE_SHARE_AGE_AND_SEX_CL, PREVALENCE_SHARE_AGE_AND_SEX_CU, PREVALENCE_SHARE_AGE_AND_SEX_M)]
    setnames(tmp1, c('PREVALENCE_SHARE_AGE_AND_SEX_CL', 'PREVALENCE_SHARE_AGE_AND_SEX_CU', 'PREVALENCE_SHARE_AGE_AND_SEX_M'), c('CL', "CU", 'M'))
  }else if(vars == 'SEX'){
    tmp1 <- unique(tmp1[, .(ROUND, COMM, SEX, PREVALENCE_SHARE_SEX_CL, PREVALENCE_SHARE_SEX_CU, PREVALENCE_SHARE_SEX_M)])
    setnames(tmp1, c('PREVALENCE_SHARE_SEX_CL', 'PREVALENCE_SHARE_SEX_CU', 'PREVALENCE_SHARE_SEX_M'), c('CL', "CU", 'M'))
  }else{
    stop()
  }
  
  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  
  # set round 14 to be the same as round 15 in inland
  tmp <- tmp1[ROUND == 15 & COMM == 'inland']
  tmp[, ROUND := 14]
  tmp1 <- rbind(tmp, tmp1)
  
  # find round label
  tmp1[, ROUND := paste0('R0', ROUND)]
  
  # merge to index round
  tmp1 <- merge(tmp1, df_round[, .(COMM, ROUND, INDEX_ROUND)], by = c('COMM', 'ROUND'))
  
  # type
  tmp1[, type := 'Share in the HIV+ census eligible individuals']
  
  return(tmp1)
}


prepare_prevalence_by_round <- function(eligible_count_round, vars, standardised.vars){
  tmp1 <- eligible_count_round[, list(M = sum(INFECTED) / sum(ELIGIBLE)), by = vars]

  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  
  tmp1[, type := 'Share in the HIV+ census eligible individuals']
}


remove_first_round <- function(tmp){
  
  rm.vars <- c('round', 'ROUND', 'ROUND_SPANYRS', 'INDEX_TIME', 'min_sample_date', 'max_sample_date')
  
  for(var in rm.vars){
    if(var %in% names(tmp)){
      tmp <- select(tmp, -var)
      
    }
  }
  
  tmp[, INDEX_ROUND:= INDEX_ROUND +1]
  tmp <- merge(tmp, df_round, by = 'INDEX_ROUND')
  
  return(tmp)
}
