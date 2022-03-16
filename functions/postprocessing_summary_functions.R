make_convergence_diagnostics_stats = function(fit, outdir)
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
  re = rstan::extract(fit)
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

summarise_var_by_age_group <- function(samples, var, df_group, df_age, transform = NULL){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[[var]]) )
  setnames(tmp1, 2:3, c('index_group','index_age'))
  
  if(!is.null(transform)){
    tmp1[, value := sapply(value, transform)]
  }
  
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group', 'index_age')]	
  tmp1 = dcast(tmp1, index_group + index_age ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  
  return(tmp1)
}

find_relative_intensity_PP <- function(samples, df_group, df_age){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group'))
  tmp1[, value := value / total_value]

  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group', 'index_age')]	
  tmp1 = dcast(tmp1, index_group + index_age ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  
  return(tmp1)
}

find_relative_intensity_PP_aggregated <- function(samples, df_group, df_age, df_age_aggregated){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'is_before_cutoff_date')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'is_before_cutoff_date'))
  tmp1[, value := value / total_value]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  tmp1 <- merge(tmp1, df_age_aggregated, by = c('age_infection.RECIPIENT', 'age_transmission.SOURCE'))
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'index_group', 'age_group_transmission.SOURCE', 'age_group_infection.RECIPIENT')]
  
  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group', 'age_group_transmission.SOURCE', 'age_group_infection.RECIPIENT')]	
  tmp1 = dcast(tmp1, index_group + age_group_transmission.SOURCE + age_group_infection.RECIPIENT ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')

  return(tmp1)
}

find_age_source_by_age_group <- function(samples, df_group, df_age){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')

  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_infection.RECIPIENT'))
  tmp1[, delta := value / total_value]
  tmp1 <- tmp1[, list(value = matrixStats::weightedMedian(age_transmission.SOURCE, delta )), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]

  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group', 'age_infection.RECIPIENT')]	
  tmp1 = dcast(tmp1, index_group + age_infection.RECIPIENT ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

find_age_source_difference_by_age_group <- function(samples, df_group, df_age){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_infection.RECIPIENT'))
  tmp1[, pi := value / total_value]
  tmp1 <- tmp1[, list(value = matrixStats::weightedMedian(age_transmission.SOURCE, pi )), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  
  tmp1[, value := value - age_infection.RECIPIENT]
  
  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group', 'age_infection.RECIPIENT')]	
  tmp1 = dcast(tmp1, index_group + age_infection.RECIPIENT ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

find_age_source_by_group <- function(samples, df_group, df_age, incidence, range_age_observed){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  
  tmp1 <- merge(tmp1, range_age_observed, by = 'index_group')
  tmp1 <- tmp1[age_infection.RECIPIENT >= min_age_infection.RECIPIENT & age_infection.RECIPIENT <= max_age_infection.RECIPIENT]
  tmp1 <- tmp1[age_transmission.SOURCE >= min_age_transmission.SOURCE & age_transmission.SOURCE <= max_age_transmission.SOURCE]
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_infection.RECIPIENT'))
  tmp1[, delta := value / total_value]
  tmp1 <- tmp1[, list(value = matrixStats::weightedMedian(age_transmission.SOURCE, delta )), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  
  # weight recipient by incidence
  di <- as.data.table(incidence)
  di <- di[ROUND == 16]
  di[, is_mf := ifelse(SEX == 'F', 1, 0)]
  di <- merge(di, range_age_observed, by = c('is_mf'), allow.cartesian=TRUE)
  di <- di[AGEYRS >= min_age_infection.RECIPIENT & AGEYRS <= max_age_infection.RECIPIENT]
  
  tmp2 <- di[, list(total_incidence = sum(INCIDENCE)), by = c('index_group')]
  tmp2 <- merge(di, tmp2, by = 'index_group')
  tmp2[, incidence_weight := INCIDENCE / total_incidence]
  setnames(tmp2, 'AGEYRS', 'age_infection.RECIPIENT')
  
  tmp1 <- merge(tmp1, tmp2, by = c('age_infection.RECIPIENT', 'index_group'))
  tmp1 <- tmp1[, list(value = sum(incidence_weight * value)), by = c('iterations', 'index_group')]
  
  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group')]	
  tmp1 = dcast(tmp1, index_group ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

find_age_recipient_by_age_group <- function(samples, df_group, df_age){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_transmission.SOURCE')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_transmission.SOURCE'))
  tmp1[, pi := value / total_value]
  tmp1 <- tmp1[, list(value = matrixStats::weightedMedian(age_infection.RECIPIENT, pi )), by = c('iterations', 'index_group', 'age_transmission.SOURCE')]
  
  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group', 'age_transmission.SOURCE')]	
  tmp1 = dcast(tmp1, index_group + age_transmission.SOURCE ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

find_age_recipient_difference_by_age_group <- function(samples, df_group, df_age){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_transmission.SOURCE')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_transmission.SOURCE'))
  tmp1[, pi := value / total_value]
  tmp1 <- tmp1[, list(value = matrixStats::weightedMedian(age_infection.RECIPIENT, pi )), by = c('iterations', 'index_group', 'age_transmission.SOURCE')]
  
  tmp1[, value := value - age_transmission.SOURCE]
  
  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group', 'age_transmission.SOURCE')]	
  tmp1 = dcast(tmp1, index_group + age_transmission.SOURCE ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

find_age_recipient_by_group <- function(samples, df_group, df_age, incidence, range_age_observed){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  tmp1 <- merge(tmp1, range_age_observed, by = 'index_group')
  tmp1 <- tmp1[age_infection.RECIPIENT >= min_age_infection.RECIPIENT & age_infection.RECIPIENT <= max_age_infection.RECIPIENT]
  tmp1 <- tmp1[age_transmission.SOURCE >= min_age_transmission.SOURCE & age_transmission.SOURCE <= max_age_transmission.SOURCE]
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_transmission.SOURCE')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_transmission.SOURCE'))
  tmp1[, pi := value / total_value]
  tmp1 <- tmp1[, list(value = sum(pi * age_infection.RECIPIENT)), by = c('iterations', 'index_group', 'age_transmission.SOURCE')]
  
  # weight recipient by incidence
  di <- as.data.table(incidence)
  di <- di[ROUND == 16]
  di[, is_mf := ifelse(SEX == 'M', 1, 0)]
  di <- merge(di, range_age_observed, by = c('is_mf'), allow.cartesian=TRUE)
  di <- di[AGEYRS >= min_age_infection.RECIPIENT & AGEYRS <= max_age_infection.RECIPIENT]
  
  tmp2 <- di[, list(total_incidence = sum(INCIDENCE)), by = c('SEX')]
  tmp2 <- merge(di, tmp2, by = 'SEX')
  tmp2[, incidence_weight := INCIDENCE / total_incidence]
  setnames(tmp2, 'AGEYRS', 'age_transmission.SOURCE')
  
  tmp1 <- merge(tmp1, tmp2, by = c('age_transmission.SOURCE', 'index_group'))
  tmp1 <- tmp1[, list(value = sum(incidence_weight * value)), by = c('iterations', 'index_group')]
  
  tmp1 = tmp1[, list(q= quantile(na.omit(value), prob=ps, na.rm = T), q_label=p_labs), 
              by=c('index_group')]	
  tmp1 = dcast(tmp1, index_group ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

find_incident_cases <- function(samples, df_group, df_age, incidence){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  tmp1 <- merge(tmp1, df_group, by = 'index_group')

  # find incident cases 
  di <- as.data.table(merge(incidence, df_round, by = 'ROUND'))
  di[, is_before_cutoff_date := ifelse(DATE_ROUND < cutoff_date, 1, 0)]
  di[, is_mf := ifelse(SEX == 'F', 1, 0)]
  setnames(di, 'AGEYRS', 'age_infection.RECIPIENT')
  di <- di[, list(INFECTIONS = runif(1, median(INFECTIONS_LB), median(INFECTIONS_UB))), by = c("is_mf", 'is_before_cutoff_date', 'age_infection.RECIPIENT')]
  tmp1 <- merge(di, tmp1, by = c('is_mf', 'is_before_cutoff_date', 'age_infection.RECIPIENT'), allow.cartesian=TRUE)

  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_infection.RECIPIENT'))
  tmp1[, pi := value / total_value]
  
  tmp1 <- tmp1[, value := INFECTIONS * pi]
  
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('index_group', 'index_age')]	
  tmp1 = dcast(tmp1, index_group + index_age ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  
  return(tmp1)
}

find_relative_incident_cases_aggregated <- function(samples, df_group, df_age, df_age_aggregated, incidence){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  # find incident cases 
  di <- as.data.table(merge(incidence, df_round, by = 'ROUND'))
  di[, is_before_cutoff_date := ifelse(DATE_ROUND < cutoff_date, 1, 0)]
  di[, is_mf := ifelse(SEX == 'F', 1, 0)]
  setnames(di, 'AGEYRS', 'age_infection.RECIPIENT')
  di <- di[, list(INFECTIONS = runif(1, median(INFECTIONS_LB), median(INFECTIONS_UB))), by = c("is_mf", 'is_before_cutoff_date', 'age_infection.RECIPIENT')]
  tmp1 <- merge(di, tmp1, by = c('is_mf', 'is_before_cutoff_date', 'age_infection.RECIPIENT'), allow.cartesian=TRUE)
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_infection.RECIPIENT'))
  tmp1[, delta := value / total_value]
  tmp1[, value := INFECTIONS * delta]
  tmp1 <- select(tmp1, - total_value)
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'is_before_cutoff_date')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'is_before_cutoff_date'))
  tmp1[, value := value / total_value]
  
  tmp1 <- merge(tmp1, df_age_aggregated, by = c('age_infection.RECIPIENT', 'age_transmission.SOURCE'))
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'index_group', 'age_group_transmission.SOURCE', 'age_group_infection.RECIPIENT')]

  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('index_group', 'age_group_transmission.SOURCE', 'age_group_infection.RECIPIENT')]	
  tmp1 = dcast(tmp1, index_group + age_group_transmission.SOURCE + age_group_infection.RECIPIENT ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}


find_incident_cases_by_age_source_group <- function(samples, df_group, df_age, incidence){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  tmp1 <- merge(tmp1, df_group, by = 'index_group')

  # find incident cases 
  di <- as.data.table(incidence)
  di <- di[ROUND == 16]
  di[, is_mf := ifelse(SEX == 'F', 1, 0)]
  setnames(di, 'AGEYRS', 'age_infection.RECIPIENT')
  tmp1 <- merge(di, tmp1, by = c('is_mf', 'age_infection.RECIPIENT'))
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_infection.RECIPIENT')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_infection.RECIPIENT'))
  tmp1[, pi := value / total_value]
  
  tmp1 <- tmp1[, value := INFECTIONS * pi]
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'index_group', 'age_transmission.SOURCE')]
  
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('index_group', 'age_transmission.SOURCE')]	
  tmp1 = dcast(tmp1, index_group + age_transmission.SOURCE ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

prepare_count_data <- function(stan_data, df_age, df_group){

  count_data <- vector(mode = 'list', length = ncol(stan_data$y))
  
  for(i in 1:ncol(stan_data$y)){
    count_data[[i]] <- data.table(count =  stan_data$y[,i])
    count_data[[i]][, index_group := i]
    count_data[[i]][, index_age := 1:nrow(count_data[[i]])]
  }
  
  count_data <- do.call('rbind', count_data)

  count_data <- merge(count_data, df_age, by = 'index_age')
  count_data <- merge(count_data, df_group, by = 'index_group')
  
  return(as.data.table(count_data))
}

find_range_age_observed <- function(tmp, df_group){
  
  tmp[, is_mf := 0]
  tmp[sex.SOURCE == 'M' & sex.RECIPIENT == 'F', is_mf := 1]
  tmp <- tmp[, list(min_age_infection.RECIPIENT = min(age_infection.RECIPIENT), 
                    max_age_infection.RECIPIENT = max(age_infection.RECIPIENT),
                    min_age_transmission.SOURCE = min(age_transmission.SOURCE), 
                    max_age_transmission.SOURCE = max(age_transmission.SOURCE)), 
             by = c('is_mf', 'date_infection_before_cutoff.RECIPIENT')]
  tmp <- tmp[, list(min_age_infection.RECIPIENT = max(min_age_infection.RECIPIENT), 
                    max_age_infection.RECIPIENT = min(max_age_infection.RECIPIENT),
                    min_age_transmission.SOURCE = max(min_age_transmission.SOURCE), 
                    max_age_transmission.SOURCE = min(max_age_transmission.SOURCE)), 
             by = c('is_mf')]
  tmp <- merge(tmp, df_group, by = 'is_mf')
}

find_sex_source <- function(samples, df_group){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  # sum across ages
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'index_group', 'index_time')]
  
  # take proportion across sex
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_time')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_time'))
  tmp1[, value := value / total_value]
  
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('index_group')]	
  tmp1 = dcast(tmp1, index_group ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}


find_sex_source_standardised <- function(samples, df_group, df_age, incidence){
  
  ps <- c(0.5, 0.025, 0.975)
  p_labs <- c('M','CL','CU')
  
  tmp1 = as.data.table( reshape2::melt(samples[['log_lambda']]) )
  setnames(tmp1, 2:3, c('index_group', 'index_age'))
  
  tmp1[, value := exp(value)]
  
  tmp1 <- merge(tmp1, df_age, by = 'index_age')
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  # find incident cases 
  di <- as.data.table(merge(incidence, df_round, by = 'ROUND'))
  di[, is_before_cutoff_date := ifelse(DATE_ROUND < cutoff_date, 1, 0)]
  di[, is_mf := ifelse(SEX == 'F', 1, 0)]
  setnames(di, 'AGEYRS', 'age_infection.RECIPIENT')
  di <- di[, list(INFECTIONS = runif(10, median(INFECTIONS_LB), median(INFECTIONS_UB)),
                  idx_draw = 1:10), by = c("is_mf", 'is_before_cutoff_date', 'age_infection.RECIPIENT')]
  tmp1 <- merge(di, tmp1, by = c('is_mf', 'is_before_cutoff_date', 'age_infection.RECIPIENT'), allow.cartesian=TRUE)
  
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_group', 'age_infection.RECIPIENT', 'idx_draw')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_group', 'age_infection.RECIPIENT', 'idx_draw'))
  tmp1[, delta := value / total_value]
  tmp1[, value := INFECTIONS * delta]
  tmp1 <- select(tmp1, - total_value)
  
  # sum across ages
  tmp1 <- tmp1[, list(value = sum(value)), by = c('iterations', 'index_group', 'index_time', 'idx_draw')]

  # take proportion across sex
  tmp2 <- tmp1[, list(total_value = sum(value)), by = c('iterations', 'index_time', 'idx_draw')]
  tmp1 <- merge(tmp1, tmp2, by = c('iterations', 'index_time', 'idx_draw'))
  tmp1[, value := value / total_value]
  
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('index_group')]	
  tmp1 = dcast(tmp1, index_group ~ q_label, value.var = "q")
  
  tmp1 <- merge(tmp1, df_group, by = 'index_group')
  
  return(tmp1)
}

