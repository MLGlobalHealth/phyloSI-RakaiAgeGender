add_stan_data_base <- function(){
  
  stan_data = list()
  
  # number of groups
  stan_data[['N_DIRECTION']] = nrow(df_direction)
  stan_data[['N_COMMUNITY']] = nrow(df_community)
  stan_data[['N_PERIOD']] = max(df_period$INDEX_TIME)
  
  # number of age 
  stan_data[['N_PER_GROUP']] = nrow(df_age)
  
  # number of age group
  stan_data[['N_AGE']] = df_age[, length(unique(AGE_INFECTION.RECIPIENT))]
  
  # number of rounds 
  tmp <- merge(df_round, df_community, by = 'COMM')[order(INDEX_COMMUNITY)]
  stan_data[['N_ROUND']] = tmp[, length(unique(INDEX_ROUND)), by = 'COMM']$V1

  # map from round to period
  stan_data[['map_round_period']] = tmp[COMM ==  df_community[order(INDEX_COMMUNITY), COMM[1]] & order(round), INDEX_TIME]
  
  # number of period
  stan_data[['N_ROUND_PER_PERIOD']] = tmp[order(INDEX_TIME) & COMM == df_community[order(INDEX_COMMUNITY), COMM[1]], length(unique(INDEX_ROUND)), by = 'INDEX_TIME']$V1
  
  return(stan_data)
}

add_phylo_data <- function(stan_data, pairs){
  
  # prepare pairs
  pairs_round <- pairs[, list(AGE_TRANSMISSION.SOURCE = floor(AGE_TRANSMISSION.SOURCE), 
                              AGE_INFECTION.RECIPIENT = floor(AGE_INFECTION.RECIPIENT), 
                              DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT, 
                              COMM.RECIPIENT = COMM.RECIPIENT, 
                              SEX.SOURCE = SEX.SOURCE)]
  
  # save count in each entry
  y = array(NA, c(stan_data[['N_PER_GROUP']], stan_data[['N_DIRECTION']], stan_data[['N_PERIOD']]))
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:stan_data[['N_PERIOD']]){
      
      # direction group
      .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
      tmp <- pairs_round[SEX.SOURCE == .SEX.SOURCE]
      
      # time group
      tmp <- tmp[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == df_period[INDEX_TIME == k, BEFORE_CUTOFF]]
      
      # count number of observation
      tmp <- tmp[, list(count = .N), by = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT')]
      tmp <- merge(df_age, tmp, 
                   by = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), all.x = T)
      tmp[is.na(count), count := 0]
      
      setkey(tmp, AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)
      
      tmp1 <- pairs[SEX.SOURCE == .SEX.SOURCE  & DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == df_period[ INDEX_TIME == k, BEFORE_CUTOFF]
                    & AGE_INFECTION.RECIPIENT >= 15 & AGE_INFECTION.RECIPIENT < 50 & AGE_TRANSMISSION.SOURCE >= 15 & AGE_TRANSMISSION.SOURCE < 50]
      stopifnot(sum(tmp$count) == nrow(tmp1))
      
      # check the order of ages is correct
      tmp <- tmp[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
      stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == tmp[, AGE_INFECTION.RECIPIENT])
      stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == tmp[, AGE_TRANSMISSION.SOURCE])
      
      cat(nrow(tmp1), 'pairs with infection', df_direction[i, LABEL_DIRECTION], 'in', df_period[INDEX_TIME == k, PERIOD], '\n')
      
      y[, i, k] = matrix(tmp$count, ncol = 1)
      
      
    }
  }
  
  
  # save stan data
  stan_data[['y']] = y
  
  # add map age source and recipient
  df_age <- df_age[order(INDEX_AGE)]
  stan_data[['map_age_source']] = df_age[, AGE_TRANSMISSION.SOURCE - min(AGE_TRANSMISSION.SOURCE) + 1]
  stan_data[['map_age_recipient']] = df_age[, AGE_INFECTION.RECIPIENT - min(AGE_INFECTION.RECIPIENT) + 1]
  
  return(stan_data)
}

add_incidence_cases <- function(stan_data, incidence_cases_round){
  
  # save count in each entry
  z = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:max(stan_data[['N_ROUND']])){
      
      .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .ROUND = df_round[INDEX_ROUND == k, ROUND]
      
      # direction group
      tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
      
      # round
      tmp <- tmp[ROUND == .ROUND]
      
      # order by age
      tmp <- tmp[order(AGEYRS)] 
      
      # sanity check
      tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT & ROUND == .ROUND]
      stopifnot(sum(tmp$INCIDENT_CASES) == sum(tmp1$INCIDENT_CASES))
      cat(sum(tmp1$INCIDENT_CASES), 'incidence cases ', df_direction[i, LABEL_DIRECTION], 'in', .ROUND, '\n')
      
      # fill
      z[, i, k] = round(tmp$INCIDENT_CASES)
      
    }
  }
  
  stan_data[['z']] = z
  
  
  return(stan_data)
  
}

add_incidence_rates <- function(stan_data, incidence_cases_round){
  
  ir = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:max(stan_data[['N_ROUND']])){
      
      .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .ROUND = df_round[INDEX_ROUND == k, ROUND]
      
      # direction group
      tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
      
      # round
      tmp <- tmp[ROUND == .ROUND]
      
      # order by age
      tmp <- tmp[order(AGEYRS)] 
      
      # sanity check
      tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT & ROUND == .ROUND]
      stopifnot(sum(tmp$INCIDENCE) == sum(tmp1$INCIDENCE))
      cat(mean(tmp1$INCIDENCE), 'incidence rates ', df_direction[i, LABEL_DIRECTION], 'in', .ROUND, '\n')
      
      # fill
      ir[, i, k] = tmp$INCIDENCE
      
    }
  }
  
  stan_data[['ir']] = ir
  
  
  return(stan_data)
  
}

add_incidence_rates_lognormal_parameters <- function(stan_data, incidence_cases_round){
  
  ir_lognorm_mean = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], max(stan_data[['N_ROUND']])))
  ir_lognorm_sd = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:max(stan_data[['N_ROUND']])){
      
      .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .ROUND = df_round[INDEX_ROUND == k, ROUND]
      
      # direction group
      tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
      
      # round
      tmp <- tmp[ROUND == .ROUND]
      
      # order by age
      tmp <- tmp[order(AGEYRS)] 
      
      # sanity check
      tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT & ROUND == .ROUND]
      stopifnot(sum(tmp$INCIDENCE) == sum(tmp1$INCIDENCE))

      # fill
      lognorm_parms <- lognorm::getParmsLognormForMedianAndUpper(median = tmp$INCIDENCE, upper = tmp$UB, sigmaFac=2)
      ir_lognorm_mean[, i, k] = lognorm_parms[, 1]
      ir_lognorm_sd[, i, k] = lognorm_parms[, 2]
    }
  }
  
  stan_data[['ir_lognorm_mean']] = ir_lognorm_mean
  stan_data[['ir_lognorm_sd']] = ir_lognorm_sd
  
  return(stan_data)
  
}

add_offset <- function(stan_data, eligible_count, df_estimated_contact_rates, 
                       use_number_susceptible_offset, use_contact_rates_prior){
  
  # add offset including the number of susceptible, the number of infected unsuppressed and potentially the contact rates
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  eligible_count_wide[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  
  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:max(stan_data[['N_ROUND']])){
      
      log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
      
      .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
      .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .ROUND <- df_round[INDEX_ROUND == k, ROUND]
      
      # add # susceptible and proportion of susceptible in recipient
      tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & ROUND == .ROUND]
      log_offset <- merge(log_offset, tmp[, .(AGEYRS, PROP_SUSCEPTIBLE, SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
      
      # add number of infected unsuppressed in source
      tmp <- eligible_count_wide[SEX == .SEX.SOURCE & ROUND == .ROUND]
      log_offset <- merge(log_offset, tmp[, .(AGEYRS, INFECTED_NON_SUPPRESSED)], by.x = 'AGE_TRANSMISSION.SOURCE', by.y = 'AGEYRS')
      
      # add contact rates
      tmp <- df_estimated_contact_rates[part.sex == .SEX.SOURCE]
      colnames(tmp) <- toupper(colnames(tmp))
      log_offset <- merge(log_offset, tmp[, .(PART.AGE, CONT.AGE, CNTCT.RATE)], 
                          by.x = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), by.y = c('PART.AGE', 'CONT.AGE'))
      
      # make log offset
      if(use_number_susceptible_offset){
        log_offset[, LOG_OFFSET := log(SUSCEPTIBLE) + log(INFECTED_NON_SUPPRESSED)]
      }else{
        log_offset[, LOG_OFFSET := log(PROP_SUSCEPTIBLE) + log(INFECTED_NON_SUPPRESSED)]
      }
      
      # potentially add contact rates
      if(use_contact_rates_prior){
        log_offset[, LOG_OFFSET := LOG_OFFSET + log(CNTCT.RATE)]
      }
      
      # check the order of ages is correct
      log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
      stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
      stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
      
      # add to array
      log_offset_array[i, k,] = log_offset[, LOG_OFFSET]
      
    }
  }
  
  stan_data[['log_offset']] = log_offset_array
  
  return(stan_data)
}

add_offset_time <- function(stan_data, df_round){
  
  # add offset including the time in years of each rounds
  
  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:max(stan_data[['N_ROUND']])){
      
      log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
      
      .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
      .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .ROUND <- df_round[INDEX_ROUND == k, ROUND]
      
      # add period in year
      log_offset[, PERIOD_SPAN := df_round[ROUND == .ROUND, ROUND_SPANYRS]]
      
      # make log offset time
      log_offset[, LOG_OFFSET_TIME :=  log(PERIOD_SPAN)]
      
      # check the order of ages is correct
      log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
      stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
      stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
      
      # add to array
      log_offset_array[i, k,] = log_offset[, LOG_OFFSET_TIME]
      
    }
  }
  
  stan_data[['log_offset_time']] = log_offset_array
  
  return(stan_data)
}

add_offset_susceptible <- function(stan_data, eligible_count){
  
  # add offset including the number of susceptible
  # that can be used to calculate incidence rates
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  
  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:max(stan_data[['N_ROUND']])){
      
      log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
      
      .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
      .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .ROUND <- df_round[INDEX_ROUND == k, ROUND]
      
      # add number of susceptible in recipient
      tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & ROUND == .ROUND]
      log_offset <- merge(log_offset, tmp[, .(AGEYRS, SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
      
      # make log offset
      log_offset[, LOG_OFFSET_SUSCEPTIBLE := log(SUSCEPTIBLE) ]
      
      # check the order of ages is correct
      log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
      stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
      stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
      
      # add to array
      log_offset_array[i, k,] = log_offset[, LOG_OFFSET_SUSCEPTIBLE]
      
    }
  }
  
  stan_data[['log_offset_susceptible']] = log_offset_array
  
  return(stan_data)
}

add_probability_sampling <- function(stan_data, proportion_sampling){
  
  # add the probability of sampling a infection event
  
  proportion_sampling <- proportion_sampling[order(SEX.RECIPIENT, COMM, BEFORE_CUTOFF, PERIOD)]
  
  log_prop_sampling_array =array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_PERIOD']], stan_data[['N_PER_GROUP']])))
  sampling_index=array(NA, c(c(stan_data[['N_PER_GROUP']], stan_data[['N_DIRECTION']], stan_data[['N_PERIOD']])))
  n_sampling_index=array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_PERIOD']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:stan_data[['N_PERIOD']]){
      
      log_prop_sampling = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
      
      .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
      .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .BEFORE_CUTOFF <- df_period[k, BEFORE_CUTOFF]
      
      # add probability of sampling recipient 
      tmp <- proportion_sampling[SEX.RECIPIENT == .SEX.RECIPIENT & BEFORE_CUTOFF == .BEFORE_CUTOFF]
      log_prop_sampling <- merge(log_prop_sampling, tmp[, .(AGEYRS.SOURCE, AGEYRS.RECIPIENT, prop_sampling)], 
                                 by.x = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), 
                                 by.y = c('AGEYRS.SOURCE', 'AGEYRS.RECIPIENT'))
      
      # make log prob sampling
      if(1){ # we will not use those entries in the likelihood anyway
        log_prop_sampling[prop_sampling == 0, prop_sampling := 0.0001]
      }
      log_prop_sampling[, LOG_PROP_SAMPLING := log(prop_sampling)]
      
      # check the order of ages is correct
      log_prop_sampling <- log_prop_sampling[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
      stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_prop_sampling[, AGE_INFECTION.RECIPIENT])
      stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_prop_sampling[, AGE_TRANSMISSION.SOURCE])
      
      # prop sampling
      log_prop_sampling_array[i, k,] = log_prop_sampling[, LOG_PROP_SAMPLING]
      
      # was the recipient sampled
      n_sampling_index[i, k] = log_prop_sampling[, sum(prop_sampling != 0.0001)]
      sampling_index[,i, k] = rep(-1, nrow(log_prop_sampling))
      
      # if sampled add the probabilties (for some sensitivity analysis we removed some phylo pairs resulting in no pairs sampleds)
      if(n_sampling_index[i, k] != 0){
        sampling_index[1:n_sampling_index[i, k], i, k] = log_prop_sampling[, which(prop_sampling != 0.0001)]
      }
    }
  }
  
  stan_data[['log_prop_sampling']] = log_prop_sampling_array
  stan_data[['n_sampling_index_y']] = n_sampling_index
  stan_data[['sampling_index_y']] = sampling_index
  stan_data[['N_OBS']] = sum(n_sampling_index)
  
  return(stan_data)
}

add_probability_detection <- function(stan_data, proportion_sampling){
  
  # add the probability of detecting a infection event
  
  pps <- unique(proportion_sampling[, .(SEX.RECIPIENT, COMM, BEFORE_CUTOFF, PERIOD, AGEYRS.RECIPIENT, prop_sampling)])
  pps <- pps[order(SEX.RECIPIENT, COMM, BEFORE_CUTOFF, PERIOD)]
  
  log_prop_sampling_array =array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_PERIOD']], stan_data[['N_AGE']])))

  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:stan_data[['N_PERIOD']]){
      
      log_prop_sampling = unique(df_age[, .(AGE_INFECTION.RECIPIENT)])
      
      .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
      .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .BEFORE_CUTOFF <- df_period[k, BEFORE_CUTOFF]
      
      # add probability of sampling recipient 
      tmp <- pps[SEX.RECIPIENT == .SEX.RECIPIENT & BEFORE_CUTOFF == .BEFORE_CUTOFF]
      log_prop_sampling <- merge(log_prop_sampling, tmp[, .(AGEYRS.RECIPIENT, prop_sampling)], 
                                 by.x = c('AGE_INFECTION.RECIPIENT'), 
                                 by.y = c('AGEYRS.RECIPIENT'))
      
      # make log prob sampling
      if(1){ # we will not use those entries in the likelihood anyway
        log_prop_sampling[prop_sampling == 0, prop_sampling := 0.0001]
      }
      log_prop_sampling[, LOG_PROP_SAMPLING := log(prop_sampling)]
      
      # check the order of ages is correct
      log_prop_sampling <- log_prop_sampling[order(AGE_INFECTION.RECIPIENT)]
      stopifnot(unique(df_age[, AGE_INFECTION.RECIPIENT]) == log_prop_sampling[, AGE_INFECTION.RECIPIENT])

      # prop sampling
      log_prop_sampling_array[i, k,] = log_prop_sampling[, LOG_PROP_SAMPLING]

    }
  }
  
  stan_data[['log_prop_detection']] = log_prop_sampling_array
  
  return(stan_data)
}

add_probability_sampling_source <- function(stan_data, proportion_unsuppressed_deepsequenced){
  
  # add the probability of sampling source
  pud <- copy(proportion_unsuppressed_deepsequenced)
  pud <- pud[, list(M = PROP_UNSUP_EVERDEEPSEQ, 
                    AGE_TRANSMISSION.SOURCE = as.numeric(gsub('(.+)-.*', '\\1', AGEGP)):as.numeric(gsub('.*-(.+)', '\\1', AGEGP))), 
             by = c('COMM', 'ROUND', 'SEX', 'AGEGP')]
  setnames(pud, c('SEX'), c('SEX.SOURCE'))
  
  pud <- pud[order(SEX.SOURCE, COMM, ROUND)]
  
  log_prop_sampling_source_array =array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_ROUND']], stan_data[['N_AGE']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(k in 1:stan_data[['N_ROUND']]){
      
      log_prop_sampling_source = unique(df_age[, .(AGE_TRANSMISSION.SOURCE)])
      
      .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
      .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
      .ROUND <- df_round[INDEX_ROUND == k, ROUND]
      
      # add probability of sampling recipient 
      tmp <- pud[SEX.SOURCE == .SEX.SOURCE & ROUND == .ROUND]
      log_prop_sampling_source <- merge(log_prop_sampling_source, tmp[, .(AGE_TRANSMISSION.SOURCE, M)], 
                                        by = c('AGE_TRANSMISSION.SOURCE'))
      
      # sanity check
      stopifnot(nrow(log_prop_sampling_source[M == 0]) == 0)
      
      # take log
      log_prop_sampling_source[, LOG_PROP_SAMPLING_SOURCE := log(M)]
      
      # check the order of ages is correct
      log_prop_sampling_source <- log_prop_sampling_source[order(AGE_TRANSMISSION.SOURCE)]
      stopifnot(unique(df_age[, AGE_TRANSMISSION.SOURCE]) == log_prop_sampling_source[, AGE_TRANSMISSION.SOURCE])
      
      # prop sampling
      log_prop_sampling_source_array[i, k,] = log_prop_sampling_source[, LOG_PROP_SAMPLING_SOURCE]
      
      
    }
  }
  
  stan_data[['log_prop_sampling_source']] = log_prop_sampling_source_array
  
  return(stan_data)
}

add_2D_splines_stan_data = function(stan_data, spline_degree = 3, n_knots_rows = 8, n_knots_columns = 8, X, Y)
{
  
  stan_data$number_rows <- length(X)
  stan_data$number_columns <- length(Y)
  
  knots_rows = X[seq(1, length(X), length.out = n_knots_rows)] 
  knots_columns = Y[seq(1, length(Y), length.out = n_knots_columns)]
  
  stan_data$num_basis_rows = length(knots_rows) + spline_degree - 1
  stan_data$num_basis_columns = length(knots_columns) + spline_degree - 1
  
  stan_data$IDX_BASIS_ROWS = 1:stan_data$num_basis_rows
  stan_data$IDX_BASIS_COLUMNS = 1:stan_data$num_basis_columns
  
  stan_data$BASIS_ROWS = bsplines(X, knots_rows, spline_degree)
  stan_data$BASIS_COLUMNS = bsplines(Y, knots_columns, spline_degree)
  
  stopifnot(all( apply(stan_data$BASIS_ROWS, 1, sum) > 0  ))
  stopifnot(all( apply(stan_data$BASIS_COLUMNS, 1, sum) > 0  ))
  
  return(stan_data)
}

bspline = function(x, k, order, intervals)
{
  
  if(order == 1){
    return(x >= intervals[k] & x < intervals[k+1])
  }
  
  w1 = 0; w2 = 0
  
  if(intervals[k] != intervals[k+order-1])
    w1 = (x - intervals[k]) / (intervals[k+order-1] - intervals[k])
  if(intervals[k+1] != intervals[k+order])
    w2 = 1 - (x - intervals[k+1]) / (intervals[k+order] - intervals[k+1])
  
  spline = w1 * bspline(x, k, order - 1, intervals) +
    w2 * bspline(x, k+1, order - 1, intervals)
  
  return(spline)
}

find_intervals = function(knots, degree, repeating = T)
{
  
  K = length(knots)
  
  intervals = vector(mode = 'double', length = 2*degree + K)
  
  # support of knots
  intervals[(degree+1):(degree+K)] = knots
  
  # extreme
  if(repeating)
  {
    intervals[1:degree] = min(knots)
    intervals[(degree+K+1):(2*degree+K)] = max(knots)
  } else {
    gamma = 0.1
    intervals[1:degree] = min(knots) - gamma*degree:1
    intervals[(degree+K+1):(2*degree+K)] = max(knots) + gamma*1:degree
  }
  
  return(intervals)
}

bsplines = function(data, knots, degree)
{
  K = length(knots)
  num_basis = K + degree - 1
  
  intervals = find_intervals(knots, degree)
  
  m = matrix(nrow = num_basis, ncol = length(data), 0)
  
  for(k in 1:num_basis)
  {
    m[k,] = bspline(data, k, degree + 1, intervals) 
  }
  
  m[num_basis,length(data)] = 1
  
  return(m)
}

add_init <- function(stan_data){
  
  #  initial values
  
  stan_init <- list()
  
  # for fit inland and fishing together
  stan_init[['log_beta_baseline']] = 0
  stan_init[['log_beta_baseline_contrast_direction']] =  0
  stan_init[['log_beta_baseline_contrast_round']] = rep(0, stan_data[['N_ROUND']]-1)
  stan_init[['log_beta_baseline_contrast_period']] = 0
  stan_init[['rho_gp1']] = array(2, dim = c(stan_data[['N_DIRECTION']]))
  stan_init[['rho_gp2']] = array(2, dim = c(stan_data[['N_DIRECTION']]))
  stan_init[['alpha_gp']] = array(1, dim = c(stan_data[['N_DIRECTION']]))
  
  return(stan_init)
}

