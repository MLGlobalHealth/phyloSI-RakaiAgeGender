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
  stan_data[['N_ROUND_INLAND']] = df_round[COMM == 'inland', length(unique(INDEX_ROUND))]
  stan_data[['N_ROUND_FISHING']] = df_round[COMM == 'fishing', length(unique(INDEX_ROUND))]
  
  # map from round to period
  stan_data[['map_round_period']] =   rbind(c(tmp[COMM ==  df_community[order(INDEX_COMMUNITY), COMM[1]] & order(round), INDEX_TIME], rep(-1,  max(stan_data[['N_ROUND']]) - stan_data[['N_ROUND']][1] )), 
                                            c(tmp[COMM ==  df_community[order(INDEX_COMMUNITY), COMM[2]] & order(round), INDEX_TIME], rep(-1,  max(stan_data[['N_ROUND']]) - stan_data[['N_ROUND']][2] )))
  
  # number of period
  stan_data[['N_ROUND_PER_PERIOD']] = rbind(tmp[order(INDEX_TIME) & COMM == df_community[order(INDEX_COMMUNITY), COMM[1]], length(unique(INDEX_ROUND)), by = 'INDEX_TIME']$V1, 
                                            tmp[order(INDEX_TIME) & COMM == df_community[order(INDEX_COMMUNITY), COMM[2]], length(unique(INDEX_ROUND)), by = 'INDEX_TIME']$V1)
  
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
  y = array(NA, c(stan_data[['N_PER_GROUP']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']]))
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:stan_data[['N_PERIOD']]){
        
        # direction group
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        tmp <- pairs_round[SEX.SOURCE == .SEX.SOURCE]
        
        # community group
        tmp <- tmp[COMM.RECIPIENT == df_community[j, COMM]]
        
        # time group
        tmp <- tmp[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == df_period[COMM == df_community[j, COMM] & INDEX_TIME == k, BEFORE_CUTOFF]]
        
        # count number of observation
        tmp <- tmp[, list(count = .N), by = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT')]
        tmp <- merge(df_age, tmp, 
                     by = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), all.x = T)
        tmp[is.na(count), count := 0]
        
        setkey(tmp, AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)
        
        tmp1 <- pairs[SEX.SOURCE == .SEX.SOURCE  & COMM.RECIPIENT == df_community[j, COMM] & DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == df_period[COMM == df_community[j, COMM] & INDEX_TIME == k, BEFORE_CUTOFF]]
        stopifnot(sum(tmp$count) == nrow(tmp1))
        
        # check the order of ages is correct
        tmp <- tmp[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == tmp[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == tmp[, AGE_TRANSMISSION.SOURCE])
        
        cat(nrow(tmp1), 'pairs with infection', df_direction[i, LABEL_DIRECTION], 'towards', df_community[j, COMM], 'in', df_period[COMM == df_community[j, COMM] & INDEX_TIME == k, PERIOD], '\n')
        
        y[, i, j, k] = matrix(tmp$count, ncol = 1)
        
        
      }
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
  z = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){ 
          # account for number of rounds in fishing that are < number of rounds in inland. Those are not used by the model
          z[, i, j, k] = -1;
          next
        }
        
        .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM = df_community[INDEX_COMMUNITY == j, COMM]
        .ROUND = df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # direction group
        tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
        
        # community group
        tmp <- tmp[COMM == .COMM]
        
        # round
        tmp <- tmp[ROUND == .ROUND]
        
        # order by age
        tmp <- tmp[order(AGEYRS)] 
        
        # sanity check
        tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT  & COMM == .COMM & ROUND == .ROUND]
        stopifnot(sum(tmp$INCIDENT_CASES) == sum(tmp1$INCIDENT_CASES))
        cat(sum(tmp1$INCIDENT_CASES), 'incidence cases ', df_direction[i, LABEL_DIRECTION], 'towards', .COMM, 'in', .ROUND, '\n')
        
        # fill
        z[, i, j, k] = round(tmp$INCIDENT_CASES)
        
      }
    }
  }
  
  stan_data[['z']] = z

  
  return(stan_data)
  
}

add_incidence_rates <- function(stan_data, incidence_cases_round){
  
  ir = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          # account for number of rounds in fishing that are < number of rounds in inland. Those are not used by the model
          ir[, i, j, k] = 0;
          next
        }
        
        .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM = df_community[INDEX_COMMUNITY == j, COMM]
        .ROUND = df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # direction group
        tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
        
        # community group
        tmp <- tmp[COMM == .COMM]
        
        # round
        tmp <- tmp[ROUND == .ROUND]
        
        # order by age
        tmp <- tmp[order(AGEYRS)] 
        
        # sanity check
        tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT  & COMM == .COMM & ROUND == .ROUND]
        stopifnot(sum(tmp$INCIDENCE) == sum(tmp1$INCIDENCE))
        cat(mean(tmp1$INCIDENCE), 'incidence rates ', df_direction[i, LABEL_DIRECTION], 'towards', .COMM, 'in', .ROUND, '\n')
        
        # fill
        ir[, i, j, k] = tmp$INCIDENCE
        
      }
    }
  }
  
  stan_data[['ir']] = ir
  
  
  return(stan_data)
  
}

add_incidence_rates_lognormal_parameters <- function(stan_data, incidence_cases_round){
  
  ir_lognorm_mean = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  ir_lognorm_sd = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          # account for number of rounds in fishing that are < number of rounds in inland. Those are not used by the model
          ir_lognorm_mean[, i, j, k] = 0
          ir_lognorm_sd[, i, j, k] = 0
          next
        }
        
        .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM = df_community[INDEX_COMMUNITY == j, COMM]
        .ROUND = df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # direction group
        tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
        
        # community group
        tmp <- tmp[COMM == .COMM]
        
        # round
        tmp <- tmp[ROUND == .ROUND]
        
        # order by age
        tmp <- tmp[order(AGEYRS)] 
        
        # sanity check
        tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT  & COMM == .COMM & ROUND == .ROUND]
        stopifnot(sum(tmp$INCIDENCE) == sum(tmp1$INCIDENCE))
        cat(mean(tmp1$INCIDENCE), 'incidence rates ', df_direction[i, LABEL_DIRECTION], 'towards', .COMM, 'in', .ROUND, '\n')
        
        # fill
        lognorm_parms <- lognorm::getParmsLognormForMedianAndUpper(median = tmp$INCIDENCE, upper = tmp$UB, sigmaFac=2)
        ir_lognorm_mean[, i, j, k] = lognorm_parms[, 1]
        ir_lognorm_sd[, i, j, k] = lognorm_parms[, 2]
      }
    }
  }
  
  stan_data[['ir_lognorm_mean']] = ir_lognorm_mean
  stan_data[['ir_lognorm_sd']] = ir_lognorm_sd
  
  return(stan_data)
  
}

add_offset <- function(stan_data, eligible_count){
  
  # add offset including the probability of susceptible and the number of infected unsuppressed
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  eligible_count_wide[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  
  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          # account for number of rounds in fishing that are < number of rounds in inland. Those are not used by the model
          log_offset_array[i, j, k,] = -1
          next
        }
        
        log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # add proportion of susceptible in recipient
        tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, PROP_SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
        
        # add number of infected unsuppressed in source
        tmp <- eligible_count_wide[SEX == .SEX.SOURCE & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, INFECTED_NON_SUPPRESSED)], by.x = 'AGE_TRANSMISSION.SOURCE', by.y = 'AGEYRS')
        
        # make log offset
        log_offset[, LOG_OFFSET := log(PROP_SUSCEPTIBLE) + log(INFECTED_NON_SUPPRESSED)]
        
        # check the order of ages is correct
        log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
        
        # add to array
        log_offset_array[i, j, k,] = log_offset[, LOG_OFFSET]
        
      }
    }
  }
  
  stan_data[['log_offset']] = log_offset_array

  return(stan_data)
}

add_offset_time <- function(stan_data, eligible_count){
  
  # add offset including the time in years of each rounds
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  
  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          # account for number of rounds in fishing that are < number of rounds in inland. Those are not used by the model
          log_offset_array[i, j, k,] = -1
          next
        }
        
        log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # add period in year
        log_offset[, PERIOD_SPAN := df_round[ROUND == .ROUND & COMM == .COMM, ROUND_SPANYRS]]
        
        # make log offset time
        log_offset[, LOG_OFFSET_TIME :=  log(PERIOD_SPAN)]
        
        # check the order of ages is correct
        log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
        
        # add to array
        log_offset_array[i, j, k,] = log_offset[, LOG_OFFSET_TIME]
        
      }
    }
  }
  
  stan_data[['log_offset_time']] = log_offset_array
  
  return(stan_data)
}

add_offset_susceptible <- function(stan_data, eligible_count){
  
  # add offset including the number of susceptible
  # that can be used to calculate incidence rates
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]

  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          # account for number of rounds in fishing that are < number of rounds in inland. Those are not used by the model
          log_offset_array[i, j, k,] = -1
          next
        }
        
        log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # add number of susceptible in recipient
        tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
        
        # make log offset
        log_offset[, LOG_OFFSET_SUSCEPTIBLE := log(SUSCEPTIBLE) ]
        
        # check the order of ages is correct
        log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
        
        # add to array
        log_offset_array[i, j, k,] = log_offset[, LOG_OFFSET_SUSCEPTIBLE]
        
      }
    }
  }
  
  stan_data[['log_offset_susceptible']] = log_offset_array
  
  return(stan_data)
}

add_probability_sampling <- function(stan_data, proportion_sampling){
  
  # add the probability of sampling a infection event
  
  proportion_sampling <- proportion_sampling[order(SEX.RECIPIENT, COMM, BEFORE_CUTOFF, PERIOD)]
  
  log_prop_sampling_array =array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']], stan_data[['N_PER_GROUP']])))
  sampling_index=array(NA, c(c(stan_data[['N_PER_GROUP']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']])))
  n_sampling_index=array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:stan_data[['N_PERIOD']]){
        
        log_prop_sampling = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .BEFORE_CUTOFF <- df_period[k, BEFORE_CUTOFF]
        
        # add probability of sampling recipient 
        tmp <- proportion_sampling[SEX.RECIPIENT == .SEX.RECIPIENT & COMM == .COMM & BEFORE_CUTOFF == .BEFORE_CUTOFF]
        log_prop_sampling <- merge(log_prop_sampling, tmp[, .(AGEYRS.SOURCE, AGEYRS.RECIPIENT, prop_sampling)], by.x = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), 
                            by.y = c('AGEYRS.SOURCE', 'AGEYRS.RECIPIENT'))
        
        # make log prob sampling
        if(1){
          log_prop_sampling[prop_sampling == 0, prop_sampling := 0.0001]
        }
        log_prop_sampling[, LOG_PROP_SAMPLING := log(prop_sampling)]
        
        # check the order of ages is correct
        log_prop_sampling <- log_prop_sampling[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_prop_sampling[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_prop_sampling[, AGE_TRANSMISSION.SOURCE])
        
        # prop sampling
        log_prop_sampling_array[i, j, k,] = log_prop_sampling[, LOG_PROP_SAMPLING]
        
        # was the recipient sampled
        n_sampling_index[i, j, k] = log_prop_sampling[, sum(prop_sampling != 0.0001)]
        sampling_index[,i, j, k] = rep(-1, nrow(log_prop_sampling))
        sampling_index[1:n_sampling_index[i, j, k], i, j, k]  = log_prop_sampling[, which(prop_sampling != 0.0001)]
      }
    }
  }
  
  stan_data[['log_prop_sampling']] = log_prop_sampling_array
  stan_data[['n_sampling_index_y']] = n_sampling_index
  stan_data[['sampling_index_y']] = sampling_index
  
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
  stan_init[['log_beta_baseline_contrast_community']] = 0
  stan_init[['log_beta_baseline_contrast_direction']] =  0
  stan_init[['log_beta_baseline_contrast_round']] = array(0, dim = c(stan_data[['N_ROUND']] - 1, stan_data[['N_DIRECTION']],  stan_data[['N_COMMUNITY']]))
  stan_init[['rho_gp1']] = array(2, dim = c(stan_data[['N_DIRECTION']]))
  stan_init[['rho_gp2']] = array(2, dim = c(stan_data[['N_DIRECTION']]))
  stan_init[['alpha_gp']] = array(1, dim = c(stan_data[['N_DIRECTION']]))
  
  return(stan_init)
}

