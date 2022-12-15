find_log_offset_by_round <- function(stan_data, eligible_count_round, df_estimated_contact_rates, 
                                     use_number_susceptible_offset, use_contact_rates_prior)
{
  # find log offset including the proportion of susceptible, number of unsuppressed and time period of each round
  # come handy when we want to specify different offset formula 
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  eligible_count_wide[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  
  
  res <- list(); index = 1
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          res[[index]]  = res[[index-1]] 
          next
        }
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- df_round[COMM == .COMM & INDEX_ROUND == k, ROUND]
        
        log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        log_offset[, INDEX_DIRECTION := df_direction[i, INDEX_DIRECTION]]
        log_offset[, INDEX_COMMUNITY := df_community[j, INDEX_COMMUNITY]]
        log_offset[, ROUND := .ROUND]
        
        # add # susceptible and proportion of susceptible in recipient
        tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, SUSCEPTIBLE, PROP_SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
        
        # add number of infected unsuppressed in source
        tmp <- eligible_count_wide[SEX == .SEX.SOURCE & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, INFECTED_NON_SUPPRESSED)], by.x = 'AGE_TRANSMISSION.SOURCE', by.y = 'AGEYRS')
        
        # add contact rates
        tmp <- df_estimated_contact_rates[part.sex == .SEX.SOURCE]
        colnames(tmp) <- toupper(colnames(tmp))
        log_offset <- merge(log_offset, tmp[, .(PART.AGE, CONT.AGE, CNTCT.RATE)], 
                            by.x = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), by.y = c('PART.AGE', 'CONT.AGE'))
        
        # add period in year
        tmp <- df_round[ROUND == .ROUND & COMM == .COMM]
        log_offset[, PERIOD_SPAN := tmp[, ROUND_SPANYRS]]
        
        # make log offset
        if(use_number_susceptible_offset){
          log_offset[, LOG_OFFSET := log(SUSCEPTIBLE) + log(INFECTED_NON_SUPPRESSED) + log(PERIOD_SPAN)]
        }else{
          log_offset[, LOG_OFFSET := log(PROP_SUSCEPTIBLE) + log(INFECTED_NON_SUPPRESSED) + log(PERIOD_SPAN)]
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
        res[[index]] = log_offset
        index = index + 1
        
      }
    }
  }
  res <- do.call('rbind', res)
  
  res[, log_CONTACT_RATES := log(CNTCT.RATE)]
  res[, log_INFECTED_NON_SUPPRESSED := log(INFECTED_NON_SUPPRESSED)]
  res[, log_PROP_SUSCEPTIBLE := log(PROP_SUSCEPTIBLE)]
  res[, log_SUSCEPTIBLE := log(SUSCEPTIBLE)]
  res[, log_PERIOD_SPAN:= log(PERIOD_SPAN)]
  
  return(res)
}

prepare_count_data <- function(stan_data){
  
  # reshape the number of source recipient pairs for figure
  # merge to the maps
  
  tmp <- as.data.table(reshape2::melt(stan_data[['y']]))
  setnames(tmp, 1:3, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_TIME'))
  tmp <- merge(tmp, df_direction, by = 'INDEX_DIRECTION')
  tmp <- merge(tmp, df_period, by = c('INDEX_TIME'))
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp <- merge(tmp, df_age, by = 'INDEX_AGE')
  
  setnames(tmp, 'value', 'count')
  
  return(as.data.table(tmp))
}

prepare_incidence_cases <- function(incidence_cases){
  
  # reshape the number of incident cases for figure
  # merge to the maps
  
  tmp <- copy(incidence_cases)
  setnames(tmp, 'AGEYRS', 'AGE_INFECTION.RECIPIENT')
  tmp[, IS_MF := as.numeric(SEX == 'F')]
  tmp <- merge(tmp, df_direction, by = 'IS_MF')
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp
}

prepare_unsuppressed_share <- function(unsuppressed_share, vars, standardised.var = NULL){
  
  # reshape the share of unsuppressed count by sex
  # merge to the maps
  
  tmp1 <- copy(unsuppressed_share)
  
  if(all(c('SEX', 'AGEYRS') %in%vars)){
    
    if(!is.null(standardised.var)){
      tmp1 <- tmp1[, .(ROUND, COMM, SEX, AGEYRS, UNSUPPRESSED_SHARE_AGE_BY_SEX_CL, UNSUPPRESSED_SHARE_AGE_BY_SEX_CU, UNSUPPRESSED_SHARE_AGE_BY_SEX_M)]
      setnames(tmp1, c('UNSUPPRESSED_SHARE_AGE_BY_SEX_CL', 'UNSUPPRESSED_SHARE_AGE_BY_SEX_CU', 'UNSUPPRESSED_SHARE_AGE_BY_SEX_M'), c('CL', "CU", 'M'))
    }else{
      tmp1 <- tmp1[, .(ROUND, COMM, SEX, AGEYRS, UNSUPPRESSED_SHARE_AGE_AND_SEX_CL, UNSUPPRESSED_SHARE_AGE_AND_SEX_CU, UNSUPPRESSED_SHARE_AGE_AND_SEX_M)]
      setnames(tmp1, c('UNSUPPRESSED_SHARE_AGE_AND_SEX_CL', 'UNSUPPRESSED_SHARE_AGE_AND_SEX_CU', 'UNSUPPRESSED_SHARE_AGE_AND_SEX_M'), c('CL', "CU", 'M'))
    }
    
  } else if(vars == 'SEX'){
    
    tmp1 <- unique(tmp1[, .(ROUND, COMM, SEX, UNSUPPRESSED_SHARE_SEX_CL, UNSUPPRESSED_SHARE_SEX_CU, UNSUPPRESSED_SHARE_SEX_M)])
    setnames(tmp1, c('UNSUPPRESSED_SHARE_SEX_CL', 'UNSUPPRESSED_SHARE_SEX_CU', 'UNSUPPRESSED_SHARE_SEX_M'), c('CL', "CU", 'M'))
  
  }else{
    stop()
  }
  
  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  
  # find round label
  tmp1[, ROUND := paste0('R0', ROUND)]
  
  # merge to index round
  tmp1 <- merge(tmp1, df_round[, .(COMM, ROUND, INDEX_ROUND, LABEL_ROUND)],  by = c('COMM', 'ROUND'))
  
  # type
  tmp1[, type := 'Share in the HIV-positive unsuppressed census eligible individuals']
  
  return(tmp1)
}

prepare_unsuppressed_median_age <- function(unsuppressed_median_age){
  tmp1 <- copy(unsuppressed_median_age)
  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  
  # find round label
  tmp1[, ROUND := paste0('R0', ROUND)]
  
  # merge to index round
  tmp1 <- merge(tmp1, df_round[, .(COMM, ROUND, INDEX_ROUND, LABEL_ROUND)],  by = c('COMM', 'ROUND'))
  return(tmp1)
}

prepare_infected_share <- function(infected_share, vars){
  
  # reshape the share of infected count by sex
  # merge to the maps
  
  tmp1 <- copy(infected_share)
  
  if(all(c('SEX', 'AGEYRS') %in%vars)){
    tmp1 <- tmp1[, .(ROUND, COMM, SEX, AGEYRS, PREVALENCE_SHARE_SEX_AND_AGE_CL, PREVALENCE_SHARE_SEX_AND_AGE_CU, PREVALENCE_SHARE_SEX_AND_AGE_M)]
    setnames(tmp1, c('PREVALENCE_SHARE_SEX_AND_AGE_CL', 'PREVALENCE_SHARE_SEX_AND_AGE_CU', 'PREVALENCE_SHARE_SEX_AND_AGE_M'), c('CL', "CU", 'M'))
  }else if(vars == 'SEX'){
    tmp1 <- unique(tmp1[, .(ROUND, COMM, SEX, PREVALENCE_SHARE_SEX_CL, PREVALENCE_SHARE_SEX_CU, PREVALENCE_SHARE_SEX_M)])
    setnames(tmp1, c('PREVALENCE_SHARE_SEX_CL', 'PREVALENCE_SHARE_SEX_CU', 'PREVALENCE_SHARE_SEX_M'), c('CL', "CU", 'M'))
  }else{
    stop()
  }
  
  tmp1[, IS_MF := as.numeric(SEX == 'M')]
  tmp1 <- merge(tmp1, df_direction, by = 'IS_MF')
  tmp1 <- merge(tmp1, df_community, by = 'COMM')
  tmp1[, ROUND := as.character(ROUND)]
  
  # find round label
  tmp1[, ROUND := paste0('R0', ROUND)]
  
  # merge to index round
  tmp1 <- merge(tmp1, df_round[, .(COMM, ROUND, INDEX_ROUND)], by = c('COMM', 'ROUND'))
  
  # type
  tmp1[, type := 'Share in the HIV-positive census eligible individuals']
  
  return(tmp1)
}


remove_first_round <- function(tmp){
  # remove first round of the contrasts (as there is no contrast on those)
  
  rm.vars <- c('round', 'ROUND', 'ROUND_SPANYRS', 'INDEX_TIME', 'min_sample_date', 'max_sample_date')
  
  for(var in rm.vars){
    if(var %in% names(tmp)){
      tmp <- select(tmp, -var)
    }
  }
  
  tmp[, INDEX_ROUND:= INDEX_ROUND +1]
  tmp <- merge(tmp, df_round, by = c('INDEX_ROUND', 'COMM'))
  
  return(tmp)
}

clean_reported_contact <- function(df_reported_contact){
  
  reported_contact <- copy(df_reported_contact)
  
  # change name of variables
  setnames(reported_contact, c('part.comm', 'round', 'cont.sex', 'cont.age.group'), c('COMM', 'ROUND', 'SEX', 'AGE_GROUP'))
  
  # create variables
  reported_contact[, LABEL_RECIPIENT := ifelse(SEX == 'Female', 'Female recipients', 'Male recipients')]
  reported_contact <- merge(reported_contact, df_direction, by = 'LABEL_RECIPIENT')
  
  # keep age of interset
  reported_contact <- reported_contact[AGE_GROUP != "[50,55)"]
  
  # reshape age group
  reported_contact[AGE_GROUP == '[15,20)', AGE_GROUP := '15-19']
  reported_contact[AGE_GROUP == '[20,25)', AGE_GROUP := '20-24']
  reported_contact[AGE_GROUP == '[25,30)', AGE_GROUP := '25-29']
  reported_contact[AGE_GROUP == '[30,35)', AGE_GROUP := '30-34']
  reported_contact[AGE_GROUP == '[35,40)', AGE_GROUP := '35-39']
  reported_contact[AGE_GROUP == '[40,45)', AGE_GROUP := '40-44']
  reported_contact[AGE_GROUP == '[45,50)', AGE_GROUP := '45-49']
  
  # find mean age group
  reported_contact[, mean_age_group := mean(c(as.numeric(gsub('(.+)-.*', '\\1', AGE_GROUP)), 
                                 as.numeric(gsub('.*-(.+)', '\\1', AGE_GROUP)))), by = 'AGE_GROUP']
  
  
  # find round label
  reported_contact <- merge(reported_contact , df_round, by = c('COMM', 'ROUND'))
  
  return(reported_contact)
}

clean_contribution_sexual_contact <- function(contribution_sexual_contact){
  
  contribution_sexual_contact <- copy(contribution_sexual_contact)
  
  # change name of variables
  setnames(contribution_sexual_contact, c('part.sex', 'part.age'), c('SEX', 'AGEYRS'))
  
  # create variables
  contribution_sexual_contact[, LABEL_SOURCE := ifelse(SEX == 'F', 'Female sources', 'Male sources')]
  
  # restrain age
  contribution_sexual_contact <- contribution_sexual_contact[AGEYRS >= 15 & AGEYRS < 50]
  
  # keep only some variable
  contribution_sexual_contact <- contribution_sexual_contact[, .(LABEL_SOURCE, AGEYRS,prop)]
  setnames(contribution_sexual_contact, 'prop', 'PROP_M')
  contribution_sexual_contact <- merge(contribution_sexual_contact, df_direction, by = 'LABEL_SOURCE')
  
  return(contribution_sexual_contact)
}

# Weighted generic quantile estimator
wquantile.generic <- function(x, probs, cdf.gen, weights = NA) {
  n <- length(x)
  if (any(is.na(weights)))
    weights <- rep(1 / n, n)
  nw <- sum(weights)^2 / sum(weights^2) # Kish's effective sample size
  
  indexes <- order(x)
  x <- x[indexes]
  weights <- weights[indexes]
  
  weights <- weights / sum(weights)
  cdf.probs <- cumsum(c(0, weights))
  
  sapply(probs, function(p) {
    cdf <- cdf.gen(nw, p)
    q <- cdf(cdf.probs)
    w <- tail(q, -1) - head(q, -1)
    sum(w * x)
  })
}

# Weighted Harrell-Davis quantile estimator
whdquantile <- function(x, probs, weights = NA) {
  cdf.gen <- function(n, p) return(function(cdf.probs) {
    pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
  })
  wquantile.generic(x, probs, cdf.gen, weights)
}

# Weighted Type 7 quantile estimator
wquantile <- function(x, probs, weights = NA) {
  cdf.gen <- function(n, p) return(function(cdf.probs) {
    h <- p * (n - 1) + 1
    u <- pmax((h - 1) / n, pmin(h / n, cdf.probs))
    u * n - h + 1
  })
  wquantile.generic(x, probs, cdf.gen, weights)
}

