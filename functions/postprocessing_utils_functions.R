find_log_offset_by_round <- function(stan_data, eligible_count_round)
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
        
        # add proportion of susceptible in recipient
        tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, PROP_SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
        
        # add number of infected unsuppressed in source
        tmp <- eligible_count_wide[SEX == .SEX.SOURCE & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, INFECTED_NON_SUPPRESSED)], by.x = 'AGE_TRANSMISSION.SOURCE', by.y = 'AGEYRS')
        
        # add period in year
        tmp <- df_round[ROUND == .ROUND & COMM == .COMM]
        log_offset[, PERIOD_SPAN := tmp[, ROUND_SPANYRS]]
        
        # make log offset
        log_offset[, LOG_OFFSET := log(PROP_SUSCEPTIBLE) + log(INFECTED_NON_SUPPRESSED) + log(PERIOD_SPAN)]
        
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
  
  res[, log_INFECTED_NON_SUPPRESSED := log(INFECTED_NON_SUPPRESSED)]
  res[, log_PROP_SUSCEPTIBLE := log(PROP_SUSCEPTIBLE)]
  res[, log_PERIOD_SPAN:= log(PERIOD_SPAN)]
  
  return(res)
}

prepare_count_data <- function(stan_data){
  
  # reshape the number of source recipient pairs for figure
  # merge to the maps
  
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
  
  # reshape the number of incident cases for figure
  # merge to the maps
  
  tmp <- copy(incidence_cases)
  setnames(tmp, 'AGEYRS', 'AGE_INFECTION.RECIPIENT')
  tmp[, IS_MF := as.numeric(SEX == 'F')]
  tmp <- merge(tmp, df_direction, by = 'IS_MF')
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp
}

prepare_unsuppressed_share <- function(unsuppressed_share, vars){
  
  # reshape the share of unsuppressed count by sex
  # merge to the maps
  
  tmp1 <- copy(unsuppressed_share)
  
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
  
  # find round label
  tmp1[, ROUND := paste0('R0', ROUND)]
  
  # merge to index round
  tmp1 <- merge(tmp1, df_round[, .(COMM, ROUND, INDEX_ROUND, LABEL_ROUND)],  by = c('COMM', 'ROUND'))
  
  # type
  tmp1[, type := 'Share in the HIV-positive unsuppressed census eligible individuals']
  
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
  setnames(reported_contact, c('part.comm', 'part.round', 'part.age', 'part.sex'), c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
  
  # create variables
  reported_contact[, LABEL_RECIPIENT := ifelse(SEX == 'F', 'Female recipients', 'Male recipients')]
  
  # find round label
  reported_contact <- merge(reported_contact , df_round, by = c('COMM', 'ROUND'))
  
  return(reported_contact)
}
