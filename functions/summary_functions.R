keep.likely.transmission.pairs <- function(dchain, threshold){
  
  dchain <- dchain[SCORE_LINKED>threshold]
  dchain[SCORE_DIR_12 <= threshold & SCORE_DIR_21 <= threshold, EST_DIR:='unclear']
  dchain[SCORE_DIR_12 > threshold, EST_DIR:='12']
  dchain[SCORE_DIR_21 > threshold, EST_DIR:='21']
  
  # find source recipient
  dchain <- dchain[EST_DIR != 'unclear']
  dchain[, `:=` (SOURCE=H1, RECIPIENT=H2)]
  dchain[EST_DIR == '21', `:=` (SOURCE=H2, RECIPIENT=H1) ]
  dchain[, `:=` (H1=NULL, H2=NULL)]
}

make.time.since.infection <- function(time.since.infection)
{
  colnames(time.since.infection)[1:3] <- tolower(colnames(time.since.infection))[1:3]
  
  setnames(time.since.infection, c('visit_dt'), c('date_collection'))
  
  time.since.infection[, date_collection := as.Date(date_collection)]
  time.since.infection <- time.since.infection[, .(study_id, date_collection, TSI_estimated_mean, TSI_estimated_min, TSI_estimated_max)]
  time.since.infection <- unique(time.since.infection)
  
  return(time.since.infection)
}

make.df.period <- function(start_observational_period, cutoff_date, stop_observational_period){
  tmp <- data.table(PERIOD = c(paste0(format(start_observational_period, '%b %Y'), '-', format(cutoff_date-31, '%b %Y')), 
                                     paste0(format(cutoff_date, '%b %Y'), '-', format(stop_observational_period, '%b %Y'))), 
             BEFORE_CUTOFF = c(T, F), 
             INDEX_TIME = 1:2, 
             PERIOD_SPAN = c(.year.diff(cutoff_date, start_observational_period), 
                             .year.diff(stop_observational_period, cutoff_date)))
  
  stopifnot(tmp[, sum(PERIOD_SPAN)] == .year.diff(stop_observational_period, start_observational_period))
  return(tmp)
}

find.time.of.infection <- function(meta, time.since.infection, use.TSI.estimate){
  
  # time.since.infection to produce age at infection and time.since.infection
  if(use.TSI.estimate)
  {
    cat("# MAKE TSI ADJUSTMENT\n")
    tmp <- unique(meta[,.(study_id, date_birth)])
    tmp <- merge(time.since.infection, tmp, by='study_id')
    tmp[, age_collection := .year.diff(date_collection, date_birth)]
    cols <- grep('TSI_estimated', colnames(tmp), value=T)
    tmp[, `:=` (age_infection_mean = age_collection - TSI_estimated_mean,
                age_infection_max = age_collection - TSI_estimated_min,
                age_infection_min = age_collection - TSI_estimated_max)]
    tmp[, (cols) := lapply(.SD, .f),.SDcols=cols]
    setnames(tmp, 'age_infection_mean', 'age_infection')
    tmp <- tmp[, .(study_id, age_infection)]
    
    meta <- merge(meta, tmp, by='study_id', all.x=T) 
    
    tmp <- unique(meta[,.(study_id, age_infection)])
    cat(tmp[, sum(!is.na(age_infection))], "TSI estimates out of", tmp[,.N], "individuals in meta data\n")
    
    cat(tmp[is.na(age_infection), .N], " HIV-positive individuals do not have an estimate for age at infection\n" )
    cat(tmp[is.na(age_infection), round(mean(study_id %in% time.since.infection$study_id)*100,2)], "% of which are included in the TSI analysis and do not have sampling date\n")
    
  }else{
    
    cat("# ASSUME INFECTION 1y PRIOR TO DIAGNOSIS\n")
    meta[, age_infection := age_first_positive - 1]
    meta[is.na(age_first_positive), age_infection := age_first_visit - 1]
    tmp <- unique(meta[, .(study_id,age_infection,round)])
    cat(tmp[age_infection <= 16, .N], "out of", tmp[, .N], "HIV-positive individuals are estimated to have been infected prior to 16 yo\n")
    
    cat(tmp[is.na(age_infection), .N], " HIV-positive individuals do not have an estimate for age at infection, from cohort:" )
    print.table(table(tmp[is.na(age_infection), round]))
  }
  
  meta[, `:=` (date_birth = as.Date(date_birth),
               date_first_visit = as.Date(date_first_visit),
               date_last_visit = as.Date(date_last_visit),
               date_first_positive = as.Date(date_first_positive))]
  
  meta[, date_infection := date_birth + 365*(age_infection), by = 'study_id']
  
  return(meta)
}

pairs.get.meta.data <- function(chain, meta, aik){
  
  # stopifnot(unique(dchain$SOURCE) %in% unique(meta_data$aid))
  # stopifnot(unique(dchain$RECIPIENT) %in% unique(meta_data$aid))
  meta <- copy(meta_data); 
  
  # find first roung 
  meta[, first_round := substr(round, 1, 4)]
  meta[round == 'neur', round := 'neuro']

  # merge by source, then recipient
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.SOURCE')
  
  tmp1 <- merge(chain, tmp, by.x = 'SOURCE', by.y = 'aid.SOURCE')
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.RECIPIENT')
  tmp1 <- merge(tmp1, tmp, by.x = 'RECIPIENT', by.y = 'aid.RECIPIENT')
  
  # find age transmission source
  tmp1[, AGE_TRANSMISSION.SOURCE := as.numeric(date_infection.RECIPIENT - date_birth.SOURCE)/365]
  
  # individuals without meta data
  missing_indiv = unique(c(chain$SOURCE, chain$RECIPIENT)[!(c(chain$SOURCE, chain$RECIPIENT) %in% meta$aid)])
  missing_indiv <- .aid2pt(missing_indiv, aik)
  cat('There are ', length(missing_indiv), 'indivs without meta data:\n' )
  missing_indiv <- grep('RK-', missing_indiv, value=T)
  cat('- ', length(missing_indiv), 'of which are in the Rakai Cohort.\n')
  cat(missing_indiv, '\n')
  
  # upper column name
  colnames(tmp1) <- toupper(colnames(tmp1))
  
  return(tmp1)
}

print.statements.about.pairs <- function(pairs){
  
  cat('\nThere is ', nrow(pairs), ' source-recipient pairs\n\n')
  
  cat(nrow(pairs[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]), ' pairs have a proxy for the age at infection of the source and recipient\n')
  cat(nrow(pairs[((SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M') | (SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F')) & (!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT))]), ' pairs are heteroxuals have a proxy for the time of infection of the source and recipient\n\n')                
  
  #   cat('\nPairs by cohort')
  #   tab <- pairs[, list(count = .N), by = c('cohort.SOURCE', 'cohort.RECIPIENT')]
  #   print_table(tab)
  
  #   cat('\nPairs enrolled in RCCS by cohort round')
  #   tab <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS', list(count = .N), by = c('cohort_round.SOURCE', 'cohort_round.RECIPIENT')]
  #   print_table(tab)
  
  cat('\nPairs by sex')
  tab <- pairs[, list(count = .N), by = c('SEX.SOURCE', 'SEX.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by community')
  tab <- pairs[, list(count = .N), by = c('COMM.SOURCE', 'COMM.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by round')
  tab <- pairs[, list(count = .N), by = c('FIRST_ROUND.SOURCE')]
  print_table(tab[order(FIRST_ROUND.SOURCE)])
  
  tab <- pairs[, list(count = .N), by = c('FIRST_ROUND.RECIPIENT')]
  print_table(tab[order(FIRST_ROUND.RECIPIENT)])
}

# Want to see how many aids in potential pairs have an associated prefix and base frequency file 
# -> can compute TSI
print.statements.about.basefreq.files <- function(chain)
{
  
  stopifnot(file.exists(file.path.phscinput) & file.exists(file.path.bflocs))
  
  # get aids of potential pairs
  aid <- chain[, unique(c(SOURCE, RECIPIENT))]
  stopifnot(all(aid %in% aik$aid))

  # Translate aid <-> PREFIXes
  tmp <- file.path(file.path.phscinput)
  tmp <- as.data.table(readRDS(tmp))
  tmp[, `:=` (aid = gsub('-fq[0-9]$','',RENAME_ID), PREFIX=gsub('_remap$','',basename(SAMPLE_ID)))]
  
  # TODO: how is it possible that there are AID in the
  # potential pairs for which we do not have any sequence identifier? 
  aik[ aid %in% aid[which(!aid %in% tmp$aid)] ]
  
  cat(sum(!aid %in% tmp$aid), "out of",length(aid),'individuals are not in the PHSC input (ie: dont know which sequence was used)\n')

  stopifnot(all(aid %in% tmp$aid))
  tmp <- unique(tmp[aid %in% aid, .(aid, PREFIX)])
 
  # Translate PREFIXES <-> bf.csv file existance 
  bflocs <- as.data.table(readRDS(file.path.bflocs)) 
  bflocs <- unique(bflocs[, list(PREFIX=PREFIX, HPC_EXISTS=!is.na(FULL))])
  
  # cat('Out of ', tmp[, length(unique(PREFIX))], ' prefixes associated with our individuals, XX are not in the \n')
  # tmp$PREFIX[which(! tmp$PREFIX %in% bflocs$PREFIX)]
  
  tmp <- merge(tmp, bflocs, by='PREFIX', all.x=T)
  tmp[is.na(HPC_EXISTS), HPC_EXISTS := FALSE]
  
  cat('PREFIXes with existing base frequency file:')
  print_table(tmp[, table(HPC_EXISTS)])
  
  cat('aids with at least one existing base frequency file:')
  tmp1 <- tmp[, list(HPC_EXISTS=any(HPC_EXISTS)), by='aid'][, table(HPC_EXISTS)]
  print_table(tmp1)
  
  if(0) # if want to send Tanya
  {
    missing_bff <- copy(tmp)
    missing_bff[, HPC_EXISTS := NULL]
    
    # if want to send Tanya:
    tmp <- missing_bff[,  .(pt_id, PREFIX)]
    name <- file.path(indir.repository, 'data/missing_bf_files_20220106.csv')
    if(!file.exists(name)){write.csv(tmp, name, row.names = F)}
    # commented out here, but all missing bf's have RCCS2 prefixes.
  }
  
  return(tmp)
}

get.age.map <- function(pairs, age_bands_reduced = 4){
  
  extended_age_length <- 0
  
  ages_source <- pairs[, {
    min_age = floor(min(c(AGE_TRANSMISSION.SOURCE,AGE_INFECTION.RECIPIENT))) - extended_age_length
    max_age = ceiling(max(c(AGE_TRANSMISSION.SOURCE,AGE_INFECTION.RECIPIENT))) + extended_age_length
    list(age = min_age:max_age)}]
  
  ages_recipient <- ages_source
  
  age_map <- data.table(expand.grid(AGE_TRANSMISSION.SOURCE = ages_source$age, 
                                    AGE_INFECTION.RECIPIENT = ages_recipient$age))
  df_age <- age_map[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]

  ages <- sort(unique(df_age$AGE_TRANSMISSION.SOURCE))
  ages <- data.table(age_infection = ages, 
                     AGE_TRANSMISSION_REDUCED.SOURCE = rep(seq(min(ages), max(ages), age_bands_reduced), each = age_bands_reduced )[1:length(ages)])
  df_age <- merge(df_age, ages, by.x = 'AGE_TRANSMISSION.SOURCE', by.y = 'age_infection')
  
  ages <- sort(unique(df_age$AGE_INFECTION.RECIPIENT))
  ages <- data.table(age_infection = ages, 
                     AGE_INFECTION_REDUCED.RECIPIENT = rep(seq(min(ages), max(ages), age_bands_reduced), each = age_bands_reduced )[1:length(ages)])
  df_age <- merge(df_age, ages, by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'age_infection')

  setkey(df_age, AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)
  
  df_age[, index_age := 1:nrow(df_age)]
  
  range_age_non_extended <<- c(min(df_age$AGE_INFECTION.RECIPIENT) + extended_age_length, 
                               max(df_age$AGE_INFECTION.RECIPIENT) - extended_age_length)
  
  colnames(df_age) <- toupper(colnames(df_age))
  
  return(df_age)
}

get.df.direction <- function(){
  df_direction <- data.table(INDEX_DIRECTION = 1:2, IS_MF = c(0, 1))
  df_direction[, LABEL_DIRECTION := ifelse(IS_MF == 1, 'Male -> Female', 'Female -> Male')]
  df_direction
}

get.df.community <- function(){
    df_community <- data.table(INDEX_COMMUNITY = 1:2, COMM = c('fishing','inland'))
    df_community[, LABEL_COMMUNITY := ifelse(COMM == 'inland', 'Inland communities', 'Fishing communities')]

  df_community
}

get.group.map <- function(stratify.by.community.recipient, df_period){
  
  df_direction <- data.table(INDEX_DIRECTION = 1:2, IS_MF = c(0, 1))
  df_direction[, LABEL_DIRECTION := ifelse(IS_MF == 1, 'Male -> Female', 'Female -> Male')]
  
  tmp <- data.table(INDEX_TIME = 1:2, BEFORE_CUTOFF = c(T, F))
  df_period <- merge(df_period, tmp, by = 'BEFORE_CUTOFF')
  
  if(stratify.by.community.recipient){
    df_community <- data.table(INDEX_COMMUNITY = 1:2, COMM = c('fishing','inland'))
    df_community[, LABEL_COMMUNITY := ifelse(COMM == 'inland', 'Inland communities', 'Fishing communities')]
  } else{
    df_community <- data.table(INDEX_COMMUNITY = 1, COMM = c('all'), LABEL_COMMUNITY = 'All communities')
  }

  df_group <- data.table(expand.grid(INDEX_DIRECTION = df_direction[, INDEX_DIRECTION],
                                     INDEX_TIME = df_period[, INDEX_TIME], 
                                     INDEX_COMMUNITY = df_community[, INDEX_COMMUNITY]))
  
  df_group[, index_group := 1:nrow(df_group)]
  
  df_group <- merge(df_group, df_direction, by = 'INDEX_DIRECTION')
  df_group <- merge(df_group, df_period, by = 'INDEX_TIME')
  df_group <- merge(df_group, df_community, by = 'INDEX_COMMUNITY')
  
  setkey(df_group, index_group)
  
  return(df_group)
}

get.age.aggregated.map <- function(age_aggregated){
  

  df_age_aggregated <- data.table(expand.grid(AGE_GROUP_INFECTION.RECIPIENT = age_aggregated, AGE_GROUP_TRANSMISSION.SOURCE = age_aggregated))
  df_age_aggregated[, age_from.RECIPIENT := gsub('(.+)-.*', '\\1', AGE_GROUP_INFECTION.RECIPIENT)]
  df_age_aggregated[, age_from.SOURCE := gsub('(.+)-.*', '\\1', AGE_GROUP_TRANSMISSION.SOURCE)]
  df_age_aggregated[, age_to.RECIPIENT := gsub('.*-(.+)', '\\1', AGE_GROUP_INFECTION.RECIPIENT)]
  df_age_aggregated[, age_to.SOURCE := gsub('.*-(.+)', '\\1', AGE_GROUP_TRANSMISSION.SOURCE)]
  
  tmp <- df_age_aggregated[, list(AGE_INFECTION.RECIPIENT = unique(age_from.RECIPIENT):unique(age_to.RECIPIENT)), by = c('AGE_GROUP_INFECTION.RECIPIENT')]
  tmp1 <- df_age_aggregated[, list(AGE_TRANSMISSION.SOURCE = unique(age_from.SOURCE):unique(age_to.SOURCE)), by = c('AGE_GROUP_TRANSMISSION.SOURCE')]
  
  df_age_aggregated <- merge(df_age_aggregated, tmp, by = 'AGE_GROUP_INFECTION.RECIPIENT', allow.cartesian=TRUE)
  df_age_aggregated <- merge(df_age_aggregated, tmp1, by = 'AGE_GROUP_TRANSMISSION.SOURCE', allow.cartesian=TRUE)
  
  df_age_aggregated[, AGE_CLASSIFICATION.SOURCE := 'Same age']
  df_age_aggregated[AGE_INFECTION.RECIPIENT < (AGE_TRANSMISSION.SOURCE - 5), AGE_CLASSIFICATION.SOURCE := 'Older']
  df_age_aggregated[AGE_INFECTION.RECIPIENT > (AGE_TRANSMISSION.SOURCE + 5), AGE_CLASSIFICATION.SOURCE := 'Younger']
  
  return(df_age_aggregated)
}

add_susceptible_infected <- function(eligible_count, proportion_prevalence){
  
  df <- copy(proportion_prevalence)
  df[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
  
  df <- merge(eligible_count[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)], df, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))

  # find infected
  df[, INFECTED := ELIGIBLE * PREVALENCE_M]
  
  # find susceptible
  df[, SUSCEPTIBLE := ELIGIBLE - INFECTED]
  
  return(df)
}

add_infected_unsuppressed <- function(eligible_count_round, proportion_unsuppressed, df_round, start_observational_period, stop_observational_period){
  
  # use round 15 for round 14 for eligible count
  di <- as.data.table(eligible_count_round)
  di15 <- di[ROUND == '15']
  di15[, ROUND := '14']
  di <- rbind(di15, di )
  
  #use round 15 for round 14
  pu <- as.data.table(proportion_unsuppressed)
  rounds_fill <- c('R014')
  pu <- pu[!ROUND %in% rounds_fill]
  for(round in rounds_fill){
    pu15 <- pu[ROUND == 'R015']
    pu15[, ROUND := round]
    pu <- rbind(pu15, pu)
  }
    
  # select variabel
  di <- di[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE, INFECTED, SUSCEPTIBLE)]
  
  # keep inside observational period
  df_round[, round := as.character(round)]
  df_round[round == '15.1', round := '15S']
  di <- merge(di, df_round, by.x = 'ROUND', by.y = 'round')
  di <- di[min_sample_date >= start_observational_period & max_sample_date <= stop_observational_period ]
  
  # find proportion of unsuppressed by round
  di[, ROUND := paste0('R0', ROUND)]
  df <- merge(di, pu, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # get infected non suppressed
  df[, INFECTED_NON_SUPPRESSED := INFECTED * M]
  df[, INFECTED_NON_SUPPRESSED_CL := INFECTED * CL]
  df[, INFECTED_NON_SUPPRESSED_CU := INFECTED * CU]
  
  # rm unecessary variable
  df <- select(df, -c('M', "CL", "CU"))
  
  return(df)
}

summarise_eligible_count_period <- function(eligible_count_round, cutoff_date, df_period){
  
  # find time intervals
  eligible_count_round[, BEFORE_CUTOFF := max_sample_date <= cutoff_date]
  
  # summarise across time periods
  df <- melt.data.table(eligible_count_round, id.vars = c('ROUND', 'SEX', 'COMM', 'AGEYRS', 'min_sample_date', 'max_sample_date', 'BEFORE_CUTOFF', 'INDEX_TIME'))
  
  # fill missing months
  df[ROUND == 'R014', max_sample_date := df_round[round == 15, min_sample_date]]
  df[ROUND == 'R015', max_sample_date := df_round[round == 16, min_sample_date]]
  df[ROUND == 'R016', max_sample_date := df_round[round == 17, min_sample_date]]
  df[ROUND == 'R017', max_sample_date := df_round[round == 18, min_sample_date]]
  
  # find length in years of each round
  df[, ROUND_SPANYRS := .year.diff(max_sample_date, min_sample_date)]
  
  # check that the lengh corresponds to the one of the period
  tmp <- unique(df[, .(ROUND, min_sample_date, max_sample_date, ROUND_SPANYRS)][order(ROUND)])
  tmp[, round := (gsub('R0(.+)', '\\1', ROUND))]
  tmp <- merge(tmp, df_round, by = 'round')
  tmp <- tmp[, list(PERIOD_SPAN_WITH_ROUND = sum(ROUND_SPANYRS)), by = 'INDEX_TIME']
  tmp <- merge(tmp, df_period, by = 'INDEX_TIME')
  stopifnot(tmp[, all(PERIOD_SPAN == PERIOD_SPAN_WITH_ROUND)])
  
  # find weight of every round in the period 
  tmp <- unique(df[, .(ROUND, ROUND_SPANYRS, BEFORE_CUTOFF)])
  tmp <- tmp[, list(WEIGHT_ROUND = ROUND_SPANYRS / sum(ROUND_SPANYRS),
                    ROUND = ROUND), by = 'BEFORE_CUTOFF']
  df <- merge(df, tmp, by = c('BEFORE_CUTOFF', 'ROUND'))
  
  # find variable over period 
  dfw <- df[, list(value = sum(value * WEIGHT_ROUND)), by = c('BEFORE_CUTOFF', 'SEX', 'COMM', 'AGEYRS', 'variable')]
  
  # plot
  if(0){
    tmp <- copy(dfw)
    tmp[, ROUND := 'average']
    tmp <- rbind(tmp, df, fill = T)
    
    tmp1 <- tmp[COMM == 'inland' & SEX == 'M' & variable %in% c('ELIGIBLE', 'SUSCEPTIBLE', 'INFECTED', 'INFECTED_NON_SUPPRESSED')]
    ggplot(tmp1, aes(x = AGEYRS, y = value))  +
      geom_line(aes(col = ROUND)) + 
      facet_grid(BEFORE_CUTOFF~variable)
  }
  
  setnames(dfw, 'value', 'count')
  
  # merge to period
  dfw <- merge(dfw, df_period, by = 'BEFORE_CUTOFF')
  
  return(dfw)
}

get_incidence_cases_round <- function(incidence, eligible_count_round, full_time_period = T){
  
  # prepare incidence
  colnames(incidence) <- toupper(colnames(incidence))
  setnames(incidence, 'AGE', 'AGEYRS')
  incidence[, COMM := 'inland']
  incidence[, SEX := substring(SEX, 1, 1)]
  incidence <- incidence[ROUND >= 14]
  
  if(grepl('R0', eligible_count_round[, ROUND[1]])){
    # add R0 in front of round index 
    incidence[, ROUND := paste0('R0', as.character(ROUND))]
  }else{
    incidence[, ROUND := as.character(ROUND)]
  }
  
  # for now set incidence in fishing to be the same as in inland
  incidence2 <- copy(incidence)
  incidence2[, COMM := 'fishing']
  incidence <- rbind(incidence, incidence2)
  
  if(0){
    ggplot(incidence, aes(x = AGEYRS)) +
      geom_line(aes(y = INCIDENCE, col = ROUND)) +
      geom_ribbon(aes(ymin = LB, ymax = UB, fill = ROUND),  alpha = 0.1) +
      labs(y = 'Incidence rate per 1 PY in inland community', x = 'Age') +
      facet_grid(SEX~COMM, label = 'label_both') +
      theme_bw() +
      theme(legend.position = 'bottom')
  }
  
  # merge to susceptible
  dir <- merge(incidence, eligible_count_round, by = c('COMM', 'AGEYRS', 'SEX', 'ROUND'))
  
  if(full_time_period){
    
    # fill missing months
    dir[ROUND == 'R014', max_sample_date := df_round[round == 15, min_sample_date]]
    dir[ROUND == 'R015', max_sample_date := df_round[round == 16, min_sample_date]]
    dir[ROUND == 'R016', max_sample_date := df_round[round == 17, min_sample_date]]
    dir[ROUND == 'R017', max_sample_date := df_round[round == 18, min_sample_date]]
    
    # find length in years of each round
    dir[, ROUND_SPANYRS := .year.diff(max_sample_date, min_sample_date)]
    
    # check that the lengh corresponds to the one of the period
    tmp <- unique(dir[, .(ROUND, min_sample_date, max_sample_date, ROUND_SPANYRS)][order(ROUND)])
    tmp[, round := (gsub('R0(.+)', '\\1', ROUND))]
    tmp <- merge(tmp, df_round, by = 'round')
    tmp <- tmp[, list(PERIOD_SPAN_WITH_ROUND = sum(ROUND_SPANYRS)), by = 'INDEX_TIME']
    tmp <- merge(tmp, df_period, by = 'INDEX_TIME')
    stopifnot(tmp[, all(PERIOD_SPAN == PERIOD_SPAN_WITH_ROUND)])
    
  } else{
    # find start and end date of rounds
    dir <- merge(dir, df_round, by.x = 'ROUND', by.y = 'round')
    
    # find length in years of each round
    dir[, ROUND_SPANYRS := .year.diff(max_sample_date, min_sample_date)]
    
  }
  
  # find incident cases
  dir[, INCIDENT_CASES:= SUSCEPTIBLE * ROUND_SPANYRS * INCIDENCE]
  dir[, INCIDENT_CASES_UB:= SUSCEPTIBLE * ROUND_SPANYRS * UB]
  dir[, INCIDENT_CASES_LB:= SUSCEPTIBLE * ROUND_SPANYRS * LB]
  
  # plot
  if(0){
    
    ggplot(dir, aes(x = AGEYRS)) +
      geom_line(aes(y = INCIDENT_CASES, col = ROUND)) +
      geom_ribbon(aes(ymin = INCIDENT_CASES_LB, ymax = INCIDENT_CASES_UB, fill = ROUND), alpha = 0.1) +
      labs(y = 'Expected number of incident cases', x = 'Age') +
      facet_grid(COMM~SEX, label = 'label_both', scale = 'free') +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    tmp <- dir[, list(INCIDENT_CASES = round(sum(INCIDENT_CASES), digits = 1),
                      INCIDENT_CASES_UB = round(sum(INCIDENT_CASES_UB), digits = 1),
                      INCIDENT_CASES_UB = round(sum(INCIDENT_CASES_UB), digits = 1)), by = c('ROUND', 'MODEL')]
    knitr::kable(subset(tmp, select = - c(MODEL)))
    
  }
  
  return(dir)
  
}


summarise_incidence_cases_period <- function(incidence_cases_round, cutoff_date, df_period){
  
  # find time intervals
  incidence_cases_round[, BEFORE_CUTOFF := max_sample_date <= cutoff_date]
  
  # summarise across time periods
  incidence_cases <- incidence_cases_round[, list(INCIDENT_CASES = sum(INCIDENT_CASES), 
                                                  INCIDENT_CASES_UB = sum(INCIDENT_CASES_UB), 
                                                  INCIDENT_CASES_LB = sum(INCIDENT_CASES_LB)), by = c('COMM', 'AGEYRS', 'SEX', 'BEFORE_CUTOFF')]
  # make period
  incidence_cases <- merge(incidence_cases, df_period, by = 'BEFORE_CUTOFF')
  
  return(incidence_cases)
}

get_proportion_sampling <- function(pairs, incidence_cases){
  # find number of pair observed
  dp <- copy(pairs)
  setnames(dp, c('SEX.RECIPIENT', 'COMM.RECIPIENT', 'AGE_INFECTION.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT'), 
           c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF'))
  dp[, AGEYRS := floor(AGEYRS)]
  dp <- dp[, list(count = .N), by = c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF')]
  
  # merge to incidence cases
  di <- merge(incidence_cases, dp, by = c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF'), all.x = T)
  di[is.na(count), count := 0]
  
  # find proportion sampling
  di[, prop_sampling := count / INCIDENT_CASES]
  
  if(0){
    tmp <- di[COMM == 'inland']
    ggplot(tmp, aes(x = AGEYRS, y = prop_sampling, col =SEX)) + 
      geom_line() + 
      facet_grid(PERIOD~COMM) 
  }
  
  di
}

prepare_unsuppressed <- function(eligible_count){
  tmp <- eligible_count[variable == 'INFECTED_NON_SUPPRESSED']
  setnames(tmp, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
  tmp[, IS_MF := as.numeric(SEX == 'M')]
  tmp <- merge(tmp, df_direction, by = 'IS_MF')
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp
}

make.df.round <- function(df_round, df_period){
  
  df_round[, INDEX_TIME := 0]
  df_round[round%in%14:15, INDEX_TIME := 1]
  df_round[round %in% 16:18, INDEX_TIME := 2]
  
  return(df_round)
}

find_log_offset_by_round <- function(stan_data, eligible_count_round, df_age, df_direction, df_community, df_period){
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  eligible_count_wide[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  
  ROUNDS <- eligible_count_wide[, unique(ROUND)]
  
  res <- list(); index = 1
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in seq_along(ROUNDS)){
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- ROUNDS[k]
        
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
        tmp <- df_round[paste0('R0', toupper(round)) == .ROUND]
        log_offset[, PERIOD_SPAN := .year.diff(tmp[, max_sample_date], tmp[, min_sample_date])]
        
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

prepare.proportion.unsuppresed <- function(proportion_unsuppressed){
  
  # use round 15 for round 16
  proportion_unsuppressed <- proportion_unsuppressed[!ROUND %in% c('R016', 'R015S')]
  proportion_unsuppressed15 <- proportion_unsuppressed[ROUND == 'R015']
  proportion_unsuppressed15[, ROUND := 'R016']
  proportion_unsuppressed <- rbind(proportion_unsuppressed, proportion_unsuppressed15)
  proportion_unsuppressed <- proportion_unsuppressed[order(ROUND)]
  
  if(0){
    # plot
    ggplot(proportion_unsuppressed, aes(x = AGEYRS)) + 
      geom_point(aes(y = PROP_NON_SUPPRESSED_EMPIRICAL, col = ROUND), alpha = 0.25) +
      geom_line(aes(y = M, col = ROUND)) + 
      # geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX), alpha = 0.5) + 
      facet_grid(COMM~SEX)+
      theme_bw() + 
      theme(legend.position='bottom')
  }
  
  return(proportion_unsuppressed)
}

find_crude_force_infection <- function(stan_data){
  
  # retrieve observed transmission events and offset
  tmp1 = as.data.table( reshape2::melt(stan_data[['y']]) )
  setnames(tmp1, 1:5, c('INDEX_AGE', 'INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'Y'))
  tmp2 = as.data.table( reshape2::melt(stan_data[['log_offset']]) )
  setnames(tmp2, 1:5, c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE', 'LOG_OFFSET'))
  tmp1 <- merge(tmp1, tmp2, by = c('INDEX_DIRECTION', 'INDEX_COMMUNITY', 'INDEX_TIME', 'INDEX_AGE'))
  tmp1 <- merge(tmp1, df_age, by = 'INDEX_AGE')
  tmp1 <- merge(tmp1, df_direction, by = 'INDEX_DIRECTION')
  tmp1 <- merge(tmp1, df_community, by = 'INDEX_COMMUNITY')
  tmp1 <- merge(tmp1, df_period, by = 'INDEX_TIME')
  
  # find offset
  tmp1[, OFFSET := exp(LOG_OFFSET)]
  
  # find crude FOI
  tmp1[, CRUDE_FOI := Y / OFFSET]
  
  
  return(tmp1)
}

