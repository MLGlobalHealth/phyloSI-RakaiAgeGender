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

make.df.period <- function(start_observational_period_inland, stop_observational_period_inland, 
                           start_observational_period_fishing, stop_observational_period_fishing, 
                           cutoff_date){
  
  tmp_inland <- data.table(PERIOD = c(paste0(format(start_observational_period_inland, '%b %Y'), '-', format(cutoff_date-31, '%b %Y')), 
                                      paste0(format(cutoff_date, '%b %Y'), '-', format(stop_observational_period_inland, '%b %Y'))), 
                           BEFORE_CUTOFF = c(T, F), 
                           INDEX_TIME = 1:2, 
                           PERIOD_SPAN = c(.year.diff(cutoff_date, start_observational_period_inland), 
                                           .year.diff(stop_observational_period_inland, cutoff_date)), 
                           COMM = 'inland')
  stopifnot(tmp_inland[, sum(PERIOD_SPAN)] == .year.diff(stop_observational_period_inland, start_observational_period_inland))
  
  tmp_fishing <- data.table(PERIOD = c(paste0(format(start_observational_period_fishing, '%b %Y'), '-', format(cutoff_date-31, '%b %Y')), 
                                      paste0(format(cutoff_date, '%b %Y'), '-', format(stop_observational_period_fishing, '%b %Y'))), 
                           BEFORE_CUTOFF = c(T, F), 
                           INDEX_TIME = 1:2, 
                           PERIOD_SPAN = c(.year.diff(cutoff_date, start_observational_period_fishing), 
                                           .year.diff(stop_observational_period_fishing, cutoff_date)), 
                           COMM = 'fishing')
  stopifnot(abs(tmp_fishing[, sum(PERIOD_SPAN)] - .year.diff(stop_observational_period_fishing, start_observational_period_fishing)) < 1e-15)
  
  tmp <- rbind(tmp_inland, tmp_fishing)
  
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
    df_community[, LABEL_COMMUNITY := factor(LABEL_COMMUNITY, levels = c('Inland communities', 'Fishing communities'))]
    
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
  df_age_aggregated[AGE_INFECTION.RECIPIENT < (AGE_TRANSMISSION.SOURCE - 2.5), AGE_CLASSIFICATION.SOURCE := 'Older']
  df_age_aggregated[AGE_INFECTION.RECIPIENT > (AGE_TRANSMISSION.SOURCE + 2.5), AGE_CLASSIFICATION.SOURCE := 'Younger']
  
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

add_infected_unsuppressed <- function(eligible_count_round, proportion_unsuppressed){
  
  # ensure that the data are data.table objects
  di <- as.data.table(eligible_count_round)
  pu <- as.data.table(proportion_unsuppressed)
  
  # use round 15 for round 14 in inland communties for eligible count
  di15 <- di[ROUND == '15' & COMM == 'inland']
  di15[, ROUND := '14']
  di <- rbind(di15, di)
  
  # use round 15 for round 14 in inland communties for proportion of unsuppressed
  pu15 <- pu[ROUND == 'R015'& COMM == 'inland']
  pu15[, ROUND := 'R014']
  pu <- rbind(pu15, pu)
  
  # use round 15 for round 15s in fishing communties for proportion of unsuppressed
  pu15 <- pu[ROUND == 'R015'& COMM == 'fishing']
  pu15[, ROUND := 'R015S']
  pu <- rbind(pu15, pu)
    
  # select variabel
  di <- di[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE, INFECTED, SUSCEPTIBLE)]
  di[, ROUND := paste0('R0', ROUND)]
  
  # find proportion of unsuppressed by round
  df <- merge(di, pu, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # get infected non suppressed
  df[, INFECTED_NON_SUPPRESSED := INFECTED * PROP_UNSUPPRESSED_M]
  df[, INFECTED_NON_SUPPRESSED_CL := INFECTED * PROP_UNSUPPRESSED_CL]
  df[, INFECTED_NON_SUPPRESSED_CU := INFECTED * PROP_UNSUPPRESSED_CU]
  
  # rm unecessary variable
  df <- select(df, -c('PROP_UNSUPPRESSED_M', "PROP_UNSUPPRESSED_CL", "PROP_UNSUPPRESSED_CU"))
  
  return(df)
}

summarise_eligible_count_period <- function(eligible_count_round, cutoff_date, df_period){
  
  # find time intervals
  eligible_count_round <- merge(eligible_count_round, df_round[, .(ROUND, ROUND_SPANYRS, INDEX_TIME, COMM)], by = c('ROUND', 'COMM'))
  eligible_count_round[, BEFORE_CUTOFF := INDEX_TIME == 1]
  
  # summarise across time periods
  df <- eligible_count_round[, .(ROUND, INDEX_TIME, SEX, COMM, AGEYRS, BEFORE_CUTOFF, ROUND_SPANYRS, 
                                 ELIGIBLE, INFECTED, SUSCEPTIBLE, PROP_UNSUPPRESSED_EMPIRICAL, INFECTED_NON_SUPPRESSED, INFECTED_NON_SUPPRESSED_CL, INFECTED_NON_SUPPRESSED_CU)]
  df <- melt.data.table(df, id.vars = c('ROUND', 'SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF', 'INDEX_TIME', 'ROUND_SPANYRS'))

  # check that the lengh corresponds to the one of the period
  tmp <- unique(df[, .(ROUND, INDEX_TIME, ROUND_SPANYRS, COMM)][order(ROUND)])
  tmp <- tmp[, list(PERIOD_SPAN_WITH_ROUND = sum(ROUND_SPANYRS)), by = c('INDEX_TIME', 'COMM')]
  tmp <- merge(tmp, df_period, by = c('INDEX_TIME', 'COMM'))
  stopifnot(tmp[, all(PERIOD_SPAN == PERIOD_SPAN_WITH_ROUND)])
  
  # find weight of every round in the period 
  tmp <- unique(df[, .(ROUND, ROUND_SPANYRS, BEFORE_CUTOFF, COMM)])
  tmp <- tmp[, list(WEIGHT_ROUND = ROUND_SPANYRS / sum(ROUND_SPANYRS),
                    ROUND = ROUND), by = c('BEFORE_CUTOFF', 'COMM')]
  df <- merge(df, tmp, by = c('BEFORE_CUTOFF', 'ROUND', 'COMM'))
  
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
  dfw <- merge(dfw, df_period, by = c('BEFORE_CUTOFF', 'COMM'))
  
  return(dfw)
}

get_incidence_cases_round <- function(incidence.inland, incidence.fishing, eligible_count_round){
  
  # prepare incidence inland
  inc.inland <- copy(incidence.inland)
  colnames(inc.inland) <- toupper(colnames(inc.inland))
  setnames(inc.inland, 'AGE', 'AGEYRS')
  inc.inland[, COMM := 'inland']
  inc.inland[, SEX := substring(SEX, 1, 1)]
  inc.inland <- inc.inland[ROUND >= 14]
  inc.inland[, ROUND := paste0('R0', as.character(ROUND))]
  inc.inland <- inc.inland[, .(SEX, ROUND, AGEYRS, INCIDENCE, LB, UB, COMM)]
  
  # prepare incidence fishing
  inc.fishing <- copy(incidence.fishing)
  colnames(inc.fishing) <- toupper(colnames(inc.fishing))
  setnames(inc.fishing, 'AGE', 'AGEYRS')
  inc.fishing[, COMM := 'fishing']
  inc.fishing[, SEX := substring(SEX, 1, 1)]
  setnames(inc.fishing, 'ROUND_LABEL', 'ROUND')
  
  # combine fishing and inland
  incidence <- rbind(inc.inland, inc.fishing)

  if(0){
    ggplot(incidence, aes(x = AGEYRS)) +
      geom_line(aes(y = INCIDENCE*100, col = ROUND)) +
      geom_ribbon(aes(ymin = LB*100, ymax = UB*100, fill = ROUND),  alpha = 0.1) +
      labs(y = 'Incidence rate per 100 PY ', x = 'Age') +
      facet_grid(COMM~SEX, label = 'label_both', scale = 'free_y') +
      theme_bw() +
      theme(legend.position = 'bottom')
  }
  
  # merge to susceptible
  dir <- merge(incidence, eligible_count_round, by = c('COMM', 'AGEYRS', 'SEX', 'ROUND'))
  
  # find length in years of each round
  dir <- merge(dir, df_round[, .(COMM, ROUND, ROUND_SPANYRS, INDEX_TIME)], by = c("COMM", 'ROUND'))
  
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
                      INCIDENT_CASES_UB = round(sum(INCIDENT_CASES_UB), digits = 1)), by = c('ROUND', 'COMM')]
    knitr::kable(tmp)
    
  }
  
  return(dir)
  
}


summarise_incidence_cases_period <- function(incidence_cases_round, cutoff_date, df_period){
  
  # summarise across time periods
  incidence_cases <- incidence_cases_round[, list(INCIDENT_CASES = sum(INCIDENT_CASES), 
                                                  INCIDENT_CASES_UB = sum(INCIDENT_CASES_UB), 
                                                  INCIDENT_CASES_LB = sum(INCIDENT_CASES_LB)), by = c('COMM', 'AGEYRS', 'SEX', 'INDEX_TIME')]
  # make period
  incidence_cases <- merge(incidence_cases, df_period, by = c('INDEX_TIME', 'COMM'))
  
  return(incidence_cases)
}

inv_logit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1-p))

get_proportion_sampling <- function(pairs, incidence_cases, outdir, 
                                    diff_proportion_sampling_source_age_mirror = F, 
                                    diff_proportion_sampling_source_age_pair = F, 
                                    smooth_proportion_sampling_probability = F){
  
  # find number of pair observed
  dp <- copy(pairs)
  setnames(dp, c('SEX.RECIPIENT', 'COMM.RECIPIENT', 'AGE_INFECTION.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT'), 
           c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF'))
  dp[, AGEYRS := floor(AGEYRS)]
  dp <- dp[, list(count = .N), by = c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF')]
  
  # merge to incidence cases
  di <- merge(incidence_cases, dp, by = c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF'), all.x = T)
  di[is.na(count), count := 0]
  
  # find empiriral proportion sampling from a source of any age to a recipient aged j
  di[, prop_sampling_empirical := count / INCIDENT_CASES]
  
  # smooth proportion of sampling
  if(smooth_proportion_sampling_probability){

    span = 0.2
    span2 = 0.2
    
    tmp <- di[BEFORE_CUTOFF == T]
    tmp <- tmp[, {
      loess.fit = loess(prop_sampling_empirical~AGEYRS, span = span, control=loess.control(surface="direct"))
      prop_sampling <- inv_logit(predict(loess.fit, newdata = AGEYRS))
      prop_sampling <- (predict(loess.fit, newdata = AGEYRS))
      list(prop_sampling_empirical = prop_sampling, AGEYRS = AGEYRS)
    },by = c('SEX', 'COMM', 'BEFORE_CUTOFF', 'PERIOD')]
    
    ds <- di[BEFORE_CUTOFF == F]
    ds <- ds[, {
      loess.fit = loess(prop_sampling_empirical~AGEYRS, span = span2, control=loess.control(surface="direct"))
      prop_sampling <- inv_logit(predict(loess.fit, newdata = AGEYRS))
      prop_sampling <- (predict(loess.fit, newdata = AGEYRS))
      list(prop_sampling_empirical = prop_sampling, AGEYRS = AGEYRS)
    },by = c('SEX', 'COMM', 'BEFORE_CUTOFF', 'PERIOD')]
    
    ds <- rbind(ds, tmp)
    ds[prop_sampling_empirical < 0, prop_sampling_empirical := 0.00]
    ds[prop_sampling_empirical > 1, prop_sampling_empirical := 1.00]
    
  }else{
    ds <- copy(di)
  }
  
  # find proportion of the porbability across age of source 
  if(diff_proportion_sampling_source_age_pair){ 

    # count pairs by age source
    dp <- copy(pairs)
    setnames(dp, c('SEX.SOURCE', 'COMM.RECIPIENT', 'AGE_INFECTION.SOURCE'), c('SEX', 'COMM', 'AGEYRS'))
    dp[, AGEYRS := floor(AGEYRS)]
    dp <- dp[, list(count = .N), by = c('SEX', 'COMM', 'AGEYRS')]
    df <- merge(dp, unique(ds[, .(AGEYRS, COMM, SEX)]), by = c('AGEYRS', 'COMM', 'SEX'), all.y = T)
    df[is.na(count), count := 0.0]
    df[, count :=  as.numeric(count)]

    # smooth count pairs
    df[, count := count + 1]
    df <- df[, {
      loess.fit = loess(count~AGEYRS, span = 0.6, control=loess.control(surface="direct"))
      count <- inv_logit(predict(loess.fit, newdata = AGEYRS))
      count <- (predict(loess.fit, newdata = AGEYRS))
      list(count = count, AGEYRS = AGEYRS)
    },by = c('SEX', 'COMM')]
    df[, min_count := min(count[count > 0]), by  = c('SEX', 'COMM')]
    df[count < 0, count := min_count]

    # find proportion of pairs by source
    df[, N_AGE_GROUP := length(AGEYRS), by  = c('SEX', 'COMM')]
    df[, proportion := count / sum(count), by  = c('SEX', 'COMM')]

    df <- as.data.table(full_join(df, data.table(BEFORE_CUTOFF = c(T,F)), by = character()))
    
    ggplot(df[COMM == 'inland'], aes(x = AGEYRS, y =proportion )) + geom_line() + 
      facet_grid(BEFORE_CUTOFF~SEX, label= 'label_both') + theme_bw() + 
      labs(x= 'Age source', y = 'Distribution of observed pairs')
    ggsave(paste0(outdir, '-data-proportion_sampling_diff_source_age_pair.png'), w = 9, h = 8)
    

  }else if(diff_proportion_sampling_source_age_mirror){
    # re-scale factor is proporitonal to the probability
    
    df <- di[, list(count = sum(count), INCIDENT_CASES = sum(INCIDENT_CASES)), 
             by = c('SEX', 'COMM', 'AGEYRS')]
    
    df[, prop_sampling_empirical := count / INCIDENT_CASES]

    # df[, prop_sampling_empirical := prop_sampling_empirical + 0.009]
    df <- df[, {
      loess.fit = loess(prop_sampling_empirical~AGEYRS, span = 0.8, control=loess.control(surface="direct"))
      prop_sampling_smooth <- inv_logit(predict(loess.fit, newdata = AGEYRS))
      prop_sampling_smooth <- (predict(loess.fit, newdata = AGEYRS))
      list(prop_sampling_empirical = prop_sampling_empirical, prop_sampling_smooth = prop_sampling_smooth, AGEYRS = AGEYRS)
    },by = c('SEX', 'COMM')]
    
    df[, min_prop_sampling_empirical := min(prop_sampling_smooth[prop_sampling_smooth > 0]), by  = c('SEX', 'COMM')]
    df[prop_sampling_smooth < 0,prop_sampling_smooth := min_prop_sampling_empirical]
    
    ggplot(df[COMM == 'inland'], aes(x = AGEYRS)) + 
      geom_line(aes(y =prop_sampling_empirical )) + 
      geom_line(aes(y =prop_sampling_smooth ), col = 'red') + 
      facet_grid(SEX~COMM, label= 'label_both', scale = 'free') + theme_bw() + 
      labs(x= 'Age recipient', y = 'Distribution of probability of observing a transmission across both periods')
    ggsave(paste0(outdir, '-data-proportion_sampling_age_aggregated.png'), w = 7, h = 6)
    
    df[, N_AGE_GROUP := length(AGEYRS), by  = c('SEX', 'COMM')]
    df[, proportion := prop_sampling_smooth / sum(prop_sampling_smooth), by  = c('SEX', 'COMM')]

    df <- as.data.table(full_join(df, data.table(BEFORE_CUTOFF = c(T,F)), by = character()))
    
    ggplot(df[COMM == 'inland' & BEFORE_CUTOFF == T], aes(x = AGEYRS, y =proportion * N_AGE_GROUP)) + geom_line() + 
      facet_grid(SEX~COMM, label= 'label_both', scale = 'free') + theme_bw() + 
      scale_y_continuous(limits = c(0, NA)) + 
      labs(x= 'Age source', y = 'Re-caling factor')
    ggsave(paste0(outdir, '-data-proportion_sampling_recaling_source_aggregated.png'), w = 7, h = 6)
    
    df <- select(df, - prop_sampling_empirical, -prop_sampling_smooth)
    
  }else{
    # re-scale factor does not change with age of source
    
    df <- copy(ds[, .(SEX, COMM, AGEYRS, BEFORE_CUTOFF, prop_sampling_empirical)])
    df[, N_AGE_GROUP := length(AGEYRS), by  = c('SEX', 'COMM', 'BEFORE_CUTOFF')]
    df[, proportion := 1 / N_AGE_GROUP]
    
    df <- select(df, - prop_sampling_empirical)
  }
  
  # find re-scale factor for source
  df[, rescale_factor := proportion * N_AGE_GROUP]
  setnames(df, c('SEX', 'AGEYRS'), c('SEX.SOURCE', 'AGEYRS.SOURCE'))
  
  # merge to recipient
  tmp <- ds[, .(SEX, COMM, AGEYRS, BEFORE_CUTOFF, prop_sampling_empirical)]
  setnames(tmp, c('SEX', 'AGEYRS'), c('SEX.RECIPIENT', 'AGEYRS.RECIPIENT'))
  df <- merge(df, tmp, by = c('BEFORE_CUTOFF', 'COMM'), allow.cartesian = T)
  df[, prop_sampling := prop_sampling_empirical * rescale_factor]
  df <- df[SEX.RECIPIENT != SEX.SOURCE]

  # check
  tmp <- df[, list(sum_prop_sampling = sum(prop_sampling)), by = c('SEX.RECIPIENT', 'AGEYRS.RECIPIENT', 'prop_sampling_empirical', 'COMM', 'BEFORE_CUTOFF', 'N_AGE_GROUP')]
  tmp[, sum_prop_sampling_empirical := prop_sampling_empirical * N_AGE_GROUP]
  stopifnot(nrow(tmp[abs(sum_prop_sampling - sum_prop_sampling_empirical) > 1e-6]) == 0)

  # merge period 
  df <- merge(df, unique(di[, .(PERIOD, BEFORE_CUTOFF, COMM)]), by = c('BEFORE_CUTOFF', 'COMM'))
  
  # warnings
  tmp <- df[prop_sampling > 1]
  if(nrow(tmp) > 0){
    cat('\n Some probabilities are greater than 1')
    cat('\n In', tmp[, unique(COMM)], 'communities at period', tmp[, unique(PERIOD)])
  }
  tmp <- df[prop_sampling < 0]
  if(nrow(tmp) > 0){
    cat('\n Some probabilities are smaller than 0')
    cat('\n In', tmp[, unique(COMM)], 'communities at period', tmp[, unique(PERIOD)])
  }
  
  if(1){ # plots
    
    tmp1 <- copy(df)
    tmp1[, Direction := 'Female -> Male']
    tmp1[SEX.RECIPIENT == 'F', Direction := 'Male -> Female']

    # age recipient plot
    # empirical probabilities
    tmp <- copy(di)
    tmp[, Direction := 'Female -> Male']
    tmp[SEX == 'F', Direction := 'Male -> Female']
    # prepare probabilities 
    tmp2 <- tmp1[, list(prop_sampling = mean(prop_sampling)), by = c('PERIOD', 'COMM', 'AGEYRS.RECIPIENT', 'Direction')]
    
    ggplot(tmp, aes(col = Direction)) + 
      geom_line(aes(x = AGEYRS, y = prop_sampling_empirical)) + 
      geom_line(data = tmp2, aes(x = AGEYRS.RECIPIENT, y = prop_sampling), linetype = 'dashed') + 
      facet_grid(PERIOD~COMM + Direction) + 
      scale_y_continuous(labels = scales::percent)  + 
      theme_bw() +
      labs(y = 'Probability of observing transmission event', x = 'Age recipient')
    ggsave(paste0(outdir, '-data-proportion_sampling_empirical_recipient_period.png'), w = 8, h = 7)
    
    # age source plot
    tmp2 <- tmp1[, list(prop_sampling = mean(prop_sampling)), by = c('PERIOD', 'COMM', 'AGEYRS.SOURCE', 'Direction')]
    ggplot(tmp2, aes(col = Direction)) + 
      geom_line(aes(x = AGEYRS.SOURCE, y = prop_sampling)) + 
      facet_grid(PERIOD~COMM + Direction) + 
      scale_y_continuous(labels = scales::percent)  + 
      theme_bw() +
      labs(y = 'Probability of observing transmission event', x = 'Age source')
    ggsave(paste0(outdir, '-data-proportion_sampling_source_period.png'), w = 8, h = 7)
    
    ggplot(tmp2, aes(col = Direction)) + 
      geom_line(aes(x = AGEYRS.SOURCE, y = prop_sampling)) + 
      facet_grid(PERIOD~COMM + Direction) + 
      scale_y_continuous(labels = scales::percent)  + 
      theme_bw() +
      labs(y = 'Probability of observing transmission event', x = 'Age source')
    
    # age source-recipient prepare pairs for plot
    dp <- pairs[, .(SEX.RECIPIENT, COMM.RECIPIENT, AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT)]
    setnames(dp, c('COMM.RECIPIENT', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT'), 
             c('COMM', 'AGEYRS.SOURCE', 'AGEYRS.RECIPIENT', 'BEFORE_CUTOFF'))
    dp[, `:=` (AGEYRS.RECIPIENT = floor(AGEYRS.RECIPIENT), AGEYRS.SOURCE = floor(AGEYRS.SOURCE))]
    dp <- merge(dp, unique(incidence_cases[, .(BEFORE_CUTOFF, PERIOD, COMM)]), by = c('BEFORE_CUTOFF', 'COMM'))
    dp[, Direction := 'Female -> Male']
    dp[SEX.RECIPIENT == 'F', Direction := 'Male -> Female']
    
    ggplot(tmp1, aes(x = AGEYRS.SOURCE, y = AGEYRS.RECIPIENT)) + 
      geom_raster(aes(fill = prop_sampling)) + 
      facet_grid(PERIOD~COMM+Direction)  +
      scale_fill_viridis_c(labels = scales::percent)  + 
      geom_point(data = dp, col = 'red') + 
      labs(x = 'Age source', y = 'Age recipient', fill = 'Probability of observing\ntransmission event') +
      scale_y_continuous(expand= c(0,0))+
      scale_x_continuous(expand= c(0,0)) +
      theme(strip.background = element_rect(colour="black", fill="white"),
            strip.text = element_text(size = rel(1)))
    ggsave(paste0(outdir, '-data-proportion_sampling_source_recipient_period.png'), w = 9, h = 8)
    
    
  }
  
  df
}

prepare_unsuppressed <- function(eligible_count){
  tmp <- eligible_count[variable == 'INFECTED_NON_SUPPRESSED']
  setnames(tmp, 'AGEYRS', 'AGE_TRANSMISSION.SOURCE')
  tmp[, IS_MF := as.numeric(SEX == 'M')]
  tmp <- merge(tmp, df_direction, by = 'IS_MF')
  tmp <- merge(tmp, df_community, by = 'COMM')
  tmp
}

make.df.round <- function(df_round_inland, df_round_fishing, df_period){
  
  #
  # for inland
  df_round_inland[, INDEX_TIME := 0]
  df_round_inland[round%in%paste0('R0', 14:15), INDEX_TIME := 1]
  df_round_inland[round %in% paste0('R0',16:18), INDEX_TIME := 2]
  
  # fill missing months
  df_round_inland[round == 'R014', max_sample_date := df_round_inland[round == 'R015', min_sample_date]]
  df_round_inland[round == 'R015', max_sample_date := df_round_inland[round == 'R016', min_sample_date]]
  df_round_inland[round == 'R016', max_sample_date := df_round_inland[round == 'R017', min_sample_date]]
  df_round_inland[round == 'R017', max_sample_date := df_round_inland[round == 'R018', min_sample_date]]
  
  #
  # for fishing
  df_round_fishing[, INDEX_TIME := 0]
  df_round_fishing[round%in%paste0('R0',c(15,'15S')), INDEX_TIME := 1]
  df_round_fishing[round %in% paste0('R0',16:18), INDEX_TIME := 2]
  
  # fill missing months
  df_round_fishing[round == 'R015', min_sample_date := start_observational_period_fishing]
  df_round_fishing[round == 'R015', max_sample_date := df_round_fishing[round == 'R015S', min_sample_date]]
  df_round_fishing[round == 'R015S', max_sample_date := cutoff_date]
  df_round_fishing[round == 'R016', min_sample_date := cutoff_date]
  df_round_fishing[round == 'R016', max_sample_date := df_round_fishing[round == 'R017', min_sample_date]]
  df_round_fishing[round == 'R017', max_sample_date := df_round_fishing[round == 'R018', min_sample_date]]
  
  #
  # combine
  df_round_inland[, COMM := 'inland']
  df_round_fishing[, COMM := 'fishing']
  df_round <- rbind(df_round_inland, df_round_fishing)
  
  # keep only round 14 to 18
  df_round <- df_round[INDEX_TIME != '0']
  df_round <- df_round[order(COMM, round)]
  
  # index 
  df_round[, INDEX_ROUND := 1:length(round), by = 'COMM']

  # find length in years of each round and check that it maches the period
  df_round[, ROUND_SPANYRS := .year.diff(max_sample_date, min_sample_date)]
  tmp <- df_round[, list(ROUND_SPANYRS = sum(ROUND_SPANYRS)), by = c('COMM', 'INDEX_TIME')]
  tmp <- merge(tmp, df_period, by = c('COMM', 'INDEX_TIME'))
  stopifnot(tmp[, all(ROUND_SPANYRS == PERIOD_SPAN)])
  
  # round in capital
  colnames(df_round) <- toupper(colnames(df_round))
  df_round[, round := gsub('R0', '', ROUND)]
  df_round[round == '15S', round := '15.1']
  df_round[, round := as.numeric(round)]
  
  # label
  df_round[, MIN_SAMPLE_DATE_LABEL := format(MIN_SAMPLE_DATE, '%b %Y')]
  df_round[, MAX_SAMPLE_DATE_LABEL := format(MAX_SAMPLE_DATE - 31, '%b %Y')]
  df_round[, LABEL_ROUND := paste0('Round ', gsub('R0', '', ROUND), '\n', MIN_SAMPLE_DATE_LABEL, '-', MAX_SAMPLE_DATE_LABEL)]
  df_round[, LABEL_ROUND := factor(LABEL_ROUND, levels = df_round[order(round), LABEL_ROUND])]
  
  return(df_round)
}

find_log_offset_by_round <- function(stan_data, eligible_count_round){
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  eligible_count_wide[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  

  res <- list(); index = 1
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:stan_data[['N_ROUND']]){
        
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



