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
  
  meta[, date_infection := date_first_visit - age_infection]
  
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
  tmp1[, age_transmission.SOURCE := as.numeric(date_infection.RECIPIENT - date_birth.SOURCE)/365]
  
  # individuals without meta data
  missing_indiv = unique(c(chain$SOURCE, chain$RECIPIENT)[!(c(chain$SOURCE, chain$RECIPIENT) %in% meta$aid)])
  missing_indiv <- .aid2pt(missing_indiv, aik)
  cat('There are ', length(missing_indiv), 'indivs without meta data:\n' )
  missing_indiv <- grep('RK-', missing_indiv, value=T)
  cat('- ', length(missing_indiv), 'of which are in the Rakai Cohort.\n')
  cat(missing_indiv, '\n')
  
  return(tmp1)
}

print.statements.about.pairs <- function(pairs){
  
  cat('\nThere is ', nrow(pairs), ' source-recipient pairs\n\n')
  
  cat(nrow(pairs[!is.na(age_transmission.SOURCE) & !is.na(age_infection.RECIPIENT)]), ' pairs have a proxy for the age at infection of the source and recipient\n')
  cat(nrow(pairs[((sex.SOURCE == 'F' & sex.RECIPIENT == 'M') | (sex.SOURCE == 'M' & sex.RECIPIENT == 'F')) & (!is.na(age_transmission.SOURCE) & !is.na(age_infection.RECIPIENT))]), ' pairs are heteroxuals have a proxy for the time of infection of the source and recipient\n\n')                
  
  #   cat('\nPairs by cohort')
  #   tab <- pairs[, list(count = .N), by = c('cohort.SOURCE', 'cohort.RECIPIENT')]
  #   print_table(tab)
  
  #   cat('\nPairs enrolled in RCCS by cohort round')
  #   tab <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS', list(count = .N), by = c('cohort_round.SOURCE', 'cohort_round.RECIPIENT')]
  #   print_table(tab)
  
  cat('\nPairs by sex')
  tab <- pairs[, list(count = .N), by = c('sex.SOURCE', 'sex.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by community')
  tab <- pairs[, list(count = .N), by = c('comm.SOURCE', 'comm.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by round')
  tab <- pairs[, list(count = .N), by = c('first_round.SOURCE')]
  print_table(tab[order(first_round.SOURCE)])
  
  tab <- pairs[, list(count = .N), by = c('first_round.RECIPIENT')]
  print_table(tab[order(first_round.RECIPIENT)])
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
  
  extended_age_length <- 5
  
  ages_source <- pairs[, {
    min_age = floor(min(c(age_transmission.SOURCE,age_infection.RECIPIENT))) - extended_age_length
    max_age = floor(max(c(age_transmission.SOURCE,age_infection.RECIPIENT))) + extended_age_length
    list(age = min_age:max_age)}]
  
  ages_recipient <- ages_source
  
  age_map <- data.table(expand.grid(age_transmission.SOURCE = ages_source$age, 
                                    age_infection.RECIPIENT = ages_recipient$age))
  df_age <- age_map[order(age_transmission.SOURCE, age_infection.RECIPIENT)]

  ages <- sort(unique(df_age$age_transmission.SOURCE))
  ages <- data.table(age_infection = ages, 
                     age_transmission_reduced.SOURCE = rep(seq(min(ages), max(ages), age_bands_reduced), each = age_bands_reduced )[1:length(ages)])
  df_age <- merge(df_age, ages, by.x = 'age_transmission.SOURCE', by.y = 'age_infection')
  
  
  ages <- sort(unique(df_age$age_infection.RECIPIENT))
  ages <- data.table(age_infection = ages, 
                     age_infection_reduced.RECIPIENT = rep(seq(min(ages), max(ages), age_bands_reduced), each = age_bands_reduced )[1:length(ages)])
  df_age <- merge(df_age, ages, by.x = 'age_infection.RECIPIENT', by.y = 'age_infection')

  setkey(df_age, age_transmission.SOURCE, age_infection.RECIPIENT)
  
  df_age[, index_age := 1:nrow(df_age)]
  
  range_age_non_extended <<- c(min(df_age$age_infection.RECIPIENT) + extended_age_length, 
                               max(df_age$age_infection.RECIPIENT) - extended_age_length)
  
  return(df_age)
}

get.group.map <- function(){
  
  df_direction <- data.table(index_direction = 1:2, is_mf = c(0, 1))
  df_direction[, label_direction := ifelse(is_mf == 1, 'Male -> Female', 'Female -> Male')]
  
  df_time <- data.table(index_time = 1:2, is_before_cutoff_date = c(1, 0))
  label_cutoff_date <- as.numeric(format(cutoff_date, '%Y'))
  label_start_observational_period <- format(start_observational_period, '%Y')
  df_time[, label_time := ifelse(is_before_cutoff_date == 1, paste0(label_start_observational_period,'-',label_cutoff_date-1), paste0(label_cutoff_date,'-2019'))]
  df_time[, label_time := factor(label_time, levels = c(paste0(label_start_observational_period,'-',label_cutoff_date-1), 
                                                        paste0(label_cutoff_date,'-2019')))]
  
  df_community <- data.table(index_community = 1:2, comm = c('fishing','inland'))
  df_community[, label_community := ifelse(comm == 'inland', 'Inland communities', 'Fishing communities')]
  
  df_group <- data.table(expand.grid(index_direction = df_direction[, index_direction],
                         index_time = df_time[, index_time], 
                         index_community = df_community[, index_community]))
  
  df_group[, index_group := 1:nrow(df_group)]
  
  df_group <- merge(df_group, df_direction, by = 'index_direction')
  df_group <- merge(df_group, df_time, by = 'index_time')
  df_group <- merge(df_group, df_community, by = 'index_community')
  
  setkey(df_group, index_group)
  
  return(df_group)
}

get.age.aggregated.map <- function(age_aggregated, incidence){
  

  df_age_aggregated <- data.table(expand.grid(age_group_infection.RECIPIENT = age_aggregated, age_group_transmission.SOURCE = age_aggregated))
  df_age_aggregated[, age_from.RECIPIENT := gsub('(.+)-.*', '\\1', age_group_infection.RECIPIENT)]
  df_age_aggregated[, age_from.SOURCE := gsub('(.+)-.*', '\\1', age_group_transmission.SOURCE)]
  df_age_aggregated[, age_to.RECIPIENT := gsub('.*-(.+)', '\\1', age_group_infection.RECIPIENT)]
  df_age_aggregated[, age_to.SOURCE := gsub('.*-(.+)', '\\1', age_group_transmission.SOURCE)]
  
  tmp <- df_age_aggregated[, list(age_infection.RECIPIENT = unique(age_from.RECIPIENT):unique(age_to.RECIPIENT)), by = c('age_group_infection.RECIPIENT')]
  tmp1 <- df_age_aggregated[, list(age_transmission.SOURCE = unique(age_from.SOURCE):unique(age_to.SOURCE)), by = c('age_group_transmission.SOURCE')]
  
  stopifnot(all(tmp1[, age_transmission.SOURCE] %in% incidence$AGEYRS))
  
  df_age_aggregated <- merge(df_age_aggregated, tmp, by = 'age_group_infection.RECIPIENT', allow.cartesian=TRUE)
  df_age_aggregated <- merge(df_age_aggregated, tmp1, by = 'age_group_transmission.SOURCE', allow.cartesian=TRUE)
  
  return(df_age_aggregated)
}

process.incidence <- function(incidence, df_round){

  di <- as.data.table(merge(incidence, df_round, by.x = 'ROUND', by.y = 'round'))
  
  di[, months_diff := .month.diff(max_sample_date, min_sample_date)]
  
  di <- di[min_sample_date >= start_observational_period]
  di[, is_before_cutoff_date := ifelse(max_sample_date < cutoff_date, 1, 0)]
  di[, is_mf := ifelse(SEX == 'F', 1, 0)]
  di[, comm := COMM]
  
  di
}


