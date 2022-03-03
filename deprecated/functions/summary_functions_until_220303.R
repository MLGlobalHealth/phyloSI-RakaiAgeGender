.read <- function(x){
  if(grepl('.csv$', x)){return(as.data.table(read.csv(x)))}
  if(grepl('.rds$|.RDS$',x)){return(as.data.table(readRDS(x)))}
}

.aid2pt <- function(x){
  if (length(unique(x)) < length(x)){stop('Error: avoid repeated entries')}
  x <- data.table(AID=x)
  x <- merge(x, anonymisation.keys, all.x=T)
  x$PT_ID
}
.pt2aid <- function(x){
  if (length(unique(x)) < length(x)){stop('Error: avoid repeated entries')}
  x <- data.table(PT_ID=x)
  x <- merge(x, anonymisation.keys, all.x=T)
  x$AID
}

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
  
  cat('There are ', nrow(dchain), ' likely source-recipient pairs')
  
  return(dchain)
}


make.date.first.positive <- function(date.first.positive.and.birthdate)
{
  # set to date format
  date.first.positive.and.birthdate[,  `:=` (visit_dt = as.Date( visit_dt ,'%d/%m/%Y'), 
                                             date_first_positive = as.Date( first_pos_dt ,'%d/%m/%Y') ) ]
  
  # Keep first visit
  tmp <- date.first.positive.and.birthdate[, list(visit_dt = min(visit_dt)), by = 'pt_id']
  date.first.positive.and.birthdate <- merge(date.first.positive.and.birthdate, tmp, by = c('pt_id', 'visit_dt'))
  setnames(date.first.positive.and.birthdate, 'visit_dt', 'date_first_visit')
  
  date.first.positive <- date.first.positive.and.birthdate[, .(pt_id, date_first_positive)]
  
  date.first.positive <- unique(date.first.positive)
  
  return(date.first.positive)
}

make.time.since.infection <- function(time.since.infection, anonymisation.keys)
{

  setnames(time.since.infection, c('visit_dt'), 'date_collection')
  time.since.infection[, date_collection := as.Date(date_collection)]
  
  time.since.infection <- merge(time.since.infection, anonymisation.keys, by='AID')
  
  time.since.infection <- time.since.infection[, .(date_collection, TSI_estimated_mean, TSI_estimated_min, TSI_estimated_max, PT_ID)]
  
  time.since.infection <- unique(time.since.infection)
  
  return(time.since.infection)
}


make.date.birth <- function(date.first.positive.and.birthdate)
{
  # set to date format
  date.first.positive.and.birthdate[,  `:=` (visit_dt = as.Date( visit_dt ,'%d/%m/%Y'), 
                                             date_first_positive = as.Date( first_pos_dt ,'%d/%m/%Y') ) ]
  
  # Keep first visit
  tmp <- date.first.positive.and.birthdate[, list(visit_dt = min(visit_dt)), by = 'pt_id']
  date.first.positive.and.birthdate <- merge(date.first.positive.and.birthdate, tmp, by = c('pt_id', 'visit_dt'))
  setnames(date.first.positive.and.birthdate, 'visit_dt', 'date_first_visit')
  
  # find birth date 
  date.first.positive.and.birthdate[, date_birth := date_first_visit - age_at_visit*365]
  
  date.birth <- date.first.positive.and.birthdate[, .(pt_id, date_birth)]
  
  date.birth <- unique(date.birth)
  
  return(date.birth)
}

get.time.collection <- function(bflocs)
{
  
  # merge estimated TSI
  tmp1 <- merge(tmp1, time.since.infection, by=c('AID', 'PT_ID'))
  tmp1[, visit_dt := as.Date(visit_dt, '%Y-%m-%d')]
  tmp1[, date_first_positive := visit_dt - as.integer(TSI_estimated_mean*365)]
  setcolorder(tmp1, c('PT_ID','AID', 'PANGEA_ID','visit_dt','date_first_positive',
                      'TSI_estimated_mean', 'TSI_estimated_min','TSI_estimated_max'))
  
  return(tmp1)
}


get.meta.data <- function(meta.rccs.1, meta.rccs.2, meta.mrc, date.first.positive, date.birth, 
                          anonymisation.keys, community.keys, use.tsi.estimates = F, time.since.infection = NULL)
{
  
  
  colnames(meta.rccs.1) <- tolower(colnames(meta.rccs.1))
  colnames(meta.rccs.2) <- tolower(colnames(meta.rccs.2))
  colnames(meta.mrc) <- tolower(colnames(meta.mrc))


  #
  # process RCCS meta-data
  #
  
  # first meta data for RCCS to get [community] and [birthdate] 
  meta.rccs.1[, pt_id := paste0('RK-', rccs_studyid)]
  meta.rccs.1[, date := as.Date(date)]
  
  tmp <- meta.rccs.1[, list(date = min(date), max_date = max(date)), by = 'pt_id'] # keep info of first visit
  meta.rccs.1 <- merge(meta.rccs.1, tmp, by = c('pt_id', 'date')) 
  setnames(meta.rccs.1, c('date', 'max_date', 'birthdate'),  c('date_first_visit', 'date_last_visit', 'date_birth'))
  
  date.range.visit <- unique(meta.rccs.1[, .(pt_id, date_first_visit, date_last_visit)]) # keep track of first and last visit
  date.birth.rccs <- unique(meta.rccs.1[, .(pt_id, date_birth)]) # keep track of birth date
  
  meta.rccs.1 <- unique( meta.rccs.1[, .(pt_id, comm_num)] ) 
  
  # second meta data for RCCS to get [sex] 
  meta.rccs.2[, visit_dt := as.Date(visit_dt)]
  meta.rccs.2[age_enrol == 'NULL', age_enrol := NA]
  
  tmp <- meta.rccs.2[, list(visit_dt = min(visit_dt), max_visit_dt = max(visit_dt)), by = 'pt_id'] # keep info of first visit
  meta.rccs.2 <- merge(meta.rccs.2, tmp, by = c('pt_id', 'visit_dt')) 
  setnames(meta.rccs.2, c('visit_dt', 'max_visit_dt'), c('date_first_visit', 'date_last_visit'))
  
  date.range.visit <- unique( rbind(date.range.visit, meta.rccs.2[, .(pt_id, date_first_visit, date_last_visit)]))# keep track of first and last visit
  meta.rccs.2[, date_birth := date_first_visit - as.numeric(age_enrol)*365]
  date.birth.rccs <- unique( rbind(date.birth.rccs, meta.rccs.2[, .(pt_id, date_birth)])) # keep track of birth date
  
  meta.rccs.2 <- unique( meta.rccs.2[, .(pt_id, sex)] )  
  
  # merge two data sources of RCCS 
  meta.rccs <- merge(meta.rccs.2, meta.rccs.1, by = 'pt_id', all.x = T)
  
  # add date first and last visit (here we are missing some without dates)
  tmp <- date.range.visit[, list(date_first_visit = min(date_first_visit), date_last_visit = max(date_last_visit)), by = 'pt_id'] # keep info of first visit
  date.range.visit <- merge(date.range.visit, tmp, by = c('pt_id', 'date_first_visit', 'date_last_visit')) 
  meta.rccs <- merge(meta.rccs, date.range.visit, by = 'pt_id', all.x=T)
  
  # add birthdate 
  date.birth.rccs <- rbind(date.birth, date.birth.rccs)[!is.na(date_birth)]
  tmp <- date.birth.rccs[, list(date_birth = min(date_birth)), by = 'pt_id'] # keep info of first visit
  date.birth.rccs <- unique( merge(date.birth.rccs, tmp, by = c('pt_id', 'date_birth')) )
  meta.rccs <- merge(meta.rccs, date.birth.rccs, by = 'pt_id', all.x = T)
  
  # add cohort round
  meta.rccs[, cohort_round := 'R15-R18']
  
  # add time first positive
  meta.rccs <- merge(meta.rccs, date.first.positive, by = 'pt_id', all.x = T)
  
  stopifnot(length(unique(meta.rccs$pt_id)) == nrow(meta.rccs))
  cat('There is ', nrow(meta.rccs), ' individuals included in the RCCS meta-data\n')
  cat('Out of them ', nrow(meta.rccs[is.na(date_first_positive)]), ' that do not have time of first positive\n')
  cat('Out of them ', nrow(meta.rccs[is.na(date_birth)]), ' that do not have birth date\n')
  
  # add Tanya's estimate time since infection 
  if(use.tsi.estimates){
    colnames(time.since.infection) <- tolower(colnames(time.since.infection))
    meta.rccs <- merge(meta.rccs, time.since.infection[, .(tsi_estimated_mean, date_collection, pt_id)], by = 'pt_id', all.x = T)
    cat('Out of them ', nrow(meta.rccs[is.na(tsi_estimated_mean)]), ' that do not have Tanya s TSI estimates\n')
  }

  
  #
  # process MRC meta-data
  #
  
  meta.mrc <- unique( meta.mrc[, .(pt_id, sex, age_enrol, visit_dt)] ) 
  meta.mrc[, visit_dt := as.Date(visit_dt, '%Y-%m-%d')]
  meta.mrc[age_enrol == 'NULL', age_enrol := NA]
  
  tmp <- meta.mrc[, list(visit_dt = min(visit_dt), max_visit_dt = max(visit_dt)), by = 'pt_id']  # keep first visit
  meta.mrc <- merge(meta.mrc, tmp, by = c('visit_dt', 'pt_id'))
  setnames(meta.mrc, c('visit_dt', 'max_visit_dt'), c('date_first_visit', 'date_last_visit'))
  
  meta.mrc[, date_birth := date_first_visit - as.numeric(age_enrol)*365] # find birth date
  
  meta.mrc <- meta.mrc[, .(pt_id, date_first_visit, date_last_visit, date_birth)]
  
  stopifnot(length(unique(meta.mrc$pt_id)) == nrow(meta.mrc))
  cat('\nThere is ', nrow(meta.mrc), ' individuals included in the MRC meta-data\n')
  cat('Out of them ', nrow(meta.mrc[is.na(date_birth)]), ' that do not have birth date\n')
  
  
  #
  # merge RCCS and MRC meta data
  #
  
  meta.rccs[, cohort := 'RCCS']
  meta.mrc[, cohort := 'MRC']
  meta <- rbind(meta.rccs, meta.mrc, fill = TRUE)
  cat('\nThere is ', nrow(meta), ' individuals included in the meta-data\n')
  
  
  #
  # Find age at first visit and age at infection
  #
  
  meta[, age_at_first_visit := as.numeric(date_first_visit - date_birth)/365]
  
  if(use.tsi.estimates){
    # Find date at infection using Tanya's time to infection estimates
    meta[, date_infection := date_collection - tsi_estimated_mean*365] 
    meta[is.na(date_collection), date_infection := date_last_visit - tsi_estimated_mean*365] 
  }
  
  if(!use.tsi.estimates){
    # Find date at infection using age at first positive, or if NA age at first visit
    meta[, date_infection := date_first_positive ] 
    meta[is.na(date_first_positive), date_infection := date_first_visit] 
  }
  
  meta[, age_infection := as.numeric(date_infection - date_birth)/365]
  
  cat('Out of them ', nrow(meta[!is.na(age_infection)]), ' have a proxy for the age at infection\n')
  
  
  #
  # last changes
  # anonymisation keys
  colnames(anonymisation.keys) <- tolower(colnames(anonymisation.keys))
  meta <- merge(meta, anonymisation.keys, by = 'pt_id')
  
  # community key
  colnames(community.keys) <- tolower(colnames(community.keys))
  community.keys[, comm := ifelse(strsplit(comm_num_a, '')[[1]][1] == 'f', 'fishing', 'island'), by = 'comm_num_a']
  meta <- merge(meta, community.keys, by.x = 'comm_num', by.y = 'comm_num_raw', all.x = T)
  
  # keep only variable of interest
  meta <- meta[, .(pt_id, aid, sex, cohort_round, cohort, comm,
                   date_birth,
                   age_infection, age_at_first_visit,
                   date_infection, date_first_positive, date_first_visit, date_last_visit)]
  
  return(meta)
}

pairs.get.meta.data <- function(dchain, meta){
  
  # stopifnot(unique(dchain$SOURCE) %in% unique(meta_data$aid))
  # stopifnot(unique(dchain$RECIPIENT) %in% unique(meta_data$aid))
  
  # merge by source, then recipient
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.SOURCE')
  tmp1 <- merge(dchain, tmp, by.x = 'SOURCE', by.y = 'aid.SOURCE')
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.RECIPIENT')
  tmp1 <- merge(tmp1, tmp, by.x = 'RECIPIENT', by.y = 'aid.RECIPIENT')
  
  # find age transmission source
  tmp1[, age_transmission.SOURCE := as.numeric(date_infection.RECIPIENT - date_birth.SOURCE)/365]
  
  # individuals without meta data
  missing_indiv = unique(c(dchain$SOURCE, dchain$RECIPIENT)[!(c(dchain$SOURCE, dchain$RECIPIENT) %in% meta$aid)])
  missing_indiv <- .aid2pt(missing_indiv)
  cat('There are ', length(missing_indiv), 'indivs without meta data:\n' )
  missing_indiv <- grep('RK-', missing_indiv, value=T)
  cat('- ', length(missing_indiv), 'of which are in the Rakai Cohort.\n')
  cat(missing_indiv, '\n')
  
  return(tmp1)
}

print.statements.about.pairs <- function(pairs, outdir){
  
  cat('\nThere is ', nrow(pairs), ' source-recipient pairs\n\n')
  
  cat(nrow(pairs[!is.na(age_transmission.SOURCE) & !is.na(age_infection.RECIPIENT)]), ' pairs have a proxy for the age at infection of the source and recipient\n')
  cat(nrow(pairs[((sex.SOURCE == 'F' & sex.RECIPIENT == 'M') | (sex.SOURCE == 'M' & sex.RECIPIENT == 'F')) & (!is.na(age_transmission.SOURCE) & !is.na(age_infection.RECIPIENT))]), ' pairs are heteroxuals have a proxy for the time of infection of the source and recipient\n\n')                
  
  cat('\nPairs by cohort')
  tab <- pairs[, list(count = .N), by = c('cohort.SOURCE', 'cohort.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs enrolled in RCCS by cohort round')
  tab <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS', list(count = .N), by = c('cohort_round.SOURCE', 'cohort_round.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by sex')
  tab <- pairs[, list(count = .N), by = c('sex.SOURCE', 'sex.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by community')
  tab <- pairs[, list(count = .N), by = c('comm.SOURCE', 'comm.RECIPIENT')]
  print_table(tab)
  
}

# Want to see how many AIDs in potential pairs have an associated prefix and base frequency file => SOMETHING
print.statements.about.basefreq.files <- function(chain)
{
  
  stopifnot(file.exists(file.path.phscinput) & file.exists(file.path.bflocs))
  
  # get aids of potential pairs
  aid <- chain[, unique(c(SOURCE, RECIPIENT))]
  stopifnot(all(aid %in% anonymisation.keys$AID))
  
  # Translate AID <-> PREFIXes
  tmp <- file.path(file.path.phscinput)
  tmp <- as.data.table(readRDS(tmp))
  tmp[, `:=` (AID = gsub('-fq[0-9]$','',RENAME_ID), PREFIX=gsub('_remap$','',basename(SAMPLE_ID)))]
  stopifnot(all(aid %in% tmp$AID))
  tmp <- unique(tmp[AID %in% aid, .(AID, PREFIX)])
  
  # Translate PREFIXES <-> bf.csv file existance 
  bflocs <- as.data.table(readRDS(file.path.bflocs)) 
  bflocs <- unique(bflocs[, list(PREFIX=PREFIX, HPC_EXISTS=!is.na(FULL))])
  
  # cat('Out of ', tmp[, length(unique(PREFIX))], ' prefixes associated with our individuals, XX are not in the \n')
  # tmp$PREFIX[which(! tmp$PREFIX %in% bflocs$PREFIX)]
  
  tmp <- merge(tmp, bflocs, by='PREFIX', all.x=T)
  tmp[is.na(HPC_EXISTS), HPC_EXISTS := FALSE]
  
  cat('PREFIXes with existing base frequency file:')
  print_table(tmp[, table(HPC_EXISTS)])
  
  cat('AIDs with at least one existing base frequency file:')
  tmp1 <- tmp[, list(HPC_EXISTS=any(HPC_EXISTS)), by='AID'][, table(HPC_EXISTS)]
  print_table(tmp1)
  
  if(0) # if want to send Tanya
  {
    missing_bff <- copy(tmp)
    missing_bff[, HPC_EXISTS := NULL]
    
    # if want to send Tanya:
    tmp <- missing_bff[,  .(PT_ID, PREFIX)]
    name <- file.path(indir.repository, 'data/missing_bf_files_20220106.csv')
    if(!file.exists(name)){write.csv(tmp, name, row.names = F)}
    # commented out here, but all missing bf's have RCCS2 prefixes.
  }
  
  return(tmp)
}

print_table <- function(table) print(knitr::kable(table))

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
  label_cutoff_date <- format(cutoff_date, '%Y')
  df_time[, label_time := paste0(ifelse(is_before_cutoff_date == 1, 'Before', 'After'), ' ', label_cutoff_date)]
  df_time[, label_time := factor(label_time, levels = paste0(c('Before', 'After'), ' ', label_cutoff_date))]
  
  df_group <- data.table(index_group = 1:(nrow(df_direction) * nrow(df_time)), 
                         index_direction = rep(df_direction$index_direction, each = nrow(df_time)),
                         index_time = rep(df_time$index_time, nrow(df_direction)))
  
  df_group <- merge(df_group, df_direction, by = 'index_direction')
  df_group <- merge(df_group, df_time, by = 'index_time')
  
  setkey(df_group, index_group)
  
  return(df_group)
}

find_incidence_rate <- function(range_age){
  incidence_female <- data.table(age_group = rep(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-49'), each = 11),
                                 year = rep(c(2000, 2001, 2003, 2004, 2005, 2007, 2009, 2010, 2012, 2014, 2015), 6),
                                 incidence = c(0.60, 0.86, 1.90, 0.99, 2.38, 0.97, 0.64, 0.52, 0.68, 0.69, 0.43,
                                               2.19, 2.03, 1.12, 0.95, 1.02, 1.69, 1.30, 2.00, 1.26, 1.55, 1.71, 
                                               1.21, 1.57, 1.41, 1.87, 1.71, 1.41, 1.30, 1.51, 0.79, 1.08, 1.46, 
                                               0.44, 0.51, 1.13, 1.37, 0.51, 1.54, 1.50, 1.13, 0.76, 0.60, 0.56, 
                                               0.84, 1.86, 0.65, 1.45, 1.21, 0.74, 1.38, 0.73, 0.69, 0.41, 0.50, 
                                               0.56, 0.37, 0.88, 0.31, 0.50, 1.01, 0.67, 0.72, 1.19, 0.45, 0.31),
                                 is_mf = 0
  )
  
  incidence_male <- data.table(age_group = rep(c('15-19', '20-24', '25-29', '30-34', '35-39', '40-49'), each = 11),
                               year = rep(c(2000, 2001, 2003, 2004, 2005, 2007, 2009, 2010, 2012, 2014, 2015), 6),
                               incidence = c(0.39, 0.50, NA, NA, NA, 0.47, 0.34, 0.31, 0.39, 0.12, NA,
                                             0.60, 1.39, 1.93, 2.01, 0.58, 0.81, 0.82, 1.36, 0.63, 0.29, 0.44,
                                             1.44, 1.14, 1.26, 1.72, 2.41, 2.15, 1.68, 0.76, 1.70, 1.51, 0.61,
                                             1.94, 0.96, 1.67, 1.83, 0.99, 0.95, 2.21, 1.49, 0.98, 0.70, 0.92,
                                             NA, 3.77, 0.78, 0.32, 0.91, 1.12, 1.03, 1.20, 0.57, 0.22, 0.33, 
                                             NA, 0.63, 0.57, 0.56, 0.81, 0.51, 1.18, 0.66, 0.35, 0.31, 0.34),
                               is_mf = 1
  )
  
  incidence <- rbind(incidence_female, incidence_male)
  
  # extend by age
  incidence[, age_from := gsub('(.+)-.*', '\\1', age_group)]
  incidence[, age_to := gsub('.*-(.+)', '\\1', age_group)]
  
  tmp <- unique(incidence[, .(age_from, age_to)])
  tmp <- tmp[, list(age = age_from:age_to), by = c('age_from', 'age_to')]
  
  incidence <- merge(incidence, tmp, by = c('age_from', 'age_to'), allow.cartesian=TRUE)
  
  # mean incidene before and after cutoff
  incidence[, date := as.Date(paste0(year, '-01-01'))]
  incidence[, is_before_cutoff_date := date < cutoff_date]
  
  if(sum(!incidence$is_before_cutoff_date) == 0){
    tmp <- incidence[year == max(year)]
    tmp[, year := as.numeric(format(cutoff_date, '%Y'))]
    tmp[, date := as.Date(paste0(year, '-01-01'))]
    tmp[, is_before_cutoff_date := date < cutoff_date]
    
    incidence <- rbind(incidence, tmp)
  }
  
  incidence <- incidence[, list(incidence = mean(na.omit(incidence))), by = c('age', 'is_mf', 'is_before_cutoff_date')]
  
  # extending to age fitted
  tmp <- unique(incidence[, .(is_mf,is_before_cutoff_date)])
  tmp[, dummy := 1]
  tmp <- merge(data.table(age = min(range_age):max(range_age), dummy = 1), 
               tmp, allow.cartesian=TRUE)
  incidence <- merge(incidence, tmp, all.y = T, by = c('age', 'is_mf', 'is_before_cutoff_date'))
  
  # filling missing incidence
  setkey(incidence, is_mf, is_before_cutoff_date, age)
  incidence[, incidence := fill_non_na(incidence), by = c('is_mf', 'is_before_cutoff_date')]
  
  
  return(incidence)
}

fill_non_na <- function(x){
  
  index <- which(is.na(x))
  index_non_na <- which(!is.na(x))
  
  if(length(index) == 0){
    return(x)
  }
  
  if(length(index_non_na) == 0){
    return(x)
  }
  
  for(i in index){
    
    index_down = index_non_na[which(index_non_na < i)]
    index_up = index_non_na[which(index_non_na > i)]
    
    if(length(index_down) > 0 & length(index_up) > 0){
      x[i] = mean(c(x[max(index_down)], x[min(index_up)]))
      next
    }
    
    if(length(index_down) > 0){
      x[i] = x[max(index_down)]
      next
    }
    
    if(length(index_up) > 0){
      x[i] = x[min(index_up)]
      next
    }
    
  }
  
  return(x)
}

print.which.NA <- function(dt,regex='SOURCE|RECIPIENT')
{
  cols <- colnames(dt); cols <- grep(regex, cols, value=T)
  .f <- function(x){sum(is.na(x))}
  tmp <- dt[, lapply(.SD, .f), .SDcols=cols]
  cols <- cols[which(tmp[1,] != 0)]
  tmp <- tmp[,
             {
               n <- names(.SD);
               cat('\n-------------------------------------- \nColumns with NA entries : # NA entries \n-------------------------------------- \n')
               lapply(seq_along(.SD),
                      FUN=function(i){ cat(n[[i]], ': ', .SD[[i]], '\n'); 0})
             }
  , .SDcols=cols]
}

resolve.duplicate.ids <- function(discarded_ids)
{
  
  discarded_ids <- discarded_ids[, lapply(.SD, function(x){paste0('RK-', x)})]
  discarded_ids[, aid_discarded := .pt2aid(discarded_ids$pt_id_discarded)]
  discarded_ids[, aid_original := .pt2aid(discarded_ids$pt_id_original)]
  
  .disc2org <- function(x){
    if (length(unique(x)) < length(x)){stop('Error: avoid repeated entries')}
    x <- data.table(pt_id_discarded=x)
    x <- merge(x, discarded_ids, all.x=T)
    x$pt_id_original
  }
  
  # substitute discarded ids with original in meta.rccs.1
  tmp <- discarded_ids[, .(pt_id_discarded, pt_id_original)]
  tmp <- tmp[, lapply(.SD, function(x){gsub('RK-', '', x)} ),]
  tmp <- merge(tmp, meta.rccs.1[RCCS_studyid %in% tmp$pt_id_discarded, ], by.y='RCCS_studyid', by.x='pt_id_discarded')
  setnames(tmp, 'pt_id_original','RCCS_studyid')
  tmp[, pt_id_discarded := NULL]

  meta.rccs.1 <- rbind(meta.rccs.1[!RCCS_studyid %in% tmp$pt_id_discarded, ],
                       tmp)
  setkey(meta.rccs.1, 'RCCS_studyid')
  
  
  # substitute discarded ids with original in meta.rccs.2
  tmp <- meta.rccs.2[pt_id %in% discarded_ids$pt_id_discarded, ]
  setnames(tmp, 'pt_id', 'pt_id_discarded')
  tmp <- merge(discarded_ids[, .(pt_id_discarded,pt_id_original)], tmp, by='pt_id_discarded')
  setnames(tmp, 'pt_id_original', 'pt_id')
  tmp[, pt_id_discarded := NULL] 
  
  meta.rccs.2 <- rbind( tmp, meta.rccs.2[!pt_id %in% discarded_ids$pt_id_discarded, ])
  setkey(meta.rccs.2, pt_id)
  
  # First positive and birthdates
  date.first.positive.and.birthdate[pt_id %in% discarded_ids$pt_id_discarded, pt_id := .disc2org(pt_id) ]
  setkey(date.first.positive.and.birthdate, 'pt_id')
  date.first.positive.and.birthdate[pt_id %in% discarded_ids$pt_id_original,]
  
  date.first.positive <- date.first.positive[! pt_id %in% discarded_ids$pt_id_discarded]
  date.birth <- date.birth[! pt_id %in% discarded_ids$pt_id_discarded]
  
  # time since infection
  stopifnot(time.since.infection[PT_ID %in% discarded_ids$pt_id_discarded, .N == 0])
  
  return(discarded_ids)
}
