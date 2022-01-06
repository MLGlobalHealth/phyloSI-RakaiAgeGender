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

make.time.first.positive <- function(time.first.positive)
{
  time.first.positive[,  `:=` (date_first_visit=as.Date( visit_dt ,'%d/%m/%Y'), date_first_positive=as.Date( first_pos_dt ,'%d/%m/%Y') ) ]
  time.first.positive[, diff := as.numeric(difftime(date_first_positive, date_first_visit, units="weeks")/52.25) ]
  time.first.positive[, age_first_positive := round(age_at_visit + diff, 1)]
  time.first.positive[, `:=` (visit_dt = NULL, first_pos_dt=NULL , diff=NULL)]
  
  # time.first.positive[, date_birth := date_first_visit - age_at_visit * 365]
  # time.first.positive[, age_first_positive := as.numeric(as.Date(date_first_positive, '%d/%m/%Y') - date_birth) / 365]

  # Check age_at_first_pos is coherent
  tmp <- time.first.positive[, max(age_first_positive)-min(age_first_positive), by='pt_id']
  tmp[,range(V1, na.rm=T),]
  time.first.positive[pt_id %in% tmp[V1 > 1, pt_id]]

  # Take median of estimated ages at first positive (given the results are fine)
  time.first.positive[, age_first_positive := round(median(age_first_positive),1) ,by='pt_id']
  time.first.positive <- unique(time.first.positive)

  return(time.first.positive)
}

get.meta.data <- function(meta.rccs.1, meta.rccs.2, meta.mrc, time.first.positive, anonymisation.keys, community.keys){
  
  colnames(meta.rccs.1) <- tolower(colnames(meta.rccs.1))
  colnames(meta.rccs.2) <- tolower(colnames(meta.rccs.2))
  colnames(meta.mrc) <- tolower(colnames(meta.mrc))
  
  #
  # process RCCS meta-data
  
  # process first meta data for RCCS
  meta.rccs.1[, pt_id := paste0('RK-', rccs_studyid)]
  tmp <- meta.rccs.1[, list(date = min(date), comm_num = min(comm_num)), by = 'pt_id'] # keep comm of first visit
  meta.rccs.1 <- unique( merge(meta.rccs.1[, .(pt_id, date, comm_num)], tmp, by = c('pt_id', 'date', 'comm_num')) )
  meta.rccs.1 <- meta.rccs.1[, .(pt_id, comm_num)]
  meta.rccs.1[, cohort_round := 'R15-R16']
  
  # process second meta data for RCCS
  meta.rccs.2 <- unique( meta.rccs.2[, .(pt_id, sex)] ) 
  
  # merge two data sources for RCCS
  meta.rccs <- merge(meta.rccs.2, meta.rccs.1, by = 'pt_id', all.x = T)
  meta.rccs[is.na(cohort_round), cohort_round := 'R17-R18']
  
  # add time first positive
  time.first.positive <- time.first.positive[pt_id %in% unique(meta.rccs$pt_id)]
  tmp <- time.first.positive[, list(date_first_visit=min(date_first_visit)), by = 'pt_id']
  time.first.positive <- merge(time.first.positive, tmp, by = c('date_first_visit', 'pt_id'))
  setnames(time.first.positive, 'age_at_visit', 'age_enrol')
  
  # merge
  meta.rccs <- merge(meta.rccs, time.first.positive[, .(pt_id, 
                                                        age_first_positive, date_first_positive, 
                                                        age_enrol, date_first_visit)], by = 'pt_id')
  meta.rccs <- unique(meta.rccs)

  stopifnot(length(unique(meta.rccs$pt_id)) == nrow(meta.rccs))
  stopifnot(nrow(meta.rccs[is.na(age_first_positive) & is.na(age_enrol)]) > 0)
  cat('There is ', nrow(meta.rccs), ' individuals included in the RCCS meta-data\n')
  
  
  #
  # process MRC meta-data
  meta.mrc <- unique( meta.mrc[, .(pt_id, sex, age_enrol, visit_dt)] ) 
  
  # keep first visit
  meta.mrc[, visit_dt := as.Date(visit_dt, '%Y-%m-%d')]
  tmp <- meta.mrc[, list(visit_dt = min(visit_dt)), by = 'pt_id']
  meta.mrc <- merge(meta.mrc, tmp, by = c('visit_dt', 'pt_id'))
  setnames(meta.mrc, 'visit_dt', 'date_first_visit')
  
  stopifnot(length(unique(meta.mrc$pt_id)) == nrow(meta.mrc))
  cat('There is ', nrow(meta.mrc), ' individuals included in the MRC meta-data\n\n')
  
  #
  # merge RCCS and MRC meta data
  meta.rccs[, cohort := 'RCCS']
  meta.mrc[, cohort := 'MRC']
  meta <- rbind(meta.rccs, meta.mrc, fill = TRUE)
  cat('There is ', nrow(meta), ' individuals included in the meta-data\n')
  
  
  # transform variables
  meta[, date_infection := date_first_positive - 365]
  meta[is.na(date_first_positive), date_infection := date_first_visit - 365]
  
  meta[, age_infection := as.numeric(age_first_positive) - 1]
  meta[is.na(age_first_positive) & !is.na(age_enrol), age_infection := as.numeric(age_enrol) - 1]
  cat('There is ', nrow(meta[!is.na(age_infection)]), ' individuals included in the meta-data with proxy for the age at infection\n')

  
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
                   age_infection,  age_first_positive, age_enrol,
                   date_infection, date_first_positive, date_first_visit)]
  meta[, birth_date := date_infection - age_infection*365]
  
  # remove very young patients (bug?)
  cat('\nExcluding very young patients')
  print_table(meta[age_infection < 11, .(aid, age_infection)])
  meta <- meta[is.na(age_infection) | age_infection >= 11]

  # rm duplicate
  meta <- unique(meta)
  
  return(meta)
}

pairs.get.meta.data <- function(dchain, meta){
  
  # stopifnot(unique(dchain$SOURCE) %in% unique(meta_data$aid))
  # stopifnot(unique(dchain$RECIPIENT) %in% unique(meta_data$aid))
  
  # merge by source
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.SOURCE')
  tmp1 <- merge(dchain, tmp, by.x = 'SOURCE', by.y = 'aid.SOURCE')
  
  # merge by recipient
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.RECIPIENT')
  tmp1 <- merge(tmp1, tmp, by.x = 'RECIPIENT', by.y = 'aid.RECIPIENT')
  
  # find age transmission source
  tmp1[, age_transmission.SOURCE := as.numeric(date_infection.RECIPIENT - birth_date.SOURCE)/365]
  
  # individuals without meta data
  missing_indiv = unique(c(dchain$SOURCE, dchain$RECIPIENT)[!(c(dchain$SOURCE, dchain$RECIPIENT) %in% meta$aid)])
  cat('There are ', length(missing_indiv), 'indivs without meta data:\n' )
  cat(missing_indiv)
  
  return(tmp1)
}

print.statements.about.pairs <- function(pairs, outdir){
  
  cat('\nThere is ', nrow(pairs), ' source-recipient pairs\n\n')
  
  cat(nrow(pairs[!is.na(age_transmission.SOURCE) & !is.na(age_infection.RECIPIENT)]), ' pairs have a proxy for the time of infection of the source and recipient\n')
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

print_table <- function(table) print(knitr::kable(table))

get.age.map <- function(pairs, age_bands_reduced = 4){
  
  ages_source <- pairs[, {
    min_age = floor(min(age_infection.RECIPIENT) - 5)
    max_age = floor(max(age_infection.RECIPIENT) + 5)
    list(age = min_age:max_age)}]
  
  ages_recipient <- pairs[, {
    min_age = floor(min(age_infection.RECIPIENT))
    max_age = floor(max(age_infection.RECIPIENT))
    list(age = min_age:max_age)}]
  
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



