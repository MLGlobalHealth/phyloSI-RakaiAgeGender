print_table <- function(table) print(knitr::kable(table))

process.quest <- function(quest)
{
  
  # one individual with reported birthyr=1882 and birthdate: 2018-08-18 , set year to 1982
  quest[study_id == 'RK-C110340', birthdat:=gsub('-18$','-82',birthdat)]
  
  # Set to date format
  cols <- c('intdate', 'birthdat')
  quest[, (cols) := lapply(.SD, .dates), .SDcols=cols, by = seq_len(nrow(quest))]
  quest[,  (cols) := lapply(.SD, function(x) {
    as.Date(x, format = '%d-%b-%Y')
  }), .SDcols = cols]
  
  # for birthdate NA, find birthdate with age at visit and visit date
  quest[is.na(birthdat), birthdat := median(intdate - 365 * ageyrs), by = 'study_id']
  
  # select columns of interest
  cols <- c('study_id', 'round', 'comm_num', 'curr_id', 'intdate', 'birthdat', 'sex', 'id')
  quest <- quest[, .SD, .SDcols=cols]
  
  cat('Range of interview date is ', format(quest[, range(intdate)]), '\n')
  cat('Range of birth date is ', format(quest[, range(birthdat)]), '\n')
  
  # NAs in each col
  print.which.NA(quest, 'quest')
  
  return(quest)
}

make_date_first_positive <- function(allhiv)
{
  
  # Check whether there are hidden NAs 
  stopifnot(allhiv[nchar(as.character(date_coll)) != 9, .N == 0])
  allhiv[nchar(as.character(firstpos_diagnosis_dt)) != 9, firstpos_diagnosis_dt := NA]
  allhiv[nchar(as.character(lastnegvd)) != 9, lastnegvd := NA]
  cols <- c('visitno','copies', 'new_copies', 'rct', 'lastnegv', 'lastnegvd', 'firstpos_diagnosis_vis')
  allhiv[, (cols):=lapply(.SD, .remove.spaces), .SDcols = cols]
  
  # Format dates
  cols <- c('date_coll', 'firstpos_diagnosis_dt')
  allhiv[, (cols) := lapply(.SD, .dates), .SDcols=cols, by = seq_len(nrow(allhiv))]
  allhiv[,  (cols) := lapply(.SD, function(x){as.Date(x, format='%d-%b-%Y')}), .SDcols=cols]
  
  tmp <- allhiv[is.na(date_coll), unique(study_id)]
  cat("There are ", length(tmp), " (person,visit)s without collection dates\n")
  cat('- all of these had visitno == "R017"\n')
  
  tmp <- allhiv[is.na(firstpos_diagnosis_dt), unique(study_id)]
  cat("There are ", length(tmp), " study_ids without date of first positive\n")
  cat("- these had visit rounds in ", allhiv[study_id %in% tmp, paste0(unique(visitno), collapse = ', ')], '\n')
  
  # rename
  setnames(allhiv, c('firstpos_diagnosis_dt', 'firstpos_diagnosis_vis'), c('date_first_positive', 'visit_first_positive'))
  
  # NAs in each col
  print.which.NA(allhiv, 'allhiv')
  
  return(allhiv)
}

process.hiv <- function(hiv)
{
  HIV <- copy(hiv)
  
  # check hidden NAs
  cols <- colnames(hiv)
  hiv[, (cols):=lapply(.SD, .remove.spaces), .SDcols = cols]
  
  # Format dates
  cols <- c('hivdate', 'firstpos_diagnosis_dt', 'lastnegvd')
  hiv[, (cols) := lapply(.SD, .dates), .SDcols=cols, by = c('study_id', 'round')]
  hiv[,  (cols) := lapply(.SD, function(x){as.Date(x, format='%d-%b-%Y')}), .SDcols=cols]
  
  if(0)
  {
    tmp <- hiv[, .(hivdate, firstpos_diagnosis_dt, lastnegvd)]; tmp <- melt(tmp)
    tmp[, range(value, na.rm = T), by='variable']
    ggplot(tmp[value >= '1980-01-01'], aes(fill=variable)) + 
      geom_histogram(aes(value)) + 
      theme_bw() + labs(x='date', y='count') 
  }
  
  tmp <- hiv[revertor == 'YES', study_id]
  
  # Last negative diagnosis not defined for X individuals with at least one negative and one positive diagnosis:
  tmp <- hiv[, any(hiv == 'N') & any(hiv == 'P'), by = 'study_id']
  tmp <- tmp[V1 == TRUE, study_id]
  cat(
    'Among ',
    length(tmp),
    ' individuals with at least one negative and one positive diagnosis:\n'
  )
  tmp <- hiv[study_id %in% tmp,]
  tmp1 <- tmp[is.na(lastnegvd)]
  cat('- ', tmp1[, .N], ' have undefined last negative visit date\n')
  tmp1 <- tmp[is.na(firstpos_diagnosis_dt)]
  cat('- ', tmp1[, .N], ' have undefined first positive visit date\n')
  
  cat('Imputing values given visits and status reported:\n\n')
  setkey(tmp, 'study_id', 'hivdate')
  tmp <- tmp[, {
    z.min <- min(which(hiv == 'P'))
    z.max <- max(which(hiv == 'N'))
    list(lastnegvd = hivdate[z.max], firstpos_diagnosis_dt = hivdate[z.min])
  }, by = study_id]
  
  hiv <- merge(tmp, hiv, by='study_id', all.y=T)
  hiv[is.na(firstpos_diagnosis_dt.x), firstpos_diagnosis_dt.x := firstpos_diagnosis_dt.y]
  hiv[is.na(lastnegvd.x), lastnegv.x := lastnegvd.y]
  hiv[, `:=` (firstpos_diagnosis_dt=firstpos_diagnosis_dt.x,
              lastnegvd=lastnegvd.x)]
  set(hiv, NULL, grep('.x$|.y$',colnames(hiv),value=T), NULL)
  
  # print statements
  tmp <- hiv[hivdate > firstpos_diagnosis_dt, ]
  cat('Test results for ',tmp[,.N], 'visits after first positive diagnosis:\n')
  print_table(tmp[, table(hiv)])
  
  tmp <- hiv[hivdate < lastnegvd]
  cat('Test results for ',tmp[,.N], 'visits prior to last negative diagnosis:\n')
  print_table(tmp[, table(hiv)])
  
  # rename
  setnames(hiv, c('firstpos_diagnosis_dt', 'firstpos_diagnosis_vis'), c('date_first_positive', 'visit_first_positive'))
  
  return(hiv)
}

compare.hiv.allhiv.firstpositivedates <- function(hiv, allhiv)
{
  hiv1 <- hiv[!is.na(date_first_positive),.(date_first_positive, visit_first_positive),by=study_id]
  allhiv1 <- allhiv[!is.na(date_first_positive),.(date_first_positive, visit_first_positive),by=study_id]
  hiv1 <- unique(hiv1)
  allhiv1 <- unique(allhiv1)
  
  # compare
  cat('\nFIRST POSITIVE ROUND:\n')
  tmp <- merge(hiv1, allhiv1, by='study_id', all.x=T, all.y=T)
  cat('hiv round == allhiv round\n')
  print_table(tmp[, table(visit_first_positive.x == visit_first_positive.y)])
  tmp[,  `:=` (visit_first_positive.y=NULL, visit_first_positive=visit_first_positive.x, visit_first_positive.x=NULL)]
  cat('hiv date == allhiv date\n')
  print_table(tmp[, table(date_first_positive.x == date_first_positive.y)])
  
  tmp1 <- tmp[date_first_positive.x != date_first_positive.y]
  if(0)
  { # based on this plots .y seems to be correct as rounds are separated through time
    ggplot(tmp1[date_first_positive.y >= '2000-01-01']) + 
      geom_histogram(aes(date_first_positive.x), fill='red') +
      geom_histogram(aes(date_first_positive.y), fill='green')
    
    ggplot(tmp1[date_first_positive.y >= '2000-01-01'], aes(fill=visit_first_positive)) + 
      geom_histogram(aes(date_first_positive.x), col='red') +
      geom_histogram(aes(date_first_positive.y), col='green')
  }
  
  # so allhiv1 should be better
  cat('When reported first positive dates are conflicting, decided to select from allhiv1\n')
  tmp[, date_first_positive := date_first_positive.y]
  tmp[is.na(date_first_positive), date_first_positive := date_first_positive.x]
  set(tmp, NULL, c('date_first_positive.x', 'date_first_positive.y'), NULL)
  
  return(tmp)
}

make.date.first.last.visit <- function(hiv)
{
  tmp <- hiv[, list(min=min(hivdate), max=max(hivdate)), by='study_id']
  stopifnot(tmp[, .lu(study_id) == .N ])
  cat(tmp[, sum(! study_id %in% tmp$study_id)], "participants do not have first and last visit dates\n")
  setnames(tmp, c('min', 'max'), c('date_first_visit', 'date_last_visit'))
  return(tmp)
}

get.meta.data <- function(quest, date.first.positive, date.first.last.visit, aik, community.keys)
{
  # Start building metadata
  meta <- copy(quest)
  meta <- unique(meta[, .(study_id, sex, birthdat, comm_num, intdate, round)])
  
  # study migrants
  setkey(meta, study_id, intdate)
  tmp <- unique(meta[,.(study_id, comm_num)])
  tmp <- tmp[, .N, by= study_id][N > 1, study_id]
  meta[, is_migrant:=ifelse(study_id %in% tmp, TRUE, FALSE)]
  meta[, is_latest:=(intdate == max(intdate)), by='study_id']
  meta <- meta[,  comm_num:=comm_num[which(is_latest)], by='study_id' ]
  
  # community key
  colnames(community.keys) <- tolower(colnames(community.keys))
  community.keys[, comm := ifelse(strsplit(as.character(comm_num_a), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'comm_num_a']
  meta <- merge(meta, community.keys[, .(comm_num_raw, comm)], by.x = 'comm_num', by.y = 'comm_num_raw', all.x = T)
  meta[, `:=` (intdate=NULL, is_latest=NULL, comm_num=NULL)]
  meta <- unique(meta)
  
  # list rounds
  tmp <- meta[, list(round = paste0(round, collapse = "_")), by ='study_id' ]
  meta <- merge(unique(select(meta, -round)), tmp, by ='study_id')
  stopifnot(meta[, .lu(study_id) == .N])
  
  # get dates of first and last visit
  meta <- merge(meta,date.first.last.visit, by='study_id')
  
  # Get anonymised ID
  meta <- merge(meta, aik, by.x='study_id', by.y='pt_id', all.x=T)
  
  # Add times first positive
  # stopifnot(meta[!study_id %in%  hiv$study_id, .N==0])
  meta <- merge(meta, date.first.positive[, .(study_id, date_first_positive)], by='study_id')
  
  # compute ages:
  meta[ , `:=` (age_first_positive = .year.diff(date_first_positive, birthdat))]
  meta[ , `:=` (age_first_visit = .year.diff(date_first_visit, birthdat))]
  
  setnames(meta, c('aid', 'birthdat'), c('aid', 'date_birth'))
  setcolorder(meta, c('study_id', 'aid','sex', 'comm', 'date_birth','round', 'age_first_positive'))
  stopifnot(meta[, .lu(study_id) == .N, ])
  
  return(meta)
}

process.meta.data <- function(raw_metadata, aik, community.keys){
  
  # make date format
  raw_metadata[, sample_date := as.Date(sample_date, format = '%Y-%m-%d')]
  raw_metadata[, firstposvd := as.Date(firstposvd, format = '%Y-%m-%d')]
  
  # remove id with no information
  raw_metadata[is.na(sample_date), sample_date := firstposvd]
  raw_metadata <- raw_metadata[!is.na(sample_date)]
  
  # study migrants
  setkey(raw_metadata, study_id, sample_date)
  tmp <- unique(raw_metadata[,.(study_id, community_number)])
  tmp <- tmp[, .N, by= study_id][N > 1, study_id]
  raw_metadata[, is_migrant:=ifelse(study_id %in% tmp, TRUE, FALSE)]
  raw_metadata[, is_latest:=(sample_date == max(sample_date)), by='study_id']
  raw_metadata <- raw_metadata[is_latest == T]
  
  # find community
  raw_metadata[, comm := 'inland']
  raw_metadata[LakeVictoria_FishingCommunity == T, comm := 'fishing']
  raw_metadata[, `:=` (community_number=NULL)]
  
  # input birthdate
  raw_metadata[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
  
  # get age first positive
  setnames(raw_metadata, c('firstposvd'), c('date_first_positive'))
  
  # find date first, last visit and round
  tmp <- raw_metadata[, list(date_first_visit = min(sample_date), 
                             date_last_visit = max(sample_date),
                             round = paste0('R0', round, collapse = "_")), by = 'study_id']
  raw_metadata <- merge(unique(select(raw_metadata, - sample_date, -round)), tmp , by = 'study_id')
  stopifnot(raw_metadata[, length(unique(study_id))] == nrow(raw_metadata))
  
  # find age at first visit and age first positive 
  raw_metadata[, age_first_visit := .year.diff(date_first_visit, date_birth)]
  raw_metadata[, age_first_positive := .year.diff(date_first_positive, date_birth)]
  
  # Get anonymised ID
  raw_metadata <- merge(raw_metadata, aik, by.x='study_id', by.y='pt_id', all.x=T)
  
  # keep variable of interest
  raw_metadata[, .(study_id, aid, sex, comm, date_birth, round, age_first_positive, is_migrant, 
                   date_first_visit, date_last_visit, date_first_positive, age_first_visit)]
  
}

process.neuro.meta.data <- function(raw_neuro_metadata, aik){
  
  # make date format
  raw_neuro_metadata[, sampleDate := as.Date(sampleDate, '%m/%d/%Y')]
  raw_neuro_metadata[, date_birth := as.Date(birthdate, '%m/%d/%Y')]
  
  # change name of variable
  setnames(raw_neuro_metadata, 'gender', 'sex')
  
  # find first and last visit
  tmp <- raw_neuro_metadata[, list(date_first_visit = min(sampleDate), 
                           date_last_visit = max(sampleDate)), by = 'study_id']
  neuro_metadata <- merge(raw_neuro_metadata, tmp, by = 'study_id')
  
  # add age at first visit
  neuro_metadata[, age_first_visit := .year.diff(date_first_visit, date_birth)]
  
  # add neuro cohort
  neuro_metadata[, round := 'neuro']
  
  # Get anonymised ID
  neuro_metadata <- merge(neuro_metadata, aik, by.x='study_id', by.y='pt_id', all.x=T)
  
  # keep variable of interest
  neuro_metadata[, .(study_id, aid, sex, date_birth, round, date_first_visit, date_last_visit, age_first_visit)]
  
  
}

