print_table <- function(table) print(knitr::kable(table))

.lu <- function(x){length(unique(x))}

.read <- function(x){
  if(grepl('.csv$', x)){return(as.data.table(read.csv(x)))}
  if(grepl('.rds$|.RDS$',x)){return(as.data.table(readRDS(x)))}
}

.aid2pt <- function(x, aik){
  if ( .lu(x) < length(x)){stop('Error: avoid repeated entries')}
  colnames(aik) <- tolower(colnames(aik))
  x <- data.table(aid=x)
  x <- merge(x, aik, all.x=T, by='aid')
  x$pt_id
}

.pt2aid <- function(x, aik){
  if (.lu(x) < length(x)){stop('Error: avoid repeated entries')}
  x <- data.table(pt_id=x)
  x <- merge(x, aik, all.x=T, by='pt_id')
  x$aid
}

.remove.spaces <- function(x) {
  x <- gsub(' ', '', x)
  x[which(x=='')] <- NA
  return(x)
}

.dates <- function(x)
{
  tmp <- as.integer(gsub('^[0-9]{2}-[A-Za-z]{3}-([0-9]{2})$', '\\1', x))
  pre <- fifelse(tmp >= 30, '19', '20')
  x <- gsub('-([0-9]{2})$', paste0('-',pre, '\\1'),x)
}

.vars.with.multiple.values <- function(dt, name)
{
  tmp <- dt[, .N, by=name][N>1, ]
  return(tmp)
}

print.which.NA <- function(dt,regex='', name='dt')
{
  cols <- colnames(dt); cols <- grep(regex, cols, value=T)
  .f <- function(x){sum(is.na(x))}
  tmp <- dt[, lapply(.SD, .f), .SDcols=cols]
  cols <- cols[which(tmp[1,] != 0)]
  cat('\nData.Table: ', name , ' | Numbers of rows: ', dt[, .N])
  tmp <- tmp[,
             {
               n <- names(.SD);
               cat('\n-------------------------------------- \nColumns with NA entries : # NA entries \n-------------------------------------- \n')
               lapply(seq_along(.SD),
                      FUN=function(i){ cat(n[[i]], ': ', .SD[[i]], '\n'); 0})
             }
             , .SDcols=cols]
}

.year.diff <- function(x, y)
{
  if (!is.Date(x)) {x <- as.Date(x, format = '%Y-%m-%d')}
  if (!is.Date(y)) {y <- as.Date(y, format = '%Y-%m-%d')}
  round(lubridate::time_length(difftime(x, y),"years"),1)
}

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

make.time.since.infection <- function(time.since.infection)
{
        colnames(time.since.infection)[1:3] <- tolower(colnames(time.since.infection))[1:3]

        setnames(time.since.infection, c('visit_dt'), c('date_collection'))

        time.since.infection[, date_collection := as.Date(date_collection)]
        time.since.infection <- time.since.infection[, .(study_id, date_collection, TSI_estimated_mean, TSI_estimated_min, TSI_estimated_max)]
        time.since.infection <- unique(time.since.infection)

        return(time.since.infection)
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
  colnames(community.keys) <- tolower(colnames(community.keys))
  community.keys[, comm := ifelse(strsplit(as.character(comm_num_a), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'comm_num_a']
  raw_metadata <- merge(raw_metadata, community.keys[, .(comm_num_raw, comm)], by.x = 'community_number', by.y = 'comm_num_raw', all.x = T)
  raw_metadata[, `:=` (community_number=NULL)]
  raw_metadata <- unique(raw_metadata)
  
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

.f <- function(x){round(x, 1)}

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
    tmp <- unique(meta[, .(study_id,age_infection)])
    cat(tmp[age_infection <= 16, .N], "out of", tmp[, .N], "HIV-positive individuals are estimated to have been infected prior to 16 yo\n")
    
    cat(tmp[is.na(age_infection), .N], " HIV-positive individuals do not have an estimate for age at infection\n" )
    
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

get.age.aggregated.map <- function(age_aggregated){
  # age_aggregated <- c('5-14', '15-24', '25-34', '35-54')
  df_age_aggregated <- data.table(expand.grid(age_group_infection.RECIPIENT = age_aggregated, age_group_transmission.SOURCE = age_aggregated))
  df_age_aggregated[, age_from.RECIPIENT := gsub('(.+)-.*', '\\1', age_group_infection.RECIPIENT)]
  df_age_aggregated[, age_from.SOURCE := gsub('(.+)-.*', '\\1', age_group_transmission.SOURCE)]
  df_age_aggregated[, age_to.RECIPIENT := gsub('.*-(.+)', '\\1', age_group_infection.RECIPIENT)]
  df_age_aggregated[, age_to.SOURCE := gsub('.*-(.+)', '\\1', age_group_transmission.SOURCE)]
  tmp <- df_age_aggregated[, list(age_infection.RECIPIENT = unique(age_from.RECIPIENT):unique(age_to.RECIPIENT)), by = c('age_group_infection.RECIPIENT')]
  tmp1 <- df_age_aggregated[, list(age_transmission.SOURCE = unique(age_from.SOURCE):unique(age_to.SOURCE)), by = c('age_group_transmission.SOURCE')]
  df_age_aggregated <- merge(df_age_aggregated, tmp, by = 'age_group_infection.RECIPIENT', allow.cartesian=TRUE)
  df_age_aggregated <- merge(df_age_aggregated, tmp1, by = 'age_group_transmission.SOURCE', allow.cartesian=TRUE)
  
  return(df_age_aggregated)
}

get.round.map <- function(){
  data.table(ROUND = 6:18, DATE_ROUND = as.Date(paste0(2004:2016,'-01-01')))
}
