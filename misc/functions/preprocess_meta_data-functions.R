print_table <- function(table) print(knitr::kable(table))

.plot.testing.history <- function(DT, col_lab='')
{

  ggplot(DT, aes(y=reorder(study_id, DUMMY)) ) + 
          geom_point(aes(x=date_coll, col=POS), pch='+') + 
          geom_point(aes(x=firstpos_diagnosis_dt), color='red') + 
          geom_point(aes(x=lastnegvd), color='blue') + 
          theme_bw() + 
          theme( axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 legend.position='bottom') +
          labs(x='Testing dates', y='Participants', col=col_lab)
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

make.date.first.positive <- function(allhiv)
{
  
  # Set empty strings to NA (after check)
  .f <- function(x) nchar(as.character(x)) != 9 
  stopifnot(allhiv[.f(date_coll) , .N == 0])
  allhiv[.f(firstpos_diagnosis_dt) , firstpos_diagnosis_dt := NA]
  allhiv[.f(lastnegvd) , lastnegvd := NA]

  cols <- c('visitno','copies', 'new_copies', 'rct', 'lastnegv', 'lastnegvd', 'firstpos_diagnosis_vis')
  allhiv[, (cols):=lapply(.SD, .remove.spaces), .SDcols = cols]
  
  # Format dates
  cols <- c('date_coll', 'firstpos_diagnosis_dt', 'lastnegvd')
  allhiv[, (cols) := lapply(.SD, .dates), .SDcols=cols, by = seq_len(nrow(allhiv))]
  allhiv[,  (cols) := lapply(.SD, function(x){as.Date(x, format='%d-%b-%Y')}), .SDcols=cols]
  
  tmp <- allhiv[is.na(date_coll), unique(study_id)]
  cat("There are ", length(tmp), " (person,visit)s without collection dates\n")
  cat('- all of these had visitno == "R017"\n')
  
  tmp <- allhiv[is.na(firstpos_diagnosis_dt), unique(study_id)]
  cat("There are ", length(tmp), " study_ids without date of first positive\n")
  cat("- these had visit rounds in ", allhiv[study_id %in% tmp, paste0(unique(visitno), collapse = ', ')], '\n')

  if(make.hiv.history.plots)
  {
          library(ggplot2)
          
          # order according to first/last ... dates
          tmp <- unique(allhiv[, .(study_id, firstpos_diagnosis_dt, lastnegvd)])
          setkey(tmp, firstpos_diagnosis_dt, lastnegvd, study_id)
          tmp[, DUMMY := 1:.N]

          tmp1 <- copy(allhiv)
          tmp1 <- merge(allhiv, tmp[, .(study_id, DUMMY)], by='study_id')

          .f <- function(cps)
          {
                  cps <- gsub('<|>|,','', cps)
                  cps[cps %like% '[0-9]']
                  cps[cps %like% '[A-z]'] <- '0'
                  cps <- as.integer(cps)
          }

          # hist(tmp1[.f(copies) <= 100, .f(copies)])
          tmp1[, POS := .f(copies) > 20  ]


          plot.testing.timeline <- .plot.testing.history(tmp1, col_lab='Viremic test')
          plot.testing.timeline.firstpos <- .plot.testing.history(tmp1[! is.na(firstpos_diagnosis_dt)], col_lab='Viremic test')

          filename <- file.path(outdir, 'testing_timeline.png')
          ggsave(plot.testing.timeline, file=filename, w=10, h=15)
          filename <- file.path(outdir, 'testing_timeline_firstpos.png')
          ggsave(plot.testing.timeline.firstpos, file=filename, w=10, h=15)

          # it seems there are some participants with positive tests predating the first pos.
          tmp1[POS==T & date_coll < firstpos_diagnosis_dt,
          {
               cat('There are ', .N, ' visits with positive test preceding the reported first positive,\n')
               cat('accounting for', uniqueN(study_id),'participants. \n')
               cat(sum(visitno %like% '^R0'), 'of these visits were part of the RCCS.\n')
          }]

          # maybe firstpositive diagnosis is not the same as first positive date.
          # it could be that previous blood samples are sequenced...
          idx <- tmp1[POS==T & date_coll < firstpos_diagnosis_dt & visitno %like% '^R0', study_id]

          cat('Update first positive date in allhiv for', length(idx),' participants...\n')
          allhiv[study_id %in% idx, firstpos_diagnosis_dt := {
                  z <- as.integer(copies) > 20
                  z <- !is.na(z) & z == TRUE
                  z <- min(c(date_coll[z], unique(firstpos_diagnosis_dt)) )
                  z
          } ,by='study_id']
  }
  
  # rename
  setnames(allhiv,
           c('firstpos_diagnosis_dt', 'firstpos_diagnosis_vis', 'lastnegvd'),
           c('date_first_positive', 'visit_first_positive', 'date_last_negative'))
  
  # NAs in each col
  print.which.NA(allhiv, 'allhiv')
  
  return(allhiv)
}

process.hiv <- function(hiv)
{
  HIV <- copy(hiv)
  # hiv <- copy(HIV)
  
  # check hidden NAs
  cols <- colnames(hiv)
  hiv[, (cols):=lapply(.SD, .remove.spaces), .SDcols = cols]
  
  # Format dates
  cols <- c('hivdate', 'firstpos_diagnosis_dt', 'lastnegvd')
  hiv[, (cols) := lapply(.SD, .dates), .SDcols=cols, by = c('study_id', 'round')]
  hiv[,  (cols) := lapply(.SD, function(x){as.Date(x, format='%d-%b-%Y')}), .SDcols=cols]
  
  if(make.hiv.history.plots)
  {
          tmp <- hiv[, .(hivdate, firstpos_diagnosis_dt, lastnegvd)]; tmp <- melt(tmp)
          tmp[, range(value, na.rm = T), by='variable']
          ggplot(tmp[value >= '1980-01-01'], aes(fill=variable)) + 
                  geom_histogram(aes(value)) + 
                  theme_bw() + labs(x='date', y='count') 

          tmp <- copy(hiv)
          setnames(tmp, c('hivdate'), c('date_coll'))
          tmp[, POS := hiv == 'P']

          tmp1 <- unique(hiv[, .(study_id, firstpos_diagnosis_dt, lastnegvd, hivdate=min(hivdate)), by='study_id'])
          setkey(tmp1, firstpos_diagnosis_dt, lastnegvd, hivdate, study_id)
          tmp1
          tmp1[, DUMMY := 1:.N] 
          # tmp1[, uniqueN(DUMMY), by='study_id'][, all(V1==1)]
          tmp <- merge(tmp, tmp1[, .(study_id, DUMMY)], by='study_id')

          tmp[, all(is.na(firstpos_diagnosis_dt)), by='study_id'][, mean(V1)]
          idx <- tmp[, uniqueN(DUMMY), by='study_id'][V1 > 1, study_id ]
          tmp[study_id %in% idx]

          plot.testing.timeline.include.negs <- .plot.testing.history(tmp, col_lab='Positive participant')
          filename <- file.path(outdir, 'testing_timeline_with_negatives.png')
          ggsave(plot.testing.timeline.include.negs, file=filename, w=10, h=15)
  }
  
  tmp <- hiv[revertor == 'YES', study_id]
  
  # Last negative diagnosis not defined for X individuals with at least one negative and one positive diagnosis:
  tmp <- hiv[, any(hiv == 'N') & any(hiv == 'P'), by = 'study_id']
  tmp <- tmp[V1 == TRUE, study_id]
  cat(
    'Among ', length(tmp),
    ' individuals with at least one negative and one positive diagnosis:\n'
  )
  tmp <- hiv[study_id %in% tmp,]
  tmp1 <- tmp[is.na(lastnegvd)]

  if(any(tmp1$study_id %in% allhiv$study_id))
          print('allhiv has some extra info on date last neg')

  cat('-\t ', tmp1[, uniqueN(study_id)], ' have undefined last negative visit date\n')

  tmp1 <- tmp1[, ! (min(hiv[hiv =='P']) < max(hiv[hiv=='N'])), by='study_id']
  cat('-\t but', tmp1[, .N], ' have a neg test following a pos\n')

  tmp1 <- tmp[is.na(firstpos_diagnosis_dt)]
  cat('-\t', tmp1[, uniqueN(study_id)], 'participants have NA first pos diag\n')

  cat('Imputing values given visits and status reported...\n\n')
  setkey(tmp, 'study_id', 'hivdate')
  tmp1 <- tmp1[, {
    z.min <- min(which(hiv == 'P'))
    z.max <- max(which(hiv == 'N'))

    list(lastnegvd = hivdate[z.max], firstpos_diagnosis_dt = hivdate[z.min])
  }, by = study_id]
  
  hiv <- merge(tmp1, hiv, by='study_id', all.y=T)
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
  # store date_last_negative to join it at the end
  dlastneg <- unique(allhiv[,.(study_id, date_last_negative)])

  hiv1 <- hiv[!is.na(date_first_positive),.(date_first_positive, visit_first_positive),by=study_id]
  allhiv1 <- allhiv[!is.na(date_first_positive),.(date_first_positive, visit_first_positive),by=study_id]
  hiv1 <- unique(hiv1)
  allhiv1 <- unique(allhiv1)
  
  # check
  .c <- function(DT) DT[, .N, by='study_id'][, all(N == 1)]
  stopifnot( .c(hiv1), .c(allhiv1) )
  
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
  tmp <- merge(tmp, dlastneg, all.x=TRUE, by='study_id')

  # add updated first positive from Joseph
  # add round to new 
  firstpos_update <- merge(firstpos_update,
                           quest[, .(study_id,
                                     visit_first_positive=round,
                                     date_first_positive=intdate)])
  cols <- c("date_first_positive", "date_last_negative")
  firstpos_update[, (cols) := lapply(.SD,as.Date) , .SDcols=cols]
  setcolorder(firstpos_update, names(tmp))

  idx <- firstpos_update[ ! study_id %in% tmp$study_id, study_id]
  cat('Adding', length(idx), "serohistories based on Joseph's update from 2022-11-28..\n")
  tmp <- rbind(firstpos_update, tmp)
  
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
  # meta <- meta[,  comm_num:=comm_num[which(is_latest)], by='study_id' ]
  
  # community key
  colnames(community.keys) <- tolower(colnames(community.keys))
  community.keys[, comm := ifelse(strsplit(as.character(comm_num_a), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'comm_num_a']
  meta <- merge(meta, community.keys[, .(comm_num_raw, comm)], by.x = 'comm_num', by.y = 'comm_num_raw', all.x = T)
  meta[, `:=` (is_latest=NULL, comm_num=NULL)]
  setnames(meta, 'intdate', 'sample_date')
  meta <- unique(meta)
  
  # list rounds
  # TODO: ? We have multiple rows for study id's living in multiple communities.
  # It is also unclear in which comm they were at a spec round.
  # Maybe do something like: comm = inland_fishing?
  # tmp <- meta[, list(round = paste0(round, collapse = "_")), by ='study_id' ]
  # meta <- merge(unique(select(meta, -round)), tmp, by ='study_id')
  #   anyDuplicated(meta, by=c('study_id', 'round'))
  #   stopifnot(meta[, .lu(study_id) == .N])
  
  # get dates of first and last visit
  meta <- merge(meta,date.first.last.visit, by='study_id')
  
  # Get anonymised ID
  meta <- merge(meta, aik, by.x='study_id', by.y='pt_id', all.x=T)
  
  # Add times first positive and last negative
  # stopifnot(meta[!study_id %in%  hiv$study_id, .N==0])
  meta <- merge(meta, date.first.positive[, .(study_id, date_last_negative, date_first_positive)], by='study_id')
  
  # compute ages:
  meta[ , `:=` (age_first_positive = .year.diff(date_first_positive, birthdat))]
  meta[ , `:=` (age_first_visit = .year.diff(date_first_visit, birthdat))]
  
  setnames(meta, c('aid', 'birthdat'), c('aid', 'date_birth'))
  setcolorder(meta, c('study_id', 'aid','sex', 'comm', 'date_birth','round', 'age_first_positive'))
  #   stopifnot(meta[, .lu(study_id) == .N, ])
  meta[, uniqueN(study_id)]
  
  return(meta)
}

process.meta.data <- function(raw_metadata, aik, community.keys){
  
  # make date format
  raw_metadata[, sample_date := as.Date(sample_date, format = '%Y-%m-%d')]
  raw_metadata[, firstposvd := as.Date(firstposvd, format = '%Y-%m-%d')]
  raw_metadata[, lastnegvd := as.Date(lastnegvd, format = '%Y-%m-%d')]
  
  # remove id with no information
  raw_metadata[is.na(sample_date), sample_date := firstposvd]
  raw_metadata <- raw_metadata[!is.na(sample_date)]
  
  # study migrants
  setkey(raw_metadata, study_id, sample_date)
  tmp <- unique(raw_metadata[,.(study_id, community_number)])
  tmp <- tmp[, .N, by= study_id][N > 1, study_id]
  raw_metadata[, is_migrant:=ifelse(study_id %in% tmp, TRUE, FALSE)]
  raw_metadata[, is_latest:=(sample_date == max(sample_date)), by='study_id']
  # raw_metadata <- raw_metadata[is_latest == T]
  
  # find community
  raw_metadata[, comm := 'inland']
  raw_metadata[LakeVictoria_FishingCommunity == T, comm := 'fishing']
  raw_metadata[, `:=` (community_number=NULL)]
  
  # input birthdate
  raw_metadata[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
  
  # get age first positive
  setnames(raw_metadata, c('firstposvd', 'lastnegvd' ), c('date_first_positive', 'date_last_negative'))
  
  # find date first, last visit 
  raw_metadata[, date_first_visit := min(sample_date), by = 'study_id']
  raw_metadata[, date_last_visit := max(sample_date), by = 'study_id']
  # stopifnot(raw_metadata[, length(unique(study_id))] == nrow(raw_metadata))
  
  # find age at first visit and age first positive 
  raw_metadata[, age_first_visit := .year.diff(date_first_visit, date_birth)]
  raw_metadata[, age_first_positive := .year.diff(date_first_positive, date_birth)]
  
  # Get anonymised ID
  raw_metadata <- merge(raw_metadata, aik, by.x='study_id', by.y='pt_id', all.x=T)
  
  # round format
  raw_metadata[, round := paste0('R0', round)]
  raw_metadata[round == 'R015.1', round := 'R015S']
  
  # keep variable of interest
  raw_metadata[, .(study_id, aid, sex, comm, date_birth, round, age_first_positive, is_migrant, 
                   date_first_visit, date_last_visit, date_last_negative, date_first_positive, age_first_visit, 
                   sample_date)]
  
}

process.neuro.meta.data <- function(raw_neuro_metadata, aik){
  
  # make date format
  raw_neuro_metadata[, sampleDate := as.Date(sampleDate, '%m/%d/%Y')]
  raw_neuro_metadata[, date_birth := as.Date(birthdate, '%m/%d/%Y')]
  
  # change name of variable
  setnames(raw_neuro_metadata, 'gender', 'sex')
  
  # find first and last visit
  raw_neuro_metadata[, date_first_visit := min(sampleDate), by = 'study_id']
  raw_neuro_metadata[, date_last_visit := max(sampleDate), by = 'study_id']
  setnames(raw_neuro_metadata, 'sampleDate', 'sample_date')

  # add age at first visit
  raw_neuro_metadata[, age_first_visit := .year.diff(date_first_visit, date_birth)]
  
  # add neuro cohort
  raw_neuro_metadata[, round := 'neuro']
  
  # Get anonymised ID
  neuro_metadata <- merge(raw_neuro_metadata, aik, by.x='study_id', by.y='pt_id', all.x=T)
  
  # keep variable of interest
  neuro_metadata <- neuro_metadata[, .(study_id, aid, sex, date_birth, round, date_first_visit, date_last_visit, age_first_visit)]

  # two rows are duplicated: remove them
  unique(neuro_metadata)
  
  
}

#### MISC ##### 
# Ask Joseph which date_first_positive and date_first_negative should we trust more. 
# Further, why isn t lastnegvd defined in hiv even when there are reported negatives?
# I wrote that after hiv <- process.hiv(hiv) line.
if(0)
{
        tmp <- hiv[hiv=='N', {
                z <- max(hivdate);
                sum(lastnegvd < z)
        }, by='study_id']
        stopifnot(tmp[,  ! any(V1>0, na.rm=T )])
        hiv[hiv=='N', lastnegvd:= max(hivdate), by='study_id']
        hiv[lastnegvd > date_first_positive,]
}
if(0)
{
        tmp <- unique(date.first.positive[, .(study_id, date_first_positive, date_last_negative)])
        tmp1 <- unique(hiv[, .(study_id, date_first_positive, date_last_negative=lastnegvd)])
        tmp[, all(study_id %in% tmp1$study_id)]
        dcomp <- merge(tmp, tmp1, by='study_id')

        # compare date_first_positive
        dcomp[, table(date_first_positive.x == date_first_positive.y, useNA='ifany'),]
        ids_firstpositive_disagree <- dcomp[date_first_positive.x != date_first_positive.y, study_id]
        dcomp[study_id %in% ids_firstpositive_disagree] # earlier first positives are not birthdates
        allhiv[study_id %in% ids_firstpositive_disagree, .(study_id, date_coll, date_first_positive) ] %>% unique
        hiv[study_id %in% ids_firstpositive_disagree, .(study_id, hivdate, hiv,  date_first_positive)] %>% unique
        colnames(hiv)

        # compare date_last_negative
        dcomp[, table(date_last_negative.x == date_last_negative.y, useNA='ifany'),]
        ids_lastnegative_disagree <- dcomp[date_last_negative.x != date_last_negative.y, study_id]
        dcomp[study_id %in% ids_lastnegative_disagree] # earlier first positives are not birthdates
        allhiv[study_id %in% ids_lastnegative_disagree, .(study_id, date_coll, date_last_negative) ] %>% unique
        hiv[study_id %in% ids_lastnegative_disagree, .(study_id, hivdate, hiv,  date_last_negative=lastnegvd)] %>% unique
        colnames(hiv)
        colnames(dcomp) <- gsub('x$','allhiv',colnames(dcomp))
        colnames(dcomp) <- gsub('y$','hiv',colnames(dcomp))
        dcomp <- dcomp[study_id %in% ids_firstpositive_disagree | study_id %in% ids_lastnegative_disagree]
        dcomp[, lastneg_disagree := study_id %in% ids_lastnegative_disagree ]
        dcomp[, firstpos_disagree := study_id %in% ids_firstpositive_disagree ]
        dcomp
        write.csv(dcomp, 'HIV_ALLHIV_firstpos_lastneg_disagree_220401.csv', row.names=F)
}

