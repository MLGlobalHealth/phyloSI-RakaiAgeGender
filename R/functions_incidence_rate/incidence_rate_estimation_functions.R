make_df_round <- function(rounds_group_1, rounds_group_2, rounds_group_3){
  
  df_round <- data.table(visit = c(rounds_group_1, rounds_group_2, rounds_group_3)) 
  
  df_round <- subset(df_round, visit != 'R015S')
  
  df_round[, round_numeric := 1:nrow(df_round)]
  df_round[, ROUND := paste0('Round: ', gsub('R0(.+)', '\\1', visit ))]
  
  df_round[, min_sample_date := as.Date(c("1999-04-06", "2000-03-20", "2001-04-03", "2002-07-15", "2003-09-26", "2005-02-15", "2006-08-30", "2008-06-17",
                                          "2010-01-18", "2011-08-10", "2013-07-08", "2015-02-23", "2016-10-03", "2018-10-03"))]
  df_round[, max_sample_date :=as.Date(c("2000-02-14", "2001-03-26", "2002-05-31", "2003-08-01", "2004-11-23", "2006-06-30", "2008-06-06", "2009-12-07",
                                         "2011-06-21", "2013-07-05", "2015-01-30", "2016-09-02", "2018-05-22", "2020-05-22"))]
  
  df_round[, midpoint_date := min_sample_date + as.numeric(max_sample_date - min_sample_date)/2]
  
  return(df_round)
}

read_hiv_data <- function(file.path.hiv.614, file.path.hiv.1518, file.path.hiv.19){
  
  # round 6 to 14
  hivstatus_vlcopies_1 <- as.data.table(read_dta(file.path.hiv.614))
  
  # round 15 to 18
  hivstatus_vlcopies_1.1518<-as.data.table(read.csv(file.path.hiv.1518))
  
  # round 19
  hivstatus_vlcopies_1.19<-as.data.table(read.csv(file.path.hiv.19))
  
  # clean round 6 to 14
  setnames(hivstatus_vlcopies_1,  'intdate','hivdate')
  hivstatus_vlcopies_1[, hivdate := as.character(hivdate)]
  hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[!round %in% rounds_group_2]
  hivstatus_vlcopies_1 <- select(hivstatus_vlcopies_1, c('study_id', 'round', 'hivdate', 'hiv'))
  
  # combine round 15 to 19
  hivstatus_vlcopies_1.1519 <- rbind(hivstatus_vlcopies_1.1518, hivstatus_vlcopies_1.19)
  hivstatus_vlcopies_1.1519 <- select(hivstatus_vlcopies_1.1519, c('study_id', 'round', 'hivdate', 'hiv'))
  
  # combine all
  hivstatus_vlcopies_1 <- rbind(hivstatus_vlcopies_1, hivstatus_vlcopies_1.1519)
  hivstatus_vlcopies_1[, round := gsub(' ', '', round)] # remove space after string
  
  # add last update on infected by Joseph
  infected_id <- c("C117824", "C119303", "K067249")
  hivstatus_vlcopies_1[study_id %in% infected_id, hiv := 'P']
  hivstatus_vlcopies_1[study_id %in% infected_id, hivstatus := T]
  
  return(hivstatus_vlcopies_1)
}

read_hiv_data_230218 <- function(file.path.hiv.68, file.path.hiv.914, file.path.hiv.1518, file.path.hiv.19){
  
  # round 6 to 8
  hivstatus_vlcopies_1 <- as.data.table(read_dta(file.path.hiv.68))
  
  # round 9 to 14
  hivstatus_vlcopies_1.914 <- as.data.table(read.csv(file.path.hiv.914))
  
  # round 15 to 18
  hivstatus_vlcopies_1.1518<-as.data.table(read.csv(file.path.hiv.1518))
  
  # round 19
  hivstatus_vlcopies_1.19<-as.data.table(read.csv(file.path.hiv.19))
  
  # clean round 6 to 8
  setnames(hivstatus_vlcopies_1,  'intdate','hivdate')
  hivstatus_vlcopies_1[, hivdate := as.character(hivdate)]
  hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[round %in% rounds_group_0]
  hivstatus_vlcopies_1 <- select(hivstatus_vlcopies_1, c('study_id', 'round', 'hivdate', 'hiv'))
  
  # clean round 9 to 14 and combine to round 6 to 8
  hivstatus_vlcopies_1.914[, hivdate := as.character(hivdate)]
  hivstatus_vlcopies_1.914 <- hivstatus_vlcopies_1.914[!round %in% rounds_group_2]
  hivstatus_vlcopies_1.914 <- select(hivstatus_vlcopies_1.914, c('study_id', 'round', 'hivdate', 'hiv'))
  hivstatus_vlcopies_1 <- rbind(hivstatus_vlcopies_1, hivstatus_vlcopies_1.914)
  
  # combine round 15 to 19
  hivstatus_vlcopies_1.1519 <- rbind(hivstatus_vlcopies_1.1518, hivstatus_vlcopies_1.19)
  hivstatus_vlcopies_1.1519 <- select(hivstatus_vlcopies_1.1519, c('study_id', 'round', 'hivdate', 'hiv'))
  
  # combine all
  hivstatus_vlcopies_1 <- rbind(hivstatus_vlcopies_1, hivstatus_vlcopies_1.1519)
  hivstatus_vlcopies_1[, round := gsub(' ', '', round)] # remove space after string
  
  return(hivstatus_vlcopies_1)
}

read_flow_data <- function(file.path.flow.614, file.path.flow.1518, file.path.flow.19){
  
  # round 6 to 14
  verif_1<-as.data.table(read_dta(file.path.flow.614))
  
  # round 15 to 18
  verif_1.1518<-as.data.table(read.csv(file.path.flow.1518))
  
  # round 19
  verif_1.19<-as.data.table(read.csv(file.path.flow.19))
  
  # clean round 6 to 14
  verif_1[, locdate := as.character(locdate)]
  verif_1[, birthdat := as.character(birthdat)]
  verif_1 <- verif_1[!round %in% rounds_group_2]
  setnames(verif_1, 'locdate', 'loc_date')
  verif_1 <- select(verif_1, c('study_id', 'round', 'loc_date', 'locdate1', 'birthdat', 'sex', 
                               'comm_num', 'ageyrs', 'resident'))
  
  # clean round 15 to 18
  verif_1.1518 <- select(verif_1.1518,-'X_merge')

  # combine  and clean round 15 to 19
  verif_1.1519 <- rbind(verif_1.1518, verif_1.19)
  verif_1.1519[, loc_date := as.character(loc_date)]
  verif_1.1519[, birthdat := as.character(birthdat)]
  verif_1.1519 <- select(verif_1.1519, c('study_id', 'round', 'loc_date', 'locdate1', 'birthdat', 'sex', 
                                         'comm_num', 'ageyrs', 'resident'))

  # combine all
  verif_1 <- rbind(verif_1, verif_1.1519)
  
  return(verif_1)
}

read_flow_data_230218 <- function(file.path.flow.68, file.path.flow.914, file.path.flow.1518, file.path.flow.19){
  
  # round 6 to 8
  verif_1<-as.data.table(read_dta(file.path.flow.68))
  
  # round 9 to 14
  verif_1.914<-as.data.table(read.csv(file.path.flow.914))
  
  # round 15 to 18
  verif_1.1518<-as.data.table(read.csv(file.path.flow.1518))
  
  # round 19
  verif_1.19<-as.data.table(read.csv(file.path.flow.19))
  
  # clean round 6 to 8
  verif_1[, locdate := as.character(locdate)]
  verif_1[, birthdat := as.character(birthdat)]
  verif_1 <- verif_1[round %in% rounds_group_0]
  setnames(verif_1, 'locdate', 'loc_date')
  verif_1 <- select(verif_1, c('study_id', 'round', 'loc_date', 'locdate1', 'birthdat', 'sex', 
                               'comm_num', 'ageyrs', 'resident'))
  
  # clean round 9 to 14
  verif_1.914[, locdate := as.character(locdate)]
  verif_1.914[, birthdat := NA_character_]
  setnames(verif_1.914, 'comm', 'comm_num')
  verif_1.914 <- verif_1.914[!round %in% rounds_group_2]
  setnames(verif_1.914, 'locdate', 'loc_date')
  verif_1.914 <- select(verif_1.914, c('study_id', 'round', 'loc_date', 'locdate1', 'birthdat', 'sex', 
                               'comm_num', 'ageyrs', 'resident'))
  
  # combine round 6 to 14
  verif_1 <- rbind(verif_1, verif_1.914)
  
  # clean round 15 to 18
  verif_1.1518 <- select(verif_1.1518,-'X_merge')
  
  # combine  and clean round 15 to 19
  verif_1.1519 <- rbind(verif_1.1518, verif_1.19)
  verif_1.1519[, loc_date := as.character(loc_date)]
  verif_1.1519[, birthdat := as.character(birthdat)]
  verif_1.1519 <- select(verif_1.1519, c('study_id', 'round', 'loc_date', 'locdate1', 'birthdat', 'sex', 
                                         'comm_num', 'ageyrs', 'resident'))
  
  # combine all
  verif_1 <- rbind(verif_1, verif_1.1519)
  
  return(verif_1)
}

format_hiv_data <- function(hivstatus_vlcopies_raw){
  
  hivstatus_vlcopies_1 <- copy(hivstatus_vlcopies_raw)
  
  #code visit dates
  hivstatus_vlcopies_1[, hivdate2 := as.Date('1900-01-01')] 
  hivstatus_vlcopies_1[round %in% rounds_group_1, hivdate2 := as.Date(hivdate)]
  hivstatus_vlcopies_1[!round %in% rounds_group_1, hivdate2 := as.Date(hivdate, format='%d-%B-%y')]
  hivstatus_vlcopies_1 <- select(hivstatus_vlcopies_1, -'hivdate')
  setnames(hivstatus_vlcopies_1, 'hivdate2', 'hivdate')
  
  # remove individual with hiv test date in 1977
  hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[!(hivdate == "1977-01-01")]
  range(hivstatus_vlcopies_1$hivdate)
  hivstatus_vlcopies_1[, table(format(hivdate, '%Y'),round)]
  
  #code HIV status
  hivstatus_vlcopies_1[, hivstatus:=NA]
  hivstatus_vlcopies_1[hiv=="N",hivstatus:=0]
  hivstatus_vlcopies_1[hiv=="P",hivstatus:=1]
  hivstatus_vlcopies_1[, table(hivstatus)]
  
  #Generate a numeric variable for round   
  hivstatus_vlcopies_1 <- merge(df_round, hivstatus_vlcopies_1, by.x = 'visit', by.y = 'round')
  setnames(hivstatus_vlcopies_1, 'round_numeric', 'round')
  hivstatus_vlcopies_1[, round := as.numeric(round)]
  hivstatus_vlcopies_1[, table(round)]
  
  #Rename study id as research id per Adam's original code
  hivstatus_vlcopies_1[, research_id:=study_id]
  
  #Sort HIV data by round 
  hivstatus_vlcopies_1 <- hivstatus_vlcopies_1[order(research_id, round)] 
  
  #Generate a variable to indicate year of hiv trst 
  hivstatus_vlcopies_1[, year := format(hivdate, format="%Y")]
  hivstatus_vlcopies_1[year < 1999, hivdate := NA]
  hivstatus_vlcopies_1[year < 1999, year := NA]
  hivstatus_vlcopies_1[, table(round,year)]
  
  # Create a visit id
  hivstatus_vlcopies_1[, visit_id:=paste(research_id, round, sep=":")]
  
  return(hivstatus_vlcopies_1)
}

format_flow_data <- function(verif_raw){
  
  verif_1 <- copy(verif_raw)
  
  #Generate a numeric variable for round   
  verif_1 <- merge(df_round, verif_1, by.x = 'visit', by.y = 'round')
  setnames(verif_1, 'round_numeric','round')
  verif_1[, round:=as.numeric(round)]
  verif_1[, table(round)]
  
  #Rename study id as research id per Adam's original code
  verif_1[, research_id:=study_id]
  
  # Create a visit id
  verif_1[, visit_id:=paste(research_id, round, sep = ":")]
  
  # find date of visit
  verif_1[,loc_date2 := as.Date('1900-01-01')] #create another variable with date format 
  verif_1[round %in% rounds_numeric_group_1, loc_date2 := as.Date(loc_date, format='%d/%m/%Y')]
  verif_1[round %in% rounds_numeric_group_2, loc_date2 := as.Date(loc_date, format='%d-%B-%y')]
  verif_1[round %in% rounds_numeric_group_3, loc_date2 := as.Date(loc_date, format='%d%b%Y')]
  verif_1 <- select(verif_1, -'loc_date')
  setnames(verif_1, 'loc_date2', 'loc_date')
  verif_1[, range(na.omit(loc_date))]
  
  ## fix obvious mistake in date of visit 
  verif_1[loc_date == '2818-10-17', loc_date := as.Date('2018-10-17')]
  verif_1[loc_date == '2220-01-21', loc_date := as.Date('2020-01-21')]
  verif_1[loc_date == '2201-08-21', loc_date := as.Date('2021-08-21')]
  verif_1[loc_date == '2048-10-16', loc_date := as.Date('2018-10-16')]
  verif_1[loc_date == '2047-02-17', loc_date := as.Date('2017-02-17')]
  verif_1[loc_date == '1971-03-27', loc_date := as.Date('2017-03-27')]
  
  # remove those for whom I cannot figure out the date of visit
  verif_1 <- verif_1[!(loc_date == '2001-04-04' & visit == 'R018')] # 2001 in round 18
  
  # check
  verif_1[, range(na.omit(loc_date))]
  verif_1[, table(format(loc_date, '%Y'), round)]
  
  # now use locdate1 for the na
  verif_1[is.na(loc_date), loc_date := as.Date(locdate1, '%d/%m/%Y')]
  
  # check
  verif_1[, range(loc_date)]
  
  # find birthdate 
  ## initiate variable for birthdays
  verif_1[, birthdat2 := as.Date('1900-01-01')] 
  
  ## use birthdate for those with non-missing birthday
  verif_1[round %in% rounds_numeric_group_1, birthdat2 := as.Date(birthdat, format='%Y-%m-%d')]
  verif_1[round %in% rounds_numeric_group_2, birthdat2 := as.Date(birthdat, format='%d-%B-%y')]
  verif_1[round %in% rounds_numeric_group_3, birthdat2 := as.Date(birthdat, format='%d%b%Y')]
  
  ## use age at visit for one individual with incorrect birthdate 
  verif_1[birthdat == '18-Aug-18', birthdat2:= as.Date(loc_date - ageyrs*365)]
  
  ## year in two digits, so have to change 20 to 19 for crazy years, i.e. 2064, to 1964
  verif_1[!is.na(birthdat2) & birthdat2 > as.Date('2006',format = '%Y'), range(birthdat2)]
  verif_1[!is.na(birthdat2) & birthdat2 > as.Date('2006',format = '%Y'), birthdat2 := as.Date(gsub('20', '19',birthdat2))]
  verif_1[!is.na(birthdat2), table(format(na.omit(birthdat2), '%Y'), round)]
  
  # check
  verif_1[, range(na.omit(birthdat2))]
  
  # use age at visit if not missing for missing birthday
  verif_1[(birthdat == '' | is.na(birthdat)) & ageyrs != 98, birthdat2 := as.Date(loc_date - ageyrs*365)]
  verif_1[, range(na.omit(birthdat2))]
  
  # check that the na are missing birthday + missing age
  verif_1[is.na(verif_1$birthdat2), table(birthdat)]
  verif_1[is.na(verif_1$birthdat2), table(ageyrs)]
  
  # remove individuals with missing birthdate
  verif_1 <- verif_1[!is.na(verif_1$birthdat2)]
  
  # for one individual born on february 29, change the date to February 28 to have valid birthdays every year
  verif_1[birthdat2 == '1984-02-29', birthdat2 := as.Date('1984-02-28')]
  verif_1[birthdat2 == '1968-02-29', birthdat2 := as.Date('1968-02-28')]
  verif_1[birthdat2 == '1988-02-29', birthdat2 := as.Date('1988-02-28')]
  
  # set name to birthdat and check no missing
  verif_1 <- select(verif_1, -'birthdat')
  setnames(verif_1, 'birthdat2', 'birthdat')
  verif_1[, range(birthdat)]
  verif_1[, sum(is.na(birthdat))]
  
  return(verif_1)
}

prepare_hiv_status <- function(hivstatus_vlcopies_1, verif_1){
  
  # keep non-na hiv test
  hivstatus_vlcopies_1_complete <- hivstatus_vlcopies_1[hiv != '']
  
  # add missing rounds between first and last round enrolled 
  hivstatus_vlcopies_1_complete <- hivstatus_vlcopies_1_complete %>%
    group_by(research_id) %>%
    complete(round = min(round):max(round)) %>%
    as.data.table()
  
  # to impute year and date at visit in missing round, use median 
  tmp <- hivstatus_vlcopies_1[, list(year.median = as.character(round(median(na.omit(as.numeric(year))))), 
                                     hivdate.median = median(na.omit(hivdate))), by = c('round')]
  hivstatus_vlcopies_1_complete <- as.data.table(merge(tmp, hivstatus_vlcopies_1_complete, by = 'round'))
  hivstatus_vlcopies_1_complete[is.na(year), year := year.median]
  hivstatus_vlcopies_1_complete[is.na(hivdate), hivdate := hivdate.median]
  
  # append sex and birthdate 
  df_birthdate = verif_1[, list(birthdate=min(birthdat), 
                                sex = min(sex)), by = 'research_id'] 
  hivstatus_vlcopies_1_complete <- merge(hivstatus_vlcopies_1_complete, df_birthdate, by=c("research_id"), all.x = T) 
  
  # append birthday for every year (remember that the analysis is stratified by age)
  df_birthdate = expand.grid("year"=as.character(seq(1999,2020,1)), research_id=unique(hivstatus_vlcopies_1$research_id)) %>%
    left_join(df_birthdate, by="research_id") %>%
    as.data.table()
  df_birthdate[, birthday := as.Date(paste0(year, '-', gsub('.*?-(.+)', '\\1', as.character(birthdate))))]
  
  tmp <- as.data.table(hivstatus_vlcopies_1_complete)[, list(min_date = min(hivdate)), by = c('round', 'research_id')]
  tmp <- tmp[order(research_id, round)]
  tmp[, number_visit := length(min_date), by = c('research_id')]
  tmp[number_visit == 1, max_date := min_date + 365, by = c('research_id')]
  tmp[number_visit > 1, max_date := c(min_date[2:length(min_date)], max(min_date) +365) , by = c('research_id')]
  tmp <- as.data.table(merge(df_birthdate, tmp, by = c('research_id'), allow.cartesian=TRUE))
  tmp <- tmp[(birthday >= min_date & birthday <= max_date)]
  tmp <- select(tmp, c('research_id', 'round', 'birthday'))
  hivstatus_vlcopies_1_fac = merge(hivstatus_vlcopies_1_complete, tmp, by=c("research_id","round"), all.x = T) 
  
  stopifnot(unique(hivstatus_vlcopies_1_complete[, .(research_id, round)])[order(research_id, round)] == 
              unique(hivstatus_vlcopies_1_fac[, .(research_id, round)])[order(research_id, round)])
  
  # extend hiv status up and down in time
  hivstatus_vlcopies_1_fac <- hivstatus_vlcopies_1_fac %>% 
    mutate(hivstatus_impute_up = hivstatus, hivstatus_impute_down = hivstatus) %>%
    fill(hivstatus_impute_down, .direction = "down") %>% 
    fill(hivstatus_impute_up, .direction = "up") %>% 
    as.data.table()
  
  # impute hiv status with negative if in between two negative visit
  # impute hiv status with poistive if in between two positive visit, else NA
  hivstatus_vlcopies_1_fac[, hivstatus_imputed := NA_integer_]
  hivstatus_vlcopies_1_fac[hivstatus_impute_down==hivstatus_impute_up, hivstatus_imputed := hivstatus_impute_down]
  hivstatus_vlcopies_1_fac <- select(hivstatus_vlcopies_1_fac, -c('hivstatus_impute_down', 'hivstatus_impute_up'))

  # re-format such that each round has a row for birthdays and for hiv test 
  hivstatus_vlcopies_1_fac <- hivstatus_vlcopies_1_fac %>%
    pivot_longer(cols=c("birthday","hivdate"), names_to = "datetype", values_to = "date") %>%
    as.data.table()
  hivstatus_vlcopies_1_fac <- unique(hivstatus_vlcopies_1_fac[!is.na(date)])
  hivstatus_vlcopies_1_fac[, round := floor(na.approx(round, na.rm=FALSE)), by = c('research_id')]
  hivstatus_vlcopies_1_fac <- hivstatus_vlcopies_1_fac[!is.na(round)]
  hivstatus_vlcopies_1_fac <- hivstatus_vlcopies_1_fac[order(research_id, date, datetype, round)]

  # start with first hiv date not first birthday (e.g., A045032)
  tmp <- hivstatus_vlcopies_1_fac[, list(first_date_birthday = datetype[1] == 'birthday', 
                                         min_date = min(date)), by = 'research_id']
  hivstatus_vlcopies_1_inc <- merge(hivstatus_vlcopies_1_fac, tmp, by = c('research_id')) 
  hivstatus_vlcopies_1_inc <- hivstatus_vlcopies_1_inc[!(first_date_birthday == T & date== min_date& datetype == 'birthday')]
  set(hivstatus_vlcopies_1_inc, NULL, 'first_date_birthday', NULL)
  hivstatus_vlcopies_1_inc <- hivstatus_vlcopies_1_inc[order(research_id, date, rev(datetype), round)]
  
  # end with last hiv date not last birthday 
  tmp <- hivstatus_vlcopies_1_inc[, list(last_date_birthday = datetype[length(datetype)] == 'birthday', 
                                         max_date = max(date)), by = 'research_id']
  hivstatus_vlcopies_1_inc <- merge(hivstatus_vlcopies_1_inc, tmp, by ='research_id')
  hivstatus_vlcopies_1_inc <- hivstatus_vlcopies_1_inc[!(last_date_birthday == T & date==max_date & datetype == 'birthday')]
  set(hivstatus_vlcopies_1_inc, NULL, 'last_date_birthday', NULL)
  hivstatus_vlcopies_1_inc <- hivstatus_vlcopies_1_inc[order(research_id, date, round)]
  
  # add number of missing visits and seroconvert status
  tmp <- hivstatus_vlcopies_1_inc[which(hivstatus_vlcopies_1_inc$datetype == 'hivdate')] 
  tmp <- tmp[, list(number_missing_visits = sum(is.na(hivstatus_imputed)), 
            seroconvert = all(0:1 %in% hivstatus)), by = 'research_id']
  hivstatus_vlcopies_1_inc <- merge(tmp, hivstatus_vlcopies_1_inc, by = 'research_id')
  
  # as data table
  hivstatus_vlcopies_1_inc <- as.data.table(hivstatus_vlcopies_1_inc)
  
  return(hivstatus_vlcopies_1_inc)
}

find_seroconvert_status <- function(hivstatus_vlcopies_1_inc){
  
  # Create dataset to indicate first neg, last neg, and first positive to generate the cohort
  firstpos <- hivstatus_vlcopies_1_inc %>% 
    group_by(research_id) %>% 
    filter(datetype == 'hivdate') %>% 
    arrange(research_id, round) %>%
    slice(match(1, hivstatus_imputed)) %>% # keep first positive test
    select(research_id, round, date) %>%
    mutate(first_pos=date, first_pos_round = round) %>% 
    select(-date,-round)
  
  lastneg <- hivstatus_vlcopies_1_inc %>% 
    group_by(research_id) %>% 
    filter(datetype == 'hivdate') %>% 
    arrange(research_id, desc(round)) %>%
    slice(match(0, hivstatus_imputed)) %>% # keep last negative test
    select(research_id, round, date) %>%
    mutate(last_neg=date, last_neg_round =round) %>% 
    select(-date,-round) %>%
    full_join(firstpos, by="research_id")
  
  status_df <- hivstatus_vlcopies_1_inc %>% 
    group_by(research_id) %>% 
    filter(datetype == 'hivdate') %>% 
    arrange(research_id, round) %>%
    slice(match(0, hivstatus_imputed)) %>%
    select(research_id, round, date) %>%
    mutate(first_neg=date) %>% select(-date,-round) %>%
    full_join(lastneg, by="research_id")
  nrow(status_df) #46512 unique id's
  
  # find seroconversion status
  status_df <- status_df %>% 
    mutate(cohortclass=ifelse(first_neg!=last_neg & is.na(first_neg)==F & is.na(first_pos)==T, "Serially negative", 
                              ifelse(first_neg==last_neg & is.na(first_neg)==F & is.na(first_pos)==T, "One test negative",
                                     ifelse(is.na(first_neg)==T & is.na(first_pos)==F, "One test positive",
                                            ifelse(is.na(first_neg)==F & is.na(first_pos)==F, "Seroconversion","No Testing")))))
  table(status_df$cohortclass)
  summary(status_df$last_neg)
  summary(status_df$first_pos)
  
  # detect first pos is always after the last neg and the last neg always after the first neg
  status_df$date_error= 0
  status_df$date_error[status_df$first_pos <= status_df$last_neg | status_df$last_neg < status_df$first_neg] = 1
  table(status_df$date_error) # dates with errors
  status_df[which(status_df$date_error == 1),]
  status_df <- subset(status_df, date_error==0)
  
  # only include those with follow-up test after a negative and no date error
  status_df <- subset(status_df, (cohortclass=="Serially negative" | cohortclass== "Seroconversion")) 
  status_df <- status_df[order(status_df$research_id, status_df$first_neg),] #sort by id and follow-up
  
  # as data table
  status_df <- as.data.table(status_df)
  
  return(status_df)
}


find_seroconvert_cohort <- function(status_df, hivstatus_vlcopies_1_inc){
  
  ##############################
  
  # FORMAT SEROCONVERT STATUS #
  
  ##############################
  
  # start
  join_df <- copy(status_df)
  
  # set seroconversion date randomly between the last negative + 1 and the first positive date
  join_df[, seroconv_date := last_neg + (first_pos - last_neg)*runif(nrow(join_df))]
  
  # set seroconversion round as the round corresponding to the year of the seroconversion for each participant
  seroconvdate <- data.table(seroconv_date = join_df$seroconv_date, research_id = join_df$research_id) 
  tmp <- as.data.table(hivstatus_vlcopies_1_inc)[datetype == 'hivdate' & seroconvert == T, list(min_date = min(date)), by = c('round', 'research_id')]
  tmp <- tmp[order(research_id, round)]
  tmp[, max_date := c(min_date[2:length(min_date)] - 1, as.Date('2020-12-31')), by = c('research_id')]
  seroconvdate <- merge(seroconvdate, tmp, by = c('research_id'))
  seroconvdate <- seroconvdate[(seroconv_date >= min_date & seroconv_date <= max_date)]
  setnames(seroconvdate, 'round', 'seroconv_round')
  stopifnot(sum(!is.na(seroconvdate$seroconv_round)) == sum(!is.na(join_df$seroconv_date))) # sanity check 
  table(seroconvdate$seroconv_round)
  join_df <- merge(join_df, seroconvdate[, .(research_id, seroconv_round)], by = c('research_id'), all.x = T)

  # sanity checks: no seroconv date for serial negative
  stopifnot(all(which(is.na(join_df$seroconv_date)) == which(join_df$cohortclass == "Serially negative") )) 
  
  # sanity checks: non-missing seroconv date for seroconverted
  stopifnot(all(which(!is.na(join_df$seroconv_date)) == which(join_df$cohortclass == "Seroconversion")))
  
  
  ##############################
  
  # CREATE SEROCONVERT COHORT #
  
  ##############################
  
  # joint seroconversion status and date to meta data 
  hiv_m <- merge(hivstatus_vlcopies_1_inc, join_df, all.y=T, by="research_id") 
  hiv_m <- hiv_m[order(hiv_m$research_id, hiv_m$date),] #sort by id and follow-up
  
  #replace HIV date as seroconversion date for those who serconvert
  hiv_m[first_pos==date & is.na(first_pos)==F & is.na(seroconv_date)==F & is.na(date)==F, date := seroconv_date]

  #replace date type as seroconversion for the date of serconversion
  hiv_m[seroconv_date==date & is.na(first_pos)==F & is.na(seroconv_date)==F & is.na(date)==F& datetype!='birthday', datetype := "Seroconversion"]
  hiv_m_d <- hiv_m[order(research_id, date)] #resort by id and follow-up
  
  length(unique(hiv_m_d$research_id)) #24558 unique individuals
  
  # replace round as seroconversion round at the date of seroconversion
  hiv_m_d[datetype== "Seroconversion", round := seroconv_round]
  
  # remove birthday if they fall on the same time as seroconversion
  hiv_m_d <- hiv_m_d[is.na(seroconv_date) | !(datetype == 'birthday' & date == seroconv_date)]
  
  # find person year by round
  seroconverter_cohort <- hiv_m_d[order(research_id, round, date)]
  seroconverter_cohort <- seroconverter_cohort[!(cohortclass == 'Seroconversion' & date > seroconv_date)]# remove all dates after seroconversion
  seroconverter_cohort <- seroconverter_cohort[, next_date:= c(date[-1], NA), by = 'research_id']
  seroconverter_cohort <- seroconverter_cohort[, hivinc:=ifelse(is.na(seroconv_date), 0, 
                                                                ifelse(next_date == seroconv_date, 1, 0)), by = 'research_id']
  seroconverter_cohort <- seroconverter_cohort[!(!is.na(next_date) & date > next_date)]
  seroconverter_cohort[, episode_start:=date]
  seroconverter_cohort[, episode_end:=next_date]
  seroconverter_cohort[, py:=as.numeric(episode_end - episode_start) / 365.25]
  
  #age based on birthday rounded down
  seroconverter_cohort[, age := floor(as.numeric(episode_start-birthdate)/365)]
  seroconverter_cohort <- seroconverter_cohort[age >= 15 & age <= 49]

  # sanity check: one participant cannot spend more than 1 year per age of life
  tmp <- seroconverter_cohort[,list(py = sum(na.omit(py))), by = c('age', 'research_id')]
  tmp[, range(py)]# errors with substracting date
  
  # sanity check: no negative py
  nrow(seroconverter_cohort[py < 0])
  
  # remove na person year
  seroconverter_cohort <- seroconverter_cohort[!is.na(py)]
  
  return(seroconverter_cohort)
}

find_incidence_rate_by_sex <- function(modelpreds.age.1218.all, seroconverter_cohort_imputation){
  
  iters <- modelpreds.age.1218.all[, sort(unique(iterations))]
  seroconverter_cohort_imputation[, Sex := ifelse(sex == 'F', 'Female', 'Male')]
  
  model_pred.list <- vector(mode = 'list', length = length(iters))
  set.seed(12)
  for(i in iters){
    
    cat('\n iterations', i)
    
    # predictions
    tmp <- modelpreds.age.1218.all[iterations == i]
    
    # merge to data
    tmp <- merge(tmp, seroconverter_cohort_imputation, by = c('round', 'Sex', 'age', 'iterations'))
    
    # find incidence rate
    tmp[, total_py := sum(pysum), by = c('round', 'Sex', 'iterations', 'iterations_within')]
    model_pred.list[[i]] <- tmp[, list(inc = sum(inc*pysum / total_py)), by = c('round', 'Sex', 'iterations', 'iterations_within')]
    
  }
  model_pred <- rbindlist(model_pred.list)
  
  return(model_pred)
}

find_incidence_rate_all <- function(modelpreds.age.1218.all, seroconverter_cohort_imputation){
  
  iters <- modelpreds.age.1218.all[, sort(unique(iterations))]
  seroconverter_cohort_imputation[, Sex := ifelse(sex == 'F', 'Female', 'Male')]
  
  model_pred.list <- vector(mode = 'list', length = length(iters))
  set.seed(12)
  for(i in iters){
    
    cat('\n iterations', i)
    
    # predictions
    tmp <- modelpreds.age.1218.all[iterations == i]
    
    # merge to data
    tmp <- merge(tmp, seroconverter_cohort_imputation, by = c('round', 'Sex', 'age', 'iterations'))
    
    # find incidence rate
    tmp[, total_py := sum(pysum), by = c('round', 'iterations', 'iterations_within')]
    model_pred.list[[i]] <- tmp[, list(inc = sum(inc*pysum / total_py)), by = c('round', 'iterations', 'iterations_within')]
    
  }
  model_pred <- rbindlist(model_pred.list)
  
  return(model_pred)
}

find_incidence_rate_by_age_group <- function(modelpreds.age.1218.all, seroconverter_cohort_imputation){
  
  iters <- modelpreds.age.1218.all[, sort(unique(iterations))]
  seroconverter_cohort_imputation[, Sex := ifelse(sex == 'F', 'Female', 'Male')]
  
  df_age_aggregated <- data.table(age = 15:49, age_group = c(rep('15-24', 10),  rep("25-34", 10), rep("35-49", 15)))
  
  model_pred.list <- vector(mode = 'list', length = length(iters))
  set.seed(12)
  for(i in iters){
    
    cat('\n iterations', i)
    
    # predictions
    tmp <- modelpreds.age.1218.all[iterations == i]
    
    # merge to data
    tmp <- merge(tmp, seroconverter_cohort_imputation, by = c('round', 'Sex', 'age', 'iterations'))
    
    # merge to age group
    tmp <- merge(tmp, df_age_aggregated, by = 'age')
    
    # find incidence rate
    tmp[, total_py := sum(pysum), by = c('round', 'age_group', 'Sex', 'iterations', 'iterations_within')]
    model_pred.list[[i]] <- tmp[, list(inc = sum(inc*pysum / total_py)), by = c('round', 'age_group', 'Sex', 'iterations', 'iterations_within')]
    
  }
  model_pred <- rbindlist(model_pred.list)
  
  return(model_pred)
}

save_statistics_incidence_cohort <- function(hivstatus_vlcopies_1_inc, status_df){
  
  stats <- list()
  comma_thousands <- function(x) format(x, big.mark=",")
  
  # number of participants
  RESEARCH_IDS <- unique(as.data.table(hivstatus_vlcopies_1_inc)$research_id)
  stats$N = comma_thousands(length(RESEARCH_IDS))
  
  # number of participants followed up in round 10 to 18
  verified_ids_r1018 <- unique(gsub('(.+)\\:.*', '\\1', verified_visit_ids_r1018))
  RESEARCH_IDS1018 <- RESEARCH_IDS[RESEARCH_IDS %in% verified_ids_r1018]
  stats$N_FOLLOWUP <- comma_thousands(as.data.table(status_df)[research_id %in% RESEARCH_IDS1018 & (cohortclass=="Serially negative" | cohortclass== "Seroconversion"), length(unique(research_id))])

  # proportion of seroconvert participants with missed visits
  df_missed_visits <- unique(select(hivstatus_vlcopies_1_inc, c('research_id', 'number_missing_visits')))
  tmp <- merge(status_df, df_missed_visits, by = 'research_id')
  tmp$number_missing_visits_status = ifelse(tmp$number_missing_visits > 1, 'More than 1 missed visits', 
                                            ifelse(tmp$number_missing_visits == 1, '1 missed visit', 'No missed visits'))
  stats$TABLE_STATUS_SEROCONVERTED = as.data.table(tmp)[research_id %in% RESEARCH_IDS & cohortclass== "Seroconversion", table(number_missing_visits_status)]
  stats$TABLE_STATUS_SEROCONVERTED_PROP <- round(stats$TABLE_STATUS_SEROCONVERTED / sum(stats$TABLE_STATUS_SEROCONVERTED) * 100, 1)
  stats$TABLE_STATUS_SEROCONVERTED <- comma_thousands(stats$TABLE_STATUS_SEROCONVERTED)
  
  # subset serially negative or seroconvert 
  part <- as.data.table(hivstatus_vlcopies_1_inc)
  part[, age := floor(as.numeric(date-birthdate)/365.25)]
  part[age < 15, age := 15]; part[age > 49, age := 49] # discrepancies in age
  part <- unique(part[, list(age = min(age), sex = unique(sex), 
                             hivstatus_imputed = min(hivstatus_imputed),
                             missed_visit = is.na(visit_id)), by = c('research_id', 'round')])
  .research_id <- as.data.table(status_df)[(cohortclass=="Serially negative" | cohortclass== "Seroconversion"), research_id]
  part <- part[research_id %in% .research_id] # keep only serially negative or seroconvert
  
  # keep rounds before and on seroconversion
  df_ne <- unique(part[, .(research_id, hivstatus_imputed, round)])
  df_ne[, index_round_after_sero := 0] # round after seroconversion
  df_ne[hivstatus_imputed == 1, index_round_after_sero := 1:length(hivstatus_imputed), by = 'research_id']
  df_ne[, index_round := 1:length(round), by = 'research_id'] # index round observed
  part <- merge(part, df_ne, by = c('research_id', 'hivstatus_imputed', 'round'))
  part2 <- part[index_round_after_sero <= 1]    # keep rounds before and on seroconversion
  
  # number of followed-up participants by age group, sex and round
  df_age_aggregated <- data.table(age = 15:49, age_group = c(rep('15-24', 10),  rep("25-34", 10), rep("35-49", 15)))
  part2 <- merge(part2, df_age_aggregated, by = 'age')
  part_all <- part2[, list(N = length(unique(research_id))), by = c('round', 'age_group', 'sex')]
  part_s <- part2[, list(N = length(unique(research_id))), by = c('round', 'sex')]
  part_r <- part2[, list(N = length(unique(research_id))), by = c('round')]
  
  part_s[, age_group := 'Total']
  part_r[, age_group := 'Total']
  part_r[, sex := 'Total']
  part_all <- do.call('rbind', list(part_all, part_s, part_r))
  part_all[, N := comma_thousands(N)]
  
  # merge
  tab <- copy(part_all)
  
  tab[, age_group := factor(age_group, levels = c('Total', unique(df_age_aggregated$age_group)))]
  tab[, sex := factor(sex, levels = c('Total', 'F', 'M'))]
  tab <- tab[order(round, sex, age_group)]
  tab <- tab[round %in%  df_round[visit >= 'R010' & visit <= 'R018', round_numeric]]
  
  stats[['table']] <- tab
  
  return(stats)
}

save_statistics_estimates <- function(df_round, 
                            seroconverter_cohort_total, seroconverter_cohort_agg, 
                            seroconverter_cohort_agg_s, seroconverter_cohort_agg_r,
                            modelpreds.agegroup.1218, modelpreds){
  
  stats <- list()
  comma_thousands <- function(x) format(x, big.mark=",")
  
  # find average duration of rounds
  tmp <- subset(df_round, round_numeric %in% c(5:9, rounds_numeric_group_2))
  stats[['AVERAGE_MONTHS_ROUND']] <- round(mean(as.numeric(c(tmp$max_sample_date - tmp$min_sample_date)) / 366 * 12))
  
  # number of incidence infection
  stats$TABLE_STATUS = comma_thousands(round(seroconverter_cohort_total[, hivinc]))
  
  # total py
  stats[['PY']] = comma_thousands(round(seroconverter_cohort_total[, .(py, CL_py, CU_py)]))
  
  # py and number of incidence by age group, sex and round
  seroconverter_cohort_agg <- seroconverter_cohort_agg[round %in% c(5:9, rounds_numeric_group_2)]
  seroconverter_cohort_agg_s <- seroconverter_cohort_agg_s[round %in% c(5:9, rounds_numeric_group_2)]
  seroconverter_cohort_agg_r <- seroconverter_cohort_agg_r[round %in% c(5:9, rounds_numeric_group_2)]
  
  seroconverter_cohort_agg_s[, age_group := 'Total']
  seroconverter_cohort_agg_r[, sex := 'Total']
  seroconverter_cohort_agg_r[, age_group := 'Total']
  
  seroconverter_cohort_agg <- do.call('rbind', list(seroconverter_cohort_agg, seroconverter_cohort_agg_s, seroconverter_cohort_agg_r))
  
  seroconverter_cohort_agg[, PY := paste0(comma_thousands(round(py, 2)), '[', comma_thousands(round(CL_py, 2)), '-', 
                                          comma_thousands(round(CU_py, 2)), ']')]
  seroconverter_cohort_agg[, INC := paste0(comma_thousands(round(hivinc, 2)), '[', comma_thousands(round(CL_hivinc, 2)), '-', 
                                           comma_thousands(round(CU_hivinc, 2)), ']')]
  seroconverter_cohort_agg[, PY := gsub('\\[', ' [', gsub(' ', '', PY))]
  seroconverter_cohort_agg[, INC := gsub('\\[', ' [', gsub(' ', '', INC))]
  seroconverter_cohort_agg <- seroconverter_cohort_agg[, .(round, sex, age_group, INC, PY)]
  
  # incidence rate by age group, sex and round
  inc_rate <- copy(modelpreds.agegroup.1218)
  inc_rate[, sex:= substr(Sex,1,1)]
  
  inc_rate_sex <- copy(modelpreds)
  inc_rate_sex[, sex:= substr(Sex,1,1)]
  inc_rate_sex[Sex == 'All', sex := 'Total']
  inc_rate_sex[, age_group := 'Total']
  setnames(inc_rate_sex, 'inc', 'incidence')
  
  inc_rate <- rbind(inc_rate, inc_rate_sex, fill=TRUE)
  
  inc_rate[, EST := paste0(comma_thousands(round(incidence*100, 2)), ' [', comma_thousands(round(lb*100, 2)), '-', 
                           comma_thousands(round(ub*100, 2)), ']')]
  inc_rate <- inc_rate[, .(round, sex, age_group, EST)]
  
  # merge
  tab <- merge(seroconverter_cohort_agg, inc_rate, by = c('round', 'sex', 'age_group'))

  .levels <- unique(seroconverter_cohort_agg$age_group)[unique(seroconverter_cohort_agg$age_group) != 'Total' ]
  tab[, age_group := factor(age_group, levels = c('Total', .levels))]
  tab[, sex := factor(sex, levels = c('Total', 'F', 'M'))]
  tab[, dummy := 1]
  tab <- tab[order(round, sex, age_group)]
  tab <- tab[, .(round, sex, age_group, dummy, INC, PY, EST)]
  
  stats[['table']] <- tab
  
  return(stats)
}

plot_data <- function(seroconverter_cohort, seroconverter_cohort_agg_s, outdir){
  
  date <- format(Sys.Date(), '%y%m%d')
  
  # find cases and py
  crude_cases_py <- subset(seroconverter_cohort) %>%
    group_by(sex, round, age, number_missing_visits_status) %>%
    summarise(pysum=py, seroconvsum=hivinc, inc=seroconvsum/pysum) %>%
    as.data.table()
  
  # expand to be sure to have all entries
  df_grid <- data.table(expand.grid(sex = crude_cases_py[, unique(sex)], 
                                    round = crude_cases_py[, unique(round)], 
                                    age = crude_cases_py[, unique(age)], 
                                    number_missing_visits_status = crude_cases_py[, unique(number_missing_visits_status)]))
  crude_cases_py <- merge(crude_cases_py, df_grid, by=c('sex', 'round','age','number_missing_visits_status'), all.y = T)
  crude_cases_py[is.na(pysum), pysum := 0]
  crude_cases_py[is.na(seroconvsum), seroconvsum := 0]
  
  # last changes in labels
  crude_cases_py <- merge(crude_cases_py, select(df_round, -visit), by.x = 'round', by.y= 'round_numeric') 
  crude_cases_py$SEX <- ifelse(crude_cases_py$sex == 'M', 'Male', 'Female')
  crude_cases_py$SEX_LABEL <- ifelse(crude_cases_py$sex == 'M', 'Men', 'Women')
  crude_cases_py.1218 <- as.data.table(crude_cases_py[which(crude_cases_py$round %in% c(5:9, rounds_numeric_group_2)), ])
  crude_cases_py.1218$ROUND <- gsub(':', '', crude_cases_py.1218$ROUND)
  
  colors <- c('#8BBCCC', '#256D85', '#06283D')
  
  # sum py accross number_missing_visits_status
  crude_cases_py_sum.1218 <- crude_cases_py.1218[, list(pysum = sum(pysum[is.na(seroconvsum)==F])), by = c('SEX', 'SEX_LABEL', 'ROUND', 'age')]
  
  # person-years male
  p1M <- ggplot(subset(crude_cases_py_sum.1218[SEX == 'Male'] )) +
    geom_bar(aes(x=age, y=pysum), stat = 'identity', fill= 'grey50') +
    facet_grid(ROUND~SEX_LABEL) +
    theme(legend.position="bottom") +
    xlab("Age")+ylab("Person-years at risk in the RCCS incidence cohort")+
    labs(fill = 'Missed preceding visit') + 
    scale_y_continuous(expand=expansion(mult = c(0, .1)),limits = c(0, 420)) +
    scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 49)) + 
    theme_bw() +
    theme(legend.position="bottom",
          # axis.text.x = element_text(angle = 30, hjust = 1),
          panel.spacing.x = unit(0.35, "cm"),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black"))
  
  # person-years female
  p1F <- ggplot(subset(crude_cases_py_sum.1218[SEX == 'Female'])) +
    geom_bar(aes(x=age, y=pysum), stat = 'identity', fill= 'grey50') +
    facet_grid(ROUND~SEX_LABEL) +
    theme(legend.position="bottom") +
    xlab("Age")+ylab("Person-years at risk in the RCCS incidence cohort")+
    labs(fill = 'Missed preceding visit') + 
    scale_y_continuous(expand=expansion(mult = c(0, .1)),limits = c(0, 420)) +
    scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 49)) + 
    theme_bw() +
    theme(legend.position="bottom",
          # axis.text.x = element_text(angle = 30, hjust = 1),
          panel.spacing.x = unit(0.35, "cm"),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black"))
  
  # incident HIV infections male
  p2M <- ggplot(subset(crude_cases_py.1218[SEX == 'Male'], is.na(seroconvsum)==F)) +
    geom_bar(aes(x=age, y=seroconvsum, fill = number_missing_visits_status ), stat = 'identity') +
    facet_grid(ROUND~SEX_LABEL) +
    theme(legend.position="bottom") +
    xlab("Age")+ylab("HIV incidence events in the RCCS incidence cohort")+
    labs(fill = 'Missed preceding visit') + 
    scale_fill_manual(values = colors) + 
    scale_y_continuous(expand=expansion(mult = c(0, .1)),limits = c(0, 8)) +
    scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 49)) + 
    scale_color_viridis_d(end = 0.85, option = 'A') +
    theme_bw() +
    theme(legend.position="bottom",
          # axis.text.x = element_text(angle = 30, hjust = 1),
          panel.spacing.x = unit(0.35, "cm"),
          strip.text = element_text(size = 10),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black"))
  
  # incident HIV infections female
  p2F <- ggplot(subset(crude_cases_py.1218[SEX == 'Female'], is.na(seroconvsum)==F)) +
    geom_bar(aes(x=age, y=seroconvsum, fill = number_missing_visits_status ), stat = 'identity') +
    facet_grid(ROUND~SEX_LABEL) +
    theme(legend.position="bottom") +
    xlab("Age")+ylab("HIV incidence events in the RCCS incidence cohort")+
    labs(fill = 'Missed preceding visit') + 
    scale_fill_manual(values = colors) + 
    scale_y_continuous(expand=expansion(mult = c(0, .1)),limits = c(0, 8)) +
    scale_x_continuous(expand = c(0,0), breaks= c(seq(15, 45, 5), 49)) + 
    scale_color_viridis_d(end = 0.85, option = 'A') +
    theme_bw() +
    theme(legend.position="bottom",
          # axis.text.x = element_text(angle = 30, hjust = 1),
          panel.spacing.x = unit(0.35, "cm"),
          strip.text.x = element_text(size = 10),
          strip.text.y = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black"))
  
  # arrange
  p <- ggarrange(p1F, p2F, p1M, p2M, nrow = 1, common.legend = T, legend='bottom', labels = c('a', 'b', 'c', 'd'))
  
  # save
  file = file.path(outdir, paste0('data-inland_', date, '.pdf'))
  ggsave(p, file= file, w = 9, h = 12)
  
  # person-years and incidence cases by Gender comparison to NEJM 
  # NEJM hiv incidence rate estimate (https://www.nejm.org/doi/suppl/10.1056/NEJMoa1702150/suppl_file/nejmoa1702150_appendix.pdf)
  nejm.f <- data.table(sex = 'F', 
                       visit = c('R011', 'R012', 'R013', 'R014', 'R015', 'R016', 'R017'), 
                       pysum = c(4425,4978,5610,5319,5587,7090,6689), 
                       seroconvsum = c(50,65,67,61,50,55,56))
  nejm.m <- data.table(sex = 'M', 
                       visit = c('R011', 'R012', 'R013', 'R014', 'R015', 'R016', 'R017'), 
                       pysum = c(3348,3791,4591,4497,4765,6069,5904), 
                       seroconvsum = c(36,40,58,44,36,32,27))
  nejm <- do.call('rbind', list(nejm.f, nejm.m))%>%
    merge(df_round, by = 'visit') %>%
    mutate(Sex = ifelse(sex == 'M', 'Male', 'Female'))
  
  # our cases and py (not stratified by age)
  crude_cases_py_noage <- copy(seroconverter_cohort_agg_sm) %>%
    subset(round >= 5)%>%
    mutate(Sex = ifelse(sex == 'M', 'Male', 'Female'), 
           number_missing_visits_status = factor(number_missing_visits_status, 
                                                 levels = rev(levels(seroconverter_cohort_agg_sm$number_missing_visits_status)))) 
  
  # plot py
  ggplot() +
    geom_bar(data=subset(crude_cases_py_noage), aes(y=py, x=round + 5, 
                                                    fill = number_missing_visits_status),size=1, stat = 'identity')+ 
    geom_point(data=subset(nejm),aes(y=pysum, x=round_numeric + 5, col = 'NEJM analysis'))+
    facet_grid(~Sex) +
    xlab("Round")+ ylab("Person-years")+
    labs(col = '', fill = '') + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_bw(base_size=16) +
    scale_color_manual(values = c('darkred', 'black' )) + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="black", fill="white"))+
    theme(legend.title=element_blank(), legend.position='bottom') + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  ggsave(file.path(outdir, paste0('data-person_years_comparison_NEJM', date, '.png')), w = 7, h = 5)
  
  # plot incidence cases
  ggplot() +
    geom_bar(data=subset(crude_cases_py_noage), aes(y=hivinc, x=round + 5, 
                                                    fill = number_missing_visits_status),size=1, stat = 'identity')+ 
    geom_point(data=subset(nejm),aes(y=seroconvsum, x=round_numeric + 5, col = 'NEJM analysis'))+
    facet_grid(~Sex) +
    xlab("Round")+ ylab("Incident cases")+
    labs(col = '', fill = '') + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_bw(base_size=16) +
    scale_color_manual(values = c('darkred', 'black' )) + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="black", fill="white"))+
    theme(legend.title=element_blank(), legend.position='bottom') + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  ggsave(file.path(outdir, paste0('data-incidentcases_comparison_NEJM', date, '.png')), w = 7, h = 5)
  
}

plot_model_fit <- function(modelpreds, modelpreds.age.1218, outdir){
  
  #date
  date <- format(Sys.Date(), '%y%m%d')
  
  # incidence trends by Gender
  ggplot() +
    geom_line(data=subset(modelpreds), aes(y=inc*100, x=round + 5),size=1)+ 
    geom_ribbon(data=subset(modelpreds),aes(ymin=lb*100, ymax=ub*100, x=round + 5), alpha = 0.3)+
    facet_grid(~Sex) +
    xlab("Round")+ ylab("Incidence rate per 100 py")+
    scale_y_continuous(expand=c(0,0)) + 
    theme_bw(base_size=16) +
    theme_bw() + 
    theme(strip.background = element_rect(colour="black", fill="white"))+
    theme(legend.title=element_blank(), legend.position='bottom') + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  ggsave(file.path(outdir, paste0('incidence_rate_sex_inland_', date, '.png')), w = 7, h = 5)
  
  ggplot() +
    geom_line(data=subset(modelpreds), aes(y=inc*100, x=round + 5,col = Sex),size=1)+ 
    geom_ribbon(data=subset(modelpreds),aes(ymin=lb*100, ymax=ub*100, x=round + 5, fill = Sex), alpha = 0.3)+
    xlab("Round")+ ylab("Incidence rate per 100 py")+
    scale_y_continuous(expand=c(0,0)) + 
    theme_bw(base_size=16) +
    theme_bw() + 
    theme(strip.background = element_rect(colour="black", fill="white"))+
    theme(legend.title=element_blank(), legend.position='bottom') + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  ggsave(file.path(outdir, paste0('incidence_rate_sex_inland2_', date, '.png')), w = 7, h = 5)
  
  # incidence trends by Gender comparison to NEJM 
  # NEJM hiv incidence rate estimate (https://www.nejm.org/doi/suppl/10.1056/NEJMoa1702150/suppl_file/nejmoa1702150_appendix.pdf)
  
  nejm.f <- data.table(Sex = 'Female', 
                       visit = c('R011', 'R012', 'R013', 'R014', 'R015', 'R016', 'R017'), 
                       inc = c(1.13, 1.31, 1.19, 1.15, 0.89, 0.78, 0.84), 
                       lb = c(0.84, 1.01, 0.93, 0.88, 0.67, 0.59, 0.64), 
                       ub = c(1.47, 1.65, 1.5, 1.46, 1.17, 1, 1.08))
  nejm.m <- data.table(Sex = 'Male', 
                       visit = c('R011', 'R012', 'R013', 'R014', 'R015', 'R016', 'R017'), 
                       inc = c(1.08, 1.06, 1.26, 0.98, 0.76, 0.53, 0.46), 
                       lb = c(0.76, 0.76, 0.97, 0.72, 0.53, 0.37, 0.31), 
                       ub = c(1.47, 1.42, 1.62, 1.3, 1.03, 0.73, 0.65))
  nejm.a <- data.table(Sex = 'All', 
                       visit = c('R011', 'R012', 'R013', 'R014', 'R015', 'R016', 'R017'), 
                       inc = c(1.11, 1.2, 1.23, 1.07, 0.83, 0.66, 0.66), 
                       lb = c(0.89, 0.98, 1.02, 0.88, 0.67, 0.53, 0.53), 
                       ub = c(1.36, 1.44, 1.45, 1.29, 1.02, 0.81, 0.81))
  nejm <- do.call('rbind', list(nejm.f, nejm.m, nejm.a))
  nejm <- merge(nejm, df_round, by = 'visit')
  
  ggplot() +
    geom_line(data=subset(modelpreds), aes(y=inc*100, x=round + 5, col = 'Our estimates'),size=1)+ 
    geom_ribbon(data=subset(modelpreds),aes(ymin=lb*100, ymax=ub*100, x=round + 5, fill = 'Our estimates'), alpha = 0.3)+
    geom_line(data=subset(nejm), aes(y=inc, x=round_numeric + 5, col = 'NEJM estimates'),size=1)+ 
    geom_ribbon(data=subset(nejm),aes(ymin=lb, ymax=ub, x=round_numeric + 5, fill = 'NEJM estimates'), alpha = 0.3)+
    facet_grid(~Sex) +
    xlab("Round")+ ylab("Incidence rate per 100 py")+
    labs(col = '', fill = '') + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_bw(base_size=16) +
    scale_color_manual(values = c('darkred', 'black')) + 
    scale_fill_manual(values = c('darkred', 'black')) + 
    theme_bw() + 
    theme(strip.background = element_rect(colour="black", fill="white"))+
    theme(legend.title=element_blank(), legend.position='bottom') + 
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  ggsave(file.path(outdir, paste0('incidence_rate_sex_inland_comparison_NEJM_', date, '.png')), w = 7, h = 5)
  
  # incidence trends by Gender and age uncertainty from estimate
  ggplot(subset(modelpreds.age.1218, model == 'model_1')) +
    geom_line(aes(x=age, y=incidence*100, col=Sex))+ 
    geom_ribbon(aes(x=age, ymin=lb*100, ymax=ub*100, fill =Sex), alpha = 0.5)+ 
    facet_grid(.~ROUND) +
    theme(legend.position="bottom") +
    ylab("ncidence rate per 100 py")+xlab("Age")+
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0), limits = c(0, 2.2)) +
    theme_bw() + 
    theme(legend.position="bottom", legend.key=element_rect(color="black"),
          strip.background = element_rect(colour="black", fill="white"),
          panel.border = element_rect(colour = "black"))
  ggsave(file.path(outdir, paste0('incidence_rate_sex_age_inland_', date, '.png')), w = 16, h = 5)

ggplot(subset(modelpreds.age.1218, model == 'model_1')) +
  geom_raster(aes(x=round+5, y=age, fill = incidence*100))+ 
  facet_grid(.~Sex) +
  scale_fill_gradientn(name="Incidence rate (per 100 py)",colours=c("white", "darkblue", "darkred")) +
  theme(legend.position="bottom") +
  labs(fill = "Crude incidence rate per 100 py")+ylab("Age")+
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  theme_bw() + 
  theme(legend.position="bottom", legend.key=element_rect(color="black"),
        strip.background = element_rect(colour="black", fill="white"),
        panel.border = element_rect(colour = "black"))
ggsave(file.path(outdir, paste0('incidence_rate_heat_sex_age_inland_', date, '.png')), w = 7, h = 8)

ggplot(subset(modelpreds.age.1218, model == 'model_1')) +
  geom_line(aes(x=age, y=incidence*100, color=as.factor(ROUND))) + 
  geom_ribbon(aes(x=age, ymin=lb*100, ymax=ub*100, fill=as.factor(ROUND)), alpha = 0.1)+
  facet_grid(~Sex) +
  theme(legend.position="bottom") +
  labs(col = "Year", fill = "Year",
       y = "Incidence rate (per 100 py)", x = "Age")+
  scale_x_continuous(expand=c(0,0)) + 
  theme_bw() +
  theme(legend.position="bottom",legend.title=element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        panel.border = element_rect(colour = "black")) 
ggsave(file.path(outdir, paste0('incidence_rate_heat_sex_age_round_inland_', date, '.png')), w = 8, h = 7)

}

find_prediction_1y <- function(modelpreds.age.1218.all, seroconverter_cohort_imputation){
  
  iters <- modelpreds.age.1218.all[, sort(unique(iterations))]
  seroconverter_cohort_imputation[, Sex := ifelse(sex == 'F', 'Female', 'Male')]

  model_pred.list <- vector(mode = 'list', length = length(iters))
  set.seed(12)
  for(i in iters){
    
    cat('\n iterations', i)
    
    tmp <- modelpreds.age.1218.all[iterations == i]
    
    # merge to data
    tmp <- merge(tmp, seroconverter_cohort_imputation, by = c('round', 'Sex', 'age', 'iterations'))
    
    # predict number of incidence infection
    tmp[, pred_hivinc := rpois(1, exp(fit)*pysum), by = c('round', 'Sex', 'age', 'iterations', 'iterations_within', 'model')]
    
    # crude incidence rate 
    tmp[, inc_crude := seroconvsum / pysum]
    
    # predict number of incidence infection
    tmp <- tmp[, list(pred_inc = median(pred_hivinc / pysum), 
                      pred_inc_lb = quantile(pred_hivinc / pysum, 0.025), 
                      pred_inc_ub = quantile(pred_hivinc / pysum, 0.975), 
                      pred_hivinc = median(pred_hivinc), 
                      pred_hivinc_lb = quantile(pred_hivinc , 0.025), 
                      pred_hivinc_ub = quantile(pred_hivinc , 0.975),
                      est_hivinc = median(inc*pysum), 
                      est_hivinc_lb= quantile(inc*pysum , 0.025), 
                      est_hivinc_ub = quantile(inc*pysum , 0.975), 
                      est_inc = median(inc), 
                      est_inc_lb= quantile(inc, 0.025), 
                      est_inc_ub = quantile(inc, 0.975)), by = c('round', 'Sex', 'age', 'iterations', 'inc_crude', 'seroconvsum', 'model', 'pysum')]
    
    model_pred.list[[i]] <- tmp
  }
  model_pred <- rbindlist(model_pred.list)
  
  # find percentage data within pred interval
  model_pred[, within.CI := seroconvsum >= pred_hivinc_lb &seroconvsum <= pred_hivinc_ub]
  
  
  return(model_pred)
}

find_statistics_prediction <- function(model_pred){
  
  stats_prediction <- list()
  comma_thousands <- function(x) format(x, big.mark=",")
  
  # prediction inside CI by sex
  model_pred_summary <- model_pred[, list(within.CI.prop =  mean(within.CI)), by = c('Sex', 'model', 'iterations')]
  model_pred_summary <- model_pred_summary[, list(q= c(quantile(within.CI.prop, prob=ps, na.rm = T), mean(within.CI.prop)),
                                                  q_label=c(p_labs, 'mean')), by = c('Sex', 'model')]	
  model_pred_summary <- dcast(model_pred_summary, Sex + model ~ q_label, value.var = "q")
  model_pred_summary[, INCI.M := paste0(comma_thousands(round(mean*100, 2)), '\\%')]
  model_pred_summary[, INCI.CI := paste0('[', comma_thousands(round(CL*100, 2)), '-', 
                                         comma_thousands(round(CU*100, 2)), ']')]
  model_pred_summary[, INCI.CI := gsub(' ', '', INCI.CI)]
  model_pred_summary <- model_pred_summary[order(Sex, model), .(Sex, model, INCI.M, INCI.CI)]
  stats_prediction$within.CI <- model_pred_summary
  
  # prediction inside CI total
  model_pred_summary <- model_pred[, list(within.CI.prop =  mean(within.CI)), by = c('model', 'iterations')]
  model_pred_summary <- model_pred_summary[, list(q= c(quantile(within.CI.prop, prob=ps, na.rm = T), mean(within.CI.prop)),
                                                  q_label=c(p_labs, 'mean')), by = c('model')]	
  model_pred_summary <- dcast(model_pred_summary, model ~ q_label, value.var = "q")
  model_pred_summary[, INCI.M := paste0(comma_thousands(round(mean*100, 2)), '\\%')]
  model_pred_summary[, INCI.CI := paste0('[', comma_thousands(round(CL*100, 2)), '-', 
                                         comma_thousands(round(CU*100, 2)), ']')]
  model_pred_summary[, INCI.CI := gsub(' ', '', INCI.CI)]
  model_pred_summary <- model_pred_summary[order(model), .( model, INCI.M, INCI.CI)]
  stats_prediction$within.CI.all <- model_pred_summary
  
  # aic
  modelaics.age.summary <- copy(modelaics.age)
  modelaics.age.summary[, AIC.M := paste0(comma_thousands(round(M)))]
  modelaics.age.summary[, AIC.CI := paste0('[', comma_thousands(round(CL)), '-', comma_thousands(round(CU)), ']')]
  modelaics.age.summary[, AIC.M := gsub(' ', '', AIC.M)]
  modelaics.age.summary[, AIC.CI := gsub(' ', '', AIC.CI)]
  modelaics.age.summary <- modelaics.age.summary[, .(Sex, model, AIC.M, AIC.CI)]
  stats_prediction$AIC <- modelaics.age.summary
  
  # mae total
  model_pred_summary <- model_pred[, list(MAE =  mean(abs(est_inc - inc_crude))), by = c('model', 'iterations')]
  model_pred_summary <- model_pred_summary[, list(q= c(quantile(MAE, prob=ps, na.rm = T), mean(MAE)),
                                                  q_label=c(p_labs, 'mean')), by = c('model')]	
  model_pred_summary <- dcast(model_pred_summary, model ~ q_label, value.var = "q")
  model_pred_summary[, MAE.M := round(mean, 4)]
  model_pred_summary[, MAE.CI := paste0('[', (round(CL, 4)), '-', (round(CU, 4)), ']')]
  model_pred_summary[, MAE.CI := gsub(' ', '', MAE.CI)]
  model_pred_summary <- model_pred_summary[order(model), .( model, MAE.M, MAE.CI)]
  stats_prediction$MAE.all <- model_pred_summary
  
  # mae by sex
  model_pred_summary <- model_pred[, list(MAE =  mean(abs(est_inc - inc_crude))), by = c('model', 'iterations', 'Sex')]
  model_pred_summary <- model_pred_summary[, list(q= c(quantile(MAE, prob=ps, na.rm = T), mean(MAE)),
                                                  q_label=c(p_labs, 'mean')), by = c('model', 'Sex')]	
  model_pred_summary <- dcast(model_pred_summary, model + Sex ~ q_label, value.var = "q")
  model_pred_summary[, MAE.M := round(mean, 4)]
  model_pred_summary[, MAE.CI := paste0('[', (round(CL, 4)), '-', (round(CU, 4)), ']')]
  model_pred_summary[, MAE.CI := gsub(' ', '', MAE.CI)]
  model_pred_summary <- model_pred_summary[order(model), .( model,Sex, MAE.M, MAE.CI)]
  stats_prediction$MAE <- model_pred_summary
  
  # last change 
  stats_prediction$MAE$MAE.M <- format(stats_prediction$MAE$MAE.M, scientific = F)
  
  return(stats_prediction)
}

plot_predicted_incidence_events <- function(model_pred, outdir){
  
  date <- format(Sys.Date(), '%y%m%d')
  
  palette_round <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(10)
  palette_round_inland <<- palette_round[c(1:6, 8:10)]
  
  # add some jitter
  set.seed(12)
  model_pred[, jitter := runif(length(seroconvsum), 0, 1), by= c('round', 'Sex', 'age', 'iterations', 'seroconvsum')]
  model_pred[, seroconvsum_jitter := seroconvsum + jitter]
  model_pred[, pred_hivinc_lb_jitter := pred_hivinc_lb + jitter]
  model_pred[, pred_hivinc_ub_jitter := pred_hivinc_ub + jitter]
  model_pred[, pred_hivinc_jitter := pred_hivinc + jitter]
  model_pred[, est_hivinc_lb_jitter := est_hivinc_lb + jitter]
  model_pred[, est_hivinc_ub_jitter := est_hivinc_ub + jitter]
  model_pred[, est_hivinc_jitter := est_hivinc + jitter]
  
  # add labels
  tmp <- merge(model_pred, df_round, by.x = 'round', by.y = 'round_numeric')
  tmp[, ROUND := gsub(':', '', ROUND)]
  tmp[, SEX_LABEL := 'Men']
  tmp[Sex == 'Female', SEX_LABEL := 'Women']
  tmp[, age_group := 'Age: 35-49']
  tmp[age < 35, age_group := 'Age: 25-34']
  tmp[age < 25, age_group := 'Age: 15-24']
  tmp <- tmp[model == 'model_1']
  
  ggplot(tmp) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    geom_errorbar(aes(x=seroconvsum_jitter, ymin=pred_hivinc_lb_jitter, ymax=pred_hivinc_ub_jitter),  color = 'grey70', 
                  width = 0, size = 0.2)+
    geom_point(aes(y=pred_hivinc_jitter, x=seroconvsum_jitter, fill=as.factor(ROUND), 
                   color=as.factor(ROUND)),
               size = 1) + 
    facet_wrap(age_group ~ SEX_LABEL, scales = 'free', ncol = 2) +
    theme(legend.position="bottom") +
    labs(col = "", fill = '',
         x = "Observed incident HIV infection", y = "Predicted incident HIV infection")+
    # scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values = palette_round_inland) + 
    scale_fill_manual(values = palette_round_inland) + 
    theme_bw() +
    theme(legend.position="bottom",
          strip.background = element_rect(colour="white", fill="white"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"), 
          legend.title = element_blank(), 
          panel.spacing.x = unit(0.75, "lines"),
          panel.spacing.y = unit(0.75, "lines")) + 
    guides(color = guide_legend(byrow = T, nrow = 2), shape = 'none')
  ggsave(file = file.path(outdir, paste0('incidence_cases_age_PPC_inland_', date, '.png')), w = 7, h = 8.5)
  
  ggplot(tmp) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    geom_errorbar(aes(x=seroconvsum_jitter, ymin=est_hivinc_lb_jitter, ymax=est_hivinc_ub_jitter),  color = 'grey70', 
                  width = 0, size = 0.2)+
    geom_point(aes(y=est_hivinc_jitter, x=seroconvsum_jitter, fill=as.factor(ROUND), 
                   color=as.factor(ROUND)),
               size = 1) + 
    facet_wrap(age_group ~ SEX_LABEL, scales = 'free', ncol = 2) +
    theme(legend.position="bottom") +
    labs(col = "", fill = '',
         x = "Observed incident HIV infection", y = "Expected incident HIV infection")+
    # scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values = palette_round_inland) + 
    scale_fill_manual(values = palette_round_inland) + 
    theme_bw() +
    theme(legend.position="bottom",
          strip.background = element_rect(colour="white", fill="white"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"), 
          legend.title = element_blank(), 
          panel.spacing.x = unit(0.75, "lines"),
          panel.spacing.y = unit(0.75, "lines")) + 
    guides(color = guide_legend(byrow = T, nrow = 2), shape = 'none')
  file = file.path(outdir, paste0('incidence_cases_age_PPC_est_inland_', date, '.png'))
  
  ggsave(file, w = 7, h = 8.5)

  
}

plot_estimate_incidence_rates <- function(model_pred, outdir){
  
  date <- format(Sys.Date(), '%y%m%d')
  
  palette_round <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(10)
  palette_round_inland <<- palette_round[c(1:6, 8:10)]

  # add some jitter
  set.seed(12)
  model_pred[, jitter := runif(length(inc_crude), 0, 0.004), by= c('round', 'Sex', 'age', 'iterations', 'inc_crude')]
  model_pred[, inc_crude_jitter := inc_crude + jitter]
  model_pred[, pred_inc_lb_jitter := pred_inc_lb + jitter]
  model_pred[, pred_inc_ub_jitter := pred_inc_ub + jitter]
  model_pred[, pred_inc_jitter := pred_inc + jitter]
  model_pred[, est_inc_lb_jitter := est_inc_lb + jitter]
  model_pred[, est_inc_ub_jitter := est_inc_ub + jitter]
  model_pred[, est_inc_jitter := est_inc + jitter]
  
  # add labels
  tmp <- merge(model_pred, df_round, by.x = 'round', by.y = 'round_numeric')
  tmp[, ROUND := gsub(':', '', ROUND)]
  tmp[, SEX_LABEL := 'Men']
  tmp[Sex == 'Female', SEX_LABEL := 'Women']
  tmp[, age_group := 'Age: 35-49']
  tmp[age < 35, age_group := 'Age: 25-34']
  tmp[age < 25, age_group := 'Age: 15-24']
  tmp <- tmp[model == 'model_1']
  
  ggplot(tmp) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    geom_errorbar(aes(x=inc_crude_jitter, ymin=pred_inc_lb_jitter, ymax=pred_inc_ub_jitter),  color = 'grey70', 
                  width = 0, size = 0.2)+
    geom_point(aes(y=pred_inc_jitter, x=inc_crude_jitter, fill=as.factor(ROUND), 
                   color=as.factor(ROUND)),
               size = 1) + 
    facet_wrap(age_group ~ SEX_LABEL, scales = 'free', ncol = 2) +
    theme(legend.position="bottom") +
    labs(col = "", fill = '',
         x = "Crude HIV incident rates", y = "Predicted HIV incident rates")+
    # scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values = palette_round_inland) + 
    scale_fill_manual(values = palette_round_inland) + 
    theme_bw() +
    theme(legend.position="bottom",
          strip.background = element_rect(colour="white", fill="white"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"), 
          legend.title = element_blank(), 
          panel.spacing.x = unit(0.75, "lines"),
          panel.spacing.y = unit(0.75, "lines")) + 
    guides(color = guide_legend(byrow = T, nrow = 2), shape = 'none')
  ggsave(file.path(outdir, paste0('incidence_rates_age_PPC_inland_', date, '.png')), w = 7, h = 8.5)
  
  ggplot(tmp) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    geom_errorbar(aes(x=inc_crude_jitter, ymin=est_inc_lb_jitter, ymax=est_inc_ub_jitter),  color = 'grey70', 
                  width = 0, size = 0.2)+
    geom_point(aes(y=est_inc_jitter, x=inc_crude_jitter, fill=as.factor(ROUND), 
                   color=as.factor(ROUND)),
               size = 1) + 
    facet_wrap(age_group ~ SEX_LABEL, scales = 'free', ncol = 2) +
    theme(legend.position="bottom") +
    labs(col = "", fill = '',
         x = "Crude HIV incident rates", y = "Estimated HIV incident rates")+
    # scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values = palette_round_inland) + 
    scale_fill_manual(values = palette_round_inland) + 
    theme_bw() +
    theme(legend.position="bottom",
          strip.background = element_rect(colour="white", fill="white"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"), 
          legend.title = element_blank(), 
          panel.spacing.x = unit(0.75, "lines"),
          panel.spacing.y = unit(0.75, "lines")) + 
    guides(color = guide_legend(byrow = T, nrow = 2), shape = 'none')
  ggsave(file.path(outdir, paste0('incidence_rates_age_PPC_est_inland_', date, '.png')), w = 7, h = 8.5)
  
}

plot_estimate_incidence_rates_summary <- function(modelpreds.age.1218, seroconverter_cohort_crude, outdir){
  
  date <- format(Sys.Date(), '%y%m%d')
  
  palette_round <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(10)
  palette_round_inland <<- palette_round[c(1:6, 8:10)]
  
  # combine
  tmp <- modelpreds.age.1218[model == 'model_1']
  tmp1<- copy(seroconverter_cohort_crude)
  tmp1[, Sex := ifelse(sex == 'F', 'Female', 'Male')]
  
  # labels
  tmp <- merge(tmp, tmp1, by = c('Sex', 'round', 'age'))
  tmp[, ROUND := gsub(':', '', ROUND)]
  tmp[, SEX_LABEL := 'Men']
  tmp[Sex == 'Female', SEX_LABEL := 'Women']
  tmp[, age_group := 'Age: 35-49']
  tmp[age < 35, age_group := 'Age: 25-34']
  tmp[age < 25, age_group := 'Age: 15-24']
  
  ggplot(tmp) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    geom_errorbar(aes(ymin=lb*100, ymax = ub*100, x=inc_crude*100), alpha = 0.3, col = 'grey50') + 
    geom_errorbarh(aes(xmin=CL*100, xmax = CU*100, y=incidence*100), alpha = 0.3, col = 'grey50') + 
    geom_point(aes(y=incidence*100, x=inc_crude*100, fill=as.factor(ROUND), color=as.factor(ROUND)),
               size = 1) + 
    facet_wrap(. ~ SEX_LABEL, ncol = 2) +
    theme(legend.position="bottom") +
    labs(col = "", fill = '',
         x = "Crude incident rates\nper 100 person-years", y = "Estimated incident rates\nper 100 person-years")+
    # scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values = palette_round_inland) + 
    scale_fill_manual(values = palette_round_inland) + 
    theme_bw() +
    theme(legend.position="bottom",
          strip.background = element_rect(colour="white", fill="white"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"), 
          legend.title = element_blank(), 
          panel.spacing.x = unit(0.75, "lines"),
          panel.spacing.y = unit(0.75, "lines")) + 
    guides(color = guide_legend(byrow = T, nrow = 2), shape = 'none')
  ggsave(file.path(outdir, paste0('incidence_rates_age_PPC_sum_inland_', date, '.pdf')), w = 5.2, h = 4.1)
}

find_incidence_estimates_loess <- function(model_pred){
  
  # use loess on the incidence rate
  modelpreds.loess.age.1218.all <- model_pred[model == 'model_1', {
    loessMod50 <- loess(sqrt(inc_crude) ~ age, span=0.7, weights = pysum)
    smoothed50 <- predict(loessMod50, new_data = age)^2
    
    list(age = age, INC_CRUDE_SMOOTH = smoothed50, 
         inc_crude = inc_crude, pysum = pysum)
  }, by = c('round', 'Sex', 'iterations')]
  nrow(modelpreds.loess.age.1218.all) #31500
  
  return(modelpreds.loess.age.1218.all)
}

plot_comparison_loess_gam <- function(modelpreds.loess.age.1218, modelpreds.age.1218, 
                                      seroconverter_cohort_imputation, outdir){
  # combine estimates
  loess <- copy(modelpreds.loess.age.1218)
  loess[, type := 'Estimates obtained with Local Polynomial Regression (Loess)']
  loess <- merge(loess, unique(modelpreds.age.1218[, .(round_label, round, ROUND)]), by = 'round_label')
  gam <- select(modelpreds.age.1218[model == 'model_1'], -c( 'model'))
  gam[, type := 'Estimates obtained with Generalized additive model (GAM)']
  tmp <- rbind(loess, gam)
  
  # crude estimates
  crude <- merge(seroconverter_cohort_imputation, df_round, by.x = 'round', by.y = 'round_numeric')
  crude <- crude[round %in% tmp[, unique(round)]]
  
  # label
  tmp[, ROUND := gsub(':', '', ROUND)]
  tmp[, SEX_LABEL := 'Men']
  tmp[Sex == 'Female', SEX_LABEL := 'Women']
  tmp[, age_group := 'Age: 35-49']
  tmp[age < 35, age_group := 'Age: 25-34']
  tmp[age < 25, age_group := 'Age: 15-24']
  
  crude[, ROUND := gsub(':', '', ROUND)]
  crude[, SEX_LABEL := 'Men']
  crude[Sex == 'Female', SEX_LABEL := 'Women']
  crude[, age_group := 'Age: 35-49']
  crude[age < 35, age_group := 'Age: 25-34']
  crude[age < 25, age_group := 'Age: 15-24']
  
  # plot
  ggplot(tmp) +
    geom_point(data = crude, aes(y=inc_crude, x=age, size = 'Crude estimates'), alpha = 0.1) + 
    geom_line(aes(y=incidence, x=age, color=type)) + 
    geom_ribbon(aes(ymin=lb, x=age,  ymax = ub, fill=type), alpha = 0.5) + 
    facet_grid(ROUND ~ SEX_LABEL) +
    theme(legend.position="bottom") +
    labs(col = "", fill = '',
         x = "Age", y = "HIV incident rates")+
    scale_size_manual(values = 0.4) + 
    # scale_x_continuous(expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.direction = 'vertical',
          strip.background = element_rect(colour="white", fill="white"),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black"), 
          legend.title = element_blank(), 
          panel.spacing.x = unit(0.75, "lines"),
          panel.spacing.y = unit(0.75, "lines")) +
    guides(size = guide_legend(override.aes = list(alpha = 1, size = 1), order = 1))
  ggsave(file.path(outdir, paste0('incidence_rates_age_smooth_comp_weighted2_inland_test-1y.png')), w = 7, h = 12)
  
}

find_statistics_prediction_loess <- function(modelpreds.loess.age.1218.all){
  
  stats_prediction_loess <- list()
  
  # mae total
  model_pred_summary <- modelpreds.loess.age.1218.all[, list(MAE =  mean(abs(INC_CRUDE_SMOOTH - inc_crude))), by = c('iterations')]
  model_pred_summary <- model_pred_summary[, list(q= c(quantile(MAE, prob=ps, na.rm = T), mean(MAE)),
                                                  q_label=c(p_labs, 'mean'))]	
  model_pred_summary <- dcast(model_pred_summary, . ~ q_label, value.var = "q")
  model_pred_summary[, MAE.M := round(mean, 4)]
  model_pred_summary[, MAE.CI := paste0('[', (round(CL, 4)), '-', (round(CU, 4)), ']')]
  model_pred_summary[, MAE.CI := gsub(' ', '', MAE.CI)]
  model_pred_summary <- model_pred_summary[, .( MAE.M, MAE.CI)]
  stats_prediction_loess$MAE.all <- model_pred_summary
  
  # mae by sex
  model_pred_summary <- modelpreds.loess.age.1218.all[, list(MAE =  mean(abs(INC_CRUDE_SMOOTH - inc_crude))), by = c( 'iterations', 'Sex')]
  model_pred_summary <- model_pred_summary[, list(q= c(quantile(MAE, prob=ps, na.rm = T), mean(MAE)),
                                                  q_label=c(p_labs, 'mean')), by = c( 'Sex')]	
  model_pred_summary <- dcast(model_pred_summary,  Sex ~ q_label, value.var = "q")
  model_pred_summary[, MAE.M := round(mean, 4)]
  model_pred_summary[, MAE.CI := paste0('[', (round(CL, 4)), '-', (round(CU, 4)), ']')]
  model_pred_summary[, MAE.CI := gsub(' ', '', MAE.CI)]
  model_pred_summary <- model_pred_summary[, .( Sex, MAE.M, MAE.CI)]
  stats_prediction_loess$MAE <- model_pred_summary
  
  return(stats_prediction_loess)
}

