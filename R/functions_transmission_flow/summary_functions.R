make.time.since.infection <- function(DT)
{
        # set iDats to Dates
        is.idate <- function(x) inherits(x, 'IDate')
        cols <- DT[, lapply(.SD, is.idate) ,] 
        cols <- names(which(unlist(cols))) 
        DT[, (cols) := lapply(.SD, as.Date) , .SDcols=cols]

        # Extract columns of interest
        cols <- grep("AID|visit_dt|pred_doi",names(DT), value=T)
        tmp <- DT[, ..cols]
        setcolorder(tmp, 'AID')

        cnd <- tmp[, .N>1, by='AID'][, any(V1)]
        if(cnd)
        {
                cat('Multiple samples per participant used:\n',
                    '\textract TSI prediction from first pos sample\n')
                setorder(tmp, AID, visit_dt)
                tmp <- tmp[, lapply(.SD, function(x) x[1] ) , by='AID']
        }

        # Get study_id's 
        tmp <- merge(aik[, .(study_id = PT_ID, AID)], tmp, all.y=TRUE)

        if('pred_doi_adjusted_mid' %in% names(tmp))
        {
                cat('Using age-deltaT adjusted TSI\n')

                names(tmp)
                tmp[, bias_adj := pred_doi_adjusted_mid - pred_doi_mid ]
                cols <- c('pred_doi_max', 'pred_doi_min')
                tmp[,  (cols) := lapply(.SD,
                              function(x) {x + bias_adj}
                              ) , .SDcols=cols]
                

                setnames(tmp, 'pred_doi_adjusted_mid', 'date_infection')
                tmp[ ,pred_doi_mid := NULL]

        }else{
                cat('Using Tanya-like TSIs\n')
                setnames(tmp, 'pred_doi_mid', 'date_infection')
        }

        tmp
}

make.df.period <- function(start_first_period_inland, stop_first_period_inland, 
                           start_second_period_inland, stop_second_period_inland, df_round)
{
  ## make map for the two time periods 
  
  tmp_inland <- data.table(PERIOD = c(paste0(format(start_first_period_inland, '%d %b %Y'), '-', format(stop_first_period_inland, '%d %b %Y')), 
                                      paste0(format(start_second_period_inland, '%d %b %Y'), '-', format(stop_second_period_inland, '%d %b %Y'))), 
                           BEFORE_CUTOFF = c(T, F), 
                           INDEX_TIME = 1:2, 
                           COMM = 'inland', 
                           MIN_PERIOD_DATE = c(start_first_period_inland, start_second_period_inland),
                           MAX_PERIOD_DATE = c(stop_first_period_inland, stop_second_period_inland))

  # find period span
  df_round[, ROUND_SPANYRS := .year.diff(MAX_SAMPLE_DATE, MIN_SAMPLE_DATE)]
  tmp <- df_round[, list(PERIOD_SPAN = sum(ROUND_SPANYRS)), by = c('COMM', 'INDEX_TIME')]
  tmp_inland <- merge(tmp, tmp_inland, by = c('COMM', 'INDEX_TIME'))
  
  # make period a factor
  tmp_inland[, PERIOD := factor(PERIOD, levels = PERIOD)]
  
  return(tmp_inland)
}

make.df.round <- function(df_round_inland)
{
  
  ## map for rounds 
  
  # define the time period corresponding to the rounds
  df_round_inland[, INDEX_TIME := 0]
  df_round_inland[round%in%paste0('R0', 10:15), INDEX_TIME := 1] 
  df_round_inland[round %in% paste0('R0',16:18), INDEX_TIME := 2]
  
  # keep only rounds that correspond to time periods
  df_round <- df_round_inland[INDEX_TIME != '0']
  df_round <- df_round[order(COMM, round)]
  
  # add index of rounds
  df_round[, INDEX_ROUND := 1:length(round), by = 'COMM']
  
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

add_susceptible_infected <- function(eligible_count_smooth, proportion_prevalence, participation, 
                                     nonparticipants.male.relative.infection, 
                                     nonparticipants.female.relative.infection)
{
  
  # prevalence 
  df <- copy(proportion_prevalence)
  df[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
  
  # merge prevalence to count of eligible population
  df <- merge(eligible_count_smooth[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)], df, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))
  
  # add participation
  df <- merge(df, participation, by = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))
  
  # find infected count
  df[SEX == 'M', INFECTED := ELIGIBLE * PARTICIPATION * PREVALENCE_M + 
       ELIGIBLE * (1-PARTICIPATION) * PREVALENCE_M * nonparticipants.male.relative.infection]
  df[SEX == 'F', INFECTED := ELIGIBLE * PARTICIPATION * PREVALENCE_M + 
       ELIGIBLE * (1-PARTICIPATION) * PREVALENCE_M * nonparticipants.female.relative.infection]
  
  # find susceptible count
  df[, SUSCEPTIBLE := ELIGIBLE - INFECTED]
  
  # do not include 15S in inland
  df <- df[!(COMM == 'inland' & ROUND == '15S')]
  
  return(df)
}

add_infected_unsuppressed <- function(eligible_count_round, treatment_cascade, participation, 
                                      nonparticipants.treated.like.participants, 
                                      nonparticipants.not.treated)
{
  
  # ensure that the data are data.table objects
  di <- as.data.table(eligible_count_round)
  pu <- as.data.table(treatment_cascade)
  pa <- as.data.table(participation)
  
  # select variable of interest
  di <- di[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE, INFECTED, SUSCEPTIBLE)]
  di[, ROUND := paste0('R0', ROUND)]
  
  # merge to the proportion of diagnosed and unsuppressed given diagnosed 
  pu[, .(ROUND, COMM, AGEYRS, SEX, 
         PROP_UNSUPPRESSED_PARTICIPANTS_M, PROP_UNSUPPRESSED_PARTICIPANTS_CL, PROP_UNSUPPRESSED_PARTICIPANTS_CU, 
         PROP_UNSUPPRESSED_NONPARTICIPANTS_M, PROP_UNSUPPRESSED_NONPARTICIPANTS_CL, PROP_UNSUPPRESSED_NONPARTICIPANTS_CU)]
  df <- merge(di, pu, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # merge to participation (i.e., proporiton of census eligible population that participated)
  pa <- pa[, .(ROUND, COMM, AGEYRS, SEX, PARTICIPATION)]
  pa[, ROUND := paste0('R0', ROUND)]
  df <- merge(df, pa, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'))
  
  # get infected un suppressed
  
  if(nonparticipants.not.treated){
    # assuming that only participants are treated and non-participants are all unsuppressed
    df[, INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS_M + INFECTED * (1-PARTICIPATION) * 1]
    df[, INFECTED_NON_SUPPRESSED_CL := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS_CL + INFECTED * (1-PARTICIPATION) * 1]
    df[, INFECTED_NON_SUPPRESSED_CU := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS_CU + INFECTED * (1-PARTICIPATION) * 1]
  }else if(nonparticipants.treated.like.participants){
    # assuming that participant and non-participant are diagnosed and treated at the same proportion
    df[, INFECTED_NON_SUPPRESSED := INFECTED * PROP_UNSUPPRESSED_PARTICIPANTS_M]
    df[, INFECTED_NON_SUPPRESSED_CL := INFECTED * PROP_UNSUPPRESSED_PARTICIPANTS_CL]
    df[, INFECTED_NON_SUPPRESSED_CU := INFECTED * PROP_UNSUPPRESSED_PARTICIPANTS_CU]
  }else{
    # assuming that participant and non-participant are treated not at the same proportion
    # and using newly registered participants to inform suppression of non-participnats
    df[, INFECTED_NON_SUPPRESSED := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS_M + INFECTED * (1-PARTICIPATION) * PROP_UNSUPPRESSED_NONPARTICIPANTS_M]
    df[, INFECTED_NON_SUPPRESSED_CL := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS_CL + INFECTED * (1-PARTICIPATION) * PROP_UNSUPPRESSED_NONPARTICIPANTS_CL]
    df[, INFECTED_NON_SUPPRESSED_CU := INFECTED * PARTICIPATION * PROP_UNSUPPRESSED_PARTICIPANTS_CU + INFECTED * (1-PARTICIPATION) * PROP_UNSUPPRESSED_NONPARTICIPANTS_CU]
    
  }

  # merge to df round
  df <- merge(df, unique(df_round[, .(COMM)]), by = c('COMM'))
  
  # rm unecessary variable
  df <- select(df, -c('PROP_UNSUPPRESSED_PARTICIPANTS_M', 'PROP_UNSUPPRESSED_PARTICIPANTS_CL', 'PROP_UNSUPPRESSED_PARTICIPANTS_CU', 
                      'PROP_UNSUPPRESSED_NONPARTICIPANTS_M', 'PROP_UNSUPPRESSED_NONPARTICIPANTS_CL', 'PROP_UNSUPPRESSED_NONPARTICIPANTS_CU',
                      'PARTICIPATION'))
  
  return(df)
}

get_incidence_cases_round <- function(incidence.inland, eligible_count_round)
{
  
  # prepare incidence rates inland
  inc.inland <- copy(incidence.inland)
  colnames(inc.inland) <- toupper(colnames(inc.inland))
  setnames(inc.inland, 'AGE', 'AGEYRS')
  setnames(inc.inland, 'ROUND_LABEL', 'ROUND')
  inc.inland[, COMM := 'inland']
  inc.inland[, SEX := substring(SEX, 1, 1)]
  inc.inland <- inc.inland[ROUND >= 10]
  inc.inland[, ROUND := paste0('R0', as.character(ROUND))]
  incidence <- inc.inland[, .(SEX, ROUND, AGEYRS, INCIDENCE, LB, UB, COMM)]
  
  if(0)
  {
    ggplot(incidence, aes(x = AGEYRS)) +
      geom_line(aes(y = INCIDENCE*100, col = ROUND)) +
      geom_ribbon(aes(ymin = LB*100, ymax = UB*100, fill = ROUND),  alpha = 0.1) +
      labs(y = 'Incidence rate per 100 PY ', x = 'Age') +
      facet_grid(COMM~SEX, label = 'label_both', scales = 'free_y') +
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
  if(0)
  {
    
    ggplot(dir, aes(x = AGEYRS)) +
      geom_line(aes(y = INCIDENT_CASES, col = ROUND)) +
      geom_ribbon(aes(ymin = INCIDENT_CASES_LB, ymax = INCIDENT_CASES_UB, fill = ROUND), alpha = 0.1) +
      labs(y = 'Expected number of incident cases', x = 'Age') +
      facet_grid(COMM~SEX, label = 'label_both', scales = 'free') +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    tmp <- dir[, list(INCIDENT_CASES = round(sum(INCIDENT_CASES), digits = 1),
                      INCIDENT_CASES_UB = round(sum(INCIDENT_CASES_UB), digits = 1),
                      INCIDENT_CASES_UB = round(sum(INCIDENT_CASES_UB), digits = 1)), by = c('ROUND', 'COMM')]
    knitr::kable(tmp)
    
  }

  # save
  file.name <- file.incidence.cases.round
  if(! file.exists(file.name) | config$overwrite.existing.files )
  {
    cat("Saving file:", file.name, '\n')
    saveRDS(dir, file = file.name)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
  
  return(dir)
  
}

summarise_incidence_cases_period <- function(incidence_cases_round, df_period)
{
  
  # sum across time periods
  incidence_cases <- incidence_cases_round[, list(INCIDENT_CASES = sum(INCIDENT_CASES), 
                                                  INCIDENT_CASES_UB = sum(INCIDENT_CASES_UB), 
                                                  INCIDENT_CASES_LB = sum(INCIDENT_CASES_LB)), by = c('COMM', 'AGEYRS', 'SEX', 'INDEX_TIME')]
  # merge to map period
  incidence_cases <- merge(incidence_cases, df_period, by = c('INDEX_TIME', 'COMM'))
  
  return(incidence_cases)
}

find.time.of.infection <- function(meta, time.since.infection, use.TSI.estimate)
{

  if(use.TSI.estimate)
  {
    cat("# MAKE TSI ADJUSTMENT\n")
    tmp <- unique(meta[,.(study_id, date_birth)])
    tmp <- merge(tmp, time.since.infection, by='study_id')

    # Compute estimated infection ages
    tmp[, age_collection := .year.diff(visit_dt, date_birth)]
    tmp[, age_infection  := .year.diff(date_infection, date_birth)]

    
    meta
    meta <- merge(meta, tmp, by=c('study_id', 'date_birth'), all.x=T) 

    # Print statements
    tmp <- unique(meta[,.(study_id, age_infection)])

    cat(tmp[, sum(!is.na(age_infection))], 
        "TSI estimates out of", tmp[,.N], "individuals in meta data\n")
    cat(tmp[is.na(age_infection), .N],
        " HIV-positive individuals do not have an estimate for age at infection\n" )
    cat(tmp[is.na(age_infection), round(mean(study_id %in% time.since.infection$study_id)*100,2)],
        "% of which are included in the TSI analysis and do not have sampling date\n")
    
  }else{
    
    cat("# ASSUME INFECTION 1y PRIOR TO DIAGNOSIS\n")
    meta[, age_infection := age_first_positive - 1]
    meta[is.na(age_first_positive), age_infection := age_first_visit - 1]
    tmp <- unique(meta[, .(study_id,age_infection,round)])
    cat(tmp[age_infection <= 16, .N], "out of", tmp[, .N], "HIV-positive individuals are estimated to have been infected prior to 16 yo\n")
    
    cat(tmp[is.na(age_infection), .N], " HIV-positive individuals do not have an estimate for age at infection, from cohort:" )
    print.table(table(tmp[is.na(age_infection), round]))
  }
  
  # Transform date format
  cols <- grep("date_", names(meta), value=TRUE)
  meta[, (cols) := lapply(.SD, as.Date) , .SDcols=cols]
  
  if(! 'date_infection' %in% names(meta))
          meta[, date_infection := date_birth + 365*(age_infection), by = 'study_id']
  
  return(meta)
}

keep.likely.transmission.pairs <- function(dchain, threshold)
{
  
  # keep pairs with a score greater than some threshold
  dchain <- dchain[SCORE_LINKED>threshold]
  
  # find direction of transmission 
  dchain[SCORE_DIR_12 <= threshold & SCORE_DIR_21 <= threshold, EST_DIR:='unclear']
  dchain[SCORE_DIR_12 > threshold, EST_DIR:='12']
  dchain[SCORE_DIR_21 > threshold, EST_DIR:='21']
  
  # find source recipient
  dchain <- dchain[EST_DIR != 'unclear']
  dchain[, `:=` (SOURCE=H1, RECIPIENT=H2)]
  dchain[EST_DIR == '21', `:=` (SOURCE=H2, RECIPIENT=H1) ]
  dchain[, `:=` (H1=NULL, H2=NULL)]
}

pairs.get.meta.data <- function(chain, metadata, aik)
{

  meta <- copy(metadata)

  # individuals without meta data
  idx <- chain[, unique(c(SOURCE, RECIPIENT)) ]
  missing_indiv <- .aid2pt(idx[! idx %in% meta$aid], aik)
  length(missing_indiv)

  cat('There are ', length(missing_indiv), 'indivs without meta data:\n' )
  missing_indiv <- grep('RK-', missing_indiv, value=T)
  cat('- ', length(missing_indiv), 'of which are in the Rakai Cohort.\n')
  cat(missing_indiv, '\n')
  
  # stopifnot(unique(dchain$SOURCE) %in% unique(meta_data$aid))
  # stopifnot(unique(dchain$RECIPIENT) %in% unique(meta_data$aid))
  
  # remove neuro individual
  meta <- meta[round != 'neuro']
  
  # merge by source, then recipient
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.SOURCE')
  tmp1 <- merge(chain, tmp, by.x = 'SOURCE', by.y = 'aid.SOURCE')
  
  tmp <- copy(meta)
  names(tmp) = paste0(names(tmp), '.RECIPIENT')
  tmp1 <- merge(tmp1, tmp, by.x = 'RECIPIENT', by.y = 'aid.RECIPIENT')
  
  # keep meta data from round the closest to the date at infection
  tmp1[, diff_date_sample_transmission.SOURCE := abs(sample_date.SOURCE - date_infection.RECIPIENT)]
  tmp1[, diff_date_sample_infection.RECIPIENT := abs(sample_date.RECIPIENT - date_infection.RECIPIENT)]
  tmp1[, select_this_date := (diff_date_sample_transmission.SOURCE == min(diff_date_sample_transmission.SOURCE) &
         diff_date_sample_infection.RECIPIENT == min(diff_date_sample_infection.RECIPIENT)), 
       by = c('SOURCE', 'RECIPIENT')]
  # In 1 case, there are 2 dates satisfying the above equality. Take any:
  tmp1 <- tmp1[select_this_date == T, lapply(.SD, `[[`, 1), by=c('SOURCE', 'RECIPIENT')]

  set(tmp1, NULL, 'select_this_date', NULL)
  set(tmp1, NULL, 'diff_date_sample_transmission.SOURCE', NULL)
  set(tmp1, NULL, 'diff_date_sample_infection.RECIPIENT', NULL)
  
  # check
  tmp2 <- chain[SOURCE %in% meta$aid & RECIPIENT %in% meta$aid]

  stopifnot(nrow(tmp1) == nrow(tmp2))

  # find age transmission source
  tmp1[, AGE_TRANSMISSION.SOURCE := .year.diff(date_infection.RECIPIENT, date_birth.SOURCE)]
  
  # upper column name
  colnames(tmp1) <- toupper(colnames(tmp1))
  
  return(tmp1)
}

find.center.of.mass.plausible.region <- function(xmin, xmax, ymin, ymax)
{
  # check
  stopifnot(xmin < xmax)
  stopifnot(ymin < ymax)
  
  # helpers
  find.center.of.mass.right.triangle.with.identity <- function( z )
  {
    x <- z[1]; y <- z[2]
    dyx <- as.integer(y - x)
    
    com <- c(x + dyx/3, y - dyx/3 )
    
    if(0) # just checking
    {
      dt <- data.table(x=x, y=y)
      ggplot(dt) + 
        geom_abline(slope=1, color='red', linetype='dotted') +
        geom_point(aes(x=x, y=y)) + 
        geom_point(aes(x=y, y=y)) + 
        geom_point(aes(x=x, y=x)) + 
        geom_point(aes(x=dom[1], y=dom[2]), pch='X') +
        theme_bw()
    }
    
    return(list(com=com, A=dyx^2/2) )
  }
  
  find.center.of.mass.rectangle.from.segment.to.identity <- function(z0, z1) 
  {
    
    if(z0[1] == z1[1] & z0[2] == z1[2])
      stop('Two points specifying side should be distinct')
    
    if(z0[2] < z0[1]) 
      stop('Segment must lie above identity y=x line')
    if(z1[2] < z1[1])
      stop('Segment must lie above identity y=x line')
    
    if(z0[1] != z1[1] & z0[2] != z1[2])
      stop('Input segment must be horizontal or vertical')
    
    x <- z1[1]
    
    
    ymin=min(z0[2], z1[2])
    ymax=max(z0[2], z1[2])
    xmin=min(z0[1], z1[1])
    xmax=max(z0[1], z1[1])
    
    # If input segment is horizontal
    if(ymax==ymin)
    {
      ymin <- xmax
      p3 <- c(xmin, ymin)
    }
    # If input segment is vertical
    if(xmin==xmax)
    {
      xmax <- ymin
      p3 <- c(xmax, ymax)
    }
    
    com <- c(xmin + .5 * as.integer(xmax - xmin), 
             ymin + .5 * as.integer(ymax - ymin))
    
    if(0)
    {
      dt <- data.table(
        ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax
      )
      
      ggplot(dt, aes(ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)) + 
        geom_rect(col='black') + 
        geom_point(aes(y=com[2], x=com[1])) +
        geom_abline(slope=1, color='red', linetype='dotted')
    }
    
    Area = as.integer(ymax - ymin) * as.integer(xmax - xmin)
    return(list(com=com, A=Area, p3=p3))
  }
  
  merge.centers.of.gravity <- function(input1, input2)
  {
    # input lists must contain Area A and Center of Gravity com
    A = input1$A + input2$A
    w1 <- input1$A/A 
    w2 <- 1 - w1
    
    # cog <- w1 * input1$com + w2 * input2$com 
    cog <- c(
      weighted.mean( x=c(input1$com[1],input2$com[1]), w=c(input1$A, input2$A) ),
      weighted.mean( x=c(input1$com[2],input2$com[2]), w=c(input1$A, input2$A) )
    )
    
    return(cog)
  }
  
  if.case.1 <- function(z=c(xmin, ymax))
    find.center.of.mass.right.triangle.with.identity(z)$com
  
  if.case.2 <- function()
  {
    z1 <- c(xmin, ymax)
    if(ymax > xmax){ # horizontal
      z2 <- c(xmax, ymax)
    }else{           # vertical
      z2 <- c(xmin, ymin)
    }
    
    out1 <- find.center.of.mass.rectangle.from.segment.to.identity(z1, z2)
    out2 <- find.center.of.mass.right.triangle.with.identity(out1$p3)
    
    merge.centers.of.gravity(out1, out2)
  }
  
  if.case.3 <- function() 
  {
    # solve using COM equations
    z <- c(xmax, ymin)
    triangle <- find.center.of.mass.right.triangle.with.identity(z)
    
    dx <- as.integer(xmax - xmin)
    dy <- as.integer(ymax - ymin)
    rect <- list(com = c(xmin + .5*dx, ymin + .5*dy), A=dx*dy)
    
    w = c(rect$A, -triangle$A)
    
    # com <- (rect$com * rect$A - triangle$com * triangle$A)/(rect$A - triangle$A)
    cog <- c(
      weighted.mean( x=c(rect$com[1], triangle$com[1]), w=w),
      weighted.mean( x=c(rect$com[2], triangle$com[2]), w=w)
    )
    
    cog
  }
  
  
  # function to transform 2d vecs into lists as wanted
  .f <- function(z)
    list(x=z[1], y=z[2])
  
  # classify geometric cases
  count <- sum( ymax > xmin, ymax > xmax, ymin > xmin, ymin > xmax)
  
  if(count == 0){
    cog <- list(NA_Date_, NA_Date_)
  }else if(count == 1){
    cog <- .f(if.case.1())
  }else if(count == 2){
    cog <- .f(if.case.2())
  }else if(count == 3){
    cog <- .f(if.case.3())
  }else if(count == 4){
    cog <- list(NA_Date_, NA_Date_)
  }
}

update.meta.pairs.after.doi.attribution <- function(path, outdir, overwrite = F)
{
  if(!file.exists(path) | overwrite)
  {
    cat('\tRunning infection dates estimation...\n\n\n')
    output <- study_infection_pairs_dates(CHAIN=chain)
    saveRDS(output, path)
  }else{
    cat('\tLoading previously estimated dates of infection...\n')
    output <- readRDS(path)
  }
  
  # extract info
  dinfrange <- output$pairs
  dage <- output$age
  
  # plot
  g <- plot.phylopair.dates.scores(dinfrange,
                                   add.dots=TRUE,
                                   daterange.vars=c('MIN','MAX'),
                                   doi.center.var='DOI',
                                   title='Midpoint/Network Attribution')
  filename=paste0(outdir, 'phylopairs_dates_scores_after_networkattribution.png')
  ggsave(filename, g, w=15, h=15, units='cm')
  
  # fix meta data
  aik2 <- copy(aik)
  names(aik2) <- tolower(names(aik2))
  
  dage[, AID := .pt2aid(study_id, aik2)]
  
  meta_data <- merge(meta_data[, -c('date_infection', 'age_infection')],
                     dage[, .(date_infection=DOI, age_infection=age_at_infection, aid=AID)],
                     by='aid', all.x=T)
  
  cols <- grep('^RECIPIENT$|^SOURCE$|SCORE|CLU', names(dinfrange), value=TRUE)
  chain <- dinfrange[, ..cols,]
  cols <- c('SOURCE', 'RECIPIENT')
  chain[, (cols) := lapply(.SD, function(x) .pt2aid(x, aik2)) , .SDcols=cols, by=cols]
  
  pairs.all <- pairs.get.meta.data(chain, meta_data, aik)
  
  list(pairs.all=pairs.all, chain=chain, meta_data=meta_data)
}

study_infection_pairs_dates <- function(CHAIN)
{
  cat('== Estimating network-coherent infection dates == \n')
  .is.mrc <- function(x)
    grepl('^MRC',x)
  
  
  .aid.2.studyid <- function(DT, by_var)
  {
    tmp <- merge(DT, aik, all.x=TRUE, by.x=by_var, by.y='AID' )
    stopifnot(tmp[, all(!is.na(PT_ID))])
    tmp[, (by_var) := NULL]
    setnames(tmp, 'PT_ID', by_var )
    tmp
  }
  
  .double.merge <- function(DT1, DT2, by_var='study_id')
  {
    cols <- setdiff(names(DT2), c('AID', 'study_id'))
    cols.SOURCE <- paste0(cols, '.SOURCE')
    cols.RECIPIENT <- paste0(cols, '.RECIPIENT')
    tmp <- merge(DT1, DT2, by.x='SOURCE', by.y=by_var, all.x=TRUE)
    setnames(tmp, cols, cols.SOURCE)
    tmp <- merge(tmp, DT2, by.x='RECIPIENT', by.y=by_var, all.x=TRUE)
    setnames(tmp, cols, cols.RECIPIENT)
    tmp
  }
  
  # Get dates and ages data
  .select <- function(cols)
    unique(meta_data[, ..cols])
  ddates <- .select(c('study_id', 'date_last_negative', 'date_infection', 'date_first_positive')) 
  dage <- .select(c('study_id', 'age_first_positive', 'date_birth')) 
  
  .check <- function(DT)
    DT[, .N, by='study_id'][, all(N == 1)]
  stopifnot(.check(dage) & .check(ddates))
  
  dage <- merge(dage, aik, all.x=TRUE, by.x='study_id', by.y='PT_ID')
  ddates <- merge(ddates, aik, all.x=TRUE, by.x='study_id', by.y='PT_ID')
  
  chain2 <- copy(CHAIN)
  
  # remove uninteresting cols
  chain2[, CNTRL_ANY := CNTRL1 | CNTRL2]
  cols <- grep('CATEGORISATION|PTY_RUN|CNTRL[0-9]', names(chain2), value=TRUE)
  cols1 <- grep('^LINK|EST_DIR',names(chain2), value=TRUE)
  chain2[,  (cols) := NULL ]
  chain2[SCORE_DIR_12 == SCORE_DIR_21, cat('- There are ', .N, 'pairs where direction of transmission is exactly equally likely\n') ]
  # if(0) chain2[SCORE_DIR_12 > SCORE_DIR_21 , unique(.SD), .SDcols=cols1]
  chain2[, (cols1) := NULL]
  chain2[, all(round(SCORE_DIR_12 + SCORE_DIR_21, 2) == 1)]
  chain2[, SCORE_DIR_SR := pmax(SCORE_DIR_12, SCORE_DIR_21)]
  cols <- grep('SCORE_NW|SCORE_DIR_[0-9][0-9]', names(chain2), value=TRUE)
  chain2[, (cols):= NULL ]
  
  # round entries for readability
  cols <- names(which( unlist(lapply(chain2, is.numeric)) ))
  chain2[, (cols) := lapply(.SD, function(x) round(x,2)) , .SDcols=cols]
  
  # Translate AIDs to study_ids
  chain2 <- .aid.2.studyid(chain2, by_var='SOURCE')
  chain2 <- .aid.2.studyid(chain2, by_var='RECIPIENT')
  
  # merge with dates
  chain2 <- .double.merge(DT1=chain2, DT2=ddates[, -'AID'])
  stopifnot(
    chain2[is.na(date_infection.RECIPIENT), all(.is.mrc(RECIPIENT))],
    chain2[is.na(date_infection.SOURCE), all(.is.mrc(SOURCE))]
  ); cat('- CHECK: All individuals in infection pairs without ddates are in the MRC\n')
  cat('\t(removing them...)\n')
  chain2 <- chain2[! (.is.mrc(SOURCE) | .is.mrc(RECIPIENT))]
  
  if(0)
  {
    plot.phylopair.dates.scores(chain2, title='Testing history')
    plot.phylopair.dates.scores(chain2, only.coherent=TRUE)
    plot.phylopair.dates.scores(chain2, only.contradict=TRUE)
    plot.phylopair.dates.scores(chain2, only.crossing=TRUE)
    plot.phylopair.dates.scores(chain2, only.contradict=TRUE, only.rect.pairs=FALSE)
  }
  
  tmp <- chain2[, date_last_negative.SOURCE >= date_first_positive.RECIPIENT]
  cat('- There are ', sum(tmp == TRUE, na.rm=TRUE), ' couples with inconsistent first positive and last negative dates...\n')
  
  # Use age to give a min bound 
  # setting date of last negative to be the minimum of 15th birthday and age first pos
  chain2 <- .double.merge(DT1=chain2, DT2=dage[, -c('AID', 'age_first_positive')])
  cols <- grep('date_birth', names(chain2), value=TRUE)
  cols1 <- gsub('birth', '15yr', cols)
  chain2[, (cols) := lapply(.SD, function(x) x+as.integer(365.25*15)) , .SDcols=cols]
  setnames(chain2, cols, cols1)
  cols <- grep('date_last', names(chain2), value=TRUE)
  chain2[, date_last_negative.SOURCE := pmax(date_last_negative.SOURCE, date_15yr.SOURCE, na.rm=T) , by='SOURCE']
  chain2[, date_last_negative.RECIPIENT := pmax(date_last_negative.RECIPIENT, date_15yr.RECIPIENT, na.rm=T) , by='RECIPIENT']
  
  # check
  # cols <- grep('date_first|date_last', names(chain2), value=TRUE)
  # chain2[, lapply(.SD,is.na) , .SDcols=cols][, lapply(.SD, sum) , .SDcols=cols]
  cols_src <- grep('SOURCE', names(chain2), value=TRUE)
  cols_rcp <- grep('RECIPIENT', names(chain2), value=TRUE)
  stopifnot(chain2[, all(date_first_positive.RECIPIENT >= date_last_negative.RECIPIENT, na.rm=T)])
  stopifnot(chain2[, all(date_first_positive.SOURCE >= date_last_negative.SOURCE, na.rm=T)])
  
  if(0) # date of infection does may not fall in range if we use "1 yr prior" DoI.
  {
    chain2[date_infection.RECIPIENT < date_last_negative.RECIPIENT, ..cols_rcp]
    chain2[date_infection.SOURCE < date_last_negative.SOURCE, ..cols_src]
    plot.phylopair.dates.scores(chain2[date_infection.SOURCE < date_last_negative.SOURCE])
    plot.phylopair.dates.scores(chain2, only.coherent=TRUE)
    plot.phylopair.dates.scores(chain2, only.contradict=TRUE)
    plot.phylopair.dates.scores(chain2, only.crossing=TRUE)
    plot.phylopair.dates.scores(chain2, only.contradict=TRUE, only.rect.pairs=FALSE)
  }
  
  cat('- Change direction of tranmsission for couple with contradicting times of infection...\n')
  
  cols <- grep('SOURCE|RECIPIENT', names(chain2), value=TRUE)
  .swap.src.rec <- function(x)
  {
    y <- copy(x)
    src <- grep('SOURCE',x)
    rcp <- grep('RECIPIENT', x)
    y[src] <- gsub('SOURCE','RECIPIENT',x[src])
    y[rcp] <- gsub('RECIPIENT','SOURCE',x[rcp])
    y
  }
  cols1 <- .swap.src.rec(cols)
  # chain2[, (uniqueN(RECIPIENT) == .N) ]
  idx <- chain2[ date_last_negative.SOURCE >= date_first_positive.RECIPIENT, SOURCE]
  idx0 <- idx[! idx %in% chain2$RECIPIENT]
  idx1 <- idx[ idx %in% chain2$RECIPIENT]
  cat('\tonly if the SOURCE does not appear as a RECIPIENT elsewhere(',length(idx1),'\n')
  chain2[ date_last_negative.SOURCE >= date_first_positive.RECIPIENT & SOURCE %in% idx0,
          (cols1) := (.SD), .SDcols=cols]
  cat('\t(exclude those others)')
  chain2 <- chain2[ date_last_negative.SOURCE < date_first_positive.RECIPIENT ]
  
  # individuals involved in multiple events:
  # SOURCEs should be infected 
  
  cat('- Attribute date of infections based on transmission network structure...\n')
  
  shrink.ranges.doi.coherently.to.network <- function(CHAIN)
  {
    cat('Shrinking date of infection dates based on network relationships...\n')
    get_range_doi_by_studyid <- function(DT)
    {
      cols <- c('ID', 'MIN', 'MAX')
      cols0 <- c('SOURCE', "date_last_negative.SOURCE", "date_first_positive.SOURCE")
      cols1 <- c('RECIPIENT', "date_last_negative.RECIPIENT", "date_first_positive.RECIPIENT")
      tmp0 <- DT[, unique(.SD) , .SDcols=cols0]
      tmp1 <- DT[, unique(.SD) , .SDcols=cols1]
      setnames(tmp0, cols0, cols)
      setnames(tmp1, cols1, cols)
      tmp0[, uniqueN(ID) == .N]; tmp1[, uniqueN(ID) == .N]
      tmp <- unique(rbind(tmp0, tmp1))
      stopifnot(tmp[, uniqueN(ID) == .N])
      return(tmp)
    }
    
    get_doi <- function(s, col)
    {
      stopifnot(length(col)==1)
      dpairs_doi[ID %in% s, ..col ][[1]]
    }
    
    shrink_max_doi_source <- function(PAIR, DATE)
    {
      # for each source, the MAX doi must preceed the MAX doi of the recipients
      tmp <- PAIR[, list(MAX2 = min(get_doi(RECIPIENT, "MAX"))) ,by=SOURCE]
      DATE <- merge( DATE, tmp, all.x=T, by.y='SOURCE', by.x='ID')
      
      if(DATE[ MAX2 < MAX, .N == 0])
      {
        DATE[, MAX2 := NULL]
        return(list(DATE=DATE, done=TRUE))
      }
      
      DATE[ MAX2 < MAX, MAX := MAX2]
      DATE[, MAX2 := NULL]
      list(DATE=DATE, done=FALSE)
    }
    
    shrink_min_doi_recipient <- function(PAIR, DATE)
    {
      # for each recipient, the MIN doi must succeed the MIN doi of the source
      tmp <- PAIR[, list(MIN2 = get_doi(SOURCE, "MIN")), by=RECIPIENT]
      DATE <- merge(DATE,tmp, all.x=T, by.y='RECIPIENT', by.x='ID')
      
      if(DATE[ MIN2 > MIN, .N == 0])
      {
        DATE[, MIN2 := NULL]
        return(list(DATE=DATE, done=TRUE))
      }
      
      DATE[ MIN2 > MIN, MIN:=MIN2]
      DATE[, MIN2 := NULL]
      list(DATE=DATE, done=FALSE)
    }
    
    dpairs <- chain2[, .(SOURCE, RECIPIENT)] 
    dpairs_doi <- get_range_doi_by_studyid(CHAIN)
    
    if(dpairs_doi[, any(is.na(MAX))])
    {
      warning('Setting unknown date first positives as max date. Need changing\n')
      tmp <- dpairs_doi[, max(MAX, na.rm=T)]
      dpairs_doi[is.na(MAX), MAX := tmp]
    }
    
    done <- FALSE; count=0
    while(done==FALSE)
    {
      count <- count + 1
      cat('\tIteration number:',count,'...\n')
      tmp <- shrink_max_doi_source(dpairs, dpairs_doi)
      dpairs_doi <- tmp$DATE
      done <- tmp$done
      tmp <- shrink_min_doi_recipient(dpairs, dpairs_doi)
      dpairs_doi <- tmp$DATE
      done <- done & tmp$done
    }
    
    return(list(P=dpairs, D=dpairs_doi))
  }
  
  tmp <- shrink.ranges.doi.coherently.to.network(chain2)
  dpairs <- copy(tmp$P)
  dpairs_doi <- copy(tmp$D)
  rm(tmp)
  setnames(dpairs_doi, 'ID', 'study_id')
  
  if(0)
  {
    tmp <- .double.merge(DT1=chain2, dpairs_doi)
    idx <- tmp[, .N, by='SOURCE'][N==2, SOURCE]
    plot.phylopair.dates.scores(tmp[SOURCE%in%idx[3]], daterange.vars=c('MIN', 'MAX'), add.dots=FALSE)
    plot.phylopair.dates.scores(tmp[SOURCE%in%idx[3]], add.dots=FALSE)
    plot.phylopair.dates.scores(tmp, daterange.vars = c('MIN', 'MAX'), add.dots=FALSE) 
    plot.phylopair.dates.scores(tmp, daterange.vars = c('MIN', 'MAX'), add.dots=FALSE)
    plot.phylopair.dates.scores(tmp,daterange.vars = c('MIN', 'MAX'), only.contradict = T, add.dots=FALSE)
    plot.phylopair.dates.scores(tmp,daterange.vars = c('MIN', 'MAX'), only.coherent=T, add.dots=FALSE)
  }
  
  
  # now need to think about how to attribute the date of infection...
  # for simple case: midpoint attribution
  # for sources with multiple recipients: generalisation averaging pfds
  # for chains of A -> B -> C ... attribute the last DoI, and then "backpropagate" with midpoint assignment
  
  assign_dates_of_infection <- function(dpairs, dpairs_doi)
  {
    # attribute "depth" of transmission chain
    depth <- 1
    dpairs[, DEPTH := depth]
    tmp <- dpairs[DEPTH == depth, unique(SOURCE)]
    while(length(tmp))
    {
      depth <- depth + 1
      dpairs[ RECIPIENT %in% tmp, DEPTH := depth]
      tmp <- dpairs[DEPTH == depth, unique(SOURCE)]
    }
    dpairs[, table(DEPTH)]
    
    traceback_depth <- function(src, rcp, arrows=F, df=T)
    {
      tmp <- dpairs[SOURCE == src & RECIPIENT == rcp]
      if(tmp$DEPTH == 1)
      {
        if(df)
          return(data.frame(SOURCE=scr, RECIPIENT=rcp, DEPTH=tmp$DEPTH))
        out <- c(src, rcp)
        if(arrows)
          out <- paste0(out, collapse=' -> ')
        return(out)
      }
      
      tmp <- dpairs[SOURCE == rcp & DEPTH == tmp$DEPTH-1, .(SOURCE, RECIPIENT)]
      
      if(df)
        return( rbind(
          data.frame(SOURCE=src, RECIPIENT=rcp, DEPTH=tmp$DEPTH),
          traceback_depth(tmp$SOURCE, tmp$RECIPIENT, df=T))
        )
      
      out <- c( rcp, traceback_depth(tmp$SOURCE, tmp$RECIPIENT) )
      if(arrows)
        out <- paste0(out, collapse=' -> ')
      return(out)
    }
    
    # get doi source | recipients bounds, and then assign midpoint for recps.
    assign.date.infection.sources <- function(src, d, precision=50)
    {
      
      # src <- dpairs[,.N, by='SOURCE'][N==2, SOURCE[1]]
      # cat(src, '\n')
      tmp <- dpairs[SOURCE == src & DEPTH == d, -"DEPTH"]
      tmp <- .double.merge(tmp, dpairs_doi)
      
      # if any source has a known date of infection, 
      # assign midpoint of plausible interval
      if(tmp[!is.na(DOI.RECIPIENT), .N])
      {
        x_max <- tmp[!is.na(DOI.RECIPIENT), min(DOI.RECIPIENT)]
        return(tmp[, mean.Date(c(MIN.SOURCE[1], x_max))])
      }
      
      # there may be some analytic formula, but lets just use numeric integration
      x_range <- tmp[, seq(from=MIN.SOURCE[1], to=MAX.SOURCE[1], length.out = precision)]
      rcp_cross_section_length <- function(x)
      {
        tmp1 <- tmp[, .(L = as.numeric(MAX.RECIPIENT - max(MIN.RECIPIENT, x))/365.25), by='RECIPIENT']
        tmp1[, prod(L)]
      }
      x_pdf <- sapply(x_range, rcp_cross_section_length)
      x_cdf <- cumsum(x_pdf)/sum(x_pdf)
      x_range[which.min(abs(x_cdf - .5))]
    }
    
    assign.date.infection.recipient <- function(rcp, d)
    {
      # rcp <- idx[1]
      # cat(rcp, '\n')
      src <- dpairs[RECIPIENT == rcp, SOURCE]
      dsrc <- dpairs_doi[study_id == src]
      drcp <- dpairs_doi[study_id == rcp]
      x_min = dsrc$DOI
      x_min = max(x_min, drcp$MIN)
      x_mid = mean.Date(c(x_min, drcp$MAX))
      x_mid
    }
    
    
    dpairs_doi[, DOI := NA_Date_]
    depth <- dpairs[, max(DEPTH)] + 1
    
    while(depth > 0 & dpairs_doi[, any(is.na(DOI))])
    {
      depth <- depth - 1
      cat('Assigning dates of infection for pairs at depth:', depth,'...\n')
      
      dpairs_doi[study_id %in% dpairs[DEPTH==depth, SOURCE] & is.na(DOI),
                 DOI := assign.date.infection.sources(study_id, depth),
                 by='study_id']
      cat(dpairs_doi[!is.na(DOI), .N], '\n')
      
      dpairs_doi[study_id %in% dpairs[DEPTH==depth, RECIPIENT] & is.na(DOI),
                 DOI := assign.date.infection.recipient(study_id, depth),
                 by='study_id']
      cat(dpairs_doi[!is.na(DOI), .N ],'\n')
      
      .check <- function(src, rcp)
      {
        tmp_s0 <- dpairs_doi[study_id == src, DOI >= MIN & DOI <= MAX]
        tmp_r0 <- dpairs_doi[study_id == rcp, DOI >= MIN & DOI <= MAX]
        tmp_s1 <- dpairs_doi[study_id == src, DOI]
        tmp_r1 <- dpairs_doi[study_id == rcp, DOI]
        tmp_s0 & tmp_r0 & (tmp_s1 <= tmp_r1)
      }
      tmp <- dpairs[, .check(SOURCE, RECIPIENT), by=c('SOURCE', 'RECIPIENT')]
      stopifnot(tmp[, all(V1==T, na.rm=T)])
    }
    dpairs_doi
  }
  
  # prepare output
  
  dpairs_doi <- assign_dates_of_infection(dpairs, dpairs_doi)
  
  dage <- merge(dage[, .(study_id, date_birth)], 
                dpairs_doi[, .(study_id, DOI)], by='study_id')
  dage[, age_at_infection := .year.diff(DOI, date_birth)]
  
  list(pairs= .double.merge(DT1=chain2, dpairs_doi), age=dage)
}

print.statements.about.pairs <- function(pairs)
{
  
  cat('\nThere is ', nrow(pairs), ' source-recipient pairs\n\n')
  
  cat(nrow(pairs[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]), ' pairs have a proxy for the age at infection of the source and recipient\n')
  cat(nrow(pairs[((SEX.SOURCE == 'F' & SEX.RECIPIENT == 'M') | (SEX.SOURCE == 'M' & SEX.RECIPIENT == 'F')) & (!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT))]), ' pairs are heteroxuals have a proxy for the time of infection of the source and recipient\n\n')                
  
  cat('\nPairs by sex')
  tab <- pairs[, list(count = .N), by = c('SEX.SOURCE', 'SEX.RECIPIENT')]
  print_table(tab)
  
  cat('\nPairs by community')
  tab <- pairs[, list(count = .N), by = c('COMM.SOURCE', 'COMM.RECIPIENT')]
  print_table(tab)
  
}

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

inv_logit <- function(x) 1 / (1 + exp(-x))
logit <- function(p) log(p / (1-p))

get_proportion_sampling <- function(pairs, incidence_cases, incidence_cases_round, outdir, byperiod = T){
  
  # proportion sampling by round
  proportion_sampling_round <- get_proportion_sampling_by_round(pairs, incidence_cases_round)
  
  # proportion sampling by period
  proportion_sampling_period <- get_proportion_sampling_by_period(pairs, incidence_cases, outdir)
  
  if(byperiod){
    return(proportion_sampling_period)
  }else{
    return(proportion_sampling_round)
  }
    
}

get_proportion_sampling_by_round <- function(pairs, incidence_cases_round)
{
  
  # find number of pairs observed
  dp <- copy(pairs)
  
  # we add pairs to the corresponding rounds (move their date of infection so that it falls within the observational period)
  tmp <- df_round[, .(ROUND, MIN_SAMPLE_DATE, MAX_SAMPLE_DATE, COMM)]
  tmp[, MIN_SAMPLE_DATE_NEXT_ROUND := c(MIN_SAMPLE_DATE[2:nrow(tmp)], MAX_SAMPLE_DATE[nrow(tmp)])]
  stopifnot(tmp[, all(MAX_SAMPLE_DATE <= MIN_SAMPLE_DATE_NEXT_ROUND)])
  
  dp <- merge(dp, tmp,  by.x = 'COMM.RECIPIENT', by.y = 'COMM', allow.cartesian=TRUE)
  dp <- dp[DATE_INFECTION.RECIPIENT >= MIN_SAMPLE_DATE & DATE_INFECTION.RECIPIENT <= MIN_SAMPLE_DATE_NEXT_ROUND]
  stopifnot(nrow(dp) == nrow(pairs))
  set(dp, NULL, 'MIN_SAMPLE_DATE', NULL)
  set(dp, NULL, 'MAX_SAMPLE_DATE', NULL)
  set(dp, NULL, 'MIN_SAMPLE_DATE_NEXT_ROUND', NULL)
  
  setnames(dp, c('SEX.RECIPIENT', 'COMM.RECIPIENT', 'AGE_INFECTION.RECIPIENT'), 
           c('SEX', 'COMM', 'AGEYRS'))
  dp[, AGEYRS := floor(AGEYRS)]
  dp <- dp[, list(count = .N), by = c('SEX', 'COMM', 'AGEYRS', 'ROUND')]
  
  # merge to incidence cases
  di <- merge(incidence_cases_round, dp, by = c('SEX', 'COMM', 'AGEYRS', 'ROUND'), all.x = T)
  di[is.na(count), count := 0]
  
  # find empiriral proportion sampling from a source of any age to a recipient aged j
  di[, prop_sampling := count / INCIDENT_CASES]
  setnames(di, c('SEX', 'AGEYRS'), c('SEX.RECIPIENT', 'AGEYRS.RECIPIENT'))
  
  # warnings
  tmp <- di[prop_sampling > 1]
  if(nrow(tmp) > 0){
    cat('\n Some probabilities are greater than 1')
    cat('\n In', tmp[, unique(COMM)], 'communities at round', tmp[, unique(ROUND)])
  }
  tmp <- di[prop_sampling < 0]
  if(nrow(tmp) > 0){
    cat('\n Some probabilities are smaller than 0')
    cat('\n In', tmp[, unique(COMM)], 'communities at round', tmp[, unique(ROUND)])
  }
  
  # save
  file.name <- file.detection.probability.round
  if(! file.exists(file.name) | config$overwrite.existing.files )
  {
    cat("Saving file:", file.name, '\n')
    saveRDS(di, file = file.name)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
  
  return(di)
}

get_proportion_sampling_by_period <- function(pairs, incidence_cases, outdir)
{
  
  # find number of pairs observed
  dp <- copy(pairs)
  setnames(dp, c('SEX.RECIPIENT', 'COMM.RECIPIENT', 'AGE_INFECTION.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT'), 
           c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF'))
  dp[, AGEYRS := floor(AGEYRS)]
  dp <- dp[, list(count = .N), by = c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF')]
  
  # merge to incidence cases
  di <- merge(incidence_cases, dp, by = c('SEX', 'COMM', 'AGEYRS', 'BEFORE_CUTOFF'), all.x = T)
  di[is.na(count), count := 0]
  
  # find empiriral proportion sampling from a source of any age to a recipient aged j
  di[, prop_sampling := count / INCIDENT_CASES]
  setnames(di, c('SEX', 'AGEYRS'), c('SEX.RECIPIENT', 'AGEYRS.RECIPIENT'))
  
  # enlarge to all sources
  tmp <- di[, .(SEX.RECIPIENT, COMM, AGEYRS.RECIPIENT, BEFORE_CUTOFF)]
  setnames(tmp, c('SEX.RECIPIENT', 'AGEYRS.RECIPIENT'), c('SEX.SOURCE', 'AGEYRS.SOURCE'))
  df <- merge(di, tmp, by = c('BEFORE_CUTOFF', 'COMM'), allow.cartesian = T)
  df <- df[SEX.RECIPIENT != SEX.SOURCE]
  
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

  
  if(dir.exists(dirname( outdir)))
  { # plots
    
    tmp1 <- copy(df)
    tmp1[, Direction := 'Female -> Male']
    tmp1[SEX.RECIPIENT == 'F', Direction := 'Male -> Female']
    tmp1[, INDEX_TIME2 := paste0('Period: ', INDEX_TIME)]
    
    # age recipient plot
    ggplot(tmp1, aes(col = Direction)) + 
      geom_line(aes(x = AGEYRS.RECIPIENT, y = prop_sampling)) + 
      facet_grid(INDEX_TIME2~COMM + Direction) + 
      scale_y_continuous(labels = scales::percent)  + 
      theme_bw() +
      labs(y = 'Probability of observing transmission event', x = 'Age recipient')
    ggsave(paste0(outdir, '-data-proportion_sampling_empirical_recipient_period.png'), w = 8, h = 7)
    
    # age source plot
    tmp2 <- tmp1[, list(prop_sampling = mean(prop_sampling)), by = c('INDEX_TIME2', 'COMM', 'AGEYRS.SOURCE', 'Direction')]
    ggplot(tmp2, aes(col = Direction)) + 
      geom_line(aes(x = AGEYRS.SOURCE, y = prop_sampling)) + 
      facet_grid(INDEX_TIME2~COMM + Direction) + 
      scale_y_continuous(labels = scales::percent)  + 
      theme_bw() +
      labs(y = 'Probability of observing transmission event', x = 'Age source')
    ggsave(paste0(outdir, '-data-proportion_sampling_source_period.png'), w = 8, h = 7)
    
    # age source-recipient prepare pairs for plot
    dp <- pairs[, .(SEX.RECIPIENT, COMM.RECIPIENT, AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT)]
    setnames(dp, c('COMM.RECIPIENT', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT', 'DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT'), 
             c('COMM', 'AGEYRS.SOURCE', 'AGEYRS.RECIPIENT', 'BEFORE_CUTOFF'))
    dp[, `:=` (AGEYRS.RECIPIENT = floor(AGEYRS.RECIPIENT), AGEYRS.SOURCE = floor(AGEYRS.SOURCE))]
    dp <- merge(dp, unique(incidence_cases[, .(BEFORE_CUTOFF, INDEX_TIME, COMM)]), by = c('BEFORE_CUTOFF', 'COMM'))
    dp[, Direction := 'Female -> Male']
    dp[SEX.RECIPIENT == 'F', Direction := 'Male -> Female']
    dp[, INDEX_TIME2 := paste0('Period: ', INDEX_TIME)]
    
    ggplot(tmp1, aes(x = AGEYRS.SOURCE, y = AGEYRS.RECIPIENT)) + 
      geom_raster(aes(fill = prop_sampling)) + 
      facet_grid(INDEX_TIME2~COMM+Direction)  +
      scale_fill_viridis_c(labels = scales::percent)  + 
      # geom_point(data = dp, col = 'red') + 
      labs(x = 'Age source', y = 'Age recipient', fill = 'Probability of observing\ntransmission event') +
      scale_y_continuous(expand= c(0,0))+
      scale_x_continuous(expand= c(0,0)) +
      theme(strip.background = element_rect(colour="black", fill="white"),
            strip.text = element_text(size = rel(1)))
    ggsave(paste0(outdir, '-data-proportion_sampling_source_recipient_period.png'), w = 12, h = 8)
  }
  
  df
}

get.age.map <- function(age_bands_reduced = 4)
{
  
  extended_age_length <- 0
  
  ages_source <- data.table(min_age = 15 - extended_age_length, max_age = 49 + extended_age_length)
  ages_source <- ages_source[, list(age = min_age:max_age)]

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

get.df.direction <- function()
{
  df_direction <- data.table(INDEX_DIRECTION = 1:2, IS_MF = c(0, 1))
  df_direction[, LABEL_DIRECTION := ifelse(IS_MF == 1, 'Male -> Female', 'Female -> Male')]
  df_direction[, LABEL_SOURCE := ifelse(IS_MF == 1, 'Male sources', 'Female sources')]
  df_direction[, LABEL_RECIPIENT := ifelse(IS_MF == 1, 'Female recipients', 'Male recipients')]
  df_direction[, LABEL_GENDER_SOURCE := ifelse(IS_MF == 1, 'Men', 'Women')]
  df_direction[, LABEL_GENDER_RECIPIENT := ifelse(IS_MF == 1, 'Women', 'Men')]
  df_direction[, LABEL_TRANSMITTING_PARTNER := ifelse(IS_MF == 1, 'Male transmitting partner', 'Female transmitting partner')]
  df_direction[, LABEL_INFECTED_PARTNER := ifelse(IS_MF == 1, 'Female infected partner', 'Male infected partner')]
  
  df_direction
}

get.df.community <- function()
{
    df_community <- data.table(INDEX_COMMUNITY = 1, COMM = c('inland'))
    df_community[, LABEL_COMMUNITY :='Inland communities']

  df_community
}

get.group.map <- function(stratify.by.community.recipient, df_period)
{
  
  df_direction <- data.table(INDEX_DIRECTION = 1:2, IS_MF = c(0, 1))
  df_direction[, LABEL_DIRECTION := ifelse(IS_MF == 1, 'Male -> Female', 'Female -> Male')]
  
  tmp <- data.table(INDEX_TIME = 1:2, BEFORE_CUTOFF = c(T, F))
  df_period <- merge(df_period, tmp, by = 'BEFORE_CUTOFF')
  
  if(stratify.by.community.recipient)
  {
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

get.age.aggregated.map <- function(age_aggregated)
{
  

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

load_incidence_rates_samples <- function(file.incidence.samples.inland){
  
  incidence_rates_round.samples <- fread(file.incidence.samples.inland)
  incidence_rates_round.samples[,COMM := 'inland']
  incidence_rates_round.samples[, SEX := substr(Sex,1,1)]
  if('ROUND' %in% names(incidence_rates_round.samples)){
    incidence_rates_round.samples[, ROUND := gsub('Round: (.+)', '\\1', ROUND)]
  }else{
    setnames(incidence_rates_round.samples, 'round_label', 'ROUND')
  }
  incidence_rates_round.samples[, ROUND := paste0('R0', ROUND)]
  setnames(incidence_rates_round.samples, 'age', 'AGEYRS')
  setnames(incidence_rates_round.samples, 'inc', 'INCIDENCE.DRAW')
  if(!'iterations_within' %in% names(incidence_rates_round.samples)){
    incidence_rates_round.samples[, iterations_within := 1]
  }
  
  # iterations: iterations over 50 data with imputed date of infection
  # iterations within: iterations within dataset of estimated incidence rate using MLE mean/sd and assuming normality
  # subsample otherwise the memory is excausted
  incidence_rates_round.samples <- incidence_rates_round.samples[iterations_within %in%1:500]
  incidence_rates_round.samples[, iterations := paste0(iterations, '-', iterations_within)]
  
  # keep var of interest
  incidence_rates_round.samples <- incidence_rates_round.samples[, .(SEX, COMM, ROUND, AGEYRS, iterations, INCIDENCE.DRAW)]

  # add rounds
  incidence_rates_round.samples <- merge(incidence_rates_round.samples, df_round, by = c('ROUND', 'COMM'))
  
  return(incidence_rates_round.samples)
}

read_treatment_cascade <- function(file.treatment.cascade.prop.participants, 
                                         file.treatment.cascade.prop.nonparticipants){
  
  # PROP_SUPPRESSED_M: Proportion of suppressed among infected
  
  # participants
  treatment_cascade_participants <- fread(file.treatment.cascade.prop.participants)
  treatment_cascade_participants <- treatment_cascade_participants[, .(AGEYRS, SEX, COMM, ROUND, 
                                                                       PROP_SUPPRESSED_M, PROP_SUPPRESSED_CL, PROP_SUPPRESSED_CU)]
  setnames(treatment_cascade_participants, 'PROP_SUPPRESSED_M', 'PROP_SUPPRESSED_PARTICIPANTS_M')
  setnames(treatment_cascade_participants, 'PROP_SUPPRESSED_CL', 'PROP_SUPPRESSED_PARTICIPANTS_CL')
  setnames(treatment_cascade_participants, 'PROP_SUPPRESSED_CU', 'PROP_SUPPRESSED_PARTICIPANTS_CU')
  treatment_cascade_participants[, PROP_UNSUPPRESSED_PARTICIPANTS_M := 1 - PROP_SUPPRESSED_PARTICIPANTS_M]
  treatment_cascade_participants[, PROP_UNSUPPRESSED_PARTICIPANTS_CL := 1 - PROP_SUPPRESSED_PARTICIPANTS_CU]
  treatment_cascade_participants[, PROP_UNSUPPRESSED_PARTICIPANTS_CU := 1 - PROP_SUPPRESSED_PARTICIPANTS_CL]
  treatment_cascade_participants <- select(treatment_cascade_participants, -c('PROP_SUPPRESSED_PARTICIPANTS_M', 
                                                                              'PROP_SUPPRESSED_PARTICIPANTS_CU', 
                                                                              'PROP_SUPPRESSED_PARTICIPANTS_CL'))

  # non-participants
  treatment_cascade_nonparticipants <- fread(file.treatment.cascade.prop.nonparticipants)
  treatment_cascade_nonparticipants <- treatment_cascade_nonparticipants[, .(AGEYRS, SEX, COMM, ROUND, 
                                                                             PROP_SUPPRESSED_M, PROP_SUPPRESSED_CL, PROP_SUPPRESSED_CU)]
  setnames(treatment_cascade_nonparticipants, 'PROP_SUPPRESSED_M', 'PROP_SUPPRESSED_NONPARTICIPANTS_M')
  setnames(treatment_cascade_nonparticipants, 'PROP_SUPPRESSED_CL', 'PROP_SUPPRESSED_NONPARTICIPANTS_CL')
  setnames(treatment_cascade_nonparticipants, 'PROP_SUPPRESSED_CU', 'PROP_SUPPRESSED_NONPARTICIPANTS_CU')
  treatment_cascade_nonparticipants[, PROP_UNSUPPRESSED_NONPARTICIPANTS_M := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_M]
  treatment_cascade_nonparticipants[, PROP_UNSUPPRESSED_NONPARTICIPANTS_CL := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_CU]
  treatment_cascade_nonparticipants[, PROP_UNSUPPRESSED_NONPARTICIPANTS_CU := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_CL]

  treatment_cascade_nonparticipants <- select(treatment_cascade_nonparticipants, -c('PROP_SUPPRESSED_NONPARTICIPANTS_M', 
                                                                                    'PROP_SUPPRESSED_NONPARTICIPANTS_CU',
                                                                                    'PROP_SUPPRESSED_NONPARTICIPANTS_CL'))
  
  # merge
  treatment_cascade <- merge(treatment_cascade_participants, treatment_cascade_nonparticipants, 
                             by = c('AGEYRS', 'SEX', 'COMM', 'ROUND'))
  
  return(treatment_cascade)
}

read_treatment_cascade_samples <- function(file.treatment.cascade.prop.participants.samples, 
                                           file.treatment.cascade.prop.nonparticipants.samples){
  
  # posterior samples of the proportion of suppressed among infected
  
  # load treatment participants
  treatment_cascade_participants <- as.data.table(readRDS(file.treatment.cascade.prop.participants.samples))
  treatment_cascade_participants <- treatment_cascade_participants[, .(AGEYRS, SEX, COMM, ROUND, iterations, PROP_SUPPRESSED_POSTERIOR_SAMPLE)]
  setnames(treatment_cascade_participants, 'PROP_SUPPRESSED_POSTERIOR_SAMPLE', 'PROP_SUPPRESSED_PARTICIPANTS')
  
  # load treatment participants
  treatment_cascade_nonparticipants <- as.data.table(readRDS(file.treatment.cascade.prop.nonparticipants.samples))
  treatment_cascade_nonparticipants <- treatment_cascade_nonparticipants[, .(AGEYRS, SEX, COMM, ROUND, iterations, PROP_SUPPRESSED_POSTERIOR_SAMPLE)]
  setnames(treatment_cascade_nonparticipants, 'PROP_SUPPRESSED_POSTERIOR_SAMPLE', 'PROP_SUPPRESSED_NONPARTICIPANTS')
  
  # merge
  treatment_cascade <- merge(treatment_cascade_participants, treatment_cascade_nonparticipants, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'))
  
  return(treatment_cascade)
}

read_pairs <- function(file.pairs){
  pairs.all <- as.data.table(readRDS(file.pairs))
  setnames(pairs.all, 'M', 'DATE_INFECTION.RECIPIENT')
  pairs.all <- select(pairs.all, -c('CL', 'IL', 'IU', 'CU', 'ROUND.M', 'DIRECTION'))
  return(pairs.all)
}
 
select.pairs.for.analysis <- function(DPAIRS,
                                      only.one.community, 
                                      use_30com_pairs, 
                                      only.transmission.after.start.observational.period, 
                                      only.transmission.before.stop.observational.period, 
                                      remove.pairs.from.rounds
                                      )
{
    
    out <- copy(DPAIRS)

    sprintfcat <- function(fmt, ...)
        cat(sprintf(fmt=fmt, ...), '\n\n')

    n <- nrow(out)
    out <- subset(out, SEX.RECIPIENT != SEX.SOURCE)
    n.out <- nrow(out)
    cat('Keep only heterosexual pairs\n')
    sprintfcat('Removing %s pairs, resulting in a total of %s pairs.', n - n.out, n.out )

    n <- nrow(out)
    out <- subset(out, COMM.SOURCE != 'neuro' & COMM.RECIPIENT != 'neuro')
    n.out <- nrow(out)
    cat('Keep only RCCS participants\n')
    sprintfcat('Removing %s pairs, resulting in a total of %s pairs.', n - n.out, n.out )

    if(!is.null(only.one.community)){
        n <- nrow(out)
        out <- subset(out, COMM.SOURCE == only.one.community & COMM.RECIPIENT == only.one.community)
        n.out <- nrow(out)
        cat('Keep only pairs in selected community\n')
        sprintfcat('Removing %s pairs, resulting in a total of %s pairs.', n - n.out, n.out )
    }
    if(use_30com_pairs){

        comm_continuously_surveyed <- c(1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 
            89, 94, 106, 107, 108, 120, 391, 602, 754)

        n <- nrow(out)
        out <- subset(out,(COMM_NUM.SOURCE %in% comm_continuously_surveyed & COMM_NUM.RECIPIENT %in% comm_continuously_surveyed))
        n.out <- nrow(out)

        cat('\nExcluding sources and recipients outside of the 30 continuously surveyed communities\n')
        sprintfcat('Removing %s pairs, resulting in a total of %s pairs.', n - n.out, n.out )
    }
    if(only.transmission.after.start.observational.period){
        n <- nrow(out)
        out <- subset(out, !(DATE_INFECTION.RECIPIENT < start_first_period_inland & COMM.RECIPIENT == 'inland'))
        n.out <- nrow(out)
        cat('\nFor inland excluding recipients infected before ', as.character(start_first_period_inland), '\n')
        sprintfcat('Removing %s pairs, resulting in a total of %s pairs.', n - n.out, n.out )
    }
    if(only.transmission.before.stop.observational.period){
        n <- nrow(out)
        out <- subset(out, !(DATE_INFECTION.RECIPIENT > stop_second_period_inland & COMM.RECIPIENT == 'inland'))
        n.out <- nrow(out)
        cat('\nFor inland excluding recipients infected after ', as.character(stop_second_period_inland), '\n')
        sprintfcat('Removing %s pairs, resulting in a total of %s pairs.', n - n.out, n.out )
    }
    if(!is.null(remove.pairs.from.rounds)){

        tmp <- df_round_inland[round %in% remove.pairs.from.rounds, list(
            min_exclusion = min(min_sample_date), 
            max_exclusion = max(max_sample_date))
        ]
        out[, DATE.COLLECTION.PAIR := max(c(DATE.COLLECTION.SOURCE, DATE.COLLECTION.RECIPIENT)), by = c('RECIPIENT', 'SOURCE')]
        n <- nrow(out)
        out <- subset(out, !(COMM.RECIPIENT == 'inland' & DATE.COLLECTION.PAIR <= tmp[, max_exclusion] & DATE.COLLECTION.PAIR >= tmp[,min_exclusion ]))
        n.out <- nrow(out)

        cat('\nExcluding pairs in inland community from round', remove.pairs.from.rounds, '\n')
        cat('\nFor inland excluding recipients infected after ', as.character(stop_second_period_inland), '\n')
        sprintfcat('Removing %s pairs, resulting in a total of %s pairs.', n - n.out, n.out )
    }

    return(out)
}

get.sample.collection.dates <- function(select_aid=NULL, get_first_visit=FALSE)
{
    # condition on whether we are using randomized or not.
    stopifnot("
        file.path.sequence.dates does not exist. 
        Make sure you are defining the variable in the script" = 
        file.exists(file.path.sequence.dates)
    )

    is_randomized <- basename(file.path.sequence.dates) %like% 'randomized'
    msg <- fifelse(is_randomized, 
        yes="Using randomized data",
        no="CARE Using confidential data")

    cat('\n(get.sample.collection.dates:', msg, ')\n')

    ddates <- readRDS(file.path.sequence.dates)

    if(!is.null(select_aid))
        ddates <- ddates[aid %in% select_aid]

    if(get_first_visit)
        ddates <- ddates[, .(date_collection=min(visit_dt)),by='aid']

    return(ddates)
}

get.sample.collection.dates.deprecated <- function(select_aid=NULL, get_first_visit=FALSE)
{
    # get collection dates 
    path.sdates.rccs <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS','200316_pangea_db_sharing_extract_rakai.csv')
    path.sdates.mrc <- file.path(indir.deepsequencedata, 'PANGEA2_MRC','200319_pangea_db_sharing_extract_mrc.csv')

    files <- c(path.sdates.rccs, path.sdates.mrc)
    cols <- c('pt_id', 'pangea_id', 'visit_dt')
    ddates <- rbindlist(lapply(files, fread, select=cols))
    ddates <- unique(ddates)
    ddates <- merge(ddates, aik, by.x='pt_id', by.y='PT_ID')
    ddates[, pt_id := NULL]
    stopifnot(ddates[, uniqueN(pangea_id)==.N])
    setnames(ddates, 'AID', 'aid')

    if(!is.null(select_aid))
        ddates <- ddates[aid %in% select_aid]

    if(get_first_visit)
        ddates <- ddates[, .(date_collection=min(visit_dt)),by='aid']
    ddates
}

set.sensitivity.indicators.from.jobname <- function(jobname)
{
  # map from analysis to jobname
  
  # central = central analysis
  
  # incloess = Using incidence rates estimated with LOESS regression
  # 30com = Using incidence rates estimated on a data subset to 28 continuously surveyed communities
  # tsinonrefined = Using non-refined infection time estimates
 
  # woR18 = Without source-recipients pairs for which the source or recipient was sequenced after round 17
  # woR1718 = Without source-recipients pairs for which the source or recipient was sequenced after round 16
  # woR161718 = Without source-recipients pairs for which the source or recipient was sequenced after round 15
  
  # seed12 = Using a bootstrap sample of the source-recipient pairs (first draw)
  # seed13 = Using a bootstrap sample of the source-recipient pairs (second draw)
  # seed14 = Using a bootstrap sample of the source-recipient pairs (third draw)
  
  # nonparttreatedaspart = Assuming the same suppression rate in non-participants as in participants
  # nonpartnottreated = Assuming that non-participants are not suppressed
  # nonpartfemalemale125moreinfectious = Assuming that prevalence in non-participants is 25% higher than in participants
  # nonpartmale125moreinfectious = Assuming that prevalence in men non-participants is 25% higher than in men participants
  # nonpartfemale125moreinfectious = Assuming that prevalence in women non-participants is 25% higher than in women participants
  
  # vl200 = Defining viral suppression as a viral load measurement below 200 copies/mL plasma blood
  
  
  token <- 0
  
  if(jobname == "incloess"){
    cat("Modifying flag: use_loess_inc_estimates...\n"  )
    use_loess_inc_estimates <<- TRUE
    token <- token + 1
  }
  
  if(jobname == "30com"){
    cat("Modifying flags: use_30com_inc_estimates ...\n"  )
    cat("Modifying flags: use_30com_pairs ...\n"  )
    use_30com_inc_estimates <<- TRUE
    use_30com_pairs <<- TRUE
    token <- token + 1
  }
  
  if(jobname == "tsinonrefined"){
    cat("Modifying flag:  use_tsi_non_refined...\n"  )
    use_tsi_non_refined <<- TRUE
    token <- token + 1
  }
  
  if(jobname %like% 'woR[0-9]+$') {
    cat("Modifying flag: remove.pairs.from.rounds...\n")
    rounds_to_remove <- gsub('woR([0-9]+)$', '\\1', jobname) %>% as.integer()
    while(rounds_to_remove > 1){
      tmp <- paste0('R0', rounds_to_remove %% 100)
      rounds_to_remove <- rounds_to_remove %/% 100
      remove.pairs.from.rounds <<- c(remove.pairs.from.rounds, tmp)
    }
    token <- token + 1
  }
  
  if(jobname %like% 'seed[0-9]+$') {
    cat("Modifying flag: pairs_replicates.seed  ...\n")
    pairs_replicates.seed <<- gsub('seed([0-9]+)$',"\\1",jobname) %>%  as.integer()
    token <- token + 1
  }
  
  if(jobname == "nonparttreatedaspart"){
    cat("Modifying flag: nonparticipants.treated.like.participants  ...\n")
    nonparticipants.treated.like.participants <<- TRUE
    token <- token + 1
  }
  
  if(jobname == "nonpartnottreated"){
    cat("Modifying flag: nonparticipants.not.treated  ...\n")
    nonparticipants.not.treated <<- TRUE
    token <- token + 1
  }
  
  if(jobname == "nonpartfemalemale125moreinfectious"){
    cat("Modifying flag: nonparticipants.male.relative.infection  ...\n")
    nonparticipants.male.relative.infection <<- 1.25
    cat("Modifying flag: nonparticipants.female.relative.infection  ...\n")
    nonparticipants.female.relative.infection <<- 1.25
    token <- token + 1
  }
  
  if(jobname == "nonpartmale125moreinfectious"){
    cat("Modifying flag: nonparticipants.male.relative.infection  ...\n")
    nonparticipants.male.relative.infection <<- 1.25
    token <- token + 1
  }
  
  if(jobname == "nonpartfemale125moreinfectious"){
    cat("Modifying flag: nonparticipants.female.relative.infection  ...\n")
    nonparticipants.female.relative.infection <<- 1.25
    token <- token + 1
  }
  
  if(jobname == "vl200"){
    cat("Modifying flag: viremic_viral_load_200ml  ...\n")
    viremic_viral_load_200ml <<- TRUE
    token <- token + 1
  }
  
  # sanity check
  if(jobname == 'central'){
    stopifnot(token == 0)
  }else{
    stopifnot(token == 1)
  }

  
}
