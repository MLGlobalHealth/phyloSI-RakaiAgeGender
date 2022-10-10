ageanalysis <- function(infile.inference=NULL,infile.prior.samples=NULL,opt=NULL,M=30,D=2,outdir){
  
  require(data.table)	
  #
  #	input args
  #
  if(is.null(opt))
  {
    opt									<- list()
    opt$adjust.sequencing.bias			<- 1
    opt$adjust.participation.bias		<- 1
    opt$migration.def.code				<- '24'
    opt$set.missing.migloc.to.inland	<- 0
    opt$set.missing.migloc.to.fishing	<- 1-opt$set.missing.migloc.to.inland		
  }
  if(is.null(infile.inference))
  {
    infile.inference	<- file.path(outdir, "RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda")
  }
  if(is.null(infile.prior.samples))
  {
    infile.prior.samples <- file.path(outdir,"samples_fit.rda")
  }
  
  cat('\ninfile.inference=',infile.inference)
  cat('\ninfile.prior.samples=',infile.prior.samples)
  
  cat('\nopt=',unlist(opt))			
  indir					<- dirname(infile.inference)	
  outfile.base			<- '~/ageanalysis/'
  load(infile.inference)
  #
  #	prepare data on observed transmission flows
  #
  #	subset to variables needed, using RTR3	
  rtr	<- copy(rtr3)
  if(opt$migration.def.code=='06')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_05YR, REC_INMIGRATE_05YR, TR_COMM_NUM_A_MIG_05YR')
    setnames(rtr, 'TR_INMIGRATE_05YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_05YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_05YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='12')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_1YR, REC_INMIGRATE_1YR, TR_COMM_NUM_A_MIG_1YR')
    setnames(rtr, 'TR_INMIGRATE_1YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_1YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_1YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='24')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_2YR, REC_INMIGRATE_2YR, TR_COMM_NUM_A_MIG_2YR')
    setnames(rtr, 'TR_INMIGRATE_2YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_2YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_2YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='36')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_3YR, REC_INMIGRATE_3YR, TR_COMM_NUM_A_MIG_3YR')
    setnames(rtr, 'TR_INMIGRATE_3YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_3YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_3YR', 'TR_COMM_NUM_A_MIG')
  }
  if(opt$migration.def.code=='48')
  {
    cat('\nSource attribution analysis based on TR_INMIGRATE_4YR, REC_INMIGRATE_4YR, TR_COMM_NUM_A_MIG_4YR')
    setnames(rtr, 'TR_INMIGRATE_4YR', 'TR_INMIGRATE')
    setnames(rtr, 'REC_INMIGRATE_4YR', 'REC_INMIGRATE')
    setnames(rtr, 'TR_COMM_NUM_A_MIG_4YR', 'TR_COMM_NUM_A_MIG')
  }
  
  rtr	<- subset(rtr, select=c('PAIRID','TR_RID','TR_COMM_NUM','TR_COMM_NUM_A','TR_COMM_NUM_A_MIG',
                              'TR_SEX','TR_BIRTHDATE','TR_COMM_TYPE','TR_INMIG_LOC','TR_INMIGRATE',
                              'REC_RID','REC_COMM_NUM','REC_COMM_NUM_A',
                              'REC_SEX','REC_BIRTHDATE','REC_COMM_TYPE','REC_INMIGRATE'))
  
  #	set unknown origin to either fishing or inland
  tmp	<- rtr[, which(TR_INMIGRATE=='inmigrant_from_unknown')]
  if(opt$set.missing.migloc.to.inland)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_inland')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'imig')
  }		
  if(opt$set.missing.migloc.to.fishing)
  {
    set(rtr, tmp, 'TR_INMIGRATE', 'inmigrant_from_fish')
    set(rtr, tmp, 'TR_COMM_NUM_A_MIG', 'fmig')
  }
  
  # add age 
  rtr[,TR_AGE_AT_MID:=2013.25-TR_BIRTHDATE]
  rtr[,REC_AGE_AT_MID:=2013.25-REC_BIRTHDATE]
  
  # impute age
  tmp	<- which(is.na(rtr$TR_AGE_AT_MID))
  set(rtr, tmp, 'TR_AGE_AT_MID', mean(rtr$TR_AGE_AT_MID[which(!is.na(rtr$TR_AGE_AT_MID))]) )
  tmp	<- which(is.na(rtr$REC_AGE_AT_MID))
  set(rtr, tmp, 'REC_AGE_AT_MID', mean(rtr$REC_AGE_AT_MID[which(!is.na(rtr$REC_AGE_AT_MID))]) )
  
  # fixup from latest surveillance data
  set(rtr, rtr[,which(TR_RID=="C036808")], 'TR_AGE_AT_MID', 39.946)	
  set(rtr, rtr[,which(REC_RID=="G036802")], 'REC_AGE_AT_MID',	44.946)	
  set(rtr, rtr[, which(REC_RID=="H103745")], 'REC_AGE_AT_MID', 20.42)	
  set(rtr, rtr[, which(REC_RID=="C121534")],'REC_AGE_AT_MID', 28.549)
  
  #	stratify age
  rtr[, TR_AGE_AT_MID_C:= as.character(cut(TR_AGE_AT_MID, breaks=c(15,16:49,50), labels=paste0(15:49,'-',16:50), right=FALSE))]
  rtr[, REC_AGE_AT_MID_C:= as.character(cut(REC_AGE_AT_MID, breaks=c(15,16:49,50), labels=paste0(15:49,'-',16:50), right=FALSE))]
  # stopifnot( nrow(subset(rtr, is.na(TR_AGE_AT_MID_C)))==0 )
  # stopifnot( nrow(subset(rtr, is.na(REC_AGE_AT_MID_C)))==0 )
  rtr <- subset(rtr, !is.na(TR_AGE_AT_MID_C))
  rtr <- subset(rtr, !is.na(REC_AGE_AT_MID_C))
  
  # define TR_COMM_TYPE_F, REC_COMM_TYPE_F (i: inland; f: fishing) 
  rtr[,TR_COMM_TYPE_F:=substr(TR_COMM_TYPE,1,1)]
  unique(rtr$TR_COMM_TYPE_F)
  rtr[substr(TR_COMM_TYPE,1,1)!='f',TR_COMM_TYPE_F:='i']
  rtr[,REC_COMM_TYPE_F:=substr(REC_COMM_TYPE,1,1)]
  unique(rtr$REC_COMM_TYPE_F)
  rtr[substr(REC_COMM_TYPE,1,1)!='f',REC_COMM_TYPE_F:='i']
  
  # define TR_COMM_TYPE_F_MIG (i: inland; f: fishing; e: external) 
  rtr[,TR_COMM_TYPE_F_MIG:=substr(TR_COMM_NUM_A_MIG,1,1)]
  unique(rtr$TR_COMM_TYPE_F_MIG)
  rtr[substr(TR_COMM_NUM_A_MIG,1,1)=='a' | substr(TR_COMM_NUM_A_MIG,1,1)=='i'|
        substr(TR_COMM_NUM_A_MIG,1,1)=='t',TR_COMM_TYPE_F_MIG:='i']
  
  #	build category to match with sampling data tables 
  rtr[, REC_SAMPLING_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_SAMPLING_CATEGORY:= paste0(TR_COMM_TYPE_F,':',TR_SEX,':',TR_AGE_AT_MID_C)]
  #	build transmission flow category 
  rtr[, REC_TRM_CATEGORY:= paste0(REC_COMM_TYPE_F,':',REC_SEX,':',REC_AGE_AT_MID_C)]
  rtr[, TR_TRM_CATEGORY:= paste0(TR_COMM_TYPE_F_MIG,':',TR_SEX,':',TR_AGE_AT_MID_C)]
  
  table(paste0(rtr$TR_COMM_TYPE_F_MIG,'-',rtr$TR_COMM_TYPE_F))
  # make all combinations of variables
  dac <- expand.grid( COMM_TYPE_F= c('i','f'),
                      SEX=  c('F','M'),
                      AGE_AT_MID_C= paste0(15:49,'-',16:50))
  
  dac <- as.data.table(dac)  				
  dac[, CATEGORY:= paste0(COMM_TYPE_F, ':', SEX, ':', AGE_AT_MID_C)]  					
  dac <- as.data.table(expand.grid(TR_SAMPLING_CATEGORY= dac$CATEGORY, REC_SAMPLING_CATEGORY= dac$CATEGORY))
  # ignore Male-Male and Female-Female combinations
  dac <- subset(dac, !(grepl('F',TR_SAMPLING_CATEGORY)&grepl('F',REC_SAMPLING_CATEGORY)) &
                  !(grepl('M',TR_SAMPLING_CATEGORY)&grepl('M',REC_SAMPLING_CATEGORY))  
  )
  # add transmission categories
  dac[, REC_TRM_CATEGORY:= REC_SAMPLING_CATEGORY]
  dac[, TR_TRM_CATEGORY:= TR_SAMPLING_CATEGORY]  
  
  # add inmigrants from external communities
  tmp <- copy(dac)
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^[f|i]','e',TR_SAMPLING_CATEGORY)]) 
  dac <- rbind(dac, tmp)  
  # add inmigrants sampled in inland and migrated from fishing
  tmp <- dac[grepl('^i',TR_SAMPLING_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^i','f',TR_SAMPLING_CATEGORY)]) 
  dac <- rbind(dac, tmp)  
  # add inmigrants sampled in fishing and migrated from inland
  tmp <- dac[grepl('^f',TR_SAMPLING_CATEGORY)]
  set(tmp, NULL, 'TR_TRM_CATEGORY', tmp[,gsub('^f','i',TR_SAMPLING_CATEGORY)]) 
  dac <- rbind(dac, tmp)   
  # remove duplicated rows 
  # TR_SAMPLING_CATEGORY f: F: 15-24: 1 and i: F: 15-24: 1 are all set to e: F: 15-24: 1
  dac <- unique(dac)
  # setnames(dac, c('TR_CATEGORY', 'REC_CATEGORY'), c('TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'))
  # 
  
  #	calculate observed number of transmissions
  #
  dobs	<- rtr[, list( TRM_OBS=length(unique(PAIRID))), by=c('TR_TRM_CATEGORY','REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY','REC_SAMPLING_CATEGORY')]
  dac[, DUMMY:= 1]
  dobs <- merge(dac, dobs, by=c('TR_TRM_CATEGORY', 'REC_TRM_CATEGORY','TR_SAMPLING_CATEGORY', 'REC_SAMPLING_CATEGORY'), all=TRUE)
  stopifnot( dobs[, !any(is.na(DUMMY))] )
  set(dobs, NULL, 'DUMMY', NULL)
  set(dobs, dobs[, which(is.na(TRM_OBS))], 'TRM_OBS', 0L)
  
  #	make PAIR_ID
  dobs[,TR_SMOOTH_CATEGORY:= as.numeric(substr(TR_TRM_CATEGORY,5,6))+0.5]
  dobs[,REC_SMOOTH_CATEGORY:= as.numeric(substr(REC_TRM_CATEGORY,5,6))+0.5]
  dobs[,OUTPUT:=paste0(substr(TR_TRM_CATEGORY,1,1),':',
                       substr(TR_TRM_CATEGORY,3,3),':',
                       substr(TR_SAMPLING_CATEGORY,1,1),':',
                       substr(REC_TRM_CATEGORY,1,1),':',
                       substr(REC_TRM_CATEGORY,3,3))]
  
  # combinations
  tmp <- subset(dobs,select = 'OUTPUT')
  tmp <- unique(tmp)
  setkey(tmp, OUTPUT)
  tmp[, OUTPUT_ID:= seq_len(nrow(tmp))]
  dobs <- merge(dobs, tmp, by='OUTPUT')
  setkey(dobs,OUTPUT_ID,TR_SMOOTH_CATEGORY,REC_SMOOTH_CATEGORY)	
  dobs[, TRM_CAT_PAIR_ID:= seq_len(nrow(dobs))]
  setkey(dobs, TRM_CAT_PAIR_ID)
  
  # sampling prior
  load(infile.prior.samples)
  set(dprior.fit, NULL, 'SAMPLING_CATEGORY', dprior.fit[, gsub('^2','e',gsub('^1','f',gsub('^0','i',SAMPLING_CATEGORY)))] )		
  
  dprior.id <- subset(dprior.fit,select = c('SAMPLING_CATEGORY','WHO'))
  setkey(dprior.id, SAMPLING_CATEGORY,WHO)
  dprior.id[,ID:= seq_len(nrow(dprior.id))]
  
  tmp <- subset(dobs, select = c('TR_SAMPLING_CATEGORY',
                                 'REC_SAMPLING_CATEGORY',
                                 'TRM_CAT_PAIR_ID'))
  
  setnames(dprior.id,colnames(dprior.id),paste0('TR_',colnames(dprior.id)))
  tmp <- merge(tmp,subset(dprior.id[TR_WHO=='TR_SAMPLING_CATEGORY',],select=c('TR_ID','TR_SAMPLING_CATEGORY')),by='TR_SAMPLING_CATEGORY',all.x = TRUE)
  setnames(dprior.id,colnames(dprior.id),gsub('TR_','REC_',colnames(dprior.id)))
  tmp <- merge(tmp,subset(dprior.id[REC_WHO=='REC_SAMPLING_CATEGORY',],select=c('REC_ID','REC_SAMPLING_CATEGORY')),by='REC_SAMPLING_CATEGORY',all.x = TRUE)
  setkey(tmp, TRM_CAT_PAIR_ID)
  xi_id <- cbind(tmp$TR_ID, tmp$REC_ID)
  
  setkey(dprior.fit,SAMPLING_CATEGORY,WHO)
  
  
  indices <- matrix(NA, M^D, D)
  mm=0;
  for (m1 in 1:M){
    for (m2 in 1:M){
      mm = mm+1
      indices[mm,] = c(m1, m2)
    }
  }
  
  B1 <- max(dobs$TR_SMOOTH_CATEGORY)
  B2 <- max(dobs$REC_SMOOTH_CATEGORY)
  L <-  matrix(rep(c(B1, B2) * 5/4,each=nrow(indices)),nrow=nrow(indices))
  sevalue <- pi * indices / (2 * L)
  ns <- nrow(dobs[OUTPUT_ID==1])
  efunc <-  do.call( rbind, lapply(1:ns, 
                                   function(k){
                                     as.vector(apply(sqrt(1/L) * sin(sevalue * (matrix(rep(c(dobs$TR_SMOOTH_CATEGORY[k],dobs$REC_SMOOTH_CATEGORY[k]),each=nrow(indices)),nrow=nrow(indices)) + L)),1, prod))} ))
  
  data.fit <- list(N=nrow(dobs),
                   N_group = max(dobs$OUTPUT_ID),
                   N_per_group = ns,
                   y=matrix(dobs$TRM_OBS,nrow=ns,ncol=max(dobs$OUTPUT_ID)),
                   D=D,
                   x=cbind(dobs$TR_SMOOTH_CATEGORY[1:ns],dobs$REC_SMOOTH_CATEGORY[1:ns]),
                   M_nD=M^2,
                   Xgp=efunc,
                   slambda=sevalue,
                   N_xi = nrow(dprior.fit),
                   shape = cbind(dprior.fit$SHAPE1,dprior.fit$SHAPE2),
                   xi_id_src = matrix(xi_id[,1],nrow=ns,ncol=max(dobs$OUTPUT_ID)),
                   xi_id_rec = matrix(xi_id[,2],nrow=ns,ncol=max(dobs$OUTPUT_ID)),
                   id_mf =   dobs[grepl(':M:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_fm =   dobs[grepl(':F:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_eh = dobs[grepl('e:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_el = dobs[grepl('e:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_hh = dobs[grepl('f:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_hl = dobs[grepl('f:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_lh = dobs[grepl('i:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_mf_h =  dobs[grepl(':M:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_mf_l =  dobs[grepl(':M:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_fm_h =  dobs[grepl(':F:',TR_TRM_CATEGORY) & grepl('f:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_fm_l =  dobs[grepl(':F:',TR_TRM_CATEGORY) & grepl('i:',REC_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_sh =  dobs[grepl('f:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_sl =  dobs[grepl('i:',TR_TRM_CATEGORY),unique(OUTPUT_ID)],
                   id_se =  dobs[grepl('e:',TR_TRM_CATEGORY),unique(OUTPUT_ID)]
  )
  
  save(data.fit,file=file.path(outdir, 'input.rda'))
  
  return(data.fit)
}

add_stan_data_base <- function(){
  
  stan_data = list()
  
  # number of groups
  stan_data[['N_DIRECTION']] = nrow(df_direction)
  stan_data[['N_COMMUNITY']] = nrow(df_community)
  stan_data[['N_PERIOD']] = max(df_period$INDEX_TIME)
  
  # number of age 
  stan_data[['N_PER_GROUP']] = nrow(df_age)
  
  # number of age group
  stan_data[['N_AGE']] = df_age[, length(unique(AGE_INFECTION.RECIPIENT))]
  
  # number of rounds 
  tmp <- merge(df_round, df_community, by = 'COMM')[order(INDEX_COMMUNITY)]
  stan_data[['N_ROUND']] = tmp[, length(unique(INDEX_ROUND)), by = 'COMM']$V1
  stan_data[['N_ROUND_INLAND']] = df_round[COMM == 'inland', length(unique(INDEX_ROUND))]
  stan_data[['N_ROUND_FISHING']] = df_round[COMM == 'fishing', length(unique(INDEX_ROUND))]
  
  # map from round to period
  stan_data[['map_round_period']] =   rbind(c(tmp[COMM ==  df_community[order(INDEX_COMMUNITY), COMM[1]] & order(round), INDEX_TIME], rep(-1,  max(stan_data[['N_ROUND']]) - stan_data[['N_ROUND']][1] )), 
                                            c(tmp[COMM ==  df_community[order(INDEX_COMMUNITY), COMM[2]] & order(round), INDEX_TIME], rep(-1,  max(stan_data[['N_ROUND']]) - stan_data[['N_ROUND']][2] )))
  
  # number of period
  stan_data[['N_ROUND_PER_PERIOD']] = rbind(tmp[order(INDEX_TIME) & COMM == df_community[order(INDEX_COMMUNITY), COMM[1]], length(unique(INDEX_ROUND)), by = 'INDEX_TIME']$V1, 
                                            tmp[order(INDEX_TIME) & COMM == df_community[order(INDEX_COMMUNITY), COMM[2]], length(unique(INDEX_ROUND)), by = 'INDEX_TIME']$V1)
  
  return(stan_data)
}

add_phylo_data <- function(stan_data, pairs){

  # prepare pairs
  pairs_round <- pairs[, list(AGE_TRANSMISSION.SOURCE = floor(AGE_TRANSMISSION.SOURCE), 
                    AGE_INFECTION.RECIPIENT = floor(AGE_INFECTION.RECIPIENT), 
                    DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT = DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT, 
                    COMM.RECIPIENT = COMM.RECIPIENT, 
                    SEX.SOURCE = SEX.SOURCE)]
  
  # save count in each entry
  y = array(NA, c(stan_data[['N_PER_GROUP']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']]))
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:stan_data[['N_PERIOD']]){
        
        # direction group
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        tmp <- pairs_round[SEX.SOURCE == .SEX.SOURCE]
        
        # community group
        tmp <- tmp[COMM.RECIPIENT == df_community[j, COMM]]
        
        # time group
        tmp <- tmp[DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == df_period[COMM == df_community[j, COMM] & INDEX_TIME == k, BEFORE_CUTOFF]]
        
        # count number of observation
        tmp <- tmp[, list(count = .N), by = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT')]
        tmp <- merge(df_age, tmp, 
                     by = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), all.x = T)
        tmp[is.na(count), count := 0]
        
        setkey(tmp, AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)
        
        tmp1 <- pairs[SEX.SOURCE == .SEX.SOURCE  & COMM.RECIPIENT == df_community[j, COMM] & DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT == df_period[COMM == df_community[j, COMM] & INDEX_TIME == k, BEFORE_CUTOFF]]
        stopifnot(sum(tmp$count) == nrow(tmp1))
        
        # check the order of ages is correct
        tmp <- tmp[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == tmp[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == tmp[, AGE_TRANSMISSION.SOURCE])
        
        cat(nrow(tmp1), 'pairs with infection', df_direction[i, LABEL_DIRECTION], 'towards', df_community[j, COMM], 'in', df_period[COMM == df_community[j, COMM] & INDEX_TIME == k, PERIOD], '\n')
        
        y[, i, j, k] = matrix(tmp$count, ncol = 1)
        
        
      }
    }
  }
    
    
  # save stan data
  stan_data[['y']] = y
  
  # add map age source and recipient
  df_age <- df_age[order(INDEX_AGE)]
  stan_data[['map_age_source']] = df_age[, AGE_TRANSMISSION.SOURCE - min(AGE_TRANSMISSION.SOURCE) + 1]
  stan_data[['map_age_recipient']] = df_age[, AGE_INFECTION.RECIPIENT - min(AGE_INFECTION.RECIPIENT) + 1]
  
  return(stan_data)
}

add_incidence_cases <- function(stan_data, incidence_cases_round){
  
  # save count in each entry
  z = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          z[, i, j, k] = -1;
          next
        }
        
        .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM = df_community[INDEX_COMMUNITY == j, COMM]
        .ROUND = df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # direction group
        tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
        
        # community group
        tmp <- tmp[COMM == .COMM]
        
        # round
        tmp <- tmp[ROUND == .ROUND]
        
        # order by age
        tmp <- tmp[order(AGEYRS)] 
        
        # sanity check
        tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT  & COMM == .COMM & ROUND == .ROUND]
        stopifnot(sum(tmp$INCIDENT_CASES) == sum(tmp1$INCIDENT_CASES))
        cat(sum(tmp1$INCIDENT_CASES), 'incidence cases ', df_direction[i, LABEL_DIRECTION], 'towards', .COMM, 'in', .ROUND, '\n')
        
        # fill
        z[, i, j, k] = round(tmp$INCIDENT_CASES)
        
      }
    }
  }
  
  stan_data[['z']] = z

  
  return(stan_data)
  
}

add_incidence_rates <- function(stan_data, incidence_cases_round){
  
  # save count in each entry
  ir = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          ir[, i, j, k] = 0;
          next
        }
        
        .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM = df_community[INDEX_COMMUNITY == j, COMM]
        .ROUND = df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # direction group
        tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
        
        # community group
        tmp <- tmp[COMM == .COMM]
        
        # round
        tmp <- tmp[ROUND == .ROUND]
        
        # order by age
        tmp <- tmp[order(AGEYRS)] 
        
        # sanity check
        tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT  & COMM == .COMM & ROUND == .ROUND]
        stopifnot(sum(tmp$INCIDENCE) == sum(tmp1$INCIDENCE))
        cat(mean(tmp1$INCIDENCE), 'incidence rates ', df_direction[i, LABEL_DIRECTION], 'towards', .COMM, 'in', .ROUND, '\n')
        
        # fill
        ir[, i, j, k] = tmp$INCIDENCE
        
      }
    }
  }
  
  stan_data[['ir']] = ir
  
  
  return(stan_data)
  
}

add_incidence_rates_lognormal_parameters <- function(stan_data, incidence_cases_round){
  
  # save count in each entry
  ir_lognorm_mean = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  ir_lognorm_sd = array(NA, c(stan_data[['N_AGE']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          ir[, i, j, k] = 0;
          next
        }
        
        .SEX.RECIPIENT = substr(gsub('.* -> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM = df_community[INDEX_COMMUNITY == j, COMM]
        .ROUND = df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # direction group
        tmp <- incidence_cases_round[SEX == .SEX.RECIPIENT]
        
        # community group
        tmp <- tmp[COMM == .COMM]
        
        # round
        tmp <- tmp[ROUND == .ROUND]
        
        # order by age
        tmp <- tmp[order(AGEYRS)] 
        
        # sanity check
        tmp1 <- incidence_cases_round[SEX == .SEX.RECIPIENT  & COMM == .COMM & ROUND == .ROUND]
        stopifnot(sum(tmp$INCIDENCE) == sum(tmp1$INCIDENCE))
        cat(mean(tmp1$INCIDENCE), 'incidence rates ', df_direction[i, LABEL_DIRECTION], 'towards', .COMM, 'in', .ROUND, '\n')
        
        # fill
        lognorm_parms <- lognorm::getParmsLognormForMedianAndUpper(median = tmp$INCIDENCE, upper = tmp$UB, sigmaFac=2)
        ir_lognorm_mean[, i, j, k] = lognorm_parms[, 1]
        ir_lognorm_sd[, i, j, k] = lognorm_parms[, 2]
      }
    }
  }
  
  stan_data[['ir_lognorm_mean']] = ir_lognorm_mean
  stan_data[['ir_lognorm_sd']] = ir_lognorm_sd
  
  return(stan_data)
  
}

add_2D_splines_stan_data = function(stan_data, spline_degree = 3, n_knots_rows = 8, n_knots_columns = 8, X, Y)
{
  
  stan_data$number_rows <- length(X)
  stan_data$number_columns <- length(Y)
  
  knots_rows = X[seq(1, length(X), length.out = n_knots_rows)] 
  knots_columns = Y[seq(1, length(Y), length.out = n_knots_columns)]
  
  stan_data$num_basis_rows = length(knots_rows) + spline_degree - 1
  stan_data$num_basis_columns = length(knots_columns) + spline_degree - 1
  
  stan_data$IDX_BASIS_ROWS = 1:stan_data$num_basis_rows
  stan_data$IDX_BASIS_COLUMNS = 1:stan_data$num_basis_columns
  
  stan_data$BASIS_ROWS = bsplines(X, knots_rows, spline_degree)
  stan_data$BASIS_COLUMNS = bsplines(Y, knots_columns, spline_degree)
  
  stopifnot(all( apply(stan_data$BASIS_ROWS, 1, sum) > 0  ))
  stopifnot(all( apply(stan_data$BASIS_COLUMNS, 1, sum) > 0  ))
  
  return(stan_data)
}

add_3D_splines_stan_data = function(stan_data, X, Y, Z,
                                    spline_degree = 3, 
                                    n_knots_1D = NULL, n_knots_2D = NULL, n_knots_3D = NULL, 
                                    knots_1D = NULL, knots_2D = NULL, knots_3D = NULL)
{
  stopifnot(length(X) == length(Y))
  stopifnot(stan_data$N_per_group == length(X) * length(Y) * length(Z))
  stopifnot(!is.null(n_knots_1D) | !is.null(knots_1D))
  stopifnot(!is.null(n_knots_2D) | !is.null(knots_2D))
  stopifnot(!is.null(n_knots_3D) | !is.null(knots_3D))
  
  stan_data$A <- length(X)
  stan_data$T <- length(Z)
  
  if(!is.null(n_knots_1D)){
    knots_1D = X[seq(1, length(X), length.out = n_knots_1D)] 
  } else{
    knots_1D = c(min(X), knots_1D, max(X))
  }
  
  if(!is.null(n_knots_2D)){
    knots_2D = Y[seq(1, length(Y), length.out = n_knots_2D)]
  } else{
    knots_2D = c(min(Y), knots_2D, max(Y))
  }
  
  if(!is.null(n_knots_3D)){
    knots_3D = Z[seq(1, length(Z), length.out = n_knots_3D)]
  } else{
    knots_3D = c(min(Z), knots_3D, max(Z))
  }
  
  stan_data$num_basis_1D = length(knots_1D) + spline_degree - 1
  stan_data$num_basis_2D = length(knots_2D) + spline_degree - 1
  stan_data$num_basis_3D = length(knots_3D) + spline_degree - 1
  
  stan_data$IDX_BASIS_1D = 1:stan_data$num_basis_1D
  stan_data$IDX_BASIS_2D = 1:stan_data$num_basis_2D
  stan_data$IDX_BASIS_3D = 1:stan_data$num_basis_3D
  
  stan_data$BASIS_1D = bsplines(X, knots_1D, spline_degree)
  stan_data$BASIS_2D = bsplines(Y, knots_2D, spline_degree)
  stan_data$BASIS_3D = bsplines(Z, knots_3D, spline_degree)
  
  stopifnot(all( apply(stan_data$BASIS_1D, 1, sum) > 0  ))
  stopifnot(all( apply(stan_data$BASIS_2D, 1, sum) > 0  ))
  stopifnot(all( apply(stan_data$BASIS_3D, 1, sum) > 0  ))

  
  return(stan_data)
}

bspline = function(x, k, order, intervals)
{
  
  if(order == 1){
    return(x >= intervals[k] & x < intervals[k+1])
  }
  
  w1 = 0; w2 = 0
  
  if(intervals[k] != intervals[k+order-1])
    w1 = (x - intervals[k]) / (intervals[k+order-1] - intervals[k])
  if(intervals[k+1] != intervals[k+order])
    w2 = 1 - (x - intervals[k+1]) / (intervals[k+order] - intervals[k+1])
  
  spline = w1 * bspline(x, k, order - 1, intervals) +
    w2 * bspline(x, k+1, order - 1, intervals)
  
  return(spline)
}

find_intervals = function(knots, degree, repeating = T)
{
  
  K = length(knots)
  
  intervals = vector(mode = 'double', length = 2*degree + K)
  
  # support of knots
  intervals[(degree+1):(degree+K)] = knots
  
  # extreme
  if(repeating)
  {
    intervals[1:degree] = min(knots)
    intervals[(degree+K+1):(2*degree+K)] = max(knots)
  } else {
    gamma = 0.1
    intervals[1:degree] = min(knots) - gamma*degree:1
    intervals[(degree+K+1):(2*degree+K)] = max(knots) + gamma*1:degree
  }
  
  return(intervals)
}

bsplines = function(data, knots, degree)
{
  K = length(knots)
  num_basis = K + degree - 1
  
  intervals = find_intervals(knots, degree)
  
  m = matrix(nrow = num_basis, ncol = length(data), 0)
  
  for(k in 1:num_basis)
  {
    m[k,] = bspline(data, k, degree + 1, intervals) 
  }
  
  m[num_basis,length(data)] = 1
  
  return(m)
}

add_informative_prior_gp_mean <- function(stan_data, df_age, file.partnership.rate, outfile.figures){
  
  # load estimated partnership.rate in misc/
  partnership.rate <- as.data.table(read.csv(file.partnership.rate))
  
  # keep heterosexuals partnerships
  partnership.rate <- partnership.rate[part.sex!=cont.sex]
  
  # format
  partnership.rate <- partnership.rate[, .(part.sex, part.age, cont.age, c)]
  setnames(partnership.rate, 'c', 'rate')
  
  tmp <- df_age[, .(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
  tmp <- as.data.table(full_join(tmp, data.table(count_per_group = apply(stan_data$y, 2, sum), 
                                                 index_group = 1:stan_data$N_group), 
                                 by = character()))
  tmp <- merge(tmp, df_group[, .(is_mf, index_group)], by = 'index_group')
  tmp[, part.sex := ifelse(is_mf == T, 'M', 'F')]
  
  tmp <- merge(tmp, partnership.rate, by.x = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT', 'part.sex'), by.y = c('part.age', 'cont.age', 'part.sex'), all.x = T)

  tmp[, total_rate := sum(rate), by = 'index_group']
  tmp[, relative_rate := rate / total_rate]
  
  tmp[, Elambda := relative_rate * count_per_group]
  tmp[, mu := log(Elambda)]
  
  stopifnot(nrow(tmp[is.na(mu)]) == 0)
  
  tmp <- tmp[order(index_group)]
  theta <- vector(mode = 'list', length = stan_data$N_group)
  for(i in 1:stan_data$N_group){
    tmp1 <- tmp[index_group == i]
    sex <- tmp1[, unique(part.sex)]
    
    mu <- as.matrix(dcast(tmp1, AGE_TRANSMISSION.SOURCE ~ AGE_INFECTION.RECIPIENT, value.var = 'mu')[,-1])
    theta[[i]] <- find_spectral_projection_gp_mean(mu,  paste0(outfile.figures, '_sex_', sex))
  }
  
  stan_data[['theta']] <- theta
  
  return(stan_data)
}

add_diagonal_prior_gp_mean <- function(stan_data, df_age, outfile.figures){
  
  tmp <- df_age[, .(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
  tmp <- as.data.table(full_join(tmp, data.table(count_per_group = apply(stan_data$y, 2, sum), index_group = 1:stan_data$N_group), 
                                 by = character()))
  
  # tmp[, rate := (max(AGE_INFECTION.RECIPIENT) - min(AGE_INFECTION.RECIPIENT))*4 +
  #       max(AGE_INFECTION.RECIPIENT) - abs( AGE_INFECTION.RECIPIENT - AGE_TRANSMISSION.SOURCE)*4 ]
  # tmp[, rate := 1 / log(2 + abs( AGE_INFECTION.RECIPIENT - AGE_TRANSMISSION.SOURCE)) ]
  tmp[, rate := max(AGE_INFECTION.RECIPIENT) - abs( AGE_INFECTION.RECIPIENT - AGE_TRANSMISSION.SOURCE) ]
  
  # tmp[AGE_INFECTION.RECIPIENT == AGE_TRANSMISSION.SOURCE, rate := 2]
  
  tmp[, total_rate := sum(rate), by = 'index_group']
  tmp[, relative_rate := rate / total_rate]
  
  tmp[, Elambda := relative_rate * count_per_group]
  tmp[, mu := log(Elambda)]

  tmp <- tmp[order(index_group)]
  theta <- vector(mode = 'list', length = stan_data$N_group)
  for(i in 1:stan_data$N_group){
    tmp1 <- tmp[index_group == i]
    
    mu <- as.matrix(dcast(tmp1, AGE_TRANSMISSION.SOURCE ~ AGE_INFECTION.RECIPIENT, value.var = 'mu')[,-1])
    theta[[i]] <- find_spectral_projection_gp_mean(mu, outfile.figures)
  }

  stan_data[['theta']] <- theta
  
  return(stan_data)
}

add_flat_prior_gp_mean <- function(stan_data, df_age, outfile.figures){
  
  tmp <- df_age[, .(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
  tmp <- as.data.table(full_join(tmp, data.table(count_per_group = apply(stan_data$y, 2, sum), index_group = 1:stan_data$N_group), 
                                 by = character()))
  
  tmp[, rate := 1 ]
  tmp[, total_rate := sum(rate), by = 'index_group']
  tmp[, relative_rate := rate / total_rate]
  
  tmp[, Elambda := relative_rate * count_per_group]
  tmp[, mu := log(Elambda)]
  
  tmp <- tmp[order(index_group)]
  theta <- vector(mode = 'list', length = stan_data$N_group)
  for(i in 1:stan_data$N_group){
    tmp1 <- tmp[index_group == i]
    
    mu <- as.matrix(dcast(tmp1, AGE_TRANSMISSION.SOURCE ~ AGE_INFECTION.RECIPIENT, value.var = 'mu')[,-1])
    theta[[i]] <- find_spectral_projection_gp_mean(mu, outfile.figures)
  }
  
  stan_data[['theta']] <- theta
  
  # theta <- matrix(nrow = stan_data[['num_basis_rows']], ncol = stan_data[['num_basis_columns']], 0)
  # stan_data[['theta']] <- rep(list(theta), stan_data[['N_group']])
  
  return(stan_data)
}

find_spectral_projection_gp_mean <- function(mu, outdir = NULL){
  
  A = t(stan_data[['BASIS_ROWS']])
  B = stan_data[['BASIS_COLUMNS']]
  theta <- MASS::ginv(A) %*% mu %*% MASS::ginv(B)

  # plot
  if(!is.null(outdir)){
    
    mu_pred <- A %*% theta %*% B
    tmp1 = as.data.table(reshape2::melt(mu_pred))
    setkey(tmp1, Var1, Var2)
    
    tmp <- copy(df_age)
    tmp[, transmission_rate := exp(tmp1$value)]
    tmp[, total_transmission_rate := sum(transmission_rate)]
    tmp[, transmission_flow := transmission_rate / total_transmission_rate]
    
    tmp1 <- tmp[, list(total_flow = sum(transmission_flow)), by = 'AGE_INFECTION.RECIPIENT']
    tmp1 <- merge(tmp, tmp1, by = 'AGE_INFECTION.RECIPIENT')
    tmp1[, delta := transmission_flow / total_flow]
    tmp1 <- tmp1[, list(Median = matrixStats::weightedMedian(AGE_TRANSMISSION.SOURCE, delta), 
                        FirstQuartile = as.numeric(modi::weighted.quantile(AGE_TRANSMISSION.SOURCE, delta, 0.25)), 
                        ThirdQuartile = as.numeric(modi::weighted.quantile(AGE_TRANSMISSION.SOURCE, delta, 0.75))), by = 'AGE_INFECTION.RECIPIENT']
    tmp1[, IQR := ThirdQuartile - FirstQuartile]
    
    ggplot(tmp, aes(x =AGE_INFECTION.RECIPIENT)) + 
      geom_raster(aes(fill = transmission_flow, y = AGE_TRANSMISSION.SOURCE)) + 
      geom_boxplot(data = tmp1, col = 'brown3', alpha= 0, width = 0.5,
        stat = "identity",
        aes(lower  = FirstQuartile,
            upper  = ThirdQuartile,
            middle = Median,
            ymin   = FirstQuartile - 1.5 * IQR, # optional
            ymax   = ThirdQuartile + 1.5 * IQR, 
            group = AGE_INFECTION.RECIPIENT) # optional
      ) +
      labs(x = 'age infection recipient', y = 'age transmission source', fill = 'prior median\ntransmission flow') + 
      scale_fill_viridis_c() + 
      scale_y_continuous(expand = c(0,0)) + 
      scale_x_continuous(expand = c(0,0)) + 
      coord_cartesian(xlim = range_age_non_extended, ylim = range_age_non_extended) +
      theme(legend.position = 'bottom') +
      geom_abline(intercept = 0, slope = 1, linetype= 'dashed', col = 'white') 
    ggsave(paste0(outdir, '-Prior_transmissionFlow.png'), w = 5, h = 5)
    
    
    #
    tmp1 <- tmp[, list(total_flow = sum(transmission_flow)), by = 'AGE_TRANSMISSION.SOURCE']
    tmp1 <- merge(tmp, tmp1, by = 'AGE_TRANSMISSION.SOURCE')
    tmp1[, delta := transmission_flow / total_flow]
    
    tmp1 <- tmp1[, list(M = matrixStats::weightedMedian(AGE_INFECTION.RECIPIENT, delta), 
                        CL = as.numeric(modi::weighted.quantile(AGE_INFECTION.RECIPIENT, delta, 0.1)), 
                        CU = as.numeric(modi::weighted.quantile(AGE_INFECTION.RECIPIENT, delta, 0.9))), by = 'AGE_TRANSMISSION.SOURCE']
    
    ggplot(tmp1, aes(x = AGE_TRANSMISSION.SOURCE)) + 
      geom_abline(intercept = 0, slope = 1, linetype= 'dashed', col = 'darkred') + 
      geom_line(aes(y = M)) + 
      geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
      labs(x = 'Age transmission source', y = 'Age infection recipient (median and 80% IQR)', color = 'prior median') + 
      theme(legend.position = 'bottom') +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), limits = range_age_non_extended) + 
      scale_x_continuous(expand = c(0,0), limits = range_age_non_extended) 
    ggsave(paste0(outdir, '-Prior_AgeInfection.png'), w = 5, h = 5)
    
    #
    tmp1 <- tmp[, list(total_flow = sum(transmission_flow)), by = 'AGE_INFECTION.RECIPIENT']
    tmp1 <- merge(tmp, tmp1, by = 'AGE_INFECTION.RECIPIENT')
    tmp1[, delta := transmission_flow / total_flow]
    
    tmp1 <- tmp1[, list(M = matrixStats::weightedMedian(AGE_TRANSMISSION.SOURCE, delta), 
                        CL = as.numeric(modi::weighted.quantile(AGE_TRANSMISSION.SOURCE, delta, 0.1)), 
                        CU = as.numeric(modi::weighted.quantile(AGE_TRANSMISSION.SOURCE, delta, 0.9))), by = 'AGE_INFECTION.RECIPIENT']
    
    ggplot(tmp1, aes(x = AGE_INFECTION.RECIPIENT)) + 
      geom_abline(intercept = 0, slope = 1, linetype= 'dashed', col = 'darkred') + 
      geom_line(aes(y = M)) + 
      geom_errorbar(aes(ymin = CL, ymax = CU), alpha = 0.5) + 
      labs(y = 'Age transmission source (median and 80% IQR)', x = 'Age infection recipient', color = 'prior median') + 
      theme(legend.position = 'bottom') +
      theme_bw() +
      scale_y_continuous(expand = c(0,0), limits = range_age_non_extended) + 
      scale_x_continuous(expand = c(0,0), limits = range_age_non_extended) 
    ggsave(paste0(outdir, '-Prior_Agetransmission.png'), w = 5, h = 5)
  }
  
  return(theta)
}

add_offset <- function(stan_data, eligible_count){
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  eligible_count_wide[, PROP_SUSCEPTIBLE := SUSCEPTIBLE / ELIGIBLE]
  
  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          log_offset_array[i, j, k,] = -1
          next
        }
        
        log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # add proportion of susceptible in recipient
        tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, PROP_SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
        
        # add number of infected unsuppressed in source
        tmp <- eligible_count_wide[SEX == .SEX.SOURCE & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, INFECTED_NON_SUPPRESSED)], by.x = 'AGE_TRANSMISSION.SOURCE', by.y = 'AGEYRS')
        
        # make log offset
        log_offset[, LOG_OFFSET := log(PROP_SUSCEPTIBLE) + log(INFECTED_NON_SUPPRESSED)]
        
        # check the order of ages is correct
        log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
        
        # add to array
        log_offset_array[i, j, k,] = log_offset[, LOG_OFFSET]
        
      }
    }
  }
  
  stan_data[['log_offset']] = log_offset_array

  return(stan_data)
}

add_offset_time <- function(stan_data, eligible_count){
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]
  
  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          log_offset_array[i, j, k,] = -1
          next
        }
        
        log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # add period in year
        log_offset[, PERIOD_SPAN := df_round[ROUND == .ROUND & COMM == .COMM, ROUND_SPANYRS]]
        
        # make log offset time
        log_offset[, LOG_OFFSET_TIME :=  log(PERIOD_SPAN)]
        
        # check the order of ages is correct
        log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
        
        # add to array
        log_offset_array[i, j, k,] = log_offset[, LOG_OFFSET_TIME]
        
      }
    }
  }
  
  stan_data[['log_offset_time']] = log_offset_array
  
  return(stan_data)
}

add_offset_susceptible <- function(stan_data, eligible_count){
  
  eligible_count_wide <- eligible_count_round[order(SEX, COMM, ROUND, AGEYRS)]

  log_offset_array = array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], max(stan_data[['N_ROUND']]), stan_data[['N_PER_GROUP']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:max(stan_data[['N_ROUND']])){
        
        if(k > stan_data[['N_ROUND']][j]){
          log_offset_array[i, j, k,] = -1
          next
        }
        
        log_offset = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .ROUND <- df_round[INDEX_ROUND == k & COMM == .COMM, ROUND]
        
        # add number of susceptible in recipient
        tmp <- eligible_count_wide[SEX == .SEX.RECIPIENT & COMM == .COMM & ROUND == .ROUND]
        log_offset <- merge(log_offset, tmp[, .(AGEYRS, SUSCEPTIBLE)], by.x = 'AGE_INFECTION.RECIPIENT', by.y = 'AGEYRS')
        
        # make log offset
        log_offset[, LOG_OFFSET_SUSCEPTIBLE := log(SUSCEPTIBLE) ]
        
        # check the order of ages is correct
        log_offset <- log_offset[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_offset[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_offset[, AGE_TRANSMISSION.SOURCE])
        
        # add to array
        log_offset_array[i, j, k,] = log_offset[, LOG_OFFSET_SUSCEPTIBLE]
        
      }
    }
  }
  
  stan_data[['log_offset_susceptible']] = log_offset_array
  
  return(stan_data)
}

add_probability_sampling <- function(stan_data, proportion_sampling){
  
  proportion_sampling <- proportion_sampling[order(SEX.RECIPIENT, COMM, BEFORE_CUTOFF, PERIOD)]
  
  log_prop_sampling_array =array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']], stan_data[['N_PER_GROUP']])))
  sampling_index=array(NA, c(c(stan_data[['N_PER_GROUP']], stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']])))
  n_sampling_index=array(NA, c(c(stan_data[['N_DIRECTION']], stan_data[['N_COMMUNITY']], stan_data[['N_PERIOD']])))
  
  for(i in 1:stan_data[['N_DIRECTION']]){
    for(j in 1:stan_data[['N_COMMUNITY']]){
      for(k in 1:stan_data[['N_PERIOD']]){
        
        log_prop_sampling = df_age[, .(AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE)]
        
        .SEX.SOURCE = substr(df_direction[i, LABEL_DIRECTION], 1, 1) 
        .SEX.RECIPIENT = substr(gsub('.*-> (.+)', '\\1', df_direction[i, LABEL_DIRECTION]), 1, 1) 
        .COMM <- df_community[j, COMM]
        .BEFORE_CUTOFF <- df_period[k, BEFORE_CUTOFF]
        
        # add probability of sampling recipient 
        tmp <- proportion_sampling[SEX.RECIPIENT == .SEX.RECIPIENT & COMM == .COMM & BEFORE_CUTOFF == .BEFORE_CUTOFF]
        log_prop_sampling <- merge(log_prop_sampling, tmp[, .(AGEYRS.SOURCE, AGEYRS.RECIPIENT, prop_sampling)], by.x = c('AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), 
                            by.y = c('AGEYRS.SOURCE', 'AGEYRS.RECIPIENT'))
        
        # make log prob sampling
        if(1){
          log_prop_sampling[prop_sampling == 0, prop_sampling := 0.0001]
        }
        log_prop_sampling[, LOG_PROP_SAMPLING := log(prop_sampling)]
        
        # check the order of ages is correct
        log_prop_sampling <- log_prop_sampling[order(AGE_TRANSMISSION.SOURCE, AGE_INFECTION.RECIPIENT)]
        stopifnot(df_age[, AGE_INFECTION.RECIPIENT] == log_prop_sampling[, AGE_INFECTION.RECIPIENT])
        stopifnot(df_age[, AGE_TRANSMISSION.SOURCE] == log_prop_sampling[, AGE_TRANSMISSION.SOURCE])
        
        # prop sampling
        log_prop_sampling_array[i, j, k,] = log_prop_sampling[, LOG_PROP_SAMPLING]
        
        # was the recipient sampled
        n_sampling_index[i, j, k] = log_prop_sampling[, sum(prop_sampling != 0.0001)]
        sampling_index[,i, j, k] = rep(-1, nrow(log_prop_sampling))
        sampling_index[1:n_sampling_index[i, j, k], i, j, k]  = log_prop_sampling[, which(prop_sampling != 0.0001)]
      }
    }
  }
  
  stan_data[['log_prop_sampling']] = log_prop_sampling_array
  stan_data[['n_sampling_index_y']] = n_sampling_index
  stan_data[['sampling_index_y']] = sampling_index
  
  return(stan_data)
}


add_init <- function(stan_data){
  stan_init <- list()
  
  # for fit inland and fishing together
  stan_init[['log_beta_baseline']] = 0
  stan_init[['log_beta_baseline_contrast_community']] = 0
  stan_init[['log_beta_baseline_contrast_direction']] =  0
  stan_init[['log_beta_baseline_contrast_round']] = array(0, dim = c(stan_data[['N_ROUND']] - 1, stan_data[['N_DIRECTION']],  stan_data[['N_COMMUNITY']]))
  stan_init[['rho_gp1']] = array(2, dim = c(stan_data[['N_DIRECTION']]))
  stan_init[['rho_gp2']] = array(2, dim = c(stan_data[['N_DIRECTION']]))
  stan_init[['alpha_gp']] = array(1, dim = c(stan_data[['N_DIRECTION']]))
  # 
  ## for fit inland and fishing independetly
  # stan_init[['log_beta_baseline']] = array(0, dim = c(stan_data[['N_COMMUNITY']]))
  # stan_init[['log_beta_baseline_contrast_direction']] =  array(0, dim = c(stan_data[['N_COMMUNITY']]))
  # stan_init[['log_beta_baseline_contrast_round']] = array(0, dim = c(stan_data[['N_ROUND']] - 1, stan_data[['N_DIRECTION']],  stan_data[['N_COMMUNITY']]))
  # stan_init[['rho_gp1']] = array(2, dim = c(stan_data[['N_DIRECTION']],stan_data[['N_COMMUNITY']]))
  # stan_init[['rho_gp2']] = array(2, dim = c(stan_data[['N_DIRECTION']],stan_data[['N_COMMUNITY']]))
  # stan_init[['alpha_gp']] = array(1, dim = c(stan_data[['N_DIRECTION']],stan_data[['N_COMMUNITY']]))
  # 
  return(stan_init)
}

