make.cluster.assignments <- function()
{
  # Get files 
  tmp <- data.table(F = list.files(file.path(indir.deepsequence_analyses, '210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd'), pattern = 'workspace.rda', full.names = TRUE))
  tmp[, PTY_RUN := as.integer(gsub('^ptyr([0-9]+)_.*', '\\1', basename(F)))]
  setkey(tmp, PTY_RUN)
  
  # Load individuals in each phyloscanner cluster
  cluster_assignments <- tmp[,{
    cat(PTY_RUN, '\n')
    load(F)
    unique(dc$host.1, dc$host.2)
  }, by = 'PTY_RUN']
  
  # Make pretty and save
  setnames(cluster_assignments, 'V1', 'AID')
  cluster_assignments <- unique(cluster_assignments)
  saveRDS(cluster_assignments, file.path(outdir, 'R1519MRC_cluster_assignments.rds'))
  return(cluster_assignments)
}

print.statements.about.potential.TNet <- function()
{
  cat('\n--------------\n\n')
  file.path.rakai1516 <- file.path(indir.deepsequencedata, 'RakaiAll_output_190327_w250_s20_p25_d50_stagetwo_rerun23_min30_conf60_phylogeography_data_with_inmigrants.rda')
  file.path.chains.data.MRC <- file.path( indir.deepsequence_analyses_MRC, 'MRCPopSample_phsc_stage2_output_newali_250_HKC_phsc', 'MRC_phscnetworks.rda')

  if (!file.exists(file.path.rakai1516))
  {
    break('Cannot find Rakai 15-16 potential networks "RakaiAll_..._with_inmigrants.rda", ask Andrea\n')
  }else{
    cat('POTENTIAL TRANSMISSION NETWORK: Comparison with R15-16 analysis:\n')
    print.statements.about.potential.TNet.rakai1516(file.path.rakai1516)
  }

  if (!file.exists(file.path.chains.data.MRC))
  {
    break('Cannot find MRC potential networks "MRC_phscnetworks.rda", ask Andrea\n')
  }else{
    cat('POTENTIAL TRANSMISSION NETWORK: Comparison with Nicholas Bossa MRC analysis:\n')
    print.statements.about.potential.TNet.MRC(file.path.chains.data.MRC)
  }

  cat('\n--------------\n\n')
}




print.statements.about.potential.TNet.rakai1516 <- function(filename)
{
  ### Rakai pairs, R15-16 from Phyloscanner's paper

  # Load
  load(filename)
  
  # get pairs 
  rtpdm <- as.data.table(rtpdm)
  # tmp <- grep('POSTERIOR',colnames(rtpdm),value=T)
  # rtpdm[, lapply(.SD, range), .SDcols = tmp]
  
  pairs <- rtpdm[, .(FEMALE_RID, MALE_RID, LINK_MF)]

  # Translate Pangea IDs in terms of AIDs
  tmp1 <- anonymisation.keys[grep('^RK', PT_ID)][, X:=NULL]
  tmp1 <- tmp1[, PT_ID := gsub('RK-', '', PT_ID)]
  if(! (all(pairs$FEMALE_RID %in% tmp1$PT_ID) & all(pairs$MALE_RID %in% tmp1$PT_ID)))
  {
    cat('Not all individuals in Rakai analysis are included in the anonymisation keys')
  }
  pairs <- merge(pairs, tmp1, by.x = 'FEMALE_RID', by.y = 'PT_ID')
  setnames(pairs, 'AID', 'AID_FEMALE')
  pairs <- merge(pairs, tmp1, by.x = 'MALE_RID', by.y = 'PT_ID')
  setnames(pairs, 'AID', 'AID_MALE')
  pairs[, `:=` (MALE_RID=NULL, FEMALE_RID=NULL)]
  
  # Now can look if all pairs are in the same clusters
  tmp <- file.path(outdir, 'R1519MRC_cluster_assignments.rds')
  if(!file.exists(tmp))
  {
    stop("Could not find R1519MRC_cluster_assignments.rds /nTry running 'make.cluster.assignments()' First")
  }
  cluster_assignments <- as.data.table(readRDS(tmp))
  cluster_assignments <- cluster_assignments[, list(PTY_RUN1 = PTY_RUN[1], PTY_RUN2 = PTY_RUN[2], PTY_RUN3=PTY_RUN[3], PTY_RUN4=PTY_RUN[4]) ,by = 'AID']
  
  # For each pair, find cluster of MALE and FEMALE
  pairs <- merge(pairs, cluster_assignments, by.x='AID_FEMALE', by.y='AID', all.x=TRUE)
  tmp <-  grep('^PTY_RUN',colnames(pairs), value=TRUE)
  setnames(pairs, tmp, paste0(tmp, '_FEMALE'))
  pairs <- merge(pairs, cluster_assignments, by.x='AID_MALE', by.y='AID', all.x=TRUE)
  tmp <-  grep('^PTY_RUN[0-9]$',colnames(pairs), value=TRUE)
  setnames(pairs, tmp, paste0(tmp, '_MALE'))
  
  # Want intersection of PTY_RUN1:4_sce and PTY_RUN1:4_rec to be non-empty!
  pairs[,same_cluster := 0 < length(intersect(c(PTY_RUN1_FEMALE, PTY_RUN2_FEMALE, PTY_RUN3_FEMALE, PTY_RUN4_FEMALE),
                                              c(PTY_RUN1_MALE,  PTY_RUN2_MALE, PTY_RUN3_MALE, PTY_RUN4_MALE))),]
  table(pairs$same_cluster)
  cat(nrow(pairs), '/', nrow(pairs[same_cluster == TRUE]),  'pairs were classified in the same possible transmission network cluster at least once.\n')
}





print.statements.about.potential.TNet.MRC <- function(filename, coverage_threshold = 0.95)
{
  cat('Using a coverage threshold of', coverage_threshold, ' and a transmission given linkage threshold of 0.6:\n')
  
  # Get likely transmission pairs from Nicholas' analysis
  env <- new.env(parent = emptyenv())
  load(filename, envir=env)
  dchain.MRC <- as.data.table(env$dchain)
  rm(env)
  
  # Only keep pairs with a given coverage threshold (I dont think direction should be determined with the same threshold right?)
  dchain.MRC <- dchain.MRC[SCORE_LINKED>coverage_threshold]
  dchain.MRC[SCORE_DIR_12 <= 0.6 & SCORE_DIR_21 <= 0.6, EST_DIR:='unclear'] # need threshold for direction of transmission given linkage
  dchain.MRC[SCORE_DIR_12 > 0.6, EST_DIR:='12'] 
  dchain.MRC[SCORE_DIR_21 > 0.6, EST_DIR:='21']
  
  # find source recipient
  cat('Out of ', nrow(dchain.MRC), 'pairs in the NBs analysis:\n',
      '- ',nrow(dchain.MRC[EST_DIR == 'unclear']), ' have unclear direction of transmission.\n')
  
  # Keep unclear. Names of SOURCE-RECIPIENT are for convenience and do not necessarily imply direction
  dchain.MRC[, `:=` (SOURCE=H1, RECIPIENT=H2)]
  dchain.MRC[EST_DIR == '21', `:=` (SOURCE=H2, RECIPIENT=H1) ]
  dchain.MRC[, `:=` (H1=NULL, H2=NULL)]
  
  # Change ID type to AID to make it compatible with our analysis
  tmp <- dchain.MRC[, unique(c(SOURCE,RECIPIENT))]
  tmp1 <- meta.mrc[,.(pangea_id, pt_id)]
  cat(' - ',sum(!tmp %in% tmp1$pangea_id), ' out of ', length(tmp), ' individuals did not have an Anonymised ID (were not included in our analysis)\n')
  
  dchain.MRC <- change.MRC.IDtype(dchain.MRC)
  pairs <- dchain.MRC[!(is.na(SOURCE) | is.na(RECIPIENT))]
  pairs <- pairs[,.(RECIPIENT, SOURCE)] 
  cat('Out of ', nrow(pairs), 'pairs with well defined AID:\n',
      '- ',nrow(dchain.MRC[EST_DIR == 'unclear']), ' have unclear direction of transmission.\n')

  # Now can look if all pairs are in the same clusters
  tmp <- file.path(outdir, 'R1519MRC_cluster_assignments.rds')
  if(!file.exists(tmp))
  {
    stop("Could not find R1519MRC_cluster_assignments.rds /nTry running 'make.cluster.assignments()' First")
  }

  cluster_assignments <- as.data.table(readRDS(tmp))
  cluster_assignments <- cluster_assignments[, list(PTY_RUN1 = PTY_RUN[1], PTY_RUN2 = PTY_RUN[2], PTY_RUN3=PTY_RUN[3], PTY_RUN4=PTY_RUN[4]) ,by = 'AID']

  # For each pair, find cluster for RECIPIENT AND SOURCE
  pairs <- merge(pairs, cluster_assignments, by.x='RECIPIENT', by.y='AID', all.x=TRUE)
  tmp <-  grep('^PTY_RUN',colnames(pairs), value=TRUE)
  setnames(pairs, tmp, paste0(tmp, '_RECIPIENT'))

  pairs <- merge(pairs, cluster_assignments, by.x='SOURCE', by.y='AID', all.x=TRUE)
  tmp <-  grep('^PTY_RUN[0-9]$',colnames(pairs), value=TRUE)
  setnames(pairs, tmp, paste0(tmp, '_SOURCE'))
  

  # Want intersection of PTY_RUN1:4_sce and PTY_RUN1:4_rec to be non-empty!
  pairs[,same_cluster := 0 < length(intersect(c(PTY_RUN1_RECIPIENT, PTY_RUN2_RECIPIENT, PTY_RUN3_RECIPIENT, PTY_RUN4_RECIPIENT),
                                                  c(PTY_RUN1_SOURCE,  PTY_RUN2_SOURCE, PTY_RUN3_SOURCE, PTY_RUN4_SOURCE))),]
  

  cat(' - ', nrow(pairs[same_cluster == TRUE]),  'pairs were classified in the same possible transmission network cluster at least once.\n')
}


change.MRC.IDtype <- function(dchain.MRC)
{
  tmp <- meta.mrc[,.(pangea_id, pt_id)]
  tmp <- merge(tmp, anonymisation.keys, by.x = 'pt_id', by.y = 'PT_ID')[, .(pangea_id, AID)]
  
  dchain.MRC <- merge(dchain.MRC, tmp, by.x= 'SOURCE', by.y='pangea_id', all.x = TRUE)
  dchain.MRC <- merge(dchain.MRC, tmp, by.x= 'RECIPIENT', by.y='pangea_id', all.x = TRUE)
  dchain.MRC <- unique(dchain.MRC) # Don t know why the above created repeated rows.
  # cat('Out of ',nrow(dchain.MRC), ' source-recipient pairs in the MRC dataset\n',
  # nrow(dchain.MRC[is.na(AID.x)]), ' sources and ', nrow(dchain.MRC[is.na(AID.y)]), ' recipients had unknown AID.\n')
  dchain.MRC[, `:=` (SOURCE=AID.x, RECIPIENT=AID.y, AID.x=NULL, AID.y=NULL) ]
  # cat('removed ', nrow(dchain.MRC[(is.na(RECIPIENT) & is.na(SOURCE))]), 'where both SOURCE and RECIPIENT had unknown AID\n')
  dchain.MRC <- dchain.MRC[!(is.na(RECIPIENT) & is.na(SOURCE))]
  return(dchain.MRC)
}





