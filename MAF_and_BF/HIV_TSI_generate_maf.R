require(data.table)

.read <- function(x){
  if(grepl('.csv$', x)){return(as.data.table(read.csv(x)))}
  if(grepl('.rds$|.RDS$',x)){return(as.data.table(readRDS(x)))}
}

outdir <- '~/Documents/Box/2021/phyloTSI'
indir.repo <- '~/git/phyloflows'
indir.deepsequence_analyses <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata/'

file.path.patstats <- file.path(indir.deepsequence_analyses, 'patstats_02_05_30_min_read_100_max_read.csv')
file.path.phscinput <- file.path(indir.deepsequence_analyses, '210120_RCCSUVRI_phscinput_runs.rds')
phyloTSI.repo <- '~/git/HIV-phyloTSI-main'
# I don t trust this one, am I sure this was produced by the prev script? TODO
# or I can use the pangea2bflocs_20220103.rds otherwise
# file.bflocs <- file.path(outdir,'BF_file_locations_RCCS12.RDS') 
file.bflocs <- file.path(outdir, 'bfloc2hpc_20220103.rds')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.path.meta.data.rccs.2 <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')


patstats <- .read(file.path.patstats)
patstats[, `:=` (X=NULL, X.1=NULL)]


##################################
# Guide
##################################

# phsc <----> PatStats
# AID           AID                     bflocs
#  L________________________>   INDEX  <----> HPC_LOC


##############################
# Get patstats AIDs
##############################

patstats.ids <- unique(patstats$host.id)
N <- length(patstats.ids)
cat( sum(grepl('CNTRL-',patstats.ids)),' out of ',length(unique(patstats.ids)),' ids in the PatSTATS file are control\n')
ids  <- unique(gsub('CNTRL-','',patstats.ids))
cat('Reffering to ', length(unique(ids)), ' unique AID (ie excluding CTRL) \n')

if(0)
{
  # CNTRL-AIDs always appear as AIDs (non CNTRL).
  gsub('CNTRL-','',grep('C',patstats.ids, value=T)) -> tmp
  tmp %in% ids
  tmp1 <- lapply(ids,
                 FUN=function(id){
                   id <- c(id, paste0('CNTRL-', id))
                   pt <- patstats[host.id %in% id]
                   patstats[, all(id %in% host.id),by=PTY_RUN][, any(V1)]
                 })
  any(unlist(tmp1))
}


##################################
# Compare PTY_RUNs with Phyloscanner input to get PREFIXes
##################################
phsc.input <- .read(file.path.phscinput)
phsc.input[, AID:=gsub('-fq[0-9]$', '', RENAME_ID)]
phsc.input[, PREFIX := gsub('_remap(.*?)$','',basename(SAMPLE_ID))]

idx <- patstats[, sort(unique(PTY_RUN))]; all(idx == phsc.input[, sort(unique(PTY_RUN))])
tmp  <- unique(patstats[, .(PTY_RUN, host.id )])
tmp1 <- unique(phsc.input[, .(PTY_RUN, RENAME_ID, PREFIX)])
tmp[, AID := gsub('^CNTRL-', '',host.id)]
tmp1[, AID := gsub('-fq[0-9]$','',RENAME_ID)]
setkey(tmp, 'PTY_RUN', 'AID'); 
setkey(tmp1, 'PTY_RUN', 'AID')

# Check that all AIDs in the Patstats are in the phsc.input ones before merging
for (i in idx){stopifnot(all(tmp[PTY_RUN == 1, AID] %in% tmp1[PTY_RUN == 1, AID]))}
tmp <- merge(tmp, unique(tmp1[, .(PTY_RUN, AID, PREFIX)]), by=c('PTY_RUN','AID'), all.x=T)

# Now, for a check. Are different prefixes associated with different AIDs in different host.ids?
tmp1 <- tmp[, length(unique(PREFIX)) ,by='host.id']
# Thing is, in PHSC, to each AID can be associated multiple prefixes...
# How do I know which one was used in the PatStats? (All of them!)
# Should I only chose the one with the latest date of collection? 

# Same bf files are used for same AID in different PTY_RUN
# This mean I can summarise the data.table ignoring PTY_RUN
tmp[, {z <- unique(PREFIX); list(PREFIX == z)} ,by='AID'][, all(V1)]
tmp[, PTY_RUN:=NULL]; tmp <- unique(tmp)

# Now, for all AIDs in the patstats file, I want to find the locations of the bf csv files

if(0)
{
  # TODO: maybe here I can check how many we have stored on the HPC.
  tmp1 <- merge(tmp[,.(AID, PREFIX)], unique(phsc.input[,.(PREFIX, SAMPLE_ID)]), by='PREFIX')
  pre <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
  tmp1[, SAMPLE_ID := file.path(pre, SAMPLE_ID)]
}


################
# Ready to get locs to obtain mafs (?)
###############

stopifnot(all(tmp$host.id %in% patstats$host.id) & all(patstats$host.id %in% tmp$host.id)) 

# Or I think I should have mappings from the HPC
bflocs <- as.data.table(readRDS(file.bflocs))
bflocs <- unique(bflocs[, .(PREFIX, FULL, HPC_LOC)])

if( bflocs[, !grepl('Documents/Box', FULL[1])] ){ bflocs[, FULL:=gsub('Documents','Documents/Box',FULL)]}

bflocs <- merge(tmp, bflocs, all.x=T, by='PREFIX')


# bflocs_tmp <- unique(bflocs[,.(host.id, HPC_LOC, FULL)])
# bflocs_tmp <- bflocs_tmp[!is.na(FULL)] 
# bflocs_tmp <- bflocs_tmp[, list(HPC_LOC=HPC_LOC[1],FULL=FULL[1]) ,by=host.id]
# bflocs_tmp <- merge(bflocs_tmp, unique(phsc.input[,.(PREFIX, PANGEA_ID)]), by=)

############
# summaries
###########
cat('Out of ', nrow(bflocs), ' PREFIXes in the PatStats file:\n',
    ' -',bflocs[, sum(file.exists(FULL))], ' have an associated bf file on the Fraser dir\n',
    ' -', bflocs[file.exists(FULL), sum( ! (grepl('^NA',HPC_LOC) | is.na(HPC_LOC) ))], ' have an associated bf file on the HPC\n',
    ' -',bflocs[, sum(!file.exists(FULL))], 'do not have an associated bf file on the Fraser dir\n')

tmp1 <- bflocs[,list(FULL=!all(is.na(FULL))),by='AID']
cat('Out of ', nrow(tmp1), ' AIDs in the PatStats file:\n',
    ' -',tmp1[, sum(!FULL)], ' do not have an associated bf file on the Fraser dir\n')

cat('Prefixes without file:\n')
bflocs[AID %in% tmp1[FULL == FALSE, AID], PREFIX]

tmp1 <- bflocs[,list(FULL=!all(is.na(FULL))),by='host.id']
PTY_incomplete <- patstats[host.id %in% tmp1[FULL == FALSE, host.id], sort(unique(PTY_RUN))]
cat(length(PTY_incomplete), 'out of ',max(patstats$PTY_RUN), 
    'PTY_RUNs do not have bf files for at least one of their AIDs\n')
# should not be an issue though



###########
# Merge with meta.data to obtain PANGEA_IDs with latest collection date
###########

# TODO: if we rerun, also use bf files for individuals that DO NOT have a collection date

tmp1 <- unique(phsc.input[ grepl('RCCS', PANGEA_ID) , list(PREFIX=PREFIX, PANGEA_ID=gsub('RCCS_','',PANGEA_ID)) ])
# meta.rccs.1 <- load(file.path(indir.deepsequence_analyses, 'RakaiPangeaMetaData_v2.rda'))
# meta.rccs.1 <- as.data.table(rccsData)
meta.rccs.2 <- as.data.table(read.csv(file.path.meta.data.rccs.2))

stopifnot(all(tmp1$PANGEA_ID %in% meta.rccs.2$pangea_id))
tmp2 <- unique(meta.rccs.2[, .(pangea_id, visit_dt)])
tmp2[, visit_dt:=as.Date(visit_dt,format='%Y-%m-%d')]

tmp1 <- merge(tmp1, tmp2, by.x='PANGEA_ID', by.y='pangea_id', all.x=T)
bflocs <- merge(bflocs, tmp1, by='PREFIX')
bflocs <- bflocs[, {z <- which.max(visit_dt);
                list(PREFIX=PREFIX[z], PANGEA_ID=PANGEA_ID[1], AID=AID[1],
                     FULL=FULL[z], HPC_LOC=HPC_LOC[z], visit_dt=visit_dt[z])
                } ,by='host.id']

stopifnot(bflocs[, .N == length(unique(host.id))])




####################################################
# Lele's function. Gets bflocs and computes the MAF
####################################################

# 1_ bflocs     2_ output name    3_ PatStats


generate_sample <- function(bflocs, outfilename=NULL, patstats=NA){
  # Recall to change FULL with HPC_LOC when running on the HPC
  
  if(is.null(outfilename)){
    date <- format(Sys.Date(),'%Y%m%d')
    outfilename <- paste0('MAF_', date, '.csv');
    outfilename <- file.path(outdir, outfilename)
  }
  
  # decide what to do with NAs and multiple bf files
  # bflocs_tmp <- unique(bflocs[,.(host.id, HPC_LOC, FULL)])
  # bflocs_tmp <- bflocs_tmp[!is.na(FULL)] 
  bflocs_tmp <- bflocs_tmp[, list(HPC_LOC=HPC_LOC[1],FULL=FULL[1]) ,by=host.id]
  # Copy this part here in summary_functions.R to get date of collection
  
  # Lele's part
  MAF_matrix<-matrix(0,ncol=10000,nrow=(nrow(bflocs_tmp)+1))
  MAF_matrix[1,]<-seq(1,10000)
  rownames(MAF_matrix)<-c('pos',bflocs_tmp$host.id)
  idx_nobf <- c()
  
  N <- nrow(bflocs_tmp)
  for(i in 1:N){
    
    cat(paste0(i, ' out of ',N,': ', bflocs_tmp$V2[i],'\n'))
    basefile <- ifelse(dir.exists('/home/andrea'), bflocs_tmp$FULL[i], bflocs_tmp$HPC_LOC[i])
    
    if (file.exists(basefile))
    {
      basefile<-read.csv(basefile,header = T)
      indexes_HXB2_pos<-which(basefile$Position.in.B.FR.83.HXB2_LAI_IIIB_BRU.K03455!='-')
      sample_HXB2_pos<-as.numeric(basefile[indexes_HXB2_pos,1])
      sample_MAFs<-apply(basefile[indexes_HXB2_pos,4:7], 1, function(x) 1-(max(x,na.rm = T)/sum(x,na.rm = T)))
      sample_MAFs<-as.numeric(gsub(NaN, 0,sample_MAFs))
      MAF_matrix[(1+i),sample_HXB2_pos]<-sample_MAFs 
    }else{
      cat('BF file does not exist')
      idx_nobf <- c(idx_nobf, i)
      MAF_matrix[(1+i),]<-NA 
    }
  }
  
  # Consider putting this outside of the function! 
  write.table(MAF_matrix,file = outfilename,quote = F,sep = ',',col.names = F) # consider col.names=T
  outfilename <- gsub('.csv$','_noNA.csv',outfilename) # Consider removing
  # tmp <- which(is.na(MAF_matrix[, 2]))
  # write.table(MAF_matrix[-tmp, ],file = outfilename,quote = F,sep = ',',col.names = F)
}

maf <- '~/Documents/Box/2021/phyloTSI/MAF_20220118.csv'
#generate_sample(maf)

# Now run Tanya's

date(as_date)
outdir <- paste0(outdir, '/outputs_',format(Sys.Date(), "%y%m%d"))

if(!dir.exists(outdir))
{
  dir.create(outdir)
  cmd <- paste0('source HIV_TSI_start.sh ',file.path.patstats,' ',maf,' ', outdir)
  system(cmd)  
}



# if(0)
{
####################################################
# Can try to run Tanya's Program
####################################################

# CAN I MOVE THIS TO START-HIV-TSI.R?

# WANT: bash script that runs Tanya s script for each individual PTY run
# This means:
#     - Get the PatStats file. 
#     - Use R to subset it according to PTY_identifier 
#     - Use subsetted PS + MAF to run Tanya's
#     - Delete subsetted PS

if(is.null(outfilename)){outfilename <- "~/Documents/2021/phyloTSI/MAF_20220104.csv"}
file.MAF <- outfilename
outdir2 <- file.path(outdir, 'outputs')

PTY_indices = 1:max(patstats$PTY_RUN) 
PTY_indices <- which(!PTY_indices %in% PTY_incomplete)
output_pty <- paste0('tsi_', PTY_indices, '.csv')
# set diff with PTY indices

PTY_idx <- 1

# Maybe we do not even need to loop? WE DO
cmd <- paste('#!/bin/sh','conda activate hivphylotsi','', sep='\n')
# make an outdirectory file if it doesn t exist
cmd <- paste0(cmd, 'mkdir -p ',outdir2,'\n')

# Subset Patstats
# First, remove first two useless columns and " symbols
# Then, filter according to PTY and remove the column
# Where do I want to store this?
tmp_pat <- file.path(outdir2, paste0('tmp_pat',PTY_idx,'.csv'))
tmp_maf <- file.path(outdir2, paste0('tmp_maf',PTY_idx,'.csv'))
file.maf <- '~/Documents/2021/phyloTSI/MAF_20220104.csv'  

cmd <- paste0(cmd, 'head -q -n 1 ',file.path.patstats,"| cut -d , -f 1,2,3 --complement | tr -d \'\"\' > ", tmp_pat, '\n')
cmd <- paste0(cmd, 'cut -d , -f 1,2 --complement ', file.path.patstats, "| tr -d \'\"\' ") 
# cmd <- paste0(cmd, "| awk -F , '$1 == ", PTY, "' | cut -d , -f 1 --complement >> ",tmp_pat," \n")
cmd <- paste0(cmd, "| awk -F , '$1 == ", PTY, "' | cut -d , -f 1 --complement | sed 's/\\.5//' >> ",tmp_pat," \n")

# there is a strange mistake in the maf, where x.positions are NOT integers, rather they end by.5 ...

cmd <- paste0(cmd, 'echo "------ Create MAF -------" \n')
cmd <- paste0(cmd, "Rscript ", file.path(indir.repo, 'MAF_and_BF', 'run_HIV_TSI.R'),
              ' -idx ', PTY_idx,
              ' -p ', tmp_pat,
              ' -m ', file.maf, '\n')
# I still need the first line on top of the above 



# Do I also need to subset the MAF? Probably not

# Run Tanya
cmd <- paste0(cmd, 'echo "------ Run Tanya -------" \n')
cmd <- paste0(cmd, paste('python',file.path(phyloTSI.repo,'HIVPhyloTSI.py'),
                         '-d',file.path(phyloTSI.repo,'Model'),
                         '-p', tmp_pat,
                         '-m', tmp_maf,
                         '-o', output_pty, '\n'))

cmd <- paste0(cmd, 'rm ', tmp_pat, '\n')
cmd <- paste0(cmd, 'rm ', tmp_maf, '\n')
cat(cmd[[1]])
}