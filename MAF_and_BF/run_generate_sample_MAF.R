# For each PatStatFile
# look at the host.id in the file and find where the corresponding base freq files live
# Then run the 'run_generate_sample_MAF.R' with the corresponding file locs
require(data.table)

indir.deepsequence_analyses <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
file.path.patstats <- '210325_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/patstats_02_05_30_min_read_100_max_read.csv'
file.path.patstats <- file.path(indir.deepsequence_analyses, file.path.patstats)
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.path.meta.data.rccs.1 <- file.path(indir.deepsequence_analyses, 'RakaiPangeaMetaData_v2.rda')

if(1)
{
  outdir <- '~/Documents/2021/phyloTSI'
  indir.repo <- '~/git/phyloflows'
  indir.deepsequencedata <- '~/Documents/ratmann_pangea_deepsequencedata/'
  indir.deepsequence_analyses <- '~/Documents/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  file.path.patstats <- file.path(indir.deepsequence_analyses, 'patstats_02_05_30_min_read_100_max_read.csv')
  file.bf.locs <- file.path(indir.repo,'MAF_and_BF','BF_file_locations_RCCS12.RDS') 
  file.anonymisation.keys <-  "~/Documents/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI//important_anonymisation_keys_210119.csv"
  file.path.meta.data.rccs.1 <- "~/Documents/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI//RakaiPangeaMetaData_v2.rda"
  file.path.meta.data.rccs.2 <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')
  file.path.PHSCruns <- file.path(indir.deepsequence_analyses, '210120_RCCSUVRI_phscinput_runs.rds')
}

# Maybe use selected samples RCCS to match AID to BF files: 



##################################
# Guide
######################qu############

#  patstats    'anon.keys'                      'meta.rccs'+PHSCruns            'bf.locs'
#    AID  <--------------->  STUDY_ID / PT_ID  <---------------  PANGEA_IDs  <-----------> BaseFrequency files 
#     L______________________________________________________________\ 
#                              phsc.input   

##################################
## Take the patstats file, and find the Base Freqs locations for all individuals in it
##################################

patstats <- as.data.table(read.csv(file.path.patstats))
tmp <- patstats[grep('CNTRL-',host.id), unique(host.id)]
tmp1 <- patstats[! grep('CNTRL-',host.id), unique(host.id)]
tmp <- gsub('CNTRL-','',tmp)
tmp[which(tmp %in% tmp1)]
tmp <- c(tmp,paste0('CNTRL-', tmp))
tmp <- patstats[host.id %in% tmp]
tmp[, ORDER:=gsub('CNTRL-','',host.id)]
setkey(tmp,tree.id, PTY_RUN, ORDER)
tmp[, ORDER := NULL]

ids <- unique(patstats$host.id)
N <- length(ids)
cat( sum(grepl('CNTRL-',ids)),' out of ',length(unique(ids)),' ids in the PatSTATS file are control\n')
patstats.ids<- unique(gsub('CNTRL-','',ids))
# length(patstats.ids) # 2983

##################################
# Study Phyloscanner run file
##################################
phsc.input <- as.data.table(readRDS(file.path.PHSCruns))
phsc.input[, AID:=gsub('-fq[0-9]$', '', RENAME_ID)]
tmp <- unique(phsc.input[, .(RENAME_ID, PANGEA_ID, AID)])
tmp[, lapply(.SD, function(x){length(unique(x))}), .SDcols=c('RENAME_ID', 'PANGEA_ID')]
tmp[, .N,by='PANGEA_ID'][N>1, if(length(PANGEA_ID) == 0){cat('phsc.input as bijection between RENAME_ID and PANGEA_ID')}]
if(all(patstats.ids %in% tmp$AID)){cat('All patstats IDs are contained in the phsc.input')}

patstats.keys <- phsc.input[AID %in% patstats.ids]
patstats.keys[, `:=` (IDCLU=NULL, CLU_SIZE=NULL, ID_TYPE=NULL, PTY_RUN=NULL, PTY_SIZE=NULL, UNIT_ID=NULL)]
# phsc.input[, length(unique(PANGEA_ID)) ,by='AID'][, table(V1)]
patstats.keys[, PANGEA_ID2 := gsub('(.*?)_','',PANGEA_ID)]
patstats.keys <- unique(patstats.keys)

tmp <- unique(patstats.keys[, .(PANGEA_ID2, AID), ])
tmp[, lapply(.SD, function(x){length(unique(x))}) , .SDcols=c('PANGEA_ID2', 'AID')]
tmp[,length(unique(PANGEA_ID2)),by='AID'][, table(V1)]

# Now load the locations on the HPC
bf.locs <- as.data.table(readRDS(file.bf.locs))
bf.locs <- bf.locs[, .(PANGEA_ID, FILE )]

tmp <- unique(patstats.keys[!grepl('MRC', PANGEA_ID)]$PANGEA_ID2)
tmp1 <- unique( bf.locs$PANGEA_ID)
cat('Out of ', patstats.keys[, length(unique(PANGEA_ID))], ' PANGEA_IDs in the PS file, ', length(tmp), ' are NOT from the MRC\n')
tmp <- tmp %in% tmp1
cat(sum(tmp),'RCCS PANGEA_IDs have a corresponding file in bf.locs csv\n')


patstats.keys[, SAMPLE_ID2 := basename(SAMPLE_ID)]
tmp <- unique(patstats.keys[, .(SAMPLE_ID, SAMPLE_ID2)])
if(nrow(tmp[, .N ,by=SAMPLE_ID2][N>1, ])==0){cat('No info lost')}

# We have all entries in Patstats keys except 27 MRC and 2 RCCS
# 2 of the MRC are in PANGEA2_MRC, can we maybe find those
tmp <- unique(patstats.keys$PANGEA_ID2)
sum(tmp %in% bf.locs$PANGEA_ID); length(tmp)
tmp <- tmp[which(! tmp %in% bf.locs$PANGEA_ID)]
patstats.keys[PANGEA_ID2 %in% tmp, table(grepl('MRC', PANGEA_ID)) ]

tmp <-  bf.locs[, all(is.na(FILE)), by=PANGEA_ID][V1==TRUE, PANGEA_ID]
tmp <- unique(bf.locs[PANGEA_ID %in% tmp, ])
tmp <- patstats.keys[PANGEA_ID2 %in% tmp$PANGEA_ID]
cat(nrow(tmp),' Base Frequency files were not found in the 2 Fraser Uploads')

# merge with bf locs
# First merge by SAMPLE_ID, then, for the remaining ones, use PANGEA_ID
tmp <- copy(bf.locs)
tmp[, SAMPLE_ID2 := basename(FILE)][, SAMPLE_ID2 := gsub('_BaseFreqs_WithHXB2.csv$','',SAMPLE_ID2)]
stopifnot(length(unique(tmp$SAMPLE_ID2)) == length(unique(tmp$FILE)))

# mean(patstats.keys$SAMPLE_ID2 %in% tmp$SAMPLE_ID2)
cat(sum(! patstats.keys$SAMPLE_ID2 %in% tmp$SAMPLE_ID2),' out of ',length(patstats.keys$SAMPLE_ID2),
    ' SAMPLE IDs do not have a matching file in bf.locs \n')
# mean(tmp$SAMPLE_ID2 %in% patstats.keys$SAMPLE_ID2)

patstats.keys <- merge(patstats.keys, tmp[,.(FILE, SAMPLE_ID2)], all.x=TRUE ,by='SAMPLE_ID2')
patstats.keys[, DUMMY := 1:.N]
patstats.keys[, grepl(PANGEA_ID2, PANGEA_ID), by=DUMMY ][, all(V1)]
patstats.keys[, `:=` (PANGEA_ID=NULL, DUMMY=NULL)]

tmp <- patstats.keys[,any(file.exists(FILE)),by=PANGEA_ID2]
cat('For how many PREFIXES can we find an existing FILE ID?\n') # redundant
cat(table(tmp$V1))

tmp1 <- tmp[V1==FALSE, unique(PANGEA_ID2)]
tmp <- merge(patstats.keys[PANGEA_ID2 %in% tmp1][, -"FILE"],
             bf.locs[PANGEA_ID %in% tmp1], by.x='PANGEA_ID2', by.y = 'PANGEA_ID', all.x = T)
tmp[,DUMMY := 1:.N]
if(tmp[,grepl(SAMPLE_ID2, FILE),by=DUMMY][, all(V1==FALSE)])
{
  cat('Out of the ',length(tmp1),' SAMPLES identifiers that do not have a matching BF file:\n')
  cat(' - ',tmp[file.exists(FILE), length(unique(PANGEA_ID2))],' have a BF file matching the PANGEA_ID (but not the rest)','\n')

  tmp[, DUMMY:= NULL]
  patstats.keys <- unique(rbind(tmp, patstats.keys[! PANGEA_ID2 %in% tmp$PANGEA_ID]))
}
# tmp[file.exists(FILE),]

patstats.keys[, DUMMY:= 1:.N]
patstats.keys[, V1 := grepl(SAMPLE_ID2, FILE), by=DUMMY]
tmp <- patstats.keys[, table(V1),by=file.exists(FILE),]$V1
cat('Out of ',sum(tmp),' PREFIXES:\n - ', tmp[3], ' had PREFIX matching bf files\n',
    '- ', tmp[2], ' had PANGEA_ID but not PREFIX matching bf file\n - ', tmp[1], ' had no associated bf files')
patstats.keys[, `:=` (DUMMY=NULL, V1=NULL)]

if(0)
{
  # look at 'desperate cases'
  tmp1 <- patstats.keys[,all(!file.exists(FILE)),by=PANGEA_ID2][V1==T, PANGEA_ID2]
  stopifnot(all(! tmp1 %in% bf.locs ))
  # patstats.keys[PANGEA_ID2 %in% tmp1]
}


if(0)
{
##################################
## Begin translation
##################################

# Check PatStats id are in the anonymisation key file
anon.keys <- as.data.table(read.csv(file.anonymisation.keys))
stopifnot(all(patstats.ids %in% anon.keys$AID))
anon.keys <- anon.keys[, list(PT_ID = gsub('^RK-','',PT_ID), AID=AID)]

# However need to translate the AID in terms of Pangea IDs.
# We need metadata from rccs
env.meta <- new.env()
load(file.path.meta.data.rccs.1, envir = env.meta)
meta.rccs.1 <- as.data.table(env.meta$rccsData)
meta.rccs.1 <- unique(meta.rccs.1[, .(Pangea.id, studyid)])
meta.rccs.1 <- meta.rccs.1[!(is.na(Pangea.id) & is.na(studyid))]
meta.rccs.1[, PT_ID := gsub('R[0-9][0-9]$|R[0-9][0-9]S$', '', studyid)]
setnames(meta.rccs.1,'Pangea.id','PANGEA_ID')
if(all(meta.rccs.1$PT_ID %in% anon.keys$PT_ID)){cat('all entries in meta.rccs.1 are in the anonymised keys')}
meta.rccs.2 <- as.data.table(read.csv(file.path.meta.data.rccs.2))
meta.rccs.2 <- unique(meta.rccs.2[, .(pt_id, pangea_id)])
meta.rccs.2[, PT_ID := gsub('RK-', '', pt_id)]
meta.rccs.2 <- meta.rccs.2[! PT_ID %in% meta.rccs.1$PT_ID,][, `:=` (pt_id = NULL, studyid = NA)]
setnames(meta.rccs.2,'pangea_id','PANGEA_ID')
if(all(meta.rccs.2$PT_ID %in% anon.keys$PT_ID)){cat('all entries in meta.rccs.2 are in the anonymised keys')}

meta.rccs <- rbind(meta.rccs.1, meta.rccs.2)

tmp <- meta.rccs[, sum(!is.na(PANGEA_ID)) ,by='PT_ID']
cat(nrow(tmp[V1 > 1]), ' out of ', meta.rccs[, length(unique(PT_ID))] ,' PT_ID have more than 1 corresponding PANGEA ID\n')
cat(nrow(tmp[V1 == 0 ]), ' out of ', meta.rccs[, length(unique(PT_ID))] ,' PT_ID have no corresponding PANGEA ID\n')
# meta.rccs[studyid %in% tmp[V1 > 1,studyid]]
cat(meta.rccs[!is.na(PANGEA_ID), length(unique(PT_ID)) ,], ' out of ', meta.rccs[,length(unique(PT_ID)),], ' PT_IDs have a corresponding PANGEA ID\n' )

# merge to get PangeaID
tmp <- patstats.ids %in% anon.keys$AID
cat(sum(tmp), 'out of ', length(tmp),'PatStats IDs are not found in the anonymised IDs \n')
anon.keys <- merge(anon.keys, meta.rccs, by.x='PT_ID', by.y='PT_ID', all.x = TRUE)

# check consistent assignment
anon.keys[, DUMMY:=1:.N]
stopifnot(anon.keys[!is.na(studyid), grepl(PT_ID, studyid) ,by=DUMMY][,all(V1)])
anon.keys[, DUMMY:=NULL]
cat(anon.keys[!is.na(studyid), length(unique(PANGEA_ID)),by=studyid][V1 > 1, .N],' studyids are associated with multiple PANGEA IDs \n')

patstats.keys <- anon.keys[AID %in% patstats.ids]
tmp <-  patstats.ids[!patstats.ids%in% anon.keys[, AID]]
cat( length(tmp),' anonymised IDs  in PatSTATS are not in the anonymised IDs list of length ', anon.keys[, length(unique(AID))],'.\n')
stopifnot(all(patstats.ids%in% anon.keys[, AID]))

# Now load the locations on the HPC
bf.locs <- as.data.table(readRDS(file.bf.locs))
bf.locs <- bf.locs[, .(PANGEA_ID, FILE )]

patstats.keys <- merge(patstats.keys, bf.locs, by.x='PANGEA_ID', by.y='PANGEA_ID', all.x=TRUE)
patstats.keys <- unique(patstats.keys)
stopifnot(all(patstats.ids %in% patstats.keys$AID))
stopifnot(all(patstats.keys$AID %in% patstats.ids))

cat('Out of ',nrow(patstats.keys), 'rows in patstats.keysy: the nuber of NA per column is:\n')
patstats.keys[, lapply(.SD, function(x){sum(is.na(x))} ) ,.SDcols=colnames(patstats.keys)]

tmp <- patstats.keys[,!all(is.na(FILE)),by='AID'][, table(V1)]
cat('Out of ',tmp[1]+tmp[2],' patstats AIDs, ',tmp[1],' do not have a corresponding BF file\n')
tmp <- patstats.keys[,!all(is.na(FILE)),by='AID'][V1 == FALSE, AID]
tmp <- anon.keys[AID %in% tmp, PT_ID] # not all are MRC
tmp1 <- grepl('MRC', tmp)
cat('- of these: ',sum(tmp1),' are from MRC and ',sum(!tmp1), ' are not')
cat(tmp[!tmp1])

# some AID have multiple associated entries
# For AIDs were only one file is not na, we take that file.
# Otherwise we take the first one
tmp <- patstats.keys[,.N,by='AID'][N>1, AID]
cat(length(unique(tmp)), 'AIDs have multiple associated entries in patstats.keys \n')
setkey(patstats.keys, AID)
tmp <- patstats.keys[AID %in% tmp]
tmp[,  `:=` (N_NA = sum(is.na(FILE)), N=length(FILE) ), by='AID']
tmp <- tmp[!(N - N_NA == 1 & is.na(FILE))]
tmp[, `:=` (N_NA = sum(is.na(FILE)), N=length(FILE)), by='AID']
cat(tmp[N - N_NA  > 1, length(unique(AID))], 'AIDs had multiple corresponding BF files\n')
cat('Only picking the first FILE for each\n')
tmp <- rbind(tmp[N-N_NA <= 1, .(AID, PANGEA_ID, PT_ID, studyid, FILE)],
             tmp[ N - N_NA > 1, list(PANGEA_ID=PANGEA_ID[1] ,PT_ID=PT_ID[1] , studyid=studyid[1], FILE=FILE[1]), by='AID'])
stopifnot(nrow(tmp) == tmp[, length(unique(AID))] ) 
patstats.keys <- rbind(patstats.keys[ !AID %in% tmp$AID, ], tmp)

# Summarise and study missing files
tmp <- patstats.keys[, all(is.na(FILE)) ,by=AID]
cat('We have FILEs for ', nrow(tmp[V1==F]) ,' entries  out ', length(patstats.ids) ,' individuals in the PatStats file\n')
tmp <- tmp[V1 == T, AID]
tmp <- patstats.keys[AID %in% tmp]
cat('Out of ', tmp[,length(unique(AID))],'AIDs without associated Base Frequency File:\n')
cat(' - ',tmp[, sum(!is.na(PANGEA_ID))], ' have an associated PANGEA_ID')
tmp <- tmp[is.na(PANGEA_ID), sum(!grepl('MRC', PT_ID))]
if(tmp==0)
{
  cat(' - All those without an associated PANGEA_ID are in the MRC.')
}else{
  cat(' - ', tmp, ' are not in the MRC')
}

}

###################################
# Guide 2
##################################
#               requires 
# HIVphyloTSI   <-------   PatStats   + MinorAlleleFreqs(MAF)
#   - MAF       <-------   BaseFrequencyFiles                                   


###################################
# Now run Lele's script
##################################

# make base frequency file
tmp_name <- gsub('-', '', Sys.Date())
tmp <- gsub('^patstats_|.csv$','',basename(file.path.patstats))
file.path.MAF <- paste0(tmp_name,'_', tmp, '_MAF.csv')
file.path.MAF <- file.path(outdir, file.path.MAF)
file.path.bfpath <- paste0(tmp_name, '_', tmp, '_bfpath.csv')
file.path.bfpath <-  file.path(outdir, file.path.bfpath)

# If more than one base freqs file exist, take the one labelled with -fq1
tmp <- patstats.keys[, sum(file.exists(FILE)), by=AID]
# tmp[, table(V1)]; patstats.keys[AID %in% tmp[V1 > 1, AID], length(unique(RENAME_ID)) == .N, by=AID ][, table(V1)]
patstats.keys[,sum(RENAME_ID == paste0(AID, '-fq1')),by=AID][,if(all(V1==1)){cat('Unique assignments of -fq1 for each AID')}]
patstats.keys2 <- patstats.keys[RENAME_ID == paste0(AID, '-fq1'), ]
setcolorder(patstats.keys, c('FILE', 'AID'))
patstats.keys2 <- patstats.keys2[, .(FILE, AID)]

# Patstats IDs included some 'CNTRL-' AIDs. Fix them back
tmp <- data.table(CNTRL = grep('CNTRL-', ids, value=T))
tmp[, AID := gsub('CNTRL-','',CNTRL)][ AID %in% ids, AIDinPS:= AID ]
patstats.keys2[AID %in% tmp$AID, AID:=paste0('CNTRL-', AID)]
tmp <- patstats.keys2[AID %in% tmp[!is.na(AIDinPS), CNTRL],]
tmp[,  AID:=gsub('CNTRL-','',AID)]
patstats.keys2 <- rbind(patstats.keys2, tmp)
stopifnot(all(patstats.keys2$AID %in% ids))
stopifnot(all(ids %in% patstats.keys2$AID))

# Save
write.table(patstats.keys2,  file.path.bfpath, row.names = FALSE, col.names = FALSE, sep=',')

patstats.keys2 <- as.data.table(read.csv(file.path.bfpath), header=FALSE)
tmp <- colnames(patstats.keys2)
tmp[grep("NA",tmp)] <- NA
patstats.keys2 <- rbind(as.data.table(t(tmp)), patstats.keys2, use.names=FALSE)
colnames(patstats.keys2) <- c('FILE', 'AID')

tmp <- gsub('^patstats_|.csv$','',basename(file.path.patstats))
file.path.bfpath2 <- paste0(tmp_name, '_', tmp, '_bfpath2.csv')
patstats.keys2 <- patstats.keys2[!is.na(FILE)]
write.table(patstats.keys2,  file.path.bfpath2, row.names = FALSE, col.names = FALSE, sep=',')


# file.path.bfpath2 <- "~/Documents/2021/MAFs/20211201_02_05_30_min_read_100_max_read_bfpath2.csv"
# write.table(patstats.keys2[!is.na(FILE)], file.path.bfpath2 , row.names = FALSE, col.names = FALSE, sep=',')

cmd <- paste0("Rscript generate_sample_MAF.R '", file.path.bfpath2, "' '", file.path.MAF, "'")
system(cmd)

if(0)
{
  # Include NA rows for individuals without associated Base Frequency file.
  # Done it in the code
  tmp <- as.data.table(read.csv(file.path.bfpath2))
  colnames(tmp) <- c('FILE', 'AID')
  tmp <- tmp[,{z <- read.csv(FILE); list(NR=nrow(z), NC=ncol(z), AID=AID)},by=FILE]
}

###################################
# And now Tanya's program
##################################


phyloTSI.dir <- '~/git/HIV-phyloTSI-main'

cmd <- paste0('#!/bin/bash \n')
cmd <- paste0(cmd, 'DIR=${pwd} \n' )
cmd <- paste0(cmd, 'conda activate hivphylotsi\n')
cmd <- paste0(cmd, 'cd ', phyloTSI.dir, '\n')
# cmd <- paste0(cmd, './',  phyloTSI.dir,'/')
cmd <- paste0(cmd, './HIVPhyloTSI.py -d Model -p ',file.path.patstats,' -m ',file.path.MAF,' -o ', 'HELLO.csv','\n')

# cmd <- paste0(cmd, 'python ')
# cmd <- paste0(cmd, 'HIVPhyloTSI.py -d Model -p ', file.path.patstats,
#               ' -m ', file.path.MAF ,' -o ',  file.path.MAF,'\n')
cmd <- paste0(cmd, 'cd $DIR\n')

cat(cmd, file='run_phyloTSI.sh')
# patstats

# problem with our produced MAF csv file. It is not the same format of the testmaf file.
# In ptcular, my position names are the Anonymised IDs while in the testmaf it s '2309502_739' and so forth

file.path.MAF <- '/home/andrea/Documents/2021/MAFs/20211129_02_05_30_min_read_100_max_read_MAF.csv'
tmp <- as.data.table(read.csv(file.path.MAF))
tmp1 <- file.path(phyloTSI.dir, 'ExampleInputs', 'testmaf.csv')
tmp1 <- as.data.table(read.csv(tmp1))
