require(data.table)
require(ggplot2)

.read <- function(x){
  if(grepl('.csv$', x)){return(as.data.table(read.csv(x)))}
  if(grepl('.rds$|.RDS$',x)){return(as.data.table(readRDS(x)))}
}

git.repo <- '~/git/phyloflows'
indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
outdir <- '~/Documents/Box/2021/phyloTSI/outputs_220118/'

preds <- list.files(outdir, pattern='^tsi', full.names = T)
preds <- data.table(FILE=preds)
preds[, PTY_RUN:=gsub("^(.*?)tsi_|.csv$", '', FILE)]
preds <- preds[, read.csv(FILE), by=PTY_RUN]
orignal_hostids <- unique(preds$host.id)
tmp <- c(1,2,grep('RF',colnames(preds)))
preds <- preds[, ..tmp]

##########################################################
# Study individuals with multiple associated predictions #
##########################################################

preds[, host.id:=gsub('CNTRL-','',host.id)]
tmp <- preds[, .N, by=host.id][N>1, host.id]
preds1 <- preds[host.id %in% tmp]
setkey(preds1, host.id, PTY_RUN)

# are all 'confidence intervals' overlapping?
# check that max(min) < min(max)
preds1[, max(RF_cc025) < min(RF_cc975), by=host.id][, all(V1)]

tmp1 <- preds1[, .(mean=mean(RF_pred_sqrt)), by=host.id][,GROUP:=rep(1:6, 50)]
# setkey(tmp1, mean); tmp1[, GROUP := ceiling((1:.N)/50)]
tmp <- merge(preds1, tmp1, by='host.id')
tmp[, N:=1:.N ,by='host.id']
tmp[, .(N=max(N)), by='host.id'][, table(N)]
tmp[, .(N=max(N)), by='host.id'][N > 4, all(paste0('CNTRL-', host.id) %in% orignal_hostids)]

if(0)
{
  gg <- ggplot(tmp, aes(x=host.id, y=RF_pred_sqrt, ymin=RF_cc025, ymax=RF_cc975, color=N)) +
    geom_point() + geom_errorbar() +
    facet_grid(GROUP~., scales = 'free')  
}

#############################################################
# Combine results
#############################################################

# Want to only output one prediction per individual.
# We get the median of the estimates 
preds[, PTY_RUN := NULL]
res <- preds[,lapply(.SD, FUN=median),by=host.id]
stopifnot(length(unique(res$host.id)) == nrow(res))
tmp1 <- c('host.id', 'RF_pred_linear', 'RF_pred_min_linear', 'RF_pred_max_linear')
tmp2 <- c('AID','TSI_estimated_mean', 'TSI_estimated_min', 'TSI_estimated_max')
setnames(res, tmp1, tmp2)

#############################################################
# Get infection dates
#############################################################

indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
file.bflocs <- file.path(dirname(outdir), 'bfloc2hpc_20220103.rds')
file.path.phscinput <- file.path(indir.deepsequence_analyses, '210120_RCCSUVRI_phscinput_runs.rds')
file.anonymisation.keys <- file.path(indir.deepsequence_analyses,'important_anonymisation_keys_210119.csv')
file.path.meta.data.rccs.2 <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '200316_pangea_db_sharing_extract_rakai.csv')


ak <-.read(file.anonymisation.keys)
ak[, X:=NULL]

bflocs <- .read(file.bflocs)
# bflocs[, FULL := gsub('/andrea/Documents/','/andrea/Documents/Box/', FULL )]
bflocs[, EXISTS:=!is.na(FULL)]
bflocs <- unique(bflocs[, .(PREFIX, EXISTS)])

phsc.input <- .read(file.path.phscinput)
phsc.input[, PREFIX:=gsub('_remap', '', basename(SAMPLE_ID))]
phsc.input <- unique(phsc.input[, .(RENAME_ID, PREFIX, PANGEA_ID)])
phsc.input[,  `:=` (PANGEA_ID=gsub('RCCS_|MRCUVRI_', '', PANGEA_ID))]

meta.rccs.2 <- .read(file.path.meta.data.rccs.2)
meta.rccs.2 <- unique(meta.rccs.2[!is.na(pangea_id), .(pt_id, pangea_id, visit_dt)])
# stopifnot(meta.rccs.1[, length(unique(birthdate)),by='RCCS_studyid'][, all(V1)])
setnames(meta.rccs.2, c('pt_id', 'pangea_id'), c('STUDY_ID', 'PANGEA_ID'))
meta.rccs.2[, visit_dt:=as.Date(visit_dt, '%Y-%m-%d')]

tmp <- merge(phsc.input, bflocs, by='PREFIX', all.x=T)
stopifnot(meta.rccs.2[, .N, by='PANGEA_ID'][, all(N==1)] & tmp[, .N, by='PANGEA_ID'][, all(N==1)] )
tmp <- merge(tmp, meta.rccs.2, by='PANGEA_ID', all.x=T)

tmp[, `:=`(AID=gsub('-fq[0-9]$', '', RENAME_ID), RENAME_ID=NULL)]
tmp1 <- tmp[EXISTS == T & !is.na(visit_dt), 
            {
              z <- which.max(visit_dt)
              list(PANGEA_ID=PANGEA_ID[z], PREFIX=PREFIX[z], STUDY_ID=STUDY_ID[z], visit_dt=visit_dt[z])
            }
            ,by='AID']

# tmp2 <- ak[AID %in% res[! AID %in% tmp1$AID, AID], grep('RK-', PT_ID, value=T)]

tmp1 <- unique(tmp1[, .(AID, STUDY_ID, PANGEA_ID, visit_dt)])
tmp <- merge(tmp1, res, by='AID', all.y=T)

# check whether it makes sense that we have some AID without PANGEA_IDs byt with predictions:
bflocs <- merge(unique(phsc.input[, .(PREFIX,RENAME_ID)]), bflocs, all.y=T, by='PREFIX' ) 
bflocs[, `:=` (AID=gsub('-fq[0-9]$', '', RENAME_ID), RENAME_ID=NULL) ]
# Ok we do not have PANGEA_IDs for individuals that do not have an associate BF file OR for MRC individuals
tmp1 <- tmp[ is.na(PANGEA_ID), unique(AID) ]
bflocs[, all(tmp1) %in% AID]
bflocs[AID %in% tmp1][, table(EXISTS)]
tmp1 <- bflocs[AID %in% tmp1 & EXISTS==TRUE, AID]
ak[AID %in% tmp1, table(grepl('RK',PT_ID))]


date <- format(Sys.Date(), '%y%m%d')
name <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', paste0('TSI_estimates_',date,'.csv'))

write.table(tmp, name, row.name = F, sep=',')

