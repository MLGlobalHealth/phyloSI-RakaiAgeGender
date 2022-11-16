library(data.table)
library(seqinr)

indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'

infile.sequence <- file.path(indir.deepsequencedata,"200422_PANGEA2_RCCSMRC_alignment.fasta")
infile.ind.rccs <- file.path(indir.deepsequencedata,'PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
infile.ind.mrc <- file.path(indir.deepsequencedata,'PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')

file.path.meta <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_20220329.RData')

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

# load meta data
load(file.path.meta)

# load data
alignment <- read.fasta(file = infile.sequence)
unique(do.call(c, lapply(alignment, unique)))
nsequence <- length(alignment)
npos <- unique(lengths(alignment))

# map alignments to studyid
dinfo <- data.table(pangea_id=names(alignment))
id.dt <- data.table(read.csv(infile.ind.rccs))
id.dt <- subset(id.dt,select = c("pt_id","pangea_id", 'cd4_count', 'visit_dt'))
id.dt[,pangea_id:=paste0('RCCS_',pangea_id)]
tmp <- data.table(read.csv(infile.ind.mrc))
tmp <- subset(tmp,select = c("pt_id","pangea_id", 'cd4_count', 'visit_dt'))
tmp[,pangea_id:=paste0('MRCUVRI_',pangea_id)]
id.dt <- rbind(id.dt,tmp)
id.dt <- unique(id.dt)
dinfo <- merge(dinfo, id.dt, by="pangea_id", all.x=T)
tmp <- dinfo[is.na(pt_id)]
dinfo <- dinfo[!is.na(pt_id)]
cat('No personal information found for ',nrow(dinfo[is.na(pt_id)]), ' sequences \n')

# remove mrc cohort
dinfo <- dinfo[!grepl('MRCUVRI', pangea_id)]

# find meta data
dm <- merge(dinfo, meta_data, by.x = c('pt_id'), by.y = c('study_id'))
colnames(dm) <- toupper(colnames(dm))

# remove neuro data
dm <- dm[ ROUND != 'neuro']

# keep meta info closer to sample date 
dm[, VISIT_DT := as.Date(VISIT_DT)]
dm[, DIFF_DATE := abs(VISIT_DT - SAMPLE_DATE), by = 'PANGEA_ID']
dm[, IS_MIN := DIFF_DATE == min(DIFF_DATE), by = 'PANGEA_ID']
dcount <- dm[IS_MIN == 1]
stopifnot(nrow(dcount) == dm[, length(unique(PANGEA_ID))])
dcount[, table(ROUND, COMM)]

# find age at visit
dcount[, AGEYRS := round(lubridate::time_length(difftime(VISIT_DT, DATE_BIRTH),"years"))] # round to match hivs
stopifnot(nrow(dcount[is.na(AGEYRS)]) == 0)
dcount <- dcount[AGEYRS > 14 & AGEYRS < 50]

# set round to 15 if inland 15S
dcount[COMM == 'inland' & ROUND == 'R015S', ROUND := 'R015']

# find characteristics sequenced id
dcount[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                     labels=c('15-24','25-34','35-49'))]
saveRDS(dcount, file.path(outdir, 'characteristics_sequenced_ind_R14_18.rds'))

# if multiple PANGEA_ID per round, keep the one the closest to visit data
dcount[, IS_MIN := DIFF_DATE == min(DIFF_DATE), by = c('PT_ID', 'ROUND')]
dcount <- dcount[IS_MIN == T]

# select round when first sequenced
dcount[, r:= as.numeric(gsub('R','',gsub('S','.1',ROUND)))]
dcount[, min.r := min(r), by = 'PT_ID']
dcount.ever <- unique(dcount[r == min.r, .(PT_ID, COMM, SEX, AGEYRS, ROUND)])
stopifnot(nrow(dcount.ever) == dcount.ever[, length(unique(PT_ID))])

# AGE GROUPS
dcount.ever[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                          labels=c('15-24','25-34','35-49'))]

############################################

# FIND NUMBER OF PARTICIPANTS EVER SEQUENCED

############################################

# unique participants by rounds  by comm, sex, agegp
sequ <- dcount.ever[, list(SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX','AGEGP', 'ROUND')]
tot1 <- dcount.ever[, list(SEX = 'Total', AGEGP = 'Total', SEQUENCE = length(unique(PT_ID))), by = c('COMM', 'ROUND')]
tot2 <- dcount.ever[, list(AGEGP = 'Total', SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX', 'ROUND')]
sequ <- rbind(sequ,tot1,tot2)
sequ[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
sequ[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
sequ[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
sequ <- sequ[order(ROUND, COMM,SEX,AGEGP),]
saveRDS(sequ, file.path(outdir, 'characteristics_ever_sequenced.rds'))

# unique participants across rounds 14-18 by comm, sex, agegp
sequ <- dcount.ever[ROUND %in% c('R014','R015','R016','R017','R018'),
                    list(SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX','AGEGP')]
tot1 <- dcount.ever[ROUND %in% c('R014','R015','R016','R017','R018'),
                    list(SEX = 'Total', AGEGP = 'Total', SEQUENCE = length(unique(PT_ID))), by = c('COMM')]
tot2 <- dcount.ever[ROUND %in% c('R014','R015','R016','R017','R018'),
                    list(AGEGP = 'Total', SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX')]
sequ <- rbind(sequ,tot1,tot2)
sequ[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
sequ[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
sequ[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
sequ <- sequ[order(COMM,SEX,AGEGP),]

saveRDS(sequ, file.path(outdir, 'characteristics_sequenced_R14_18.rds'))