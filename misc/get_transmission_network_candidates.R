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

# keep visit dt closer to sample date
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

# find count sequenced
sequ <- dcount[, list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Total'), by = c('COMM', 'ROUND')]
sequ <- rbind(sequ, dcount[SEX == 'F', list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Female'), by = c('COMM', 'ROUND')])
sequ <- rbind(sequ, dcount[SEX == 'M', list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Male'), by = c('COMM', 'ROUND')])
sequ <- rbind(sequ, dcount[SEX == 'F' & AGEYRS < 25, list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
sequ <- rbind(sequ, dcount[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
sequ <- rbind(sequ, dcount[SEX == 'F' & AGEYRS > 34, list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
sequ <- rbind(sequ, dcount[SEX == 'M' & AGEYRS < 25, list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
sequ <- rbind(sequ, dcount[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
sequ <- rbind(sequ, dcount[SEX == 'M' & AGEYRS > 34, list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

sequ[, TYPE := factor(TYPE, levels = c('Total', 'Female', 'Female, 15-24', "Female, 25-34", "Female, 35-49", 
                                      "Male",  "Male, 15-24", "Male, 25-34", "Male, 35-49"))]
sequ <- sequ[order(COMM, ROUND, TYPE)]
saveRDS(sequ, file.path(outdir, 'characteristics_sequenced.rds'))

sequ <- dcount[SEX == 'F', list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Female'), by = c('COMM')]
sequ <- rbind(sequ, dcount[SEX == 'M', list(SEQUENCE = length(unique(PT_ID)),  TYPE = 'Male'), by = c('COMM')])
saveRDS(sequ, file.path(outdir, 'characteristics_sequenced_brief.rds'))

# plot
tmp <- dcount[, list(SEQUENCE = length(unique(PT_ID))), by = c('COMM', 'ROUND', 'SEX', 'AGEYRS')]
tmp[, ROUND := gsub('R0(.+)','\\1', ROUND)]
tmp[, ROUND_LABEL := paste0('ROUND:', ROUND)]
tmp <- tmp[!(ROUND == '15S' & COMM == 'inland')]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

p <- ggplot(tmp, aes(x = AGEYRS)) +
  geom_bar(aes(y = SEQUENCE), stat = 'identity', fill = 'grey60') +
  labs(y = 'Count of HIV-positive participants with virus sequenced', x = 'Age') +
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)))
p
ggsave(p, file = file.path(outdir, 'Participants_sequenced_age.png'), w = 8, h = 9)

