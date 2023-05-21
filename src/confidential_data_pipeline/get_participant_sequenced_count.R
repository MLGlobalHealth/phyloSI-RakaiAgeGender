library(data.table)
library(seqinr)
library(dplyr)

gitdir <- here()
source(file.path(gitdir, "config.R"))

c(  infile.sequence,
    infile.ind.rccs,
    infile.ind.mrc ,
    infile.seq.criteria ,
    file.path.metadata,
    file.path.hiv,
    file.path.quest, 
    file.path.neuro.metadata, 
    file.community.keys) |> file.exists() |> all() |> stopifnot()

# load files
community.keys <- fread(file.community.keys)
community.keys[, comm := fifelse(COMM_NUM_A %like% 'f', 'fishing', 'inland')]

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 14:18)),
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(14, 15, '15S', 16:18))))


###############################################

# GET META DATA

###############################################

#
# Meta data
#

meta_data <- fread(file.path.metadata) #additional meta_data

# find age
meta_data[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
meta_data[, AGEYRS := round(lubridate::time_length(difftime(sample_date, date_birth),"years"))]
meta_data[is.na(AGEYRS), AGEYRS := round(lubridate::time_length(difftime(firstposvd, date_birth),"years"))]
meta_data[is.na(AGEYRS)]

# find community
meta_data[, COMM := 'inland']
meta_data[LakeVictoria_FishingCommunity == 'yes', COMM := 'fishing']

# find hiv status
meta_data[, HIV := ifelse(is.na(firstposvd), 'N', 'P')]

# set sample date to hiv test if missing
meta_data[is.na(sample_date) & !is.na(round), sample_date := firstposvd]

# keep variable of interest
meta_data[, round := paste0('R0', round)]
colnames(meta_data) <- toupper(colnames(meta_data))
meta_data <- meta_data[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV, SAMPLE_DATE)]

# set 15.1 to be 15S
meta_data[ROUND == 'R015.1', ROUND := 'R015S']


#
# Quest
#

quest <- as.data.table(read.csv(file.path.quest))

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate)]

# find  community
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# sample date
setnames(rinc, 'INTDATE', 'SAMPLE_DATE')

# format date
rinc[, SAMPLE_DATE := as.character(SAMPLE_DATE)]

# add meta data from Kate
tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], rinc[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
tmp[, SAMPLE_DATE := as.character(SAMPLE_DATE)]

rinc <- rbind(rinc[, .(STUDY_ID, SEX, ROUND, COMM, SAMPLE_DATE, AGEYRS)], 
              tmp[, .(STUDY_ID, SEX, ROUND, COMM, SAMPLE_DATE, AGEYRS)])


#
# HIV
#

hiv <- as.data.table(read.csv(file.path.hiv))
hiv[, round := gsub(' ', '', round)] # remove space in string

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))

# add meta data from Joseph
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# add meta data from Kate
tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], hivs[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
tmp[, SAMPLE_DATE := as.character(SAMPLE_DATE)]
hivs <- rbind(hivs[, .(STUDY_ID, SEX, ROUND, COMM, SAMPLE_DATE, AGEYRS, HIV)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, SAMPLE_DATE, AGEYRS, HIV)])


###############################################

# GET SEQUENCES DATA

###############################################

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

# remove neuro data
neuro.metadata <- as.data.table(read.csv(file.path.neuro.metadata))
dinfo <- dinfo[!pt_id %in% neuro.metadata[, paste0('RK-', studyid)]]

# find meta data
hivs[, pt_id := paste0('RK-', STUDY_ID)]
dm <- merge(dinfo, hivs, by = c('pt_id'))
colnames(dm) <- toupper(colnames(dm))
stopifnot(dm[, length(unique(PT_ID))] == dinfo[, length(unique(pt_id))])

# keep sequences which meet minimum criteria
load(infile.seq.criteria)
dct[,PANGEA_ID:=paste0('RCCS_',PANGEA_ID)]
dm <- merge(dm,dct,by='PANGEA_ID',all.x=T)
dm <- subset(dm,V1=='TRUE')

# keep meta info closer to sample date
dm[, VISIT_DT := as.Date(VISIT_DT)]
dm[, SAMPLE_DATE := as.Date(SAMPLE_DATE)]
dm[, DIFF_DATE := abs(VISIT_DT - SAMPLE_DATE), by = 'PANGEA_ID']
dm[, IS_MIN := DIFF_DATE == min(na.omit(DIFF_DATE)), by = 'PANGEA_ID']
dcount <- dm[IS_MIN == 1]
stopifnot(nrow(dcount) == dm[, length(unique(PANGEA_ID))])
dcount[, table(ROUND, COMM)]

dcount[HIV == 'N'] ## negative but sequenceD?

# set round to 15 if inland 15S
dcount[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'pt_id']
dcount[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
dcount <- dcount[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]
set(dcount, NULL, 'PARTICIPATED_TO_ROUND_RO15', NULL)

# keep round of interest
dcount <- merge(dcount, df_round, by = c('COMM', 'ROUND'))

# create age groups
dcount[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,49),include.lowest=T,right=F,
                     labels=c('15-24','25-34','35-49'))]

# keep necessary variable
dcount <- dcount[, .(COMM, ROUND, PANGEA_ID, PT_ID, STUDY_ID, SEX, AGEYRS, AGEGP, DIFF_DATE)]

# save sequenced id
file.name <- file.characteristics_sequenced_ind_R14_18
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(dcount, file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}


############################################

# FIND UNIQUE NUMBER OF PARTICIPANTS SEQUENCED

############################################


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

# unique participants by rounds  by comm, sex, agegp
sequ <- dcount.ever[, list(SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX','AGEGP', 'ROUND')]
tot1 <- dcount.ever[, list(SEX = 'Total', AGEGP = 'Total', SEQUENCE = length(unique(PT_ID))), by = c('COMM', 'ROUND')]
tot2 <- dcount.ever[, list(AGEGP = 'Total', SEQUENCE = length(unique(PT_ID))), by = c('COMM','SEX', 'ROUND')]
sequ <- rbind(sequ,tot1,tot2)
sequ[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
sequ[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
sequ[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
sequ <- sequ[order(ROUND, COMM,SEX,AGEGP),]

# save
file.name <- file.characteristics_ever_sequenced
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(sequ, file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}

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

# save
file.name <- file.characteristics_sequenced_R14_18
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(sequ, file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}



