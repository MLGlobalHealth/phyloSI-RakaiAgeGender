library(data.table)
indir.deepsequencedata <- '~/OneDrive - Imperial College London/PANGEA/ratmann_pangea_deepsequencedata/live'

file.path.hiv <- file.path(file.path.hiv)
file.path.quest <- file.path(file.path.quest)
file.community.keys <- file.path(file.community.keys)

# load files
hiv <- as.data.table(read.csv(file.path.hiv))
quest <- as.data.table(read.csv(file.path.quest))
community.keys <- as.data.table(read.csv(file.community.keys))

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate)]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# find participant
part <- rinc[, list(PARTICIPANT = .N,  TYPE = 'Total'), by = c('COMM', 'ROUND')]
part <- rbind(part, rinc[SEX == 'F', list(PARTICIPANT = .N,  TYPE = 'Female'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M', list(PARTICIPANT = .N,  TYPE = 'Male'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'F' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'F' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rinc[SEX == 'M' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 10:18)), 
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(15, '15S', 16:18))))
part <- merge(part, df_round, by = c('COMM', 'ROUND'))


########################################

# GET HIV-POSITIVE AMONG PARTICIPANTS #

########################################

# load datasets 
hiv <- as.data.table(read.csv(file.path.hiv))
hiv[, round := gsub(' ', '', round)] # remove space in string

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# find HIV prevalence rate for participant
rprev <- hivs[, list(COUNT = sum(HIV == 'P')), by = c('ROUND', 'SEX', 'COMM', 'AGEYRS')]

# get hiv table
hivp <- rprev[, list(HIV = sum(COUNT),  TYPE = 'Total'), by = c('COMM', 'ROUND')]
hivp <- rbind(hivp, rprev[SEX == 'F', list(HIV = sum(COUNT),  TYPE = 'Female'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M', list(HIV = sum(COUNT),  TYPE = 'Male'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'F' & AGEYRS < 25, list(HIV = sum(COUNT),  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(HIV = sum(COUNT),  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'F' & AGEYRS > 34, list(HIV = sum(COUNT),  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M' & AGEYRS < 25, list(HIV = sum(COUNT),  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(HIV = sum(COUNT),  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
hivp <- rbind(hivp, rprev[SEX == 'M' & AGEYRS > 34, list(HIV = sum(COUNT),  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
hivp <- merge(hivp, df_round, by = c('COMM', 'ROUND'))


######################################################

# GET HIV-POSITIVE AND ART NAIVE AMONG PARTICIPANTS #

######################################################

# keep variable of interest
sart <- quest[, .(ageyrs, round, study_id, sex, comm_num, arvmed)]

# find hiv status
sart <- merge(sart, hiv, by = c('round', 'study_id'))

# find  community
sart <- merge(sart, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(sart) <- toupper(colnames(sart))

# restric age
#sart <- sart[AGEYRS > 14 & AGEYRS < 50]

# keep HIV positive
sart <- sart[HIV == 'P']

# get ART status
sart[, ART := ARVMED ==1]
sart[is.na(ARVMED), ART := F]

# find participant
sartp <- sart[, list(SELF_REPORTED_ART = sum(ART == F),  TYPE = 'Total'), by = c('COMM', 'ROUND')]
sartp <- rbind(sartp, sart[SEX == 'F', list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M', list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'F' & AGEYRS < 25, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'F' & AGEYRS > 34, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M' & AGEYRS < 25, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
sartp <- rbind(sartp, sart[SEX == 'M' & AGEYRS > 34, list(SELF_REPORTED_ART = sum(ART== F),  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
sartp <- merge(sartp, df_round, by = c('COMM', 'ROUND'))

########################

# MAKE TABLE

########################

tab <- merge(part, hivp, by = c('TYPE', 'COMM', 'ROUND'))
tab <- merge(tab, sartp, by = c('TYPE', 'COMM', 'ROUND'))

# tab <- tab[COMM =='inland' & TYPE != 'Total' & TYPE != 'Female' & TYPE != 'Male' ]
# tab[,SEX:= sapply(strsplit(TYPE,','), function(x)x[1])]
# tab[,AGE:= sapply(strsplit(TYPE,','), function(x)x[2])]

########################

# HIV SEQUENCES

########################
# Col 3: at least one HIV deep sequence obtained with Illumina HiSeq / Bonsall protocol.
# Col 4: if not as in the former, at least one HIV deep sequence obtained with Illumina HiSeq / Gall protocol. 
# Col 5: if not as in the former, at least one HIV deep sequence obtained with Illumina MiSeq / Gall protocol. 
# map individuals and pangea ids 
psmap <-  data.table(read.csv(path.sdates.rccs))
psmap <- unique(subset(psmap, select = c('pt_id','pangea_id')))


# load sequences
ds <- data.table(readRDS(file.path(indir.deepsequencedata,'PANGEA2_RCCS/200422_PANGEA2_RCCS_selected_samples.rds')))
ds[grep('PANGEA1_hiseq',REMAP_BAM), TYPE := 'hiseq']
ds[grep('PANGEA1_miseq',REMAP_BAM), TYPE := 'miseq']
ds[is.na(TYPE),TYPE := 'bonsall']
table(ds$TYPE)
ds <- unique(subset(ds, select = c('PANGEA_ID','TYPE')))
ds <- merge(ds, psmap, by.x='PANGEA_ID', by.y = 'pangea_id', all.x = T)


# classify type of sequences
dcl <- ds[,list(bonsall = any(TYPE == 'bonsall')), by='pt_id']
tmp <- ds[,list(hiseq = any(TYPE == 'hiseq')), by='pt_id']
dcl <- merge(dcl, tmp, by='pt_id', all = T)
tmp <- ds[,list(miseq = any(TYPE == 'miseq')), by='pt_id']
dcl <- merge(dcl, tmp, by='pt_id', all = T)
dcl <- dcl[!is.na(pt_id)]
dcl[bonsall == T, hiseq := NA]
dcl[hiseq == T | bonsall == T, miseq := NA]
dcl[,pt_id := gsub('RK-','', pt_id)]
dcl[bonsall == T]
dcl[hiseq == T]
dcl[miseq == T]

nrow(dcl[bonsall == T]) + nrow(dcl[hiseq == T]) + nrow(dcl[miseq == T])
nrow(dcl)
# map to the demographics of hiv positive individuals
dclq <- merge(sart, dcl, by.x = 'STUDY_ID', by.y = 'pt_id', all.x = T)

load(infile.seq.criteria)
dcts <- dct[,list(HREAD=any(V1)),by='pt_id']
dcts[,pt_id:= gsub('RK-','',pt_id)]

# map to the demographics of hiv positive individuals
dclq <- merge(dclq, dcts, by.x = 'STUDY_ID', by.y = 'pt_id', all.x = T)
dclq <- dclq[AGEYRS > 14 & AGEYRS < 50]

# remove neuro data
neuro.metadata <- as.data.table(read.csv(file.path.neuro.metadata))
dclq <- dclq[!STUDY_ID %in% neuro.metadata[,(studyid)]]

# set round to 15 if inland 15S
dclq[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
dclq[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
dclq <- dclq[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]
set(dclq, NULL, 'PARTICIPATED_TO_ROUND_RO15', NULL)

# remove any individuals negative at R15
tmp <- subset(hivs,select=c('STUDY_ID','ROUND','HIV'))
setnames(tmp,'HIV','HIV_check')
dclq <- merge(dclq,tmp,by=c('STUDY_ID','ROUND'),all.x=T)
dclq <- subset(dclq,HIV_check=='P' | is.na(HIV_check))

dclq <- merge(dclq,df_round,by=c('COMM','ROUND'))

# summarize by age, sex, round
dclqs <- dclq[, list(bonsall = sum(bonsall & HREAD, na.rm = T),
                     hiseq = sum(hiseq & HREAD, na.rm = T),
                     miseq = sum(miseq & HREAD, na.rm = T),
                     TYPE = 'Total'), by = c('COMM', 'ROUND')]

dclqs[ROUND == 'R010',]
dclqs <- rbind(dclqs, dclq[SEX == 'F', list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                            hiseq = sum(hiseq & HREAD, na.rm = T),
                                            miseq = sum(miseq & HREAD, na.rm = T),
                                            TYPE = 'Female'), by = c('COMM', 'ROUND')])
dclqs <- rbind(dclqs, dclq[SEX == 'M', list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                            hiseq = sum(hiseq & HREAD, na.rm = T),
                                            miseq = sum(miseq & HREAD, na.rm = T),
                                            TYPE = 'Male'), by = c('COMM', 'ROUND')])
dclqs <- rbind(dclqs, dclq[SEX == 'F' & AGEYRS < 25, list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                                          hiseq = sum(hiseq & HREAD, na.rm = T),
                                                          miseq = sum(miseq & HREAD, na.rm = T),
                                                          TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
dclqs <- rbind(dclqs, dclq[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                                                        hiseq = sum(hiseq & HREAD, na.rm = T),
                                                                        miseq = sum(miseq & HREAD, na.rm = T),
                                                                        TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
dclqs <- rbind(dclqs, dclq[SEX == 'F' & AGEYRS > 34, list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                                          hiseq = sum(hiseq & HREAD, na.rm = T),
                                                          miseq = sum(miseq & HREAD, na.rm = T),
                                                          TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
dclqs <- rbind(dclqs, dclq[SEX == 'M' & AGEYRS < 25, list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                                          hiseq = sum(hiseq & HREAD, na.rm = T),
                                                          miseq = sum(miseq & HREAD, na.rm = T),
                                                          TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
dclqs <- rbind(dclqs, dclq[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                                                        hiseq = sum(hiseq & HREAD, na.rm = T),
                                                                        miseq = sum(miseq & HREAD, na.rm = T),
                                                                        TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
dclqs <- rbind(dclqs, dclq[SEX == 'M' & AGEYRS > 34, list(bonsall = sum(bonsall & HREAD, na.rm = T),
                                                          hiseq = sum(hiseq & HREAD, na.rm = T),
                                                          miseq = sum(miseq & HREAD, na.rm = T),
                                                          TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])
dclqs

# bind to table
tab <- merge(tab, dclqs, by = c('COMM', 'ROUND', 'TYPE'))
tab[, HREAD_COUNT := bonsall + hiseq + miseq]

tab <- tab[COMM=='inland']
setkey(tab, ROUND, TYPE)
tab[,TYPE:=factor(TYPE, levels = c('Total','Female','Female, 15-24',
                                   'Female, 25-34','Female, 35-49',
                                   'Male','Male, 15-24', 'Male, 25-34',
                                   'Male, 35-49'))]
tab[,PERC:=round(HREAD_COUNT/HIV * 100,2)]
# tab[,h_nh_read_count/(bonsall + hiseq + miseq)]
deepsq1 <- tab[ROUND %in% paste0('R0',10:14)]
deepsq2 <- tab[ROUND %in% paste0('R0',15:18)]
saveRDS(deepsq1, file = file.path(outdir,'sequence_tab1_update.rds'))
saveRDS(deepsq2, file = file.path(outdir,'sequence_tab2_update.rds'))

