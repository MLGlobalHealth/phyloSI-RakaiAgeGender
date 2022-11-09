library(data.table)

# paths
indir.repository <-'~/git/phyloflows'
indir.deepsequence.data <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live'
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live'

# file
path.tests <- file.path(indir.deepsequence.data, 'RCCS_R15_R20',"all_participants_hivstatus_vl_220729.csv")
file.path.quest <- file.path(indir.deepsequence.data, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_220909.csv')

# tuning
VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 1000 # WHO standards


#################

# PREPARE DATA #

#################

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND %in% c(15:18, 15.5)]
# dall <- dall[ROUND == round]

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# keep within census eligible age
DT <- subset(dall, AGEYRS <= 50)


#################################

# KEEP INDIVIDUALS SEEN FOR THE FIRST TIME  
# THAT ARE THE CLOSEST TO NON-PARTICIPANTS

#################################

# keep variable of interest
quest <- as.data.table(read.csv(file.path.quest))
rinc <- quest[, .(round, study_id)]

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# find index of round
rinc <- rinc[order(STUDY_ID, ROUND)]
rinc[, INDEX_ROUND := 1:length(ROUND), by = 'STUDY_ID']

# format round as in DT
rinc[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
rinc[ROUND == '15S', ROUND := '15.1']
rinc[, ROUND := as.numeric(ROUND)]

# merge
DT <- merge(DT, rinc, by= c('STUDY_ID', 'ROUND'))

# keep participants seen for the first time 
DT <- DT[INDEX_ROUND == 1]

# percent of hiv positive with viral laod measurements 
virperc <- DT[HIV_STATUS == 1 & FC == 'inland']
virperc <- virperc[, paste0(round(mean(!is.na(VL_COPIES))*100, 2)), by = 'ROUND'][order(ROUND)]
print(virperc)


#################################

# FIND INDIVIDUALS WITH VIREMIC VIRA; LOADS

#################################

# remove HIV+ individuals with missing VLs  
DT <- subset(DT, HIV_STATUS==0 | HIV_AND_VL==1)

# set ARVMED to 0 for HIV-
set(DT, DT[, which(ARVMED==1 & HIV_STATUS==0)], 'ARVMED', 0) 

# define VLC as VL_COPIES for HIV+ and as 0 for HIV-
set(DT, NULL, 'VLC', DT$VL_COPIES)
set(DT, DT[,which(HIV_STATUS==0)], 'VLC', 0)

# define detectable VL as VLD and undetectable VL as VLU (machine-undetectable)
set(DT, NULL, 'VLU', DT[, as.integer(VLC<VL_DETECTABLE)])
set(DT, NULL, 'VLD', DT[, as.integer(VLC>=VL_DETECTABLE)])

# define suppressed VL as VLS and unsuppressed as VLNS (according to WHO criteria)	
set(DT, NULL, 'VLS', DT[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
set(DT, NULL, 'VLNS', DT[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])

# find HIV+ and detectable
set(DT, NULL, 'HIV_AND_VLD', DT[, as.integer(VLD==1 & HIV_AND_VL==1)])

# reset undetectable to VLC 0
set(DT, DT[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
setkey(DT, ROUND, FC, SEX, AGEYRS)


#################################

# AGGREGATE BY ROUND, SEX, COMM AND AGE  #

#################################

# get count for every categories
tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
tmp1 <- DT[, sort(unique(ROUND))]
vla <- as.data.table(expand.grid(ROUND=tmp1,
                                 FC=c('fishing','inland'),
                                 SEX=c('M','F'),
                                 AGEYRS=tmp))
vla <- vla[, {		
  z <- which(DT$ROUND==ROUND & DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
  list(N          = length(z), # number of participants
       HIV_N      = sum(DT$HIV_STATUS[z]==1), # number of HIV+
       VLNS_N     = sum(DT$VLNS[z]==1), # number of unsuppressed from viral load
       ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z])) # number of unsuppressed from self-reporting
  )				
}, by=names(vla)]

setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
vla[, SEX:= as.integer(SEX_LABEL=='M')]
vla[, AGE:= AGE_LABEL-14L]
vla[, ROW_ID:= seq_len(nrow(vla))]


##########################################

# FIND UNSUPPRESSED PROPORTION CRUDE ESTIMATE #

##########################################

# find empirical proportions
vla[, NONVLNS := HIV_N-VLNS_N]
vla[, EMPIRICAL_NONVLNS_IN_HIV := NONVLNS / HIV_N, by = c('ROUND', 'LOC', 'SEX', 'AGE')]# proportion of suppressed
vla[, EMPIRICAL_VLNS_IN_HIV := 1 - EMPIRICAL_NONVLNS_IN_HIV]# proportion of unsuppressed

##########################################

# SAVE DE-IDENTIFIED DATA #

##########################################

write.csv(vla, file.path(indir.repository, 'data', 'aggregated_newlyregistered_count_unsuppressed.csv'), row.names = F)

