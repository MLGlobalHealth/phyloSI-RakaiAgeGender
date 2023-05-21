library(data.table)
library(here)

gitdir <- here()
source(file.path(gitdir, "config.R"))

c(  path.tests) |> file.exists() |> all() |> stopifnot()

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

# percent of hiv positive with viral laod measurements 
virperc <- dall[HIV_STATUS == 1 & COMM == 'inland']
virperc <- virperc[, paste0(round(mean(!is.na(HIV_VL))*100, 2)), by = 'ROUND']
print(virperc)

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# keep within census eligible age
DT <- subset(dall, AGEYRS <= 50)

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

# SAVE DE-INDENTIFIED DATA #

##########################################

file.name <- path.count.unsupp
if( !file.exists(file.name))
{
  cat('\n Careful: This data should already exist exist in ', file.name  )
  cat('\n check that your Zenodo path is correctly specified in config.R ' )
  cat('\nIf you wish to proceed, and save this file anyway run the commented line below')
  #  write.csv(vla, file , row.names = F)
}else{
  cat('\n Output file', file.name,'already exists.\n')
}

