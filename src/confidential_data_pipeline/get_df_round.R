library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library("haven")
library(here)

gitdir <- here()
source(file.path(gitdir, "config.R"))

c(  file.community.keys,
    file.path.hiv,
    file.path.quest ,
    file.path.hiv.614 ,
    file.path.flow.614 ) |> file.exists() |> all() |> stopifnot()

# load files
community.keys <- fread(file.community.keys)


################################

# COMBINE DATASETS ACROSS MULTIPLE ROUNDS

################################

#
# Quest

# load datasets round 14 only
flow.14<-as.data.table(read_dta(file.path.flow.614))
flow.14 <- flow.14[, .(round, study_id, ageyrs, sex, comm_num, locdate)]
setnames(flow.14, 'locdate', 'intdate')
flow.14 <- flow.14[!round %in% paste0('R0', 15:18)]
flow.14[, intdate := as.Date(intdate, format = '%d/%m/%Y')]

# load datasets ROUND 15 TO 18
quest <- as.data.table(read.csv(file.path.quest.1518))
quest<- quest[, .(round, study_id, ageyrs, sex, comm_num, intdate)]
quest[, intdate := as.Date(intdate, format = '%d-%B-%y')]
quest <- rbind(flow.14, quest)


#
# HIV

# load datasets round 14 only
hiv.14<-as.data.table(read_dta(file.path.hiv.614))
hiv.14 <- hiv.14[, .(study_id, round, hiv, intdate)]
setnames(hiv.14, 'intdate', 'hivdate')
hiv.14 <- hiv.14[!round %in% paste0('R0', 15:18)]
hiv.14[, hivdate := as.Date(hivdate)]

# load datasets ROUND 15 TO 18
hiv <- as.data.table(read.csv(file.path.hiv.1518))
hiv <- hiv[, .(study_id, round, hiv, hivdate)]
hiv[, hivdate := as.Date(hivdate, format = '%d-%B-%y')]
hiv <- rbind(hiv.14, hiv)
hiv[, round := gsub(' ', '', round)] # remove space in string


#
# Find start and end data round
#

df <- merge(hiv, quest, by = c('study_id', 'round'))
cat('Interview date = hiv test date for ', round(nrow(df[intdate == hivdate]) / nrow(df) * 100, 2), '% of the participants')

df <- df[intdate == hivdate] # remove participants who did not had their interview at the same time as the hiv test 

# fishing
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
df_round_fishing <- df[comm_num %in% community.keys[comm =='fishing', COMM_NUM_RAW], 
                       list(min_sample_date = min(intdate), max_sample_date = max(intdate)), by='round']
df_round_fishing <- df_round_fishing[order(round)]
df_round_fishing[, COMM := 'fishing']
df_round_fishing

# inland
df_round_inland <- df[!comm_num %in% community.keys[comm =='fishing', COMM_NUM_RAW], 
                      list(min_sample_date = min(intdate), max_sample_date = max(intdate)), by='round']
# exclude round 15s
df_round_inland <- df_round_inland[round != 'R015S']
df_round_inland <- df_round_inland[order(round)]
df_round_inland[, COMM := 'inland']
df_round_inland

if(0){
  library(ggplot2)
  
  tmp <- rbind(df_round_inland, df_round_fishing)
  ggplot(tmp, aes(y = as.factor(round))) + 
    geom_errorbarh(aes(xmin =min_sample_date , xmax =  max_sample_date, col = as.factor(round))) + 
    facet_grid(COMM~.)
  
}

#
# SAVE ROUND TIMELINE
#

file.name=file.path.round.timeline
if( !file.exists(file.name))
{
  cat('\n Careful: This data should already exist exist in ', file.name  )
  cat('\n check that your Zenodo path is correctly specified in config.R ' )
  cat('\nIf you wish to proceed, and save this file anyway run the commented line below')
  #   save(df_round_inland, df_round_fishing, file=filename, row.names=F)
}else{
  cat('\n Output file', file.name,'already exists.\n')
}
