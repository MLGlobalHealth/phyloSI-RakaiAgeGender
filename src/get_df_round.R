library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library("haven")

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

# round 15 to 18
file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'HIV_R15_R18_VOIs_220129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'quest_R15_R18_VoIs_220129.csv')

# round 14
file.path.hiv.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'hivincidence_1.dta')
file.path.flow.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'verif_1.dta')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))

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
quest <- as.data.table(read.csv(file.path.quest))
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
hiv <- as.data.table(read.csv(file.path.hiv))
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

save(df_round_inland, df_round_fishing, file = file.path(indir.repository, 'data', 'RCCS_round_timeline_220905.RData'), row.names = F)


