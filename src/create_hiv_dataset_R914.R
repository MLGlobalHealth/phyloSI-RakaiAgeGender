library("data.table")
library("dplyr")

# set up path to data
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'

# file to data
file.quest_R09_R14 = file.path(indir.deepsequencedata, 'RCCS_R9_R14/quest_R09_R14.csv')

# load file
quest_R09_R14 <- as.data.table(read.csv(file.quest_R09_R14))

# select variable from hiv dataset
hiv_R09_R14 <- select(quest_R09_R14, c('study_id', 'round', 'hivdate', 'int_date', 'hiv', 'revertor', 'id', 'lastnegv', 'lastnegvd', 
                                       'firstpos_diagnosis_vis', 'firstpos_diagnosis_dt'))
hiv_R09_R14 <- unique(hiv_R09_R14)

# find interview date  
hiv_R09_R14[, int_date := as.Date(int_date, format = '%d/%m/%Y')]
hiv_R09_R14[hivdate == '', hivdate := as.character(int_date)]
hiv_R09_R14 <- select(hiv_R09_R14, -c('int_date'))

# save
file.hiv_R09_R14 = file.path(indir.deepsequencedata, 'RCCS_R9_R14/HIV_R09_R14.csv')
write.csv(hiv_R09_R14, file = file.hiv_R09_R14, row.names = F)

