# update from Joseph in February 2023 giving us quest data for all communities (not only communities continuously surveyed)
# the hiv data were embedded in the quest data. So this file exctract the necessary information from the quest data
# and create a standalone hiv data

library("data.table")
library("dplyr")
library(here)

gitdir <- here()
source(file.path(gitdir, "config.R"))

file.exists(c(
    file.quest_R09_R14 ))  |> all() |> stopifnot()

# load file
quest_R09_R14 <- fread(file.quest_R09_R14)

# select variable from hiv dataset
cols <- c('study_id', 'round', 'hivdate', 'int_date', 'hiv', 'revertor', 'id', 'lastnegv', 'lastnegvd', 'firstpos_diagnosis_vis', 'firstpos_diagnosis_dt')
hiv_R09_R14 <- select(quest_R09_R14, cols)
hiv_R09_R14 <- unique(hiv_R09_R14)

# find interview date  
hiv_R09_R14[, int_date := as.Date(int_date, format = '%d/%m/%Y')]
hiv_R09_R14[hivdate == '', hivdate := as.character(int_date)]
hiv_R09_R14 <- select(hiv_R09_R14, -c('int_date'))

# save
if(! file.exists(file.hiv_R09_R14))
{
    cat('\n Saving ', file.hiv_R09_R14 ,'...\n')
    write.csv(hiv_R09_R14, file = file.hiv_R09_R14, row.names = F)
}else{
    cat('\n Output file ', file.hiv_R09_R14 ,'already exists\n')
}
cat('\n Done \n')



