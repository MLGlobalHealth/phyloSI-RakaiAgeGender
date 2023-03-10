library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library("haven")
library(here)

usr <- Sys.info()['user']
if(usr == 'andrea')
{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live'
}

indir.repository <- here::here()
indir.deepsequencedata.r151r18 <- file.path(indir.deepsequencedata, "RCCS_R15_R18")

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')

# round 15 to 18 
file.path.hiv <- file.path(indir.deepsequencedata.r151r18, 'HIV_R15_R18_VOIs_220129.csv')

# round 15 to 18 without 16
file.path.quest <- file.path(indir.deepsequencedata.r151r18, 'quest_R15_R18_VoIs_220129.csv')

# round 16
file.path.quest.16 <- file.path(indir.deepsequencedata.r151r18, 'quest_R15_R19_VoIs_Dec072022.csv')

# round < 14
file.path.hiv.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'hivincidence_1.dta')
file.path.quest.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'quest_1.dta')

# Latest update from Joseph concerning dates of infection
file.path.update.first.positive <- file.path(indir.deepsequencedata.r151r18, '221128_requested_updated_serohistory.csv')

file.exists(c(
    indir.deepsequencedata.r151r18,
    file.community.keys ,
    path.stan ,
    file.path.hiv ,
    file.path.quest ,
    file.path.quest.16 ,
    file.path.hiv.614 ,
    file.path.quest.614 ,
    file.path.update.first.positive))  |> all() |> stopifnot()

# load files
community.keys <- fread(file.community.keys)


################################

# QUEST

################################


# load datasets ROUND 15 TO 18 without 16
cols <- c('round', 'study_id', 'ageyrs', 'sex', 'comm_num', 'intdate', 'birthdat', 'arvmed', 'cuarvmed')
quest <- fread(file.path.quest, select=cols)
quest[, intdate := as.Date(intdate, format = '%d-%B-%y')]

# load datasets rond 16 and combine
quest.16 <- as.data.table(read.csv(file.path.quest.16))
quest.16<- quest.16[, .(round, study_id, ageyrs, sex, comm_num, intdate, birthdat, arvmed, cuarvmed)]
quest.16[, intdate := as.Date(intdate, format = '%d-%B-%y')]
quest <- rbind(quest[round != 'R016'], quest.16[round == 'R016'])

# find age in quest (some old age were not correct taking month for year (e.g., H109070))
quest[, birthdat2 := birthdat]
quest[, birthdat2 := as.Date(birthdat2, format = '%d-%B-%y')]

## use age at visit for one individual with incorrect birthdate 
quest[birthdat == '18-Aug-18', birthdat2 := intdate - ageyrs*365]

## year in two digits, so have to change 20 to 19 for crazy years, i.e. 2064, to 1964
quest[!is.na(birthdat2) & birthdat2 > as.Date('2006',format = '%Y'), range(birthdat2)] 
quest[!is.na(birthdat2) & birthdat2 > as.Date('2006',format = '%Y'), birthdat2 := birthdat2 - 100*365]
quest[!is.na(birthdat2) , range(birthdat2)]

# for non-missing birthdate use to detrmine age
quest[!is.na(birthdat2), ageyrs := floor(as.numeric(intdate - birthdat2) / 365)]
quest[, range(ageyrs)]

# remove unecessary column
set(quest, NULL, c('birthdat2', 'birthdat'), NULL)

# load datasets round 10-14 
quest.14<-as.data.table(read_dta(file.path.quest.614))
quest.14 <- quest.14[, .(round, study_id, ageyrs, sex, comm_num, intdate, arvmed, cuarvmed)]
quest.14 <- quest.14[!round %in% paste0('R0', 15:18)]
quest.14[, intdate := as.Date(intdate)]

# merge to 10-14 
quest <- rbind(quest.14, quest)

# `save
file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_221208.csv')
write.csv(quest, file = file.name, row.names = F)



################################

# HIV

################################

# load datasets round 10-14 
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

# add last update on infected by joseph
hiv.update <- as.data.table(read.csv(file.path.update.first.positive))
hiv[study_id %in% hiv.update[, gsub('RK-(.+)', '\\1', study_id)], hiv := 'P'] # set to positive hiv test after first positive test

# `save
file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'HIV_R6_R18_221129.csv')
if(! file.exists(file.name))
{
    cat("\n Saving output file", file.name, "\n")
    fwrite(hiv, file = file.name, row.names = F)
}else{
    cat("\n Output file", file.name, "already exists\n")
}
cat("\n Done \n")
