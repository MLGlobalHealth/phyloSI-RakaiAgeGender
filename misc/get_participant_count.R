library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

file.census.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220830.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

# round 15 to 18
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'quest_R15_R18_VoIs_220129.csv')

# round 14
file.path.quest.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'quest_1.dta')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))
ncen <- as.data.table(read.csv(file.census.count))

################################

# COMBINE DATASETS ACROSS MULTIPLE ROUNDS

################################

# load datasets round 14 only
quest.14<-as.data.table(read_dta(file.path.quest.614))
quest.14 <- quest.14[, .(round, study_id, ageyrs, sex, comm_num, intdate)]
quest.14 <- quest.14[!round %in% paste0('R0', 15:18)]
quest.14[, intdate := as.Date(intdate)]

# load datasets ROUND 15 TO 18
quest <- as.data.table(read.csv(file.path.quest))
quest<- quest[, .(round, study_id, ageyrs, sex, comm_num, intdate)]
quest[, intdate := as.Date(intdate, format = '%d-%B-%y')]
quest <- rbind(quest.14, quest)


################################

# FIND PARTICIPATION

################################

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate)]

# Set to date format
rin[, intdate := as.Date(intdate, format = '%d-%b-%y')]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# GET COUNT
rinc <- rinc[, list(PARTICIPANT  = .N), by = c('AGEYRS', 'ROUND', 'SEX', 'COMM')]

# GET PROPORTION OF PARITCIPATION
tmp <- select(ncen, c('AGEYRS', 'ROUND', 'SEX', 'COMM', 'ELIGIBLE_NOT_SMOOTH'))
rinc[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
rpr <- merge(rinc,tmp , by =  c('AGEYRS', 'ROUND', 'SEX', 'COMM'))
rpr[, PARTICIPATION := PARTICIPANT / ELIGIBLE_NOT_SMOOTH]

# PLOT
tmp <- copy(rpr)
tmp[, ROUND_LABEL := paste0('ROUND:', ROUND)]
tmp <- tmp[!(ROUND == '15S' & COMM == 'inland')]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

p <- ggplot(tmp[!ROUND %in% c("06", "07", "08", "09", "10", "11")], aes(x = AGEYRS)) +
  geom_bar(aes(y = PARTICIPANT), stat = 'identity', fill = 'grey60') +
  labs(y = 'Participants', x = 'Age') +
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)))
p
ggsave(p, file = file.path(outdir, 'Participants.png'), w = 8, h = 9)

p <- ggplot(tmp[!ROUND %in% c("06", "07", "08", "09", "10", "11")], aes(x = AGEYRS)) +
  geom_line(aes(y = PARTICIPATION)) +
  labs(y = 'Participation among census eligible individuals', x = 'Age') +
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) 
p
ggsave(p, file = file.path(outdir, 'Participation.png'), w = 8, h = 9)

## save
tmp <- rpr[, list(MEAN = paste0(round(mean(PARTICIPATION)*100))), by = c('SEX', 'COMM')]
saveRDS(tmp, file.path(outdir, 'Participation.rds'))







