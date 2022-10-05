library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)


indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

file.census.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220830.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_220909.csv')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))
ncen <- as.data.table(read.csv(file.census.count))
quest <- as.data.table(read.csv(file.path.quest))


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

# GET PARTICIPANT
rinc <- rinc[, list(PARTICIPANT  = .N), by = c('AGEYRS', 'ROUND', 'SEX', 'COMM')]

# GET PARTICIPANT SMOOTH with loess smooth
rinc <- rinc[order(COMM, ROUND, SEX, AGEYRS)]
rinc <- rinc[, {
  loessMod50 <- loess(PARTICIPANT ~ AGEYRS, span=0.5)
  smoothed50 <- predict(loessMod50, new_data = AGEYRSPREDICT) 

  list(AGEYRS = AGEYRS, PARTICIPANT_SMOOTH = smoothed50, PARTICIPANT = PARTICIPANT)
}, by = c('COMM', 'SEX', 'ROUND')]

# GET PROPORTION OF PARITCIPATION
tmp <- select(ncen, c('AGEYRS', 'ROUND', 'SEX', 'COMM', 'ELIGIBLE_NOT_SMOOTH', 'ELIGIBLE'))
rinc[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
rpr <- merge(rinc,tmp , by =  c('AGEYRS', 'ROUND', 'SEX', 'COMM'))
rpr[, PARTICIPATION := PARTICIPANT / ELIGIBLE_NOT_SMOOTH]
rpr[, PARTICIPATION_SMOOTH := PARTICIPANT_SMOOTH / ELIGIBLE]
rpr[PARTICIPATION > 1, PARTICIPATION := 1]
rpr[PARTICIPATION_SMOOTH > 1, PARTICIPATION_SMOOTH := 1]

# PLOT
tmp <- copy(rpr)
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp <- tmp[!(ROUND == '15S' & COMM == 'inland')]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

p <- ggplot(tmp[!ROUND %in% c("06", "07", "08", "09")], aes(x = AGEYRS)) +
  geom_bar(aes(y = PARTICIPANT), stat = 'identity', fill = 'grey60') +
  labs(y = 'Participants', x = 'Age') +
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)))
p
ggsave(p, file = file.path(outdir, 'Participants.png'), w = 8, h = 9)

p <- ggplot(tmp[!ROUND %in% c("06", "07", "08", "09")], aes(x = AGEYRS)) +
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

# sum over age
tmp1 <- tmp[, list(PARTICIPANT = sum(PARTICIPANT), ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'SEX_LABEL', 'COMM_LABEL', 'COMM')]
tmp1[, NON_PARTICIPANT := ELIGIBLE - PARTICIPANT]
tmp1 <- melt.data.table(tmp1, id.vars = c('SEX_LABEL', 'COMM_LABEL', 'ROUND', 'COMM'))
tmp1 <- tmp1[variable != 'ELIGIBLE']
tmp1 <- tmp1[!ROUND %in% c("06", "07", "08", "09")]
tmp1[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp1[, VARIABLE_LABEL := ifelse(variable == 'PARTICIPANT', 'Participant', 'Non-participant')]
tmp1[, COMM_LABEL := gsub('\n', '', COMM_LABEL)]
df_round <- data.table(ROUND = tmp1[, unique(ROUND)])
df_round[, ROUND_INDEX := 1:length(ROUND)]
tmp1 <- merge(tmp1, df_round, by = 'ROUND')

p <- ggplot(NULL) +
  geom_bar_pattern(data = tmp1[SEX_LABEL == 'Male'], aes(x = as.numeric(ROUND_INDEX) - 0.21, y = value, fill = 'Male', alpha = VARIABLE_LABEL, pattern = COMM_LABEL), stat = 'identity', width=0.4, 
                   color = 'grey60',
                   pattern_fill = "white",
                   pattern_colour = 'grey50',
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_spacing = 0.035,
                   pattern_key_scale_factor = 0.6) +
  geom_bar_pattern(data = tmp1[SEX_LABEL == 'Female'], aes(x = as.numeric(ROUND_INDEX) + 0.21, y = value, fill = 'Female', alpha = VARIABLE_LABEL, pattern = COMM_LABEL), stat = 'identity', width=0.4, 
                   color = 'grey60',
                    pattern_fill = "white",
                   pattern_colour = 'grey50',
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_spacing = 0.035,
                   pattern_key_scale_factor = 0.6) +
  labs(y = 'Census eligible population') +
  scale_fill_manual(values = c('Female' = 'lightpink', 'Male' = 'lightblue'))+
  scale_pattern_manual(values = c('Inland communities' = "none", 'Fishing communities' = "stripe")) +
  scale_alpha_manual(values = c('Participant' = 1, 'Non-participant' = 0.5)) + 
  theme_bw() + 
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)), 
        panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust =1), 
        axis.title.x = element_blank(), 
        legend.title = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, .05))) + 
  scale_x_continuous(breaks = tmp1[, unique(ROUND_INDEX)], label =  tmp1[, unique(ROUND_LABEL)] )  + 
  guides(alpha = guide_legend(byrow = T, nrow = 2, order = 1, override.aes = list(pattern = "none", col = 'white')),
         fill = guide_legend(byrow = T, nrow = 2, order = 2, override.aes = list(pattern = "none", col = 'white')),
         pattern = guide_legend(byrow = T, nrow = 2, override.aes = list(fill = "white", col = 'black'), order = 3))
ggsave(p, file = file.path(outdir, 'Participants_aggregated_age.png'), w = 5, h = 5)

## save
tmp <- rpr[!ROUND %in% c("06", "07", "08", "09"), list(MEAN = paste0(round(mean(PARTICIPATION)*100))), by = c('SEX', 'COMM')]
saveRDS(tmp, file.path(outdir, 'Participation.rds'))

file <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'RCCS_participation_220915.csv')
tmp <- rpr[, .(AGEYRS, ROUND, SEX, COMM, PARTICIPANT, PARTICIPATION_SMOOTH)]
setnames(tmp, 'PARTICIPATION_SMOOTH', 'PARTICIPATION')
write.csv(tmp, file = file, row.names = F)





