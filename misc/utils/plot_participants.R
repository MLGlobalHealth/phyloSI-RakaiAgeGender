library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(haven)
library(gridExtra)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

file.seq.count <- file.path(outdir, 'characteristics_sequenced.rds')

file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_220909.csv')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 12:18)), 
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(15, '15S', 16:18))))
tmp <- unique(df_round[, .(ROUND)])[order(ROUND)]
tmp[, INDEX_ROUND := 1:nrow(tmp)]
df_round <- merge(df_round, tmp, by = 'ROUND')
df_round[, ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]


#################################

# GET PARTICIPANTS #

#################################

# load datasets 
quest <- as.data.table(read.csv(file.path.quest))

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate)]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# find count by sex/community/round
rpar <- rinc[, list(PARTICIPANT = .N), by = c('SEX', 'COMM', 'ROUND')]

# keep round of interest
rpar <- merge(rpar, df_round, by = c('COMM', 'ROUND'))

# labels
rpar[, SEX_LABEL := 'Female']
rpar[SEX == 'M', SEX_LABEL := 'Male']
rpar[, COMM_LABEL := 'Inland communities']
rpar[COMM == 'fishing', COMM_LABEL := 'Fishing communities']
rpar[, COMM_LABEL := factor(COMM_LABEL, levels = c('Inland communities', 'Fishing communities'))]



########################################################################

# GET PARTICIPANTS CONSIDERED FOR TRANSMISSION NETWORK RECONSTRUCTION

########################################################################

# load seq count
sequ <- as.data.table(readRDS(file.seq.count))

# keep round of interest
sequ <- merge(sequ, df_round, by = c('COMM', 'ROUND'))

# find count by sex/community/round
sequ <- sequ[TYPE %in% c('Female', "Male")]
sequ[, SEX := substr(TYPE,1,1)]

# labels
sequ[, SEX_LABEL := 'Female']
sequ[SEX == 'M', SEX_LABEL := 'Male']
sequ[, COMM_LABEL := 'Inland communities']
sequ[COMM == 'fishing', COMM_LABEL := 'Fishing communities']
sequ[, COMM_LABEL := factor(COMM_LABEL, levels = c('Inland communities', 'Fishing communities'))]



#########

# PLOT

#########

# participant
p1 <- ggplot(rpar, aes(x = INDEX_ROUND, group = interaction(SEX_LABEL, COMM_LABEL))) + 
  geom_line(aes(y = PARTICIPANT), position=position_dodge(width = 0.3), alpha = 0.5) + 
  geom_point(aes(y = PARTICIPANT, col = SEX_LABEL, shape = COMM_LABEL), size = 2, position=position_dodge(width = 0.3)) + 
  scale_color_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
  scale_shape_manual(values=c('Fishing communities' = 17, 'Inland communities' = 19)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_blank(), 
        legend.position = c(0.17,0.805), 
        legend.key.height = unit(0.44, 'cm'), 
        legend.background=element_blank(),
        # legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.spacing.y= unit(0.00001, 'cm')) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),limits = c(0, NA)) + 
  labs(y = 'Participants') + 
  scale_x_continuous(labels = df_round[order(ROUND), unique(ROUND_LABEL)], breaks = df_round[order(ROUND), unique(INDEX_ROUND)])+
  guides(color = guide_legend(order = 1))

# seqiuence
p2 <- ggplot(sequ, aes(x = INDEX_ROUND, group = interaction(SEX_LABEL, COMM_LABEL))) + 
  geom_line(aes(y = SEQUENCE),  position=position_dodge(width = 0.3), alpha = 0.5) + 
  geom_point(aes(y = SEQUENCE, col = SEX_LABEL, shape = COMM_LABEL), size = 2, position=position_dodge(width = 0.3)) + 
  scale_color_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
  scale_shape_manual(values=c('Fishing communities' = 17, 'Inland communities' = 19)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        axis.title.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        legend.position = 'none', 
        # legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.spacing.y= unit(0.00001, 'cm')) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),limits = c(0, NA)) + 
  labs(y = 'Participants with virus sequenced') + 
  scale_x_continuous(labels = df_round[order(ROUND), unique(ROUND_LABEL)], breaks = df_round[order(ROUND), unique(INDEX_ROUND)])+
  guides(color = guide_legend(order = 1)) + 
  coord_cartesian(xlim = df_round[order(ROUND), range(INDEX_ROUND)]) 

# combine
p <- grid.arrange(p1, p2, layout_matrix = rbind(c(1, 1), c(NA, 2)), widths = c(0.01, 0.99), heights = c(0.465, 0.535))
ggsave(p, file = file.path(outdir, 'participant_sequenced.png'), w = 5.2, h= 5.2)






