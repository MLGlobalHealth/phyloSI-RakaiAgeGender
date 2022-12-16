library(data.table)
library(ggplot2)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')
file.community.keys.aggregated <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_id_index.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_221208.csv')

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'census_eligible_count_by_gender_loc_age')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))
quest <- as.data.table(read.csv(file.path.quest))
community.aggregated <- as.data.table(read.csv(file.community.keys.aggregated))

# find community location
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
quest <- merge(quest, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# subset to inland
inlcom <- quest[comm == 'inland']

# find communities by round
inlcom <- inlcom[, list(COMM_NUM = sort(unique(comm_num))), by= 'round']
inlcom <- inlcom[round != 'R015S']

# find communities new id
community.aggregated[grepl('i-', comm_id)][, table(comm_id)]
community.aggregated[comm_id=='i-20']
inlcom <- merge(inlcom, community.aggregated, by.x = 'COMM_NUM', by.y = 'comm_num')
community.aggregated[!comm_num %in% inlcom[, COMM_NUM]]

# find number of community by round
numcom <- inlcom[, list(N_COMM = length(unique(comm_id))), by = 'round'][order(round)]
numcom.total <- inlcom[, list(N_COMM = length(unique(comm_id)))]
numcom.unique <- length(Reduce(intersect, inlcom[, .(list(unique(comm_id))), round]$V1))

# save
saveRDS(list(numcom, numcom.total, numcom.unique), file.path(outdir, 'number_communities_surveyed.rds'))

# plot
df <- data.table(expand.grid(COMM_NUM = inlcom[, sort(unique(COMM_NUM))], 
                             round = inlcom[, sort(unique(round))]))
tmp <- merge(df, inlcom, by = c('COMM_NUM', 'round'), all.x = T)
tmp[, SURVEYED :=T]
tmp[is.na(comm_id), SURVEYED :=F]
ggplot(tmp, aes(x = round, y = as.factor(COMM_NUM))) + 
  geom_raster(aes(fill = SURVEYED)) + 
  theme_bw() + 
  labs(x = 'Round', y = 'Community index')
ggsave('~/Downloads/communities_surveyed.png', w = 6, h = 5)

df <- data.table(expand.grid(comm_id = inlcom[, sort(unique(comm_id))], 
                             round = inlcom[, sort(unique(round))]))
tmp <- merge(df, inlcom, by = c('comm_id', 'round'), all.x = T)
tmp[, SURVEYED :=T]
tmp[is.na(COMM_NUM), SURVEYED :=F]
ggplot(tmp, aes(x = round, y = as.factor(comm_id))) + 
  geom_raster(aes(fill = SURVEYED)) + 
  theme_bw() + 
  labs(x = 'Round', y = 'Community index')
ggsave('~/Downloads/communities_grouped_surveyed.png', w = 6, h = 5)
