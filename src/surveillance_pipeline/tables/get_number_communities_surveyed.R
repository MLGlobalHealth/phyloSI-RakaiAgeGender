library(data.table)
library(ggplot2)

usr <- Sys.info()[['user']]
if(usr!='andrea')
{
    indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
    indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
}else{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepsequence_analyses <- '/home/andrea/HPC/project/ratmann_deepseq_analyses/live'
}

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')
file.community.keys.aggregated <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_id_index.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_221208.csv')

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'census_eligible_count_by_gender_loc_age')

# load files
community.keys <- fread(file.community.keys)
quest <- fread(file.path.quest)
community.aggregated <- fread(file.community.keys.aggregated)

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
community.aggregated

inlcom[, list(N_COMM = uniqueN())]

# find number of community by round
# do merge close communities
numcom <- inlcom[, list(N_COMM = uniqueN(comm_id)), by = 'round'][order(round)]
numcom.total <- inlcom[, list(N_COMM = length(unique(comm_id)))]
numcom.unique <- length(Reduce(intersect, inlcom[, .(list(unique(comm_id))), round]$V1))
# save
saveRDS(list(numcom, numcom.total, numcom.unique), file.path(outdir, 'number_communities_surveyed.rds'))

# Do not merge close communities
numcom2 <- inlcom[, list(N_COMM = uniqueN(COMM_NUM)), by = 'round'][order(round)]
numcom.total2 <- inlcom[, list(N_COMM = uniqueN(COMM_NUM))]
numcom.unique2 <- length(Reduce(intersect, inlcom[, .(list(unique(COMM_NUM))), round]$V1))
# save
saveRDS(list(numcom2, numcom.total2, numcom.unique2), file.path(outdir, 'number_communities_surveyed_nomerging.rds'))



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

# make table
df <- data.table(expand.grid(COMM_NUM = inlcom[, sort(unique(COMM_NUM))],
                             round = inlcom[, sort(unique(round))]))
tmp <- merge(df, inlcom, by = c('COMM_NUM', 'round'), all.x = T)
tmp[, SURVEYED :=1]
tmp[is.na(comm_id), SURVEYED :=0]
tmp2 <- unique(subset(inlcom,select=c('COMM_NUM','comm_id')))
tmp <- merge(subset(tmp,select=c('COMM_NUM','round','SURVEYED')),tmp2,by=c('COMM_NUM'),all=T)
# just keep rounds 10-18
dr <- data.table(round=paste0('R0',seq(10,18,1)))
tmp <- merge(tmp,dr,by='round')
# aggregate the communities together
tmp <- dcast(tmp, comm_id~round,value.var = 'SURVEYED',fun.aggregate=sum)
tmp <- melt(tmp, id.vars=c('comm_id'),variable.name='round')
# recode
tmp[value==0, SURVEYED:= 'No']
tmp[value>0, SURVEYED:= 'Yes']
set(tmp,NULL,'value',NULL)
tmp <- dcast(tmp, comm_id~round,value.var = 'SURVEYED')

saveRDS(tmp, file.path(outdir, 'communities_surveyed_by_round.rds'))
