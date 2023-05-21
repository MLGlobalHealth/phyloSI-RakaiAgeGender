library(data.table)
library(ggplot2)

# directory of the repository
gitdir <- here()

# load file paths
source(file.path(gitdir, 'config.R'))

# outdir directory for stan fit
outdir <- file.path("../phyloSI-RakaiAgeGender-outputs","census_eligible_count_by_gender_loc_age")
if(usr == 'melodiemonod'){
  outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'census_eligible_count_by_gender_loc_age')
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

file.exists(c(
  file.community.keys ,
  file.community.keys.aggregated,
  file.path.quest ))  |> all() |> stopifnot()

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

inlcom[, list(N_COMM = uniqueN(COMM_NUM))]

# find number of community by round
# do merge close communities
numcom <- inlcom[, list(N_COMM = uniqueN(comm_id)), by = 'round'][order(round)]
numcom.total <- inlcom[, list(N_COMM = length(unique(comm_id)))]
numcom.unique <- length(Reduce(intersect, inlcom[, .(list(unique(comm_id))), round]$V1))

# save
file.name <- file.path(outdir, 'number_communities_surveyed.rds')
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
  saveRDS(list(numcom, numcom.total, numcom.unique), file.name)
}else{
  cat("File:", file.name, "already exists...\n")
}


# Do not merge close communities
numcom2 <- inlcom[, list(N_COMM = uniqueN(COMM_NUM)), by = 'round'][order(round)]
numcom.total2 <- inlcom[, list(N_COMM = uniqueN(COMM_NUM))]
numcom.unique2 <- length(Reduce(intersect, inlcom[, .(list(unique(COMM_NUM))), round]$V1))

# save
file.name <-  file.path(outdir, 'number_communities_surveyed_nomerging.rds')
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
  saveRDS(list(numcom2, numcom.total2, numcom.unique2), file.name)
}else{
  cat("File:", file.name, "already exists...\n")
}

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
# ggsave('~/Downloads/communities_surveyed.png', w = 6, h = 5)

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

# save
file.name <- file.path(outdir, 'communities_surveyed_by_round.rds')
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
  saveRDS(tmp, file.name)
}else{
  cat("File:", file.name, "already exists...\n")
}

