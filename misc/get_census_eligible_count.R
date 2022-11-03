library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library("haven")

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'census_eligible_count_by_gender_loc_age')
file.community.keys <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

# round 14
file.path.flow.614 <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'verif_1.dta')

# round 15 to 18
file.path.flow <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'FlowR15_R18_VoIs_220129.csv')


# load files
flow <- as.data.table(read.csv(file.path.flow))
flow.14<-as.data.table(read_dta(file.path.flow.614))
community.keys <- as.data.table(read.csv(file.community.keys))


#
# combine flow across rounds
# 

# up until round 14
flow.14 <- select(flow.14, c('comm_num', 'locate1', 'locate2', 'resident', 'ageyrs', 'sex', 'round'))
flow.14 <- flow.14[!round %in% paste0('R0', 15:18)]

# round 15 to 18
flow <- select(flow, c('comm_num', 'locate1', 'locate2', 'resident', 'ageyrs', 'sex', 'round'))

# merge
flow <- rbind(flow, flow.14)

#
# Find census eligible count
#

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
flow <- merge(flow, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# Code for ineligibility
flow[, reason_ineligible := NA_character_]

flow[locate1==10 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==2 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==13 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==5 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==10, reason_ineligible := "Out_migrated"]

flow[locate1==7 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==2 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==6 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==3 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==5 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==17 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==88, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==2, reason_ineligible := "Already_seen"]

flow[locate1==11 & locate2==8, reason_ineligible := "Dead"]
flow[locate2==11, reason_ineligible := "Dead"]

flow[resident==0, reason_ineligible := "not_resident"]

flow[ageyrs<15 | ageyrs > 49, reason_ineligible := "Not_within_eligible_age_range"]

flow[is.na(reason_ineligible), reason_ineligible := 'none']

# find count eligible
re <- flow[, list(count = .N), by = c('reason_ineligible', 'round', 'comm', 'ageyrs', 'sex')]
re <- dcast.data.table(re, round + comm + ageyrs + sex ~ reason_ineligible, value.var = 'count')
re[is.na(re)] = 0
re[, ELIGIBLE := round(none + Out_migrated / 2)]
re <- re[ELIGIBLE != 0]

# additional variable
colnames(re) <- toupper(colnames(re))
re[, ROUND := substring(ROUND, 3)]

# find index sex and comm
re <- re[order(ROUND, SEX, COMM, AGEYRS)]
re[, SEX_INDEX := ifelse(SEX == 'M', 1, 0)]
re[, COMM_INDEX := ifelse(COMM == 'fishing', 1, 0)]
  
# find smooth count with loess smooth
rounds <- unique(re$ROUND)
AGEYRSPREDICT <- re[, sort(unique(AGEYRS))]

ncen <- vector(mode = 'list', length = length(rounds))
for(i in seq_along(rounds)){
  
  round <- rounds[i]
  DT <- copy(re[ROUND == round] )
  DT <- DT[order(ROUND, COMM, SEX, AGEYRS)]
  
  # loess
  ncen.by.age <- DT[, {
    loessMod25 <- loess(ELIGIBLE ~ AGEYRS, span=0.25)
    loessMod50 <- loess(ELIGIBLE ~ AGEYRS, span=0.5)
    loessMod75 <- loess(ELIGIBLE ~ AGEYRS, span=0.75)
    
    smoothed25 <- predict(loessMod25, new_data = AGEYRSPREDICT) 
    smoothed50 <- predict(loessMod50, new_data = AGEYRSPREDICT) 
    smoothed75 <- predict(loessMod75, new_data = AGEYRSPREDICT) 
    
    list(AGEYRS = AGEYRSPREDICT, ELIGIBLE_SMOOTH.25 = smoothed25, 
         ELIGIBLE_SMOOTH.50 = smoothed50, ELIGIBLE_SMOOTH.75 = smoothed75, ELIGIBLE = ELIGIBLE)
  }, by = c('COMM', 'SEX')]
  
  ncen.by.age <- merge(ncen.by.age, DT, by=c('SEX','COMM', 'AGEYRS', 'ELIGIBLE'))
  
  # keep
  ncen[[i]] <- ncen.by.age
}
ncen <- do.call('rbind', ncen)

if(0){
  
  tmp <- ncen[, .(COMM, SEX, AGEYRS, ROUND, ELIGIBLE, ELIGIBLE_SMOOTH.25, ELIGIBLE_SMOOTH.50, ELIGIBLE_SMOOTH.75)]
  tmp <- melt.data.table(tmp, id.vars = c('COMM', 'SEX', 'AGEYRS', 'ELIGIBLE', 'ROUND'))
  tmp <- tmp[ROUND != '15S']
  df_label <- tmp[, list(diff = round(abs(sum(ELIGIBLE) - sum(value))), 
                         ylevel = 900 - 3*as.numeric(gsub('.*\\.(.+)', '\\1', variable))), by =  c('COMM', 'SEX', 'variable', 'ROUND')]
  p <- ggplot(tmp, aes(x = AGEYRS)) +
    geom_bar(data = unique(tmp[, .(COMM, SEX, AGEYRS, ELIGIBLE, ROUND)]), aes(y = ELIGIBLE), stat = 'identity', alpha = 0.5) +
    geom_line(aes(y = value, col = variable)) +
    labs(y = 'Census eligible individuals', x = 'Age') +
    facet_grid(ROUND~SEX+COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom') + 
    geom_label(data = df_label, aes(x = 49, y = ylevel, label=diff, col = variable), size = 3, label.size = NA)
  ggsave(p, file = file.path(outdir, 'Smooth_census_eligible_count_all_round.png'), w = 10, h = 10)
  
  tmp <- tmp[ROUND == '15']
  df_label <- df_label[ROUND == '15']
  p <- ggplot(tmp, aes(x = AGEYRS)) +
    geom_bar(data = unique(tmp[, .(COMM, SEX, AGEYRS, ELIGIBLE, ROUND)]), aes(y = ELIGIBLE), stat = 'identity', alpha = 0.5) +
    geom_line(aes(y = value, col = variable)) +
    labs(y = 'Census eligible individuals', x = 'Age') +
    facet_grid(COMM~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom') + 
    geom_label(data = df_label, aes(x = 49, y = ylevel, label=diff, col = variable), size = 3, label.size = NA)
  ggsave(p, file = file.path(outdir, 'Smooth_census_eligible_count_round15.png'), w = 6, h = 6)
}

# choose smoothing 50
ncen[, ELIGIBLE_SMOOTH := ELIGIBLE_SMOOTH.50]
ncen <- select(ncen, -c('ELIGIBLE_SMOOTH.25', 'ELIGIBLE_SMOOTH.50', 'ELIGIBLE_SMOOTH.75'))

# table and plot
tmp <- ncen[, list(count = sum(ELIGIBLE), COUNT_SMOOTH = sum(ELIGIBLE_SMOOTH)), by = c('ROUND', 'SEX', 'COMM')]
tmp[, mean(abs(count - COUNT_SMOOTH))]
knitr::kable(tmp[order(COMM,SEX,ROUND)])

if(1){
  tmp <- copy(ncen)
  tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
  tmp <- tmp[!(ROUND == '15S' & COMM == 'inland')]
  tmp[, SEX_LABEL := 'Female']
  tmp[SEX== 'M', SEX_LABEL := 'Male']
  tmp[, COMM_LABEL := 'Fishing\n communities']
  tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
  
  p <- ggplot(tmp[!ROUND %in% c("06", "07", "08", "09")], aes(x = AGEYRS)) +
    geom_bar(aes(y = ELIGIBLE), stat = 'identity', fill = 'grey70') +
    geom_line(aes(y = ELIGIBLE_SMOOTH), col = 'darkred', alpha  = 0.6) +
    labs(y = 'Census eligible individuals', x = 'Age') +
    facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)))
  p
  ggsave(p, file = file.path(outdir, 'CensusEligibleIndividuals.png'), w = 8, h = 10)

}

ncen[, ELIGIBLE_NOT_SMOOTH := ELIGIBLE]
ncen[, ELIGIBLE := ELIGIBLE_SMOOTH]
ncen <- select(ncen, -'ELIGIBLE_SMOOTH')

# save
write.csv(ncen, file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220830.csv'), row.names = F)
