library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'census_eligible_count_by_gender_loc_age')

file.path.flow <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'FlowR15_R18_VoIs_220129.csv')
file.community.keys <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

# load files
flow <- as.data.table(read.csv(file.path.flow))
community.keys <- as.data.table(read.csv(file.community.keys))


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

# choose smoothing 15
ncen[, ELIGIBLE_SMOOTH := ELIGIBLE_SMOOTH.50]
ncen <- select(ncen, -c('ELIGIBLE_SMOOTH.25', 'ELIGIBLE_SMOOTH.50', 'ELIGIBLE_SMOOTH.75'))

# table and plot
tmp <- ncen[, list(count = sum(ELIGIBLE), COUNT_SMOOTH = sum(ELIGIBLE_SMOOTH)), by = c('ROUND', 'SEX', 'COMM')]
tmp[, mean(abs(count - COUNT_SMOOTH))]
knitr::kable(tmp[order(COMM,SEX,ROUND)])

if(1){
  
  p <- ggplot(ncen, aes(x = AGEYRS)) +
    geom_bar(aes(y = ELIGIBLE_SMOOTH), stat = 'identity') +
    geom_line(aes(y = ELIGIBLE), col = 'darkred', alpha  = 0.6) +
    labs(y = 'Smooth count census eligible individuals', x = 'Age') +
    facet_grid(ROUND~SEX + COMM, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(p, file = file.path(outdir, 'CensusEligibleIndividuals.png'), w = 10, h = 10)

}

ncen[, ELIGIBLE := ELIGIBLE_SMOOTH]
ncen <- select(ncen, -'ELIGIBLE_SMOOTH')

# save
write.csv(ncen, file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220807.csv'), row.names = F)



# # remove round 15s
# re <- re[ROUND != '15S']
# # find smooth count census eligible with gp
# if(0){
#   
#   # load a 1D smooth gp model
#   # file.stan.model <- file.path(indir.repository, 'misc', 'stan_models', 'poisson_gp.stan')
#   file.stan.model <- file.path(indir.repository, 'misc', 'stan_models', 'negative_binomial_gp.stan')
#   model = rstan::stan_model(file.stan.model)
#   
#   rounds <- unique(re$ROUND)
#   AGEYRS <- re[, sort(unique(AGEYRS))]
#   
#   #find smooth count
#   for(round in rounds){
#     # round = 18
#     
#     DT <- re[ROUND == round]
#     DT <- DT[order(SEX_INDEX, COMM_INDEX,AGEYRS)]
#     
#     stan.data <- list(X = AGEYRS, 
#                       N = length(AGEYRS), 
#                       COUNT_00 = DT[SEX_INDEX == 0 & COMM_INDEX == 0, ELIGIBLE], 
#                       COUNT_01 = DT[SEX_INDEX == 0 & COMM_INDEX == 1, ELIGIBLE],
#                       COUNT_10 = DT[SEX_INDEX == 1 & COMM_INDEX == 0, ELIGIBLE], 
#                       COUNT_11 = DT[SEX_INDEX == 1 & COMM_INDEX == 1, ELIGIBLE])
#     stan.data[['rho_hyper_par_1']] <- 3.172883 
#     stan.data[['rho_hyper_par_2']] <- 8.696633
#     stan.data[['TOTAL_COUNT_00']] = sum(stan.data[['COUNT_00']])
#     stan.data[['TOTAL_COUNT_01']] = sum(stan.data[['COUNT_01']])
#     stan.data[['TOTAL_COUNT_10']] = sum(stan.data[['COUNT_10']])
#     stan.data[['TOTAL_COUNT_11']] = sum(stan.data[['COUNT_11']])
#     
#     options(mc.cores = parallel::detectCores())
#     rstan_options(auto_write = TRUE)
#     fit <- sampling(model, data=stan.data, iter=10e3, warmup=5e2, chains=1, 
#                     control = list(max_treedepth = 15, adapt_delta= 0.999))
#     filename <- paste0( '220807_census_count_gp_stan_round',round,'.rds') 
#     # poisson 220805 with total, 220806 withouttotal
#     # neg binomial 220807  withouttotal
#     saveRDS(fit, file=file.path(outdir,filename))
#     
#   }
#   
#   # load results 
#   ncen <- vector(mode = 'list', length = length(rounds))
#   for(i in seq_along(rounds)){
#     
#     round <- rounds[i]
#     DT <- copy(re[ROUND == round] )
#     
#     # load samples
#     filename <- paste0( '220806_census_count_gp_stan_round',round,'.rds')
#     fit <- readRDS(file.path(outdir,filename))
#     summ <- summary(fit)$summary
#     summ[which(summ[,9]< 100),]
#     samples <- rstan::extract(fit)
#     
#     #	summarise
#     ps <- c(0.025,0.5,0.975)
#     qlab <- c('CL','M','CU')
#     tmp <- cbind(apply(samples$lambda_00, 2, quantile, probs=ps),
#                  apply(samples$lambda_10, 2, quantile, probs=ps),
#                  apply(samples$lambda_01, 2, quantile, probs=ps),
#                  apply(samples$lambda_11, 2, quantile, probs=ps))
#     tmp <- round(tmp)
#     rownames(tmp) <- qlab
#     tmp <- as.data.table(reshape2::melt(tmp))
#     ncen.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
#     tmp <- as.data.table(expand.grid(AGEYRS=AGEYRS, SEX_INDEX=c(0,1), COMM_INDEX=c(0,1)))
#     ncen.by.age <- cbind(tmp, ncen.by.age) 
#     ncen.by.age <- merge(DT, ncen.by.age, by=c('SEX_INDEX','COMM_INDEX', 'AGEYRS'))
#     
#     # set names
#     setnames(ncen.by.age, c('M', "CL", 'CU'), c('ELIGIBLE_SMOOTH', 'ELIGIBLE_SMOOTH_CL', 'ELIGIBLE_SMOOTH_CU'))
#     
#     if(0){
#       summ <- summary(fit)$summary
#       summ[which(summ[, 9] < 100),]
#     }
#     
#     if(0){
#       ggplot(ncen.by.age, aes(x = AGEYRS)) +
#         geom_line(aes(y = ELIGIBLE_SMOOTH)) +
#         geom_ribbon(aes(ymin = ELIGIBLE_SMOOTH_CL, ymax = ELIGIBLE_SMOOTH_CU), alpha = 0.5) +
#         geom_point(aes(y = ELIGIBLE), col = 'darkred') +
#         labs(y = 'Census eligible individuals', x = 'Age') +
#         facet_grid(COMM~SEX, label = 'label_both') +
#         theme_bw() +
#         theme(legend.position = 'bottom')
#       
#       tmp <- ncen.by.age[, list(count = sum(ELIGIBLE), COUNT_SMOOTH = sum(ELIGIBLE_SMOOTH)), by = c('SEX', 'COMM')]
#       knitr::kable(tmp[order(SEX,COMM)])
#     }
#     
#     # keep
#     ncen[[i]] <- ncen.by.age
#   }
#   ncen <- do.call('rbind', ncen)
#   
# }
