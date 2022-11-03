library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library("haven")

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'HIV_R6_R18_220909.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_220909.csv')

# load files
community.keys <- as.data.table(read.csv(file.community.keys))
quest <- as.data.table(read.csv(file.path.quest))
hiv <- as.data.table(read.csv(file.path.hiv))


#################################

# HIV TESTS USING HIV DATA SET #

#################################

if(0){ # check percentage with hiv tests
  
  for(Round in hiv[, sort(unique(round))]){
    hiv_n <- hiv[round == Round, length(unique(study_id))]
    participant_n <- quest[round == gsub(' ', '', Round), length(unique(study_id))]
    cat('There is ', participant_n, 'participants in round', Round, ', ')
    cat(hiv_n, 'of them have an hiv test result (', round(hiv_n  / participant_n, 4), '%)\n')
  }
  
}

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

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# find HIV prevalence rate for participant
rprev <- hivs[, list(COUNT = sum(HIV == 'P'),
                     TOTAL_COUNT = length(HIV)), by = c('ROUND', 'SEX', 'COMM', 'AGEYRS')]

# plot
if(1){
  
  # tmp <- merge(rprev, df_round[, .(ROUND, LABEL_ROUND, COMM)], by = c('ROUND', 'COMM'))
  tmp <- copy(rprev)
  tmp[, Negative := TOTAL_COUNT - COUNT] 
  setnames(tmp, 'COUNT', 'Positive')
  tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'COMM', 'SEX', 'AGEYRS', 'TOTAL_COUNT')) #,'LABEL_ROUND'
  tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
  tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
  tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
  tmp[, SEX_LABEL := 'Female']
  tmp[SEX== 'M', SEX_LABEL := 'Male']
  tmp[, COMM_LABEL := 'Fishing\n communities']
  tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
  
  # plot
  p <- ggplot(tmp[!ROUND %in% c("06", "07", "08", "09")], aes(x = AGEYRS, y = value)) +
    geom_bar(aes(fill = variable), stat = 'identity') + 
    labs(x = 'Age', y = 'Count participants', fill = 'HIV status') +
    facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1))) + 
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  p
  ggsave(p, file=file.path(outdir, paste0('count_participants_by_gender_loc_age.png')), w=8, h=10)
  
}



########################

# FIND SMOOTH PREVENCE #

########################

setnames(rprev, c('COMM','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
rprev[, LOC:= as.integer(LOC_LABEL=='fishing')]
rprev[, SEX:= as.integer(SEX_LABEL=='M')]
rprev[, AGE:= AGE_LABEL-14L]
rprev[, ROW_ID:= seq_len(nrow(rprev))]

# find empirical proportions
rprev[, EMPIRICAL_PREVALENCE := COUNT / TOTAL_COUNT, by = c('ROUND', 'LOC', 'SEX', 'AGE')]# prevalence

# find smooth proportion
for(round in c('R010', 'R011', 'R012', 'R013', 'R014', "R015", "R016", "R017", "R018", "R015S")){
  round <- 'R010'
  
  DT <- copy(rprev[ROUND == round] )
  
  stopifnot(length(round) == 1)
  cat('Fitting stan model for round ', round, '\n')
  
  # account for unobserved entries
  tmp <- data.table(expand.grid(LOC = c(0,1), SEX = c(0,1), AGE_LABEL = rprev[, sort(unique(AGE_LABEL))]))
  DT <- merge(DT, tmp, by = c('LOC', 'SEX', 'AGE_LABEL'), all.y = T)
  DT[is.na(COUNT), COUNT := 0]
  DT[is.na(TOTAL_COUNT), TOTAL_COUNT := 0]
  DT <- DT[order(SEX, LOC, AGE_LABEL)]
  
  # predicts age 
  x_predict <- seq(rprev[, min(AGE_LABEL)], rprev[, max(AGE_LABEL)+1], 0.5)
  
  # make stan data
  stan.data <- list(
    x_predict = x_predict,
    y_observed_00 = DT[SEX==0 & LOC==0, COUNT],
    y_observed_10 = DT[SEX==1 & LOC==0, COUNT],
    y_observed_01 = DT[SEX==0 & LOC==1, COUNT],
    y_observed_11 = DT[SEX==1 & LOC==1, COUNT],
    total_observed_00 = DT[SEX==0 & LOC==0, TOTAL_COUNT],
    total_observed_10 = DT[SEX==1 & LOC==0, TOTAL_COUNT],
    total_observed_01 = DT[SEX==0 & LOC==1, TOTAL_COUNT],
    total_observed_11 = DT[SEX==1 & LOC==1, TOTAL_COUNT],
    alpha_hyper_par_00 = 2,
    alpha_hyper_par_10 = 2,
    alpha_hyper_par_01 = 2,
    alpha_hyper_par_11 = 2
  )
  stan.data$N_predict <- length(stan.data$x_predict)
  stan.data$observed_idx <- which(stan.data$x_predict%%1==0.5)
  stan.data$N_observed <- length(stan.data$observed_idx)
  stan.data$rho_hyper_par_00 <- diff(range(stan.data$x_predict))/3
  stan.data$rho_hyper_par_10 <- diff(range(stan.data$x_predict))/3
  stan.data$rho_hyper_par_01 <- diff(range(stan.data$x_predict))/3
  stan.data$rho_hyper_par_11 <- diff(range(stan.data$x_predict))/3
  
  # load stan model
  stan.model <- stan_model(path.stan, model_name='gp_all')	
  
  # run and save model
  fit <- sampling(stan.model, data=stan.data, iter=10e3, warmup=5e2, chains=1, control = list(max_treedepth= 15, adapt_delta= 0.999))
  filename <- paste0('hivprevalence_gp_stanfit_round',gsub('R0', '', round),'_220909.rds')
  saveRDS(fit, file=file.path(outdir,filename))
  # fit <- readRDS(file.path(outdir,filename))
  
}

# load results 
rounds <- c(10:15, '15S', 16:18)
nsinf <- vector(mode = 'list', length = length(rounds))
nsinf.samples <- vector(mode = 'list', length = length(rounds))
for(i in seq_along(rounds)){
  
  round <- rounds[i]
  DT <- copy(rprev[ROUND == paste0('R0', round)] )
  
  # predicts age 
  x_predict <- seq(rprev[, min(AGE_LABEL)], rprev[, max(AGE_LABEL)+1], 0.5)
  
  # load samples
  filename <- paste0('hivprevalence_gp_stanfit_round',round,'_220909.rds')
  fit <- readRDS(file.path(outdir,filename))
  re <- rstan::extract(fit)
  
  #
  #	summarise prevalence by sex and age
  ps <- c(0.025,0.5,0.975)
  qlab <- c('CL','M','CU')
  tmp <- as.data.table(reshape2::melt(re$p_predict_00))
  tmp[, `:=` (SEX = 0, LOC = 0)]
  tmp1 <- as.data.table(reshape2::melt(re$p_predict_10))
  tmp1[, `:=` (SEX = 1, LOC = 0)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(re$p_predict_01))
  tmp1[, `:=` (SEX = 0, LOC = 1)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(re$p_predict_11))
  tmp1[, `:=` (SEX = 1, LOC = 1)]
  tmp <- rbind(tmp, tmp1)

  tmp[, AGE_LABEL := x_predict[Var2]]

  nsinf.by.age = tmp[, list(q= quantile(value, prob=ps, na.rm = T), q_label=qlab), by=c('SEX', 'LOC', 'AGE_LABEL')]
  nsinf.by.age = as.data.table(reshape2::dcast(nsinf.by.age, ... ~ q_label, value.var = "q"))

  # merge 
  nsinf.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_PREVALENCE)), nsinf.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  nsinf.samples.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_PREVALENCE)), tmp, by=c('SEX','LOC', 'AGE_LABEL'))

  # load change of var name
  set(nsinf.by.age, NULL, 'SEX', NULL)
  set(nsinf.by.age, NULL, 'LOC', NULL)
  setnames(nsinf.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PREVALENCE_M', 'PREVALENCE_CL', 'PREVALENCE_CU'))
  nsinf.by.age[, ROUND := paste0('R0', round)]

  # load change of var name
  set(nsinf.samples.by.age, NULL, 'SEX', NULL)
  set(nsinf.samples.by.age, NULL, 'LOC', NULL)
  set(nsinf.samples.by.age, NULL, 'Var2', NULL)
  setnames(nsinf.samples.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'value'),
           c('COMM', 'SEX', 'AGEYRS', 'PREVALENCE_POSTERIOR_SAMPLE'))
  nsinf.samples.by.age[, ROUND := paste0('R0', round)]

  # keep
  nsinf[[i]] <- nsinf.by.age
  nsinf.samples[[i]] <- nsinf.samples.by.age
}

nsinf <- do.call('rbind', nsinf)
nsinf.samples <- do.call('rbind', nsinf.samples)

# check
stopifnot(nrow(nsinf[COMM == 'inland']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] *nsinf[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(nsinf[COMM == 'fishing']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] *nsinf[COMM == 'fishing', length(unique(ROUND))])


#########

# PLOT #

#########

tmp <- copy(nsinf)
tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

ggplot(tmp, aes(x = AGEYRS)) + 
  geom_point(aes(y = EMPIRICAL_PREVALENCE, col = SEX_LABEL), alpha = 0.5, shape = 16) + 
  geom_line(aes(y = PREVALENCE_M, col = SEX_LABEL)) + 
  geom_ribbon(aes(ymin = PREVALENCE_CL, ymax = PREVALENCE_CU, fill = SEX_LABEL), alpha = 0.7) + 
  facet_grid(ROUND_LABEL~COMM_LABEL) +
  theme_bw() +
  scale_color_manual(values = c('Male'='royalblue3','Female'='deeppink')) +
  scale_fill_manual(values = c('Male'='lightblue3','Female'='lightpink1')) +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)),
        panel.spacing.y = unit(0.7, "lines")) + 
  labs(x = 'Age', y = 'HIV prevalence in RCCS participants', col = '', fill= '')+ 
  scale_y_continuous(labels = scales::percent, limits= c(0,1), expand = c(0,0))
ggsave(file=file.path(outdir, paste0('smooth_prevalence_221031.png')),  w = 8, h = 10)


ggplot(tmp, aes(x = AGEYRS)) + 
  # geom_point(aes(y = EMPIRICAL_PREVALENCE), alpha = 0.5, col = 'darkred') + 
  geom_line(aes(y = PREVALENCE_M, col = ROUND_LABEL)) + 
  # geom_ribbon(aes(ymin = PREVALENCE_CL, ymax = PREVALENCE_CU), alpha = 0.5) + 
  facet_grid(COMM_LABEL~SEX_LABEL) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) + 
  labs(x = 'Age', y = 'Prevalence among participants')

#########

# SAVE #

#########


file.name <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('RCCS_prevalence_estimates_220811.csv'))
write.csv(nsinf, file = file.name, row.names = F)

file.name <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('RCCS_prevalence_posterior_sample_220818.csv'))
write.csv(nsinf.samples, file = file.name, row.names = F)

