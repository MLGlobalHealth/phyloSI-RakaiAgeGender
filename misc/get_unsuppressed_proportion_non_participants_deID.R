################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(Hmisc)
library(rstan)

# paths
indir.repository <-'~/git/phyloflows'
indir.repository <- '~/Imperial/phyloflows'

indir.deepsequence.data <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live'
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live'
outdir <- file.path(indir.deepsequence.analyses, 'PANGEA2_RCCS', 'vl_suppofinfected_by_gender_loc_age')

# file
path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')
path.data <- file.path(indir.repository, 'data', 'unsuppressed_proportion_non_participants.csv')

# Load data
vla <- as.data.table( read.csv(path.data) )

##########################################

# PLOT #

##########################################

if(1){
  tmp <- vla[, .(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL, HIV_N, VLNS_N)]
  tmp[, `Non viremic` := HIV_N - VLNS_N] 
  setnames(tmp, 'VLNS_N', 'Viremic')
  tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'HIV_N'))
  setnames(tmp, 'LOC_LABEL', 'COMM')
  setnames(tmp, 'SEX_LABEL', 'SEX')
  tmp[, ROUND := as.character(ROUND)]
  tmp[ROUND == '15.5', ROUND := '15S']
  tmp <- tmp[ROUND != '15S']
  tmp[, ROUND_LABEL := paste0('ROUND: ', ROUND)]
  tmp <- tmp[!(ROUND == '15S')]
  # tmp <- tmp[!(ROUND == '15')]
  tmp[, SEX_LABEL := 'Female']
  tmp[SEX== 'M', SEX_LABEL := 'Male']
  tmp[, COMM_LABEL := 'Fishing\n communities']
  tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
  
  # plot
  p <- ggplot(tmp, aes(x = AGE_LABEL, y = value)) +
    geom_bar(aes(fill = variable), stat = 'identity') + 
    labs(x = 'Age', y = 'Count HIV-positive newly registered participants', fill = 'Viral load') +
    facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1)))
  ggsave(p, file=file.path(outdir, paste0('count_unsuppressed_by_gender_loc_age_newlyregistered_221101.png')), w=9, h=8)
  
}


##########################################

# FIND UNSUPPRESSED PROPORTION SMOOTH ESTIMATE #

##########################################

# find smooth proportion
for(round in 15:18){
  # round <- 15
  # round <- 16
  # round <- 17
  # round <- 18
  
  DT <- copy(vla[ROUND == round] )
  stopifnot(length(round) == 1)
  cat('Fitting stan model for round ', round, '\n')
  
  # predicts age 
  x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
  
  # make stan data
  stan.data <- list(
    x_predict = x_predict,
    y_observed_00 = DT[SEX==0 & LOC==0, HIV_N-VLNS_N],
    y_observed_10 = DT[SEX==1 & LOC==0, HIV_N-VLNS_N],
    y_observed_01 = DT[SEX==0 & LOC==1, HIV_N-VLNS_N],
    y_observed_11 = DT[SEX==1 & LOC==1, HIV_N-VLNS_N],
    total_observed_00 = DT[SEX==0 & LOC==0, HIV_N],
    total_observed_10 = DT[SEX==1 & LOC==0, HIV_N],
    total_observed_01 = DT[SEX==0 & LOC==1, HIV_N],
    total_observed_11 = DT[SEX==1 & LOC==1, HIV_N],
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
  filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'_vl_', VIREMIC_VIRAL_LOAD, '_newlyregistered.rds')
  saveRDS(fit, file=file.path(outdir,filename))
  # fit <- readRDS(file.path(outdir,filename))
  
}


###############

# LOAD RESULTS

###############

rounds <- 15:18
nsinf <- vector(mode = 'list', length = length(rounds))
nsinf.samples <-  vector(mode = 'list', length = length(rounds))
for(i in seq_along(rounds)){
  
  round <- rounds[i]
  DT <- copy(vla[ROUND == round] )
  
  # predicts age 
  x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
  
  # load samples
  filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'_vl_', VIREMIC_VIRAL_LOAD, '_newlyregistered.rds')
  fit <- readRDS(file.path(outdir,filename))
  re <- rstan::extract(fit)
  
  #	summarise proportion of unsuppressed by sex and age
  ps <- c(0.025,0.5,0.975)
  qlab <- c('CL','M','CU')
  tmp <- as.data.table(reshape2::melt(1-re$p_predict_00))
  tmp[, `:=` (SEX = 0, LOC = 0)]
  tmp1 <- as.data.table(reshape2::melt(1-re$p_predict_10))
  tmp1[, `:=` (SEX = 1, LOC = 0)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(1-re$p_predict_01))
  tmp1[, `:=` (SEX = 0, LOC = 1)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(1-re$p_predict_11))
  tmp1[, `:=` (SEX = 1, LOC = 1)]
  tmp <- rbind(tmp, tmp1)
  
  tmp[, AGE_LABEL := x_predict[Var2]]
  
  nsinf.by.age = tmp[, list(q= quantile(value, prob=ps, na.rm = T), q_label=qlab), by=c('SEX', 'LOC', 'AGE_LABEL')]
  nsinf.by.age = as.data.table(reshape2::dcast(nsinf.by.age, ... ~ q_label, value.var = "q"))
  
  nsinf.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_VLNS_IN_HIV)), nsinf.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  nsinf.samples.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_VLNS_IN_HIV)), tmp, by=c('SEX','LOC', 'AGE_LABEL'))
  
  # change of var name
  set(nsinf.by.age, NULL, 'SEX', NULL)
  set(nsinf.by.age, NULL, 'LOC', NULL)
  setnames(nsinf.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'EMPIRICAL_VLNS_IN_HIV', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PROP_UNSUPPRESSED_EMPIRICAL',
             'PROP_UNSUPPRESSED_M', 'PROP_UNSUPPRESSED_CL', 'PROP_UNSUPPRESSED_CU'))
  nsinf.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nsinf.samples.by.age, NULL, 'SEX', NULL)
  set(nsinf.samples.by.age, NULL, 'LOC', NULL)
  set(nsinf.samples.by.age, NULL, 'Var2', NULL)
  setnames(nsinf.samples.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'EMPIRICAL_VLNS_IN_HIV', 'value'),
           c('COMM', 'SEX', 'AGEYRS', 'PROP_UNSUPPRESSED_EMPIRICAL', 'PROP_UNSUPPRESSED_POSTERIOR_SAMPLE'))
  nsinf.samples.by.age[, ROUND := paste0('R0', round)]
  
  # keep
  nsinf[[i]] <- nsinf.by.age
  nsinf.samples[[i]] <- nsinf.samples.by.age
}
nsinf <- do.call('rbind', nsinf)
nsinf.samples <- do.call('rbind', nsinf.samples)


#########

# PLOT #

#########

tmp <- copy(nsinf)
tmp[, ROUND_LABEL := paste0('ROUND: ', gsub('R0(.+)', '\\1', ROUND))]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

# plot
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_line(aes(y = 1-PROP_UNSUPPRESSED_M)) + 
  geom_ribbon(aes(ymin = 1-PROP_UNSUPPRESSED_CL, ymax = 1-PROP_UNSUPPRESSED_CU), alpha = 0.5) + 
  geom_point(aes(y = 1-PROP_UNSUPPRESSED_EMPIRICAL), alpha = 0.5, col = 'darkred') + 
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  labs(x = 'Age', y = 'Proportion of HIV-positive newly registered participants with suppressed viral load', fill = '') +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) + 
  scale_y_continuous(labels = scales::percent, limits= c(0,1))
ggsave(file=file.path(outdir, paste0('smooth_unsuppressed_proportion_newlyregistered.png')), w=9, h=8)


#########

# SAVE #

#########

file.name <- file.path(indir.deepsequence.data, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_vl_', VIREMIC_VIRAL_LOAD, '_newlyregistered_221101.csv'))
write.csv(nsinf, file = file.name, row.names = F)

file.name <- file.path(indir.deepsequence.data, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_posterior_samples_vl_', VIREMIC_VIRAL_LOAD, '_newlyregistered_221101.csv'))
write.csv(nsinf.samples, file = file.name, row.names = F)
