################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(Hmisc)
library(rstan)

# paths
indir.repository <-'~/git/phyloflows/misc'
indir.deepsequence.data <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live'
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live'
outdir <- file.path(indir.deepsequence.analyses, 'PANGEA2_RCCS', 'vl_suppofinfected_by_gender_loc_age')

# file
path.stan <- file.path(indir.repository, 'stan_models', 'binomial_gp.stan')
path.tests <- file.path(indir.deepsequence.data, 'RCCS_R15_R20',"all_participants_hivstatus_vl_220729.csv")

# functions
# source( file.path(indir.repository, 'functions', 'get_nonsuppressed_proportion_from_vl_functions.R') )

# tuning
VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 1000 # WHO standards


#################

# PREPARE DATE #

#################

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND %in% c(15:18, 15.5)]
# dall <- dall[ROUND == round]

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# keep within census eligible age
DT <- subset(dall, AGEYRS <= 50)

# remove HIV+ individuals with missing VLs  
DT <- subset(DT, HIV_STATUS==0 | HIV_AND_VL==1)

# set ARVMED to 0 for HIV-
set(DT, DT[, which(ARVMED==1 & HIV_STATUS==0)], 'ARVMED', 0) 

# define VLC as VL_COPIES for HIV+ and as 0 for HIV-
set(DT, NULL, 'VLC', DT$VL_COPIES)
set(DT, DT[,which(HIV_STATUS==0)], 'VLC', 0)

# define detectable VL as VLD and undetectable VL as VLU (machine-undetectable)
set(DT, NULL, 'VLU', DT[, as.integer(VLC<VL_DETECTABLE)])
set(DT, NULL, 'VLD', DT[, as.integer(VLC>=VL_DETECTABLE)])

# define suppressed VL as VLS and unsuppressed as VLNS (according to WHO criteria)	
set(DT, NULL, 'VLS', DT[, as.integer(VLC<VIREMIC_VIRAL_LOAD)])
set(DT, NULL, 'VLNS', DT[, as.integer(VLC>=VIREMIC_VIRAL_LOAD)])

# find HIV+ and detectable
set(DT, NULL, 'HIV_AND_VLD', DT[, as.integer(VLD==1 & HIV_AND_VL==1)])

# reset undetectable to VLC 0
set(DT, DT[, which(HIV_AND_VL==1 & VLU==1)], 'VLC', 0)
setkey(DT, ROUND, FC, SEX, AGEYRS)

# get count for every categories
tmp <- seq.int(min(DT$AGEYRS), max(DT$AGEYRS))
tmp1 <- DT[, sort(unique(ROUND))]
vla <- as.data.table(expand.grid(ROUND=tmp1,
                                 FC=c('fishing','inland'),
                                 SEX=c('M','F'),
                                 AGEYRS=tmp))
vla <- vla[, {		
  z <- which(DT$ROUND==ROUND & DT$FC==FC & DT$SEX==SEX & DT$AGEYRS==AGEYRS)	
  list(N          = length(z), # number of participants
       HIV_N      = sum(DT$HIV_STATUS[z]==1), # number of HIV+
       VLNS_N     = sum(DT$VLNS[z]==1), # number of unsuppressed from viral load
       ARV_N      = sum(DT$ARVMED[z]==0 & DT$HIV_STATUS[z]==1 & !is.na(DT$ARVMED[z])) # number of unsuppressed from self-reporting
  )				
}, by=names(vla)]

setnames(vla, c('FC','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
vla[, LOC:= as.integer(LOC_LABEL=='fishing')]
vla[, SEX:= as.integer(SEX_LABEL=='M')]
vla[, AGE:= AGE_LABEL-14L]
vla[, ROW_ID:= seq_len(nrow(vla))]



##########################################

# FIND UNSUPPRESSED PROPORTION ESTIMATE #

##########################################

# find empirical proportions
vla[, NONVLNS := HIV_N-VLNS_N]
vla[, EMPIRICAL_NONVLNS_IN_HIV := NONVLNS / HIV_N, by = c('ROUND', 'LOC', 'SEX', 'AGE')]# proportion of suppressed
vla[, EMPIRICAL_VLNS_IN_HIV := 1 - EMPIRICAL_NONVLNS_IN_HIV]# proportion of unsuppressed

if(0){
  tmp <- vla[, .(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL, HIV_N, VLNS_N)]
  tmp[, Suppressed := HIV_N - VLNS_N] 
  setnames(tmp, 'VLNS_N', 'Unsuppressed')
  tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'HIV_N'))
  
  # plot
  p <- ggplot(tmp, aes(x = AGE_LABEL, y = value)) +
    geom_bar(aes(fill = variable), stat = 'identity') + 
    labs(x = 'age at visit (years)', y = 'Count HIV+ individuals') +
    theme_bw() + 
    facet_grid(ROUND~SEX_LABEL+LOC_LABEL) 
  ggsave(p, file=file.path(outdir, paste0('count_infected_by_gender_loc_age.png')), w=10, h=9)
  
}

# find smooth proportion
for(round in 15:18){
  round <- 15
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
  filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'_vl_', VIREMIC_VIRAL_LOAD, '.rds')
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
  filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'_vl_', VIREMIC_VIRAL_LOAD, '.rds')
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

ggplot(nsinf, aes(x = AGEYRS)) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU), alpha = 0.5) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL), alpha = 0.5, col = 'darkred') + 
  facet_grid(ROUND~COMM+SEX) + 
  theme_bw()

ggplot(nsinf, aes(x = AGEYRS)) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M, col = ROUND)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU, fill = ROUND), alpha = 0.5) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL, col = ROUND), alpha = 0.5) + 
  facet_grid(COMM~SEX) + 
  theme_bw()

#########

# SAVE #

#########

file.name <- file.path(indir.deepsequence.data, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_vl_', VIREMIC_VIRAL_LOAD, '_220803.csv'))
write.csv(nsinf, file = file.name, row.names = F)

file.name <- file.path(indir.deepsequence.data, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_posterior_samples_vl_', VIREMIC_VIRAL_LOAD, '_220818.csv'))
write.csv(nsinf.samples, file = file.name, row.names = F)
