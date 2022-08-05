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
path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')
path.tests <- file.path(indir.deepsequence.data, 'RCCS_R15_R20',"all_participants_hivstatus_vl_220729.csv")

# functions
# source( file.path(indir.repository, 'functions', 'get_nonsuppressed_proportion_from_vl_functions.R') )

# tuning
VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 1000 # WHO standards

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND %in% 15:18]
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

# find empirical proportions
vla[, HIV_AND_VLNS := HIV_N-VLNS_N]
vla[, EMPIRICAL_VLNS_IN_HIV := HIV_AND_VLNS / sum(HIV_N), by = c('ROUND', 'LOC', 'SEX', 'AGE')]# proportion of suppressed
vla[, EMPIRICAL_NONVLNS_IN_HIV := 1 - EMPIRICAL_VLNS_IN_HIV]# proportion of unsuppressed

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



if(0){
  # gp

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
    
  # load results 
  rounds <- 15:18
  nsinf <- vector(mode = 'list', length = length(rounds))
  for(i in seq_along(rounds)){
    
    round <- rounds[i]
    DT <- copy(vla[ROUND == round] )
    
    # load samples
    filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'_vl_', VIREMIC_VIRAL_LOAD, '.rds')
    fit <- readRDS(file.path(outdir,filename))
    re <- rstan::extract(fit)
    
    #	summarise
    ps <- c(0.025,0.5,0.975)
    qlab <- c('CL','M','CU')
    tmp <- cbind(apply(1 - re$p_predict_00, 2, quantile, probs=ps),
                 apply(1 - re$p_predict_10, 2, quantile, probs=ps),
                 apply(1 - re$p_predict_01, 2, quantile, probs=ps),
                 apply(1 - re$p_predict_11, 2, quantile, probs=ps))
    rownames(tmp) <- qlab
    tmp <- as.data.table(reshape2::melt(tmp))
    nsinf.by.age <- dcast.data.table(tmp, Var2~Var1, value.var='value')
    tmp <- as.data.table(expand.grid(AGE_LABEL=x_predict, SEX=c(0,1), LOC=c(0,1)))
    nsinf.by.age <- cbind(tmp, nsinf.by.age) 
    nsinf.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_NONVLNS_IN_HIV)), nsinf.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
    
    # load change of var name
    set(nsinf.by.age, NULL, 'SEX', NULL)
    set(nsinf.by.age, NULL, 'LOC', NULL)
    set(nsinf.by.age, NULL, 'Var2', NULL)
    setnames(nsinf.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'EMPIRICAL_NONVLNS_IN_HIV'), c('COMM', 'SEX', 'AGEYRS', 'PROP_NON_SUPPRESSED_EMPIRICAL'))
    nsinf.by.age[, ROUND := paste0('R0', round)]
    
    # keep
    nsinf[[i]] <- nsinf.by.age
  }
  nsinf <- do.call('rbind', nsinf)
  file.name <- file.path(indir.deepsequence.data, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_vl_', VIREMIC_VIRAL_LOAD, '_220803.csv'))
  write.csv(nsinf, file = file.name, row.names = F)
  
  
}

if(0){
  # loess Regression
  DT <- copy(vla)
  DT <- DT[order(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL)]
  
  logit <- function(p) log(p / (1 - p))
  inv_logit <- function(x) 1 / (1 + exp(-x))
  
  DT <- DT[, {
    .EMPIRICAL_NONVLNS_IN_HIV <- EMPIRICAL_NONVLNS_IN_HIV[!is.nan(EMPIRICAL_NONVLNS_IN_HIV)]
    .EMPIRICAL_NONVLNS_IN_HIV[.EMPIRICAL_NONVLNS_IN_HIV == 0] = 0.001
    .EMPIRICAL_NONVLNS_IN_HIV[.EMPIRICAL_NONVLNS_IN_HIV == 1] = 0.999
    .LOGIT_EMPIRICAL_NONVLNS_IN_HIV = logit(.EMPIRICAL_NONVLNS_IN_HIV)
    .AGE <- AGE[!is.nan(EMPIRICAL_NONVLNS_IN_HIV)]
    loessMod50 <- loess(.LOGIT_EMPIRICAL_NONVLNS_IN_HIV ~ .AGE,  span=.5, control=loess.control(surface="direct"))
    loessMod75 <- loess(.LOGIT_EMPIRICAL_NONVLNS_IN_HIV ~ .AGE,  span=.75, control=loess.control(surface="direct"))
    loessMod90 <- loess(.LOGIT_EMPIRICAL_NONVLNS_IN_HIV ~ .AGE,  span=.9, control=loess.control(surface="direct"))
    
    smoothed50 <- inv_logit(predict(loessMod50, newdata = AGE) )
    smoothed75 <- inv_logit(predict(loessMod75, newdata = AGE) )
    smoothed90 <- inv_logit(predict(loessMod90, newdata = AGE) )
    
    list(AGE_LABEL = AGE_LABEL, N = N, HIV_N = HIV_N, VLNS_N = VLNS_N, ARV_N = ARV_N, AGE = AGE, 
         HIV_AND_VLNS = HIV_AND_VLNS, ROW_ID = ROW_ID, EMPIRICAL_VLNS_IN_HIV = EMPIRICAL_VLNS_IN_HIV, 
         EMPIRICAL_NONVLNS_IN_HIV = EMPIRICAL_NONVLNS_IN_HIV, 
         SMOOTH.50 = smoothed50, SMOOTH.75 = smoothed75, SMOOTH.90 = smoothed90)
  }, by = c('ROUND', 'LOC_LABEL', 'SEX_LABEL')]
  
  # load change of var name
  setnames(DT, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'EMPIRICAL_NONVLNS_IN_HIV'), c('COMM', 'SEX', 'AGEYRS', 'PROP_NON_SUPPRESSED_EMPIRICAL'))
  DT[, ROUND := paste0('R0', ROUND)]
  file.name <- file.path(indir.deepsequence.data, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_loess_vl_', VIREMIC_VIRAL_LOAD, '_220803.csv'))
  write.csv(DT, file = file.name, row.names = F)
  
}

