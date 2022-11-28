library(data.table)
library(ggplot2)
library(Hmisc)
library(rstan)

# directory of the repository
indir.repository <-'~/git/phyloflows'

# outdir directory for stan fit
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live'
outdir <- file.path(indir.deepsequence.analyses, 'PANGEA2_RCCS', 'vl_suppofinfected_by_gender_loc_age')

# files
path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')
path.data <- file.path(indir.repository, 'data', 'aggregated_newlyregistered_count_unsuppressed_vl200.csv')

# Load count of newly registered participants with unsuppressed viral loads
vla <- as.data.table( read.csv(path.data) )


##########

# PLOT #

##########

if(1){
  tmp <- vla[, .(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL, HIV_N, VLNS_N)]
  tmp[, `Non viremic` := HIV_N - VLNS_N] 
  setnames(tmp, 'VLNS_N', 'Viremic')
  tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'HIV_N'))
  tmp[, variable := factor(variable, levels= c('Non viremic', 'Viremic'))]
  setnames(tmp, 'LOC_LABEL', 'COMM')
  setnames(tmp, 'SEX_LABEL', 'SEX')
  tmp[, ROUND := as.character(ROUND)]
  tmp[ROUND == '15.5', ROUND := '15S']
  tmp <- tmp[ROUND != '15S']
  tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
  tmp <- tmp[!(ROUND == '15S')]
  # tmp <- tmp[!(ROUND == '15')]
  tmp[, SEX_LABEL := 'Women']
  tmp[SEX== 'M', SEX_LABEL := 'Men']
  tmp[, COMM_LABEL := 'Fishing\n communities']
  tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
  tmp <- tmp[AGE_LABEL > 14 & AGE_LABEL < 50]
  
  # plot
  p <- ggplot(tmp[COMM == 'inland'], aes(x = AGE_LABEL, y = value)) +
    geom_bar(aes(fill = variable), stat = 'identity') + 
    labs(x = 'Age', y = 'Count newly registered HIV-positive participants', fill = '') +
    facet_grid(ROUND_LABEL~SEX_LABEL) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1))) + 
    scale_fill_manual(values = c('#9F73AB', '#432C7A'), 
                      labels = c('Viral load <= 200 copies/mL', 'Viral load > 200 copies/mL')) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_x_continuous(expand = c(0,0))
  ggsave(p, file=file.path(outdir, paste0('count_unsuppressed_by_gender_loc_age_newlyregistered_vl200_221121.pdf')), w=7, h=5.2)
  
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
  filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'_vl_200_newlyregistered.rds')
  saveRDS(fit, file=file.path(outdir,filename))
  # fit <- readRDS(file.path(outdir,filename))
  
}


###############

# LOAD RESULTS

###############

rounds <- 15:18
nsinf <- vector(mode = 'list', length = length(rounds))
nsinf.samples <-  vector(mode = 'list', length = length(rounds))
nspred <- vector(mode = 'list', length = length(rounds))
convergence.list <- vector(mode = 'list', length = length(rounds))
for(i in seq_along(rounds)){
  
  round <- rounds[i]
  DT <- copy(vla[ROUND == round] )
  
  # age to predict
  x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL)+1], 0.5)
  
  # load samples
  filename <- paste0( '220729f_notsuppAmongInfected_gp_stan_round',round,'_vl_200_newlyregistered.rds')
  fit <- readRDS(file.path(outdir,filename))
  re <- rstan::extract(fit)
  
  #
  # Find rhat and neff
  #
  
  sum_fit <- summary(fit)
  neff <- na.omit(sum_fit$summary[, 9])
  rhat <- na.omit(sum_fit$summary[, 10])
  convergence <- data.table(ROUND = round, neff = neff, rhat = rhat)
  
  #
  #	summarise estimated unsuppressed by sex and age
  #
  
  # extract estimates
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
  
  #
  #	summarise predicted unsuppressed by sex and age
  #
  
  # merge to total count 
  tmp <- merge(tmp, DT[, .(SEX, LOC, AGE_LABEL, HIV_N)], by = c('SEX', 'LOC', 'AGE_LABEL'), all.x = T)
  
  # predict count and predict unsuppressed
  tmp[!is.na(HIV_N), COUNT_PREDICT := rbinom(1, HIV_N, value), by = c('SEX', 'LOC', 'AGE_LABEL', 'iterations')]
  tmp[, UNSUPPRESSED_PREDICT := COUNT_PREDICT / HIV_N]
  
  # summarise
  nspred.by.age = tmp[, list(q= quantile(UNSUPPRESSED_PREDICT, prob=ps, na.rm = T), q_label=qlab), by=c('SEX', 'LOC', 'AGE_LABEL')]
  nspred.by.age = as.data.table(reshape2::dcast(nspred.by.age, ... ~ q_label, value.var = "q"))
  
  # sub-sample the last 9500 iterations
  it <- data.table(iterations = tmp[, sort(unique(iterations))])
  it[, iterations_rev := max(iterations):1]
  tmp <- merge(it, tmp, by = 'iterations')
  tmp <- tmp[iterations_rev %in% 1:9500]
  tmp[, iterations := iterations - min(iterations)+ 1]
  set(tmp, NULL, 'iterations_rev', NULL)
  
  #
  # POSTPROCESING
  #
  
  # merge to data
  nsinf.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_VLNS_IN_HIV)), nsinf.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  nsinf.samples.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_VLNS_IN_HIV)), tmp, by=c('SEX','LOC', 'AGE_LABEL'))
  nspred.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_VLNS_IN_HIV)), nspred.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  
  # change of var name
  set(nsinf.by.age, NULL, 'SEX', NULL)
  set(nsinf.by.age, NULL, 'LOC', NULL)
  setnames(nsinf.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'EMPIRICAL_VLNS_IN_HIV', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PROP_UNSUPPRESSED_EMPIRICAL',
             'PROP_UNSUPPRESSED_M', 'PROP_UNSUPPRESSED_CL', 'PROP_UNSUPPRESSED_CU'))
  nsinf.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nspred.by.age, NULL, 'SEX', NULL)
  set(nspred.by.age, NULL, 'LOC', NULL)
  setnames(nspred.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PROP_UNSUPPRESSED_M', 'PROP_UNSUPPRESSED_CL', 'PROP_UNSUPPRESSED_CU'))
  nspred.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nsinf.samples.by.age, NULL, 'SEX', NULL)
  set(nsinf.samples.by.age, NULL, 'LOC', NULL)
  set(nsinf.samples.by.age, NULL, 'HIV_N', NULL)
  set(nsinf.samples.by.age, NULL, 'Var2', NULL)
  set(nsinf.samples.by.age, NULL, 'UNSUPPRESSED_PREDICT', NULL)
  set(nsinf.samples.by.age, NULL, 'COUNT_PREDICT', NULL)
  setnames(nsinf.samples.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'EMPIRICAL_VLNS_IN_HIV', 'value'),
           c('COMM', 'SEX', 'AGEYRS', 'PROP_UNSUPPRESSED_EMPIRICAL', 'PROP_UNSUPPRESSED_POSTERIOR_SAMPLE'))
  nsinf.samples.by.age[, ROUND := paste0('R0', round)]
  
  # keep
  nsinf[[i]] <- nsinf.by.age
  nsinf.samples[[i]] <- nsinf.samples.by.age
  nspred[[i]] <- nspred.by.age
  convergence.list[[i]] <- convergence
}
nsinf <- do.call('rbind', nsinf)
nsinf.samples <- do.call('rbind', nsinf.samples)
nspred <- do.call('rbind', nspred)
convergence <- do.call('rbind', convergence.list)


###########################

# STATISTICS FOR PAPER #

###########################

# get proportion of predicted art use inside credible interval
stats <- list()
tmp <- nspred[COMM == 'inland' & !is.na(EMPIRICAL_VLNS_IN_HIV)]
tmp[, within.CI := EMPIRICAL_VLNS_IN_HIV >= PROP_UNSUPPRESSED_CL & EMPIRICAL_VLNS_IN_HIV <= PROP_UNSUPPRESSED_CU]
stats[['within.CI']] <- tmp[, paste0(round(mean(within.CI)*100, 2))]

# get lowest rhat and lowest neff
stats[['min_neff']]  = convergence[, round(min(neff))]
stats[['max_rhat']] = convergence[, round(max(rhat), 4)]


#########

# SAVE #

#########

file.name <- file.path(indir.repository, 'fit', paste0('RCCS_nonsuppressed_proportion_posterior_samples_vl_200_newlyregistered_221121.rds'))
saveRDS(nsinf.samples, file = file.name)

file.name <- file.path(outdir, paste0('RCCS_nonsuppressed_proportion_model_fit_newlyregistered_vl200_221121.RDS'))
saveRDS(stats, file = file.name)
