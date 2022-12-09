library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library("haven")

# directory of the repository
indir.repository <- '~/git/phyloflows'

# outdir directory for stan fit
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'suppofinfected_by_gender_loc_age')

# files
path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')
data.path <- file.path(indir.repository, 'data', 'aggregated_participants_count_art_coverage.csv')

# find count of participants who reported art use
rart <- as.data.table(read.csv(data.path))


#################################

# PLOT  #

#################################

# plot
if(1){
  
  tmp <- copy(rart)
  tmp[, `Do not use` := TOTAL_COUNT - COUNT] 
  setnames(tmp, 'COUNT', 'Use')
  tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'COMM', 'SEX', 'AGEYRS', 'TOTAL_COUNT'))
  tmp[, variable := factor(variable, levels = c('Use', 'Do not use'))]
  tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
  tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
  tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
  tmp[, SEX_LABEL := 'Women']
  tmp[SEX== 'M', SEX_LABEL := 'Men']
  tmp[, COMM_LABEL := 'Fishing\n communities']
  tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
  tmp <- tmp[AGEYRS > 14 & AGEYRS < 50]
  
  # plot
  p <- ggplot(tmp[!ROUND %in% c("06", "07", "08", "09", "10") & COMM == 'inland'], aes(x = AGEYRS, y = value)) +
    geom_bar(aes(fill = variable), stat = 'identity') + 
    labs(x = 'Age', y = 'Count HIV-positive participants', fill = '') +
    facet_grid(ROUND_LABEL~SEX_LABEL) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_x_continuous(expand = c(0,0))+ 
    scale_fill_manual(values = c('#90B77D', '#425F57'), 
                      labels = c('Reported ART use', 'Did not report ART use')) 
  ggsave(p, file=file.path(outdir, paste0('count_selfreportedart_by_gender_loc_age_221208.pdf')), w=7, h=9)
  
}

########################

# FIND CRUDE PROPORTION #

########################

setnames(rart, c('COMM','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
rart[, LOC:= as.integer(LOC_LABEL=='fishing')]
rart[, SEX:= as.integer(SEX_LABEL=='M')]
rart[, AGE:= AGE_LABEL-14L]
rart[, ROW_ID:= seq_len(nrow(rart))]

# find empirical proportions
rart[, PROP_ART_COVERAGE_EMPIRICAL := COUNT / TOTAL_COUNT, by = c('ROUND', 'LOC', 'SEX', 'AGE')]


########################

# FIND SMOOTH PROPORTION #

########################

# find smooth proportion
for(round in c('R010', 'R011', 'R012', 'R013', 'R014', "R015", "R015S", 'R016', 'R017', 'R018')){
  # round <- 'R018'
  
  DT <- copy(rart[ROUND == round] )
  DT <- DT[order(SEX, LOC, AGE_LABEL)]
  
  stopifnot(length(round) == 1)
  cat('Fitting stan model for round ', round, '\n')
  
  # account for unobserved entries
  tmp <- data.table(expand.grid(LOC = c(0,1), SEX = c(0,1), AGE_LABEL = rart[, sort(unique(AGE_LABEL))])) 
  DT <- merge(tmp, DT, by = c('LOC', 'SEX', 'AGE_LABEL'), all.x = T)
  DT[is.na(COUNT), COUNT := 0]
  DT[is.na(TOTAL_COUNT), TOTAL_COUNT := 0]
  DT <- DT[order(SEX, LOC, AGE_LABEL)]
  
  # predicts age 
  x_predict <- seq(rart[, min(AGE_LABEL)], rart[, max(AGE_LABEL)+1], 0.5)
  
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
  filename <- paste0('art_gp_stanfit_round',gsub('R0', '', round),'_221208.rds')
  saveRDS(fit, file=file.path(outdir,filename))
  # fit <- readRDS(file.path(outdir,filename))
}

# load results 
rounds <- c(10:15, '16', '17', '18')
nsinf <- vector(mode = 'list', length = length(rounds))
nsinf.samples <- vector(mode = 'list', length = length(rounds))
nspred <- vector(mode = 'list', length = length(rounds))
convergence.list <- vector(mode = 'list', length = length(rounds))
for(i in seq_along(rounds)){
  
  round <- rounds[i]
  DT <- copy(rart[ROUND == paste0('R0', round)] )
  
  # account for unobserved age entries but not loc
  tmp <- data.table(expand.grid(LOC = DT[, unique(LOC)], SEX = c(0,1), AGE_LABEL = rart[, sort(unique(AGE_LABEL))], ROUND = DT[, unique(ROUND)])) 
  tmp <- merge(tmp, unique(rart[LOC %in% DT[, unique(LOC)], .(LOC, SEX, SEX_LABEL, LOC_LABEL)]), by = c('LOC', 'SEX'), all.x = T)
  DT <- merge(tmp, DT, by = c('LOC', 'SEX', 'AGE_LABEL', 'SEX_LABEL', 'LOC_LABEL', 'ROUND'), all.x = T)
  DT[is.na(COUNT), COUNT := 0]
  DT[is.na(TOTAL_COUNT), TOTAL_COUNT := 0]
  DT <- DT[order(SEX, LOC, AGE_LABEL)]
  
  # age to predict
  x_predict <- seq(rart[, min(AGE_LABEL)], rart[, max(AGE_LABEL)+1], 0.5)
  
  # load samples
  if(round == '15'){ # change after add of 15s
    filename <- paste0('art_gp_stanfit_round',round,'_221116.rds')
  }else if(as.numeric(round) >= 16){# change after joseph update
    filename <- paste0('art_gp_stanfit_round',round,'_221208.rds')
  }else{
    filename <- paste0('art_gp_stanfit_round',round,'_221101.rds')
  }
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
  #	summarise estimated art use by sex and age
  #
  
  # extract estimates
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
  
  #
  #	summarise predicted art use by sex and age
  #
  
  # merge to total count 
  tmp <- merge(tmp, DT[, .(SEX, LOC, AGE_LABEL, TOTAL_COUNT)], by = c('SEX', 'LOC', 'AGE_LABEL'), all.x = T)
  
  # predict count and predict art use
  tmp[!is.na(TOTAL_COUNT), COUNT_PREDICT := rbinom(1, TOTAL_COUNT, value), by = c('SEX', 'LOC', 'AGE_LABEL', 'iterations')]
  tmp[, ART_USE_PREDICT := COUNT_PREDICT / TOTAL_COUNT]
  
  # summarise
  nspred.by.age = tmp[, list(q= quantile(ART_USE_PREDICT, prob=ps, na.rm = T), q_label=qlab), by=c('SEX', 'LOC', 'AGE_LABEL')]
  nspred.by.age = as.data.table(reshape2::dcast(nspred.by.age, ... ~ q_label, value.var = "q"))
  
  
  #
  # POSTPROCESING
  #
  
  # merge to data
  nsinf.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, PROP_ART_COVERAGE_EMPIRICAL)), nsinf.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  nsinf.samples.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, PROP_ART_COVERAGE_EMPIRICAL)), tmp, by=c('SEX','LOC', 'AGE_LABEL'))
  nspred.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, PROP_ART_COVERAGE_EMPIRICAL)), nspred.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  
  # load change of var name
  set(nsinf.by.age, NULL, 'SEX', NULL)
  set(nsinf.by.age, NULL, 'LOC', NULL)
  setnames(nsinf.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PROP_ART_COVERAGE_M', 'PROP_ART_COVERAGE_CL', 'PROP_ART_COVERAGE_CU'))
  nsinf.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nspred.by.age, NULL, 'SEX', NULL)
  set(nspred.by.age, NULL, 'LOC', NULL)
  setnames(nspred.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PROP_ART_COVERAGE_M', 'PROP_ART_COVERAGE_CL', 'PROP_ART_COVERAGE_CU'))
  nspred.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nsinf.samples.by.age, NULL, 'SEX', NULL)
  set(nsinf.samples.by.age, NULL, 'LOC', NULL)
  set(nsinf.samples.by.age, NULL, 'TOTAL_COUNT', NULL)
  set(nsinf.samples.by.age, NULL, 'COUNT_PREDICT', NULL)
  set(nsinf.samples.by.age, NULL, 'ART_USE_PREDICT', NULL)
  set(nsinf.samples.by.age, NULL, 'Var2', NULL)
  setnames(nsinf.samples.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'value'),
           c('COMM', 'SEX', 'AGEYRS', 'PROP_ART_COVERAGE_POSTERIOR_SAMPLE'))
  
  nsinf.samples.by.age[, ROUND := paste0('R0', round)]
  
  # keep
  nsinf[[i]] <- nsinf.by.age
  nspred[[i]] <- nspred.by.age
  nsinf.samples[[i]] <- nsinf.samples.by.age
  convergence.list[[i]] <- convergence
}
nsinf <- do.call('rbind', nsinf)
nsinf.samples <- do.call('rbind', nsinf.samples)
nspred <- do.call('rbind', nspred)
convergence <- do.call('rbind', convergence.list)

# check all entries are complete
stopifnot(nrow(nsinf[COMM == 'inland']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] * nsinf[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(nsinf[COMM == 'fishing']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] * nsinf[COMM == 'fishing', length(unique(ROUND))])


###########################

# STATISTICS FOR PAPER #

###########################

# get proportion of predicted art use inside credible interval
stats <- list()
tmp <- nspred[COMM == 'inland' & !is.na(PROP_ART_COVERAGE_EMPIRICAL)]
tmp[, within.CI := PROP_ART_COVERAGE_EMPIRICAL >= PROP_ART_COVERAGE_CL & PROP_ART_COVERAGE_EMPIRICAL <= PROP_ART_COVERAGE_CU]
stats[['within.CI']] <- tmp[, paste0(round(mean(within.CI)*100, 2))]

# get lowest rhat and lowest neff
stats[['min_neff']]  = convergence[, round(min(neff))]
stats[['max_rhat']] = convergence[, round(max(rhat), 4)]


#########

# SAVE #

#########


file.name <- file.path(indir.repository, 'fit', paste0('RCCS_art_posterior_samples_221208.rds'))
saveRDS(nsinf.samples, file = file.name)

file.name <- file.path(outdir, paste0('RCCS_art_model_fit_221208.RDS'))
saveRDS(stats, file = file.name)
