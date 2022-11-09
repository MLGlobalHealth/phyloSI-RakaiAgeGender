library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library("haven")

# directory to repository
indir.repository <- "~/git/phyloflows"

# outdir to save stan fit
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'prevalence_by_gender_loc_age')

# files
path.data <- file.path(indir.repository, 'data', 'aggregated_count_hiv_positive.csv')
path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')

# Load count of participants by hiv status
rprev <- as.data.table( read.csv(path.data) )


#################################

# PLOT #

#################################

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

# FIND EMPIRICAL PREVENCE #

########################

setnames(rprev, c('COMM','SEX','AGEYRS'), c('LOC_LABEL','SEX_LABEL','AGE_LABEL'))
rprev[, LOC:= as.integer(LOC_LABEL=='fishing')]
rprev[, SEX:= as.integer(SEX_LABEL=='M')]
rprev[, AGE:= AGE_LABEL-14L]
rprev[, ROW_ID:= seq_len(nrow(rprev))]

# find empirical proportions
rprev[, EMPIRICAL_PREVALENCE := COUNT / TOTAL_COUNT, by = c('ROUND', 'LOC', 'SEX', 'AGE')]# prevalence


########################

# FIND SMOOTH PREVENCE #

########################

# find smooth proportion
for(round in c('R010', 'R011', 'R012', 'R013', 'R014', "R015", "R016", "R017", "R018", "R015S")){
  # round <- 'R010'
  
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
rounds <- c(10:15, 16:18)
nsinf <- vector(mode = 'list', length = length(rounds))
nspred <- vector(mode = 'list', length = length(rounds))
nsinf.samples <- vector(mode = 'list', length = length(rounds))
convergence.list <- vector(mode = 'list', length = length(rounds))
for(i in seq_along(rounds)){
  
  round <- rounds[i]
  DT <- copy(rprev[ROUND == paste0('R0', round)] )
  
  # age to predict
  x_predict <- seq(rprev[, min(AGE_LABEL)], rprev[, max(AGE_LABEL)+1], 0.5)
  
  # load samples
  filename <- paste0('hivprevalence_gp_stanfit_round',round,'_220909.rds')
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
  #	summarise estimated prevalence by sex and age
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
  set(tmp, NULL, 'Var2', NULL)
  
  # summarise
  nsinf.by.age = tmp[, list(q= quantile(value, prob=ps, na.rm = T), q_label=qlab), by=c('SEX', 'LOC', 'AGE_LABEL')]
  nsinf.by.age = as.data.table(reshape2::dcast(nsinf.by.age, ... ~ q_label, value.var = "q"))
  
  
  #
  #	summarise predicted prevalence by sex and age
  #
  
  # merge to total count 
  tmp <- merge(tmp, DT[, .(SEX, LOC, AGE_LABEL, TOTAL_COUNT)], by = c('SEX', 'LOC', 'AGE_LABEL'), all.x = T)
  
  # predict count and predict prevalence
  tmp[!is.na(TOTAL_COUNT), COUNT_PREDICT := rbinom(1, TOTAL_COUNT, value), by = c('SEX', 'LOC', 'AGE_LABEL', 'iterations')]
  tmp[, PREVALENCE_PREDICT := COUNT_PREDICT / TOTAL_COUNT]

  # summarise
  nspred.by.age = tmp[, list(q= quantile(PREVALENCE_PREDICT, prob=ps, na.rm = T), q_label=qlab), by=c('SEX', 'LOC', 'AGE_LABEL')]
  nspred.by.age = as.data.table(reshape2::dcast(nspred.by.age, ... ~ q_label, value.var = "q"))
  
  
  #
  # POSTPROCESING
  #
  
  # merge to data
  nsinf.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_PREVALENCE)), nsinf.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  nsinf.samples.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_PREVALENCE)), tmp, by=c('SEX','LOC', 'AGE_LABEL'))
  nspred.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, EMPIRICAL_PREVALENCE)), nspred.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  
  # load change of var name
  set(nsinf.by.age, NULL, 'SEX', NULL)
  set(nsinf.by.age, NULL, 'LOC', NULL)
  setnames(nsinf.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PREVALENCE_M', 'PREVALENCE_CL', 'PREVALENCE_CU'))
  nsinf.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nspred.by.age, NULL, 'SEX', NULL)
  set(nspred.by.age, NULL, 'LOC', NULL)
  setnames(nspred.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PREVALENCE_M', 'PREVALENCE_CL', 'PREVALENCE_CU'))
  nspred.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nsinf.samples.by.age, NULL, 'SEX', NULL)
  set(nsinf.samples.by.age, NULL, 'LOC', NULL)
  set(nsinf.samples.by.age, NULL, 'COUNT_PREDICT', NULL)
  set(nsinf.samples.by.age, NULL, 'TOTAL_COUNT', NULL)
  set(nsinf.samples.by.age, NULL, 'PREVALENCE_PREDICT', NULL)
  setnames(nsinf.samples.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'value'),
           c('COMM', 'SEX', 'AGEYRS', 'PREVALENCE_POSTERIOR_SAMPLE'))
  nsinf.samples.by.age[, ROUND := paste0('R0', round)]
  
  # keep
  nsinf[[i]] <- nsinf.by.age
  nspred[[i]] <- nspred.by.age
  nsinf.samples[[i]] <- nsinf.samples.by.age
  convergence.list[[i]] <- convergence
}
nsinf <- do.call('rbind', nsinf)
nspred <- do.call('rbind', nspred)
nsinf.samples <- do.call('rbind', nsinf.samples)
convergence <- do.call('rbind', convergence.list)

# check that all entries are complete
stopifnot(nrow(nsinf[COMM == 'inland']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] *nsinf[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(nsinf[COMM == 'fishing']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] *nsinf[COMM == 'fishing', length(unique(ROUND))])


##############################

# PLOT PREVALENCE #

##############################

# PREDICTED PREVALENCE
tmp <- copy(nspred)
tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp[, SEX_LABEL := 'Women']
tmp[SEX== 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

ggplot(tmp[COMM == 'inland'], aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PREVALENCE_CL, ymax = PREVALENCE_CU, fill = SEX_LABEL), alpha = 0.5) + 
  geom_line(aes(y = PREVALENCE_M, col = SEX_LABEL, linetype = 'Prediction')) + 
  geom_point(aes(y = EMPIRICAL_PREVALENCE, col = SEX_LABEL, shape = 'Data'), alpha = 0.7) + 
  facet_grid(ROUND_LABEL~.) +
  theme_bw() +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) +
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)),
        panel.spacing.y = unit(0.7, "lines")) + 
  labs(x = 'Age', y = 'HIV prevalence in RCCS participants', col = '', fill= '', 
       shape = '', linetype = '')+ 
  scale_y_continuous(labels = scales::percent, limits= c(0,1), expand = c(0,0)) +
  scale_x_continuous( expand = c(0,0)) +
  guides(shape = guide_legend(order=1), linetype = guide_legend(order=2), 
         color = guide_legend(order=3), fill = guide_legend(order=3))
ggsave(file=file.path(outdir, paste0('smooth_predicted_prevalence_221101.png')),  w = 5, h = 10)


# ESTIMATED PREVALENCE
tmp <- copy(nsinf)
tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp[, SEX_LABEL := 'Women']
tmp[SEX== 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

ggplot(tmp[COMM == 'inland'], aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PREVALENCE_CL, ymax = PREVALENCE_CU, fill = SEX_LABEL), alpha = 0.7) + 
  geom_line(aes(y = PREVALENCE_M, col = SEX_LABEL, linetype = 'Fit')) + 
  geom_point(aes(y = EMPIRICAL_PREVALENCE, col = SEX_LABEL, shape = 'Data'), alpha = 0.7) + 
  facet_wrap(~ROUND_LABEL) +
  theme_bw() +
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) +
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)),
        panel.spacing.y = unit(0.7, "lines")) + 
  labs(x = 'Age', y = 'HIV prevalence in RCCS participants', col = '', fill= '', 
       shape = '', linetype = '')+ 
  scale_y_continuous(labels = scales::percent, limits= c(0,1), expand = c(0,0)) +
  scale_x_continuous( expand = c(0,0))+ 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))
ggsave(file=file.path(outdir, paste0('smooth_estimated_prevalence_221101.pdf')),  w = 7, h = 7)


###########################

# STATISTICS FOR PAPER #

###########################

# get proportion of predicted prevalence inside credible interval
stats <- list()
tmp <- nspred[COMM == 'inland' & !is.na(EMPIRICAL_PREVALENCE)]
tmp[, within.CI := EMPIRICAL_PREVALENCE >= PREVALENCE_CL & EMPIRICAL_PREVALENCE <= PREVALENCE_CU]
stats[['within.CI']] <- tmp[, paste0(round(mean(within.CI)*100, 2))]

# get lowest rhat and lowest neff
stats[['min_neff']]  = convergence[, round(min(neff))]
stats[['max_rhat']] = convergence[, round(max(rhat), 4)]


#########

# SAVE #

#########


file.name <- file.path(indir.repository, 'fit', paste0('RCCS_prevalence_estimates_220811.csv'))
write.csv(nsinf, file = file.name, row.names = F)

file.name <- file.path(indir.repository, 'fit', paste0('RCCS_prevalence_posterior_sample_220818.rds'))
saveRDS(nsinf.samples, file = file.name)

file.name <- file.path(outdir, paste0('RCCS_prevalence_model_fit_convergence_221101.RDS'))
saveRDS(stats, file = file.name)
