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

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'suppofinfected_by_gender_loc_age')

file.community.keys <- file.path(indir.deepsequence_analyses,'PANGEA2_RCCS1519_UVRI', 'community_names.csv')

path.stan <- file.path(indir.repository, 'misc', 'stan_models', 'binomial_gp.stan')

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'HIV_R6_R18_220909.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', 'Quest_R6_R18_220909.csv')
path.tests <- file.path(indir.deepsequencedata, 'RCCS_R15_R20',"all_participants_hivstatus_vl_220729.csv")

# load files
community.keys <- as.data.table(read.csv(file.community.keys))
quest <- as.data.table(read.csv(file.path.quest))
hiv <- as.data.table(read.csv(file.path.hiv))



#################################

# FIND SELF-REPORTED ART USE  #

#################################

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate, arvmed, cuarvmed)]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# find index of round
rinc <- rinc[order(STUDY_ID, ROUND)]
rinc[, INDEX_ROUND := 1:length(ROUND), by = 'STUDY_ID']

# restric age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# keep HIV positive
rprev <- hivs[HIV == 'P']

# get ART status
rprev[, ART := ARVMED ==1]
rprev[is.na(ARVMED), ART := F]
rprev[, table(ROUND)]


#################################

# ADD VIRAL LOAD DATA  #

#################################

# for round with suppressed set art to true if indiv is suppressed

# tuning
VL_DETECTABLE = 400
VIREMIC_VIRAL_LOAD = 1000 # WHO standards

# Load data: exclude round 20 as incomplete
dall <- fread(path.tests)
dall <- dall[ROUND %in% c(15:18, 15.5)]

# rename variables according to Oli's old script + remove 1 unknown sex
setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
dall <- dall[! SEX=='']

# remove HIV+ individuals with missing VLs  
DT <- subset(dall, HIV_STATUS==0 | HIV_AND_VL==1)

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

# keep within census eligible age
DT <- subset(DT, AGEYRS <= 50)

# keep infected
DT <- DT[HIV_STATUS ==1]

# get ART status
DT[, ART := ARVMED ==1]
DT[is.na(ARVMED), ART := F]

# merge to self-reported data
tmp <- DT[, .(STUDY_ID, ROUND, SEX, FC, VLNS, AGEYRS)]
tmp[, ROUND := paste0('R0', ROUND)]
setnames(tmp, 'AGEYRS', 'AGEYRS2')
rprev <- merge(rprev, tmp, by.x = c('STUDY_ID', 'ROUND', 'SEX', 'COMM'), by.y = c('STUDY_ID', 'ROUND', 'SEX', 'FC'), all.x = T, all.y = T)

# set ageyrs to the viral load data if available
rprev[!is.na(AGEYRS2), AGEYRS := AGEYRS2]
set(rprev, NULL, 'AGEYRS2', NULL)

# set art to true if viremic viral load
rprev[VLNS == 0, ART := T]
rprev <- rprev[!is.na(ART)]

# remove na vlns for round 16 onwards otherwise it leads to % art < % suppressed
nrow(rprev[ROUND == 'R016' & is.na(VLNS)]) / nrow(rprev[ROUND == 'R016' & !is.na(VLNS)])
nrow(rprev[ROUND == 'R017' & is.na(VLNS)]) / nrow(rprev[ROUND == 'R017' & !is.na(VLNS)])
nrow(rprev[ROUND == 'R018' & is.na(VLNS)]) / nrow(rprev[ROUND == 'R018' & !is.na(VLNS)])
rprev <- rprev[!(ROUND %in% c('R016', 'R017', 'R018') & is.na(VLNS))]


#################################

# KEEP INDIVIDUALS SEEN FOR THE FIRST TIME  
# THAT ARE THE CLOSEST TO NON-PARTICIPANTS

#################################

rprev <- rprev[INDEX_ROUND == 1]


#################################

# AGGREGATE BY ROUND, SEX, COMM AND AGE  #

#################################

# find self reported under art for participant
rart <- rprev[, list(COUNT = sum(ART == T), TOTAL_COUNT = length(ART)), by = c('ROUND', 'SEX', 'COMM', 'AGEYRS')]


#################################

# PLOT  #

#################################

# plot
if(1){
  
  tmp <- copy(rart)
  tmp[, `Do not use` := TOTAL_COUNT - COUNT] 
  setnames(tmp, 'COUNT', 'Use')
  tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'COMM', 'SEX', 'AGEYRS', 'TOTAL_COUNT'))
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
    labs(x = 'Age', y = 'Count newly registered HIV-positive newly registered participants', fill = 'Self-reported ART use') +
    facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
    theme_bw() +
    theme(legend.position = 'bottom', 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = rel(1))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  ggsave(p, file=file.path(outdir, paste0('count_selfreportedart_by_gender_loc_age_newlyregistered_221101.png')), w=8, h=12)
  
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
  round <- 'R018'
  
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
  filename <- paste0('art_gp_stanfit_round',gsub('R0', '', round),'_newlyregistered_221101.rds')
  saveRDS(fit, file=file.path(outdir,filename))
  # fit <- readRDS(file.path(outdir,filename))
}

# load results 
rounds <- c(10:15, '15S', '16', '17', '18')
nsinf <- vector(mode = 'list', length = length(rounds))
nsinf.samples <- vector(mode = 'list', length = length(rounds))
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
  
  # predicts age 
  x_predict <- seq(rart[, min(AGE_LABEL)], rart[, max(AGE_LABEL)+1], 0.5)
  
  # load samples
  filename <- paste0('art_gp_stanfit_round',round,'_newlyregistered_221101.rds')
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
  nsinf.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, PROP_ART_COVERAGE_EMPIRICAL)), nsinf.by.age, by=c('SEX','LOC', 'AGE_LABEL'))
  nsinf.samples.by.age <- merge(subset(DT, select=c(SEX,SEX_LABEL,LOC,LOC_LABEL,AGE_LABEL, PROP_ART_COVERAGE_EMPIRICAL)), tmp, by=c('SEX','LOC', 'AGE_LABEL'))
  
  # load change of var name
  set(nsinf.by.age, NULL, 'SEX', NULL)
  set(nsinf.by.age, NULL, 'LOC', NULL)
  setnames(nsinf.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'M', "CL", "CU"), 
           c('COMM', 'SEX', 'AGEYRS', 'PROP_ART_COVERAGE_M', 'PROP_ART_COVERAGE_CL', 'PROP_ART_COVERAGE_CU'))
  nsinf.by.age[, ROUND := paste0('R0', round)]
  
  # load change of var name
  set(nsinf.samples.by.age, NULL, 'SEX', NULL)
  set(nsinf.samples.by.age, NULL, 'LOC', NULL)
  set(nsinf.samples.by.age, NULL, 'Var2', NULL)
  setnames(nsinf.samples.by.age, c('LOC_LABEL', 'SEX_LABEL', 'AGE_LABEL', 'value'),
           c('COMM', 'SEX', 'AGEYRS', 'PROP_ART_COVERAGE_POSTERIOR_SAMPLE'))
  
  nsinf.samples.by.age[, ROUND := paste0('R0', round)]
  
  # keep
  nsinf[[i]] <- nsinf.by.age
  nsinf.samples[[i]] <- nsinf.samples.by.age
}
nsinf <- do.call('rbind', nsinf)
nsinf.samples <- do.call('rbind', nsinf.samples)

# check all entries are complete
stopifnot(nrow(nsinf[COMM == 'inland']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] * nsinf[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(nsinf[COMM == 'fishing']) == nsinf[, length(unique(AGEYRS))] * nsinf[, length(unique(SEX))] * nsinf[COMM == 'fishing', length(unique(ROUND))])


#########

# PLOT #

#########

tmp <- copy(nsinf)
tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
tmp <- tmp[!(ROUND == 'R015')]
tmp[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
ggplot(tmp, aes(x = AGEYRS)) +
  geom_point(aes(y = PROP_ART_COVERAGE_EMPIRICAL), alpha = 0.5, col = 'darkred') +
  geom_line(aes(y = PROP_ART_COVERAGE_M)) +
  geom_ribbon(aes(ymin = PROP_ART_COVERAGE_CL, ymax = PROP_ART_COVERAGE_CU), alpha = 0.5) +
  facet_grid(ROUND_LABEL~COMM_LABEL + SEX_LABEL) +
  theme_bw() +
  theme(legend.position = 'bottom',
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) +
  labs(x = 'Age', y = 'Self-reported ART coverage among newly registered participants') + 
  scale_y_continuous(labels = scales::percent, limits= c(0,1))
ggsave(file=file.path(outdir, paste0('smooth_artcoverage_newlyregistered_221101.png')), w=8, h=8)



#########

# SAVE #

#########

file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_art_estimates_newlyregistered_221101.csv'))
write.csv(nsinf, file = file.name, row.names = F)

file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_art_posterior_samples_newlyregistered_221101.csv'))
write.csv(nsinf.samples, file = file.name, row.names = F)
