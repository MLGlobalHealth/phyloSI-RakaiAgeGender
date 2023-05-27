#################################################################################################
# Rakai data analysis for EMOD inputs
# August 12, 2022
# Adam Akullian/Kate Grabowski/Melodie Monod;
#################################################################################################

library("ggplot2")
library("data.table")
library("dplyr")
library("haven")
library("tidyr")
library("tidyverse")
library("mgcv")
library("zoo")
library("remotes")
library("metR")
library("lubridate")
library(ggpubr)
library(here)
library(optparse)

# set to directory
gitdir <- here()

# source paths 
source(file.path(gitdir, "config.R"))

# make sure file exists
c(  file.path.seroconverter_cohort,
    file.path.seroconverter_cohort.30,
    file.incidence.fits,
    file.incidence.30com.fits) |> file.exists() |> all() |> stopifnot()

# sensitivity analysis?
restrict_to_30_comms <- F

# path where to store results
if(!restrict_to_30_comms){
  outdir <- output.dir.incidence.estimation
}else{
  outdir <- output.dir.incidence.estimation.30comms
}

if(!dir.exists(outdir)) dir.create(outdir, recursive = T)

# load functions
source(file.path(gitdir.R.incid, 'incidence_rate_estimation_functions.R'))

# load data 
if(!restrict_to_30_comms){
  file <- file.path.seroconverter_cohort
}else{
  file <- file.path.seroconverter_cohort.30
}
seroconverter_cohort.list <- readRDS(file)

# load fits
if(!restrict_to_30_comms){
  file.name <- file.incidence.fits
}else{
  file.name <- file.incidence.30com.fits
}
load(file.name)

# utils
rounds_group_1 <- c("R006","R007", "R008", "R009", "R010", "R011", "R012", "R013", "R014")
rounds_group_2 <- c("R015", "R015S", "R016", "R017", "R018")
rounds_group_3 <- c('R019')

# make df round
df_round <- make_df_round(rounds_group_1, rounds_group_2, rounds_group_3)

rounds_numeric_group_1 <- df_round[visit %in% rounds_group_1, round_numeric]
rounds_numeric_group_2 <- df_round[visit %in% rounds_group_2, round_numeric]
rounds_numeric_group_3 <- df_round[visit == rounds_group_3, round_numeric]


###############################################################################################

# SUMMARISE PERSON YEARS, INCIDENCE INFECTION AND INCIDENCE ESTIMATES ACROSS IMPUTED DATA SETS

###############################################################################################

ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

#
# Summarise seroconverter_cohort
#

seroconverter_cohort.all <- rbindlist(seroconverter_cohort.list)
rm(seroconverter_cohort.list)

if(!'age_group' %in% names(seroconverter_cohort.all)){
  df_age_aggregated <- data.table(age = 15:49, age_group = c(rep('15-24', 10),  rep("25-34", 10), rep("35-49", 15)))
  seroconverter_cohort.all <- merge(seroconverter_cohort.all, df_age_aggregated, by = 'age')
}

# pysum and seroconvsum by round, sex and 1-year age band for every iterations
seroconverter_cohort_imputation <- seroconverter_cohort.all[, list(pysum=sum(na.omit(py)), 
                                                                   seroconvsum=sum(hivinc)), 
                                                            by = c('sex', 'round', 'age', 'iterations')]
seroconverter_cohort_imputation[, inc_crude := seroconvsum / pysum]

# crude incidence by round, sex and 1-year age band 
seroconverter_cohort_crude <- seroconverter_cohort_imputation[, list(q= quantile(inc_crude, prob=ps, na.rm = T),q_label=p_labs), by=c('sex', 'round', 'age')]	
seroconverter_cohort_crude <- dcast(seroconverter_cohort_crude, sex + round + age ~ q_label, value.var = "q")
setnames(seroconverter_cohort_crude, c('M'), c('inc_crude'))

# classify number of missing visits
seroconverter_cohort.all$number_missing_visits_status = ifelse(seroconverter_cohort.all$number_missing_visits > 1, '>1', 
                                                               ifelse(seroconverter_cohort.all$number_missing_visits == 1, '1', '0'))
seroconverter_cohort.all$number_missing_visits_status <- factor(seroconverter_cohort.all$number_missing_visits_status, levels = c('0', '1', '>1'))

# by round, sex, 1-year age band and number_missing_visits_status
seroconverter_cohort <- seroconverter_cohort.all[py > 0, list(seroconvsum = sum(na.omit(hivinc)), pysum = sum(na.omit(py))), by=c('sex', 'round', 'age', 'number_missing_visits_status', 'iterations')]
tmp1 = seroconverter_cohort[, list(q= quantile(seroconvsum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_hivinc')), by=c('sex', 'round', 'age', 'number_missing_visits_status')]	
tmp1 = dcast(tmp1, sex + round + age + number_missing_visits_status ~ q_label, value.var = "q")
seroconverter_cohort <- seroconverter_cohort[, list(q= quantile(pysum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_py')), by=c('sex', 'round', 'age', 'number_missing_visits_status')]	
seroconverter_cohort <- dcast(seroconverter_cohort, sex + round + age + number_missing_visits_status ~ q_label, value.var = "q")
seroconverter_cohort <- merge(seroconverter_cohort, tmp1, by = c('sex', 'round', 'age', 'number_missing_visits_status')) 
setnames(seroconverter_cohort, c('M_hivinc', 'M_py'), c('hivinc', 'py'))

# by sex, age groups and round
seroconverter_cohort_agg <- seroconverter_cohort.all[py > 0, list(seroconvsum = sum(na.omit(hivinc)), pysum = sum(na.omit(py))), by=c('sex', 'round', 'age_group', 'iterations')]
tmp1 = seroconverter_cohort_agg[, list(q= quantile(seroconvsum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_hivinc')), by=c('sex', 'round', 'age_group')]	
tmp1 = dcast(tmp1, sex + round + age_group ~ q_label, value.var = "q")
seroconverter_cohort_agg <- seroconverter_cohort_agg[, list(q= quantile(pysum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_py')), by=c('sex', 'round', 'age_group')]	
seroconverter_cohort_agg <- dcast(seroconverter_cohort_agg, sex + round + age_group ~ q_label, value.var = "q")
seroconverter_cohort_agg <- merge(seroconverter_cohort_agg, tmp1, by = c('sex', 'round', 'age_group')) 
setnames(seroconverter_cohort_agg, c('M_hivinc', 'M_py'), c('hivinc', 'py'))

# by sex and round
seroconverter_cohort_agg_s <- seroconverter_cohort.all[py > 0, list(seroconvsum = sum(na.omit(hivinc)), pysum = sum(na.omit(py))), by=c('sex', 'round', 'iterations')]
tmp1 = seroconverter_cohort_agg_s[, list(q= quantile(seroconvsum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_hivinc')), by=c('sex', 'round')]	
tmp1 = dcast(tmp1, sex + round ~ q_label, value.var = "q")
seroconverter_cohort_agg_s <- seroconverter_cohort_agg_s[, list(q= quantile(pysum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_py')), by=c('sex', 'round')]	
seroconverter_cohort_agg_s <- dcast(seroconverter_cohort_agg_s, sex + round ~ q_label, value.var = "q")
seroconverter_cohort_agg_s <- merge(seroconverter_cohort_agg_s, tmp1, by = c('sex', 'round')) 
setnames(seroconverter_cohort_agg_s, c('M_hivinc', 'M_py'), c('hivinc', 'py'))

# by sex and round and missing visits
seroconverter_cohort_agg_sm <- seroconverter_cohort.all[py > 0, list(seroconvsum = sum(na.omit(hivinc)), pysum = sum(na.omit(py))), by=c('sex', 'round', 'number_missing_visits_status', 'iterations')]
tmp1 = seroconverter_cohort_agg_sm[, list(q= quantile(seroconvsum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_hivinc')), by=c('sex', 'round', 'number_missing_visits_status')]	
tmp1 = dcast(tmp1, sex + round + number_missing_visits_status ~ q_label, value.var = "q")
seroconverter_cohort_agg_sm <- seroconverter_cohort_agg_sm[, list(q= quantile(pysum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_py')), by=c('sex', 'round', 'number_missing_visits_status')]	
seroconverter_cohort_agg_sm <- dcast(seroconverter_cohort_agg_sm, sex + round + number_missing_visits_status ~ q_label, value.var = "q")
seroconverter_cohort_agg_sm <- merge(seroconverter_cohort_agg_sm, tmp1, by = c('sex', 'round', 'number_missing_visits_status')) 
setnames(seroconverter_cohort_agg_sm, c('M_hivinc', 'M_py'), c('hivinc', 'py'))

# by round
seroconverter_cohort_agg_r <- seroconverter_cohort.all[py > 0, list(seroconvsum = sum(na.omit(hivinc)), pysum = sum(na.omit(py))), by=c('round', 'iterations')]
tmp1 = seroconverter_cohort_agg_r[, list(q= quantile(seroconvsum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_hivinc')), by=c('round')]	
tmp1 = dcast(tmp1, round ~ q_label, value.var = "q")
seroconverter_cohort_agg_r <- seroconverter_cohort_agg_r[, list(q= quantile(pysum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_py')), by=c('round')]	
seroconverter_cohort_agg_r <- dcast(seroconverter_cohort_agg_r, round ~ q_label, value.var = "q")
seroconverter_cohort_agg_r <- merge(seroconverter_cohort_agg_r, tmp1, by = c('round')) 
setnames(seroconverter_cohort_agg_r, c('M_hivinc', 'M_py'), c('hivinc', 'py'))

# in total
seroconverter_cohort_total <- seroconverter_cohort.all[round %in% c(5:9, rounds_numeric_group_2) & py > 0, list(seroconvsum = sum(na.omit(hivinc)), pysum = sum(na.omit(py))), by=c('iterations')]
tmp1 = seroconverter_cohort_total[, list(q= quantile(seroconvsum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_hivinc'))]	
tmp1 = dcast(tmp1, . ~ q_label, value.var = "q")
seroconverter_cohort_total <- seroconverter_cohort_total[, list(q= quantile(pysum, prob=ps, na.rm = T),q_label=paste0(p_labs, '_py'))]	
seroconverter_cohort_total <- dcast(seroconverter_cohort_total, . ~ q_label, value.var = "q")
seroconverter_cohort_total <- cbind(seroconverter_cohort_total, tmp1) 
setnames(seroconverter_cohort_total, c('M_hivinc', 'M_py'), c('hivinc', 'py'))


#
# summarise incidence 
#

set.seed(12)
modelpreds.age.1218.all <- rbindlist(modelpreds.age.1218.list)
modelpreds.age.1218.all <- modelpreds.age.1218.all[, list(inc = exp(rnorm(1000, fit, se.fit)), 
                                                          iterations_within = 1:1000), by=c('Sex', 'model', 'round', 'round_label', 'ROUND', 'age', 'iterations', 'fit')]

# by round
modelpreds_all <- find_incidence_rate_all(modelpreds.age.1218.all[model == 'model_1'], seroconverter_cohort_imputation)
modelpreds_all = modelpreds_all[, list(q= quantile(inc, prob=ps, na.rm = T),
                                       q_label=paste0(p_labs, '_inc')), by=c('round')]	
modelpreds_all = dcast(modelpreds_all, round ~ q_label, value.var = "q")
setnames(modelpreds_all, c('M_inc', 'CL_inc', 'CU_inc'), c('inc', 'lb', 'ub'))
modelpreds_all[, Sex := 'All']

# by sex and round
modelpreds <- find_incidence_rate_by_sex(modelpreds.age.1218.all[model == 'model_1'], seroconverter_cohort_imputation)
modelpreds = modelpreds[, list(q= quantile(inc, prob=ps, na.rm = T),
                               q_label=paste0(p_labs, '_inc')), by=c('Sex', 'round')]	
modelpreds = dcast(modelpreds, Sex + round ~ q_label, value.var = "q")
setnames(modelpreds, c('M_inc', 'CL_inc', 'CU_inc'), c('inc', 'lb', 'ub'))
modelpreds <- rbind(modelpreds_all, modelpreds)
rm(modelpreds_all)

# by sex, age group and round
modelpreds.agegroup.1218 <- find_incidence_rate_by_age_group(modelpreds.age.1218.all[model == 'model_1'], seroconverter_cohort_imputation)
modelpreds.agegroup.1218 = modelpreds.agegroup.1218[, list(q= quantile(inc, prob=ps, na.rm = T),
                                                           q_label=paste0(p_labs, '_inc')), by=c('Sex', 'round', 'age_group')]	
modelpreds.agegroup.1218 = dcast(modelpreds.agegroup.1218, Sex + round + age_group ~ q_label, value.var = "q")
setnames(modelpreds.agegroup.1218, c('M_inc', 'CL_inc', 'CU_inc'), c('incidence', 'lb', 'ub'))

# by sex, age and round
modelpreds.age.1218 = modelpreds.age.1218.all[, list(q= quantile(inc, prob=ps, na.rm = T),
                                                     q_label=paste0(p_labs, '_inc')), by=c('Sex', 'model', 'round', 'round_label', 'ROUND', 'age')]	
modelpreds.age.1218 = dcast(modelpreds.age.1218, Sex + model + round + round_label + ROUND + age ~ q_label, value.var = "q")
setnames(modelpreds.age.1218, c('M_inc', 'CL_inc', 'CU_inc'), c('incidence', 'lb', 'ub'))
rm(modelpreds.age.1218.list)


#
# Summarise aic 
#

modelaics.age <- rbindlist(modelaics.age.list)
modelaics.age <- modelaics.age[, list(q= quantile(aic, prob=ps, na.rm = T),q_label=p_labs), by = c('Sex', 'model')]	
modelaics.age <- dcast(modelaics.age, Sex + model ~ q_label, value.var = "q")
rm(modelaics.age.list)


##########################

#  SAVE STATISTICS FOR PAPER

##########################

stats <- save_statistics_estimates(df_round, 
                                   seroconverter_cohort_total, seroconverter_cohort_agg, 
                                   seroconverter_cohort_agg_s, seroconverter_cohort_agg_r,
                                   modelpreds.agegroup.1218, modelpreds)

#########

#  PLOT 

#########

#
# DATA
# 

plot_data(seroconverter_cohort, seroconverter_cohort_agg_s, outdir)


#
# ESTIMATES
# 

plot_model_fit(modelpreds, modelpreds.age.1218, outdir)


#
#  MODEL FIT
# 

# generate predictions by 1-year age band
model_pred <- find_prediction_1y(modelpreds.age.1218.all, seroconverter_cohort_imputation)

# save statistics on prediction
stats_prediction <- find_statistics_prediction(model_pred)

# plot for all iterations
plot_predicted_incidence_events(model_pred, outdir)
plot_estimate_incidence_rates(model_pred, outdir)

# plot summary across iterations
plot_estimate_incidence_rates_summary(modelpreds.age.1218, seroconverter_cohort_crude, outdir)


###########################

#  SENSITIVITY ANALYSIS 

###########################

# finding smooth incidence rate using loess for every iteration
modelpreds.loess.age.1218.all <- find_incidence_estimates_loess(model_pred)

# save statistics on prediction
stats_prediction_loess <- find_statistics_prediction_loess(modelpreds.loess.age.1218.all)

# summarise 
modelpreds.loess.age.1218.all <- merge(modelpreds.loess.age.1218.all, select(df_round, -visit), by.x = 'round', by.y= 'round_numeric') 
modelpreds.loess.age.1218.all[, round_label := as.numeric(gsub('Round: (.+)', '\\1', ROUND))]
modelpreds.loess.age.1218 = modelpreds.loess.age.1218.all[, list(q= quantile(INC_CRUDE_SMOOTH, prob=ps, na.rm = T),
                                                                 q_label=paste0(p_labs, '_inc')), by=c('Sex', 'round_label', 'age')]	
modelpreds.loess.age.1218 = dcast(modelpreds.loess.age.1218, Sex + round_label + age ~ q_label, value.var = "q")
setnames(modelpreds.loess.age.1218, c('M_inc', 'CL_inc', 'CU_inc'), c('incidence', 'lb', 'ub'))

# plot comparison smooth incidence rate using loess and gam
plot_comparison_loess_gam(
  modelpreds.loess.age.1218,
  modelpreds.age.1218, 
  seroconverter_cohort_imputation,
  outdir)

#
#
# SAVE
#
#

if(!restrict_to_30_comms){
  file.name	<- file.incidence.inland
}else{
  file.name	<- file.incidence.30com.inland
}
tmp <- select(modelpreds.age.1218[model == 'model_1'], c('Sex', 'round_label', 'age', 'incidence', 'lb', 'ub'))
if(! file.exists(file.name))
{
  cat("Saving file:", file.name, '\n')
  write.csv(tmp, file = file.name, row.names = F)
}else{
  cat("File:", file.name, "already exists...\n")
}

if(!restrict_to_30_comms){
  file.name	<- file.incidence.samples.inland
}else{
  file.name	<- file.incidence.30com.samples.inland
}
tmp <- modelpreds.age.1218.all[model == 'model_1', .(Sex, round_label, age, iterations, fit, inc, iterations_within)]
if(! file.exists(file.name))
{
  cat("Saving file:", file.name, '\n')
  write.csv(tmp, file = file.name, row.names = F)
}else{
  cat("File:", file.name, "already exists...\n")
}

if(!restrict_to_30_comms){
  
  file.name	<- file.incidence.loess.inland
  tmp <- select(modelpreds.loess.age.1218, c('Sex', 'round_label', 'age', 'incidence', 'lb', 'ub'))
  if(! file.exists(file.name))
  {
    cat("Saving file:", file.name, '\n')
    write.csv(tmp, file = file.name, row.names = F)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
  
  file.name	<- file.incidence.loess.samples.inland
  tmp <- modelpreds.loess.age.1218.all[, .(Sex, round_label, age, iterations, INC_CRUDE_SMOOTH)]
  setnames(tmp, 'INC_CRUDE_SMOOTH', 'inc')
  if(! file.exists(file.name))
  {
    cat("Saving file:", file.name, '\n')
    write.csv(tmp, file = file.name, row.names = F)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
  
  file.name <- file.path(outdir, "Rakai_inc_model_fit_inland_221107.csv")
  if(! file.exists(file.name))
  {
    cat("Saving file:", file.name, '\n')
    write.csv(model_pred, file = file.name, row.names = F)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
  
  # For paper
  file.name <- file.path(outdir, 'incidence_inland_estimates_for_paper_221129.RDS')
  if(! file.exists(file.name))
  {
    cat("Saving file:", file.name, '\n')
    saveRDS(stats, file.name)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
  
  file.name <- file.path(outdir, 'incidence_inland_prediction_for_paper_221107.RDS')
  if(! file.exists(file.name))
  {
    cat("Saving file:", file.name, '\n')
    saveRDS(stats_prediction, file.name)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
  
  file.name <- file.path(outdir, 'incidence_inland_prediction_loess_for_paper_221116.RDS')
  if(! file.exists(file.name))
  {
    cat("Saving file:", file.name, '\n')
    saveRDS(stats_prediction_loess, file.name)
  }else{
    cat("File:", file.name, "already exists...\n")
  }
}





