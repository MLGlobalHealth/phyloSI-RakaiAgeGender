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
    file.path.seroconverter_cohort.30) |> file.exists() |> all() |> stopifnot()

# sensitivity analysis?
restrict_to_30_comms <- F

# load functions
source(file.path(gitdir.functions, 'functions_incidence_rate', 'incidence_rate_estimation_functions.R'))

# utils
rounds_group_1 <- c("R006","R007", "R008", "R009", "R010", "R011", "R012", "R013", "R014")
rounds_group_2 <- c("R015", "R015S", "R016", "R017", "R018")
rounds_group_3 <- c('R019')

# make df round
df_round <- make_df_round(rounds_group_1, rounds_group_2, rounds_group_3)

rounds_numeric_group_1 <- df_round[visit %in% rounds_group_1, round_numeric]
rounds_numeric_group_2 <- df_round[visit %in% rounds_group_2, round_numeric]
rounds_numeric_group_3 <- df_round[visit == rounds_group_3, round_numeric]

# load data 
if(!restrict_to_30_comms){
  file <- file.path.seroconverter_cohort
}else{
  file <- file.path.seroconverter_cohort.30
}
seroconverter_cohort.list <- readRDS(file)


###############################

# FIND SEROCONVERSION DATE #

###############################

N <- 50
seed = 12

modelpreds.age.1218.list = modelaics.age.list = vector(mode = 'list', length = N)

set.seed(seed)
for(i in 1:N){
  
  cat('\niteration', i)
  
  ####################################
  
  # GENERATE RANDOM DATE OF INFECTION
  
  ####################################
  
  seroconverter_cohort <- as.data.table(seroconverter_cohort.list[[i]])
  
  
  ############
  
  # MODEL FIT
  
  ############
  
  summary(seroconverter_cohort$hivinc)
  
  # MODEL 1
  gamfit.m.age.int <- gam(hivinc ~ s(round, bs="gp") + s(age, bs="gp")+ ti(round, age)  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="M" & py>0)) 
  gamfit.f.age.int <- gam(hivinc ~ s(round, bs="gp") + s(age, bs="gp")+ ti(round, age)  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="F" & py>0))
  
  # MODEL 2
  gamfit.m.age.2 <- gam(hivinc ~ s(round, age, bs="gp")+ ti(round, age)  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="M" & py>0)) 
  gamfit.f.age.2 <- gam(hivinc ~ s(round, age, bs="gp")+ ti(round, age)  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="F" & py>0))
  
  # MODEL 3
  gamfit.m.age.3 <- gam(hivinc ~ s(round, bs="gp") + s(age, bs="gp")  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="M" & py>0)) 
  gamfit.f.age.3 <- gam(hivinc ~ s(round, bs="gp") + s(age, bs="gp")  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="F" & py>0))
  
  # MODEL 3
  gamfit.m.age.4 <- gam(hivinc ~ s(round, age, bs="gp")  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="M" & py>0)) 
  gamfit.f.age.4 <- gam(hivinc ~ s(round, age, bs="gp")  + offset(log(py)), family = poisson, data = subset(seroconverter_cohort,sex=="F" & py>0))
  
  # prediction
  gg.age <- expand.grid(round=df_round$round_numeric, py=1, age=seq(15,49,1))
  
  modelpreds.age <- rbind(cbind(Sex="Male", model="model_1",gg.age, as.data.frame(predict(gamfit.m.age.int, gg.age, type="link", se.fit=TRUE))),
                          cbind(Sex="Female", model="model_1",gg.age, as.data.frame(predict(gamfit.f.age.int, gg.age, type="link", se.fit=TRUE))),
                          cbind(Sex="Male", model="model_2",gg.age, as.data.frame(predict(gamfit.m.age.2, gg.age, type="link", se.fit=TRUE))),
                          cbind(Sex="Female", model="model_2",gg.age, as.data.frame(predict(gamfit.m.age.2, gg.age, type="link", se.fit=TRUE))),
                          cbind(Sex="Male", model="model_3",gg.age, as.data.frame(predict(gamfit.m.age.3, gg.age, type="link", se.fit=TRUE))),
                          cbind(Sex="Female", model="model_3",gg.age, as.data.frame(predict(gamfit.m.age.3, gg.age, type="link", se.fit=TRUE))),
                          cbind(Sex="Male", model="model_4",gg.age, as.data.frame(predict(gamfit.m.age.4, gg.age, type="link", se.fit=TRUE))),
                          cbind(Sex="Female", model="model_4",gg.age, as.data.frame(predict(gamfit.m.age.4, gg.age, type="link", se.fit=TRUE)))
  )
  
  modelpreds.age <- modelpreds.age %>% 
    mutate(incidence=exp(fit), se.fit=se.fit, lb=exp(fit - 1.96*(se.fit)), ub=exp(fit + 1.96*se.fit)) %>%
    group_by(age, Sex, model) %>%
    merge(select(df_round, -visit), by.x = 'round', by.y= 'round_numeric') %>%
    mutate(round_label = as.numeric(gsub('Round: (.+)', '\\1', ROUND)))
  modelpreds.age.1218 <- modelpreds.age[which(modelpreds.age$round_label %in% 10:18), ]
  
  # aic
  modelaics.age <- rbind(cbind(Sex="Male", model="model_1",aic = gamfit.m.age.int$aic),
                         cbind(Sex="Female", model="model_1",aic = gamfit.f.age.int$aic),
                         cbind(Sex="Male", model="model_2",aic = gamfit.m.age.2$aic),
                         cbind(Sex="Female", model="model_2",aic = gamfit.f.age.2$aic),
                         cbind(Sex="Male", model="model_3",aic = gamfit.m.age.3$aic),
                         cbind(Sex="Female", model="model_3",aic = gamfit.f.age.3$aic),
                         cbind(Sex="Male", model="model_4",aic = gamfit.m.age.4$aic),
                         cbind(Sex="Female", model="model_4",aic = gamfit.f.age.4$aic)
  )%>%
    as.data.frame()
  modelaics.age$aic <- as.numeric( modelaics.age$aic)
  

  # prepare
  modelpreds.age.1218$iterations <- i
  modelaics.age$iterations <- i
  
  modelpreds.age.1218.list[[i]] <- modelpreds.age.1218
  modelaics.age.list[[i]] <- modelaics.age
  
}

#################

# SAVE

#################

if(!restrict_to_30_comms){
  file.name <- file.incidence.fits
}else{
  file.name <- file.incidence.30com.fits
}
if(! file.exists(file.name))
{
    cat("Saving file:", file.name, '\n')
    save(modelpreds.age.1218.list, seroconverter_cohort.list, modelaics.age.list, file = file.name)
}else{
    cat("File:", file.name, "already exists...\n")
}
