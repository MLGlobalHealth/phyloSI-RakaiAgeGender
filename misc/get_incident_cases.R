library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'preliminary', 'IncidentCases')

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'HIV_R15_R18_VOIs_220129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'quest_R15_R18_VoIs_220129.csv')
file.path.flow <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'FlowR15_R18_VoIs_220129.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.incidence	<- file.path(indir.deepsequencedata, 'RCCS_R15_R18', "Rakai_incpredictions_220524.csv")
file.df_round <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_20220329.RData')

hiv <- as.data.table(read.csv(file.path.hiv))
quest <- as.data.table(read.csv(file.path.quest))
flow <- as.data.table(read.csv(file.path.flow))
community.keys <- as.data.table(read.csv(file.community.keys))
incidence <- as.data.table(read.csv(file.incidence))
# 
load(file.df_round)

source(file.path(indir.repository, 'functions', 'utils.R'))
source(file.path(indir.repository, 'misc', 'functions', 'preprocess_meta_data-functions.R'))

file.stan.model <- file.path(indir.repository, 'misc', 'stan_models', 'one_dimensional_gp_binomial.stan')
model = rstan::stan_model(file.stan.model)
ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')


###############################

# CENSUS ELIGIBLE INDIVIDUALS #

###############################

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
flow <- merge(flow, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# Code for ineligibility

flow[, reason_ineligible := NA_character_]

flow[locate1==10 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==2 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==13 & locate2==8, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==5 & locate2==10, reason_ineligible := "Out_migrated"]
flow[locate1==3 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==13, reason_ineligible := "Out_migrated"]
flow[locate1==6 & locate2==10, reason_ineligible := "Out_migrated"]

flow[locate1==7 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==2 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==6 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==3 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==5 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==7, reason_ineligible := "Already_seen"]
flow[locate1==17 & locate2==8, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==88, reason_ineligible := "Already_seen"]
flow[locate1==7 & locate2==2, reason_ineligible := "Already_seen"]

flow[locate1==11 & locate2==8, reason_ineligible := "Dead"]
flow[locate2==11, reason_ineligible := "Dead"]

flow[resident==0, reason_ineligible := "not_resident"]

flow[ageyrs<15 | ageyrs > 49, reason_ineligible := "Not_within_eligible_age_range"]

flow[is.na(reason_ineligible), reason_ineligible := 'none']

# find count eligible
re <- flow[, list(count = .N), by = c('reason_ineligible', 'round', 'comm', 'ageyrs', 'sex')]
re <- dcast.data.table(re, round + comm + ageyrs + sex ~ reason_ineligible, value.var = 'count')
re[is.na(re)] = 0
re[, ELIGIBLE := round(none + Out_migrated / 2)]
re <- re[ELIGIBLE != 0]

# additional variable
colnames(re) <- toupper(colnames(re))
re[, ROUND := substring(ROUND, 3)]

# save
write.csv(re, file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220411.csv'), row.names = F)

# find count eligible comm num
re <- flow[, list(count = .N), by = c('reason_ineligible', 'round', 'comm_num', 'ageyrs', 'sex')]
re <- dcast.data.table(re, round + comm_num + ageyrs + sex ~ reason_ineligible, value.var = 'count')
re[is.na(re)] = 0
re[, ELIGIBLE := round(none + Out_migrated / 2)]
re <- re[ELIGIBLE != 0]

# additional variable
colnames(re) <- toupper(colnames(re))
re[, ROUND := substring(ROUND, 3)]

# save
write.csv(re, file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220620.csv'), row.names = F)


# table and plot
tmp <- re[, list(count = sum(ELIGIBLE)), by = c('ROUND', 'COMM')]
knitr::kable(tmp[order(COMM,ROUND)])

if(1){
  ggplot(re[COMM == 'inland'], aes(x = AGEYRS, y = ELIGIBLE)) +
    geom_bar(stat = 'identity', position = "dodge") +
    labs(y = 'Census eligible individuals in round 16', x = 'Age') +
    facet_grid(ROUND~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-CensusEligibleIndividuals.png'), w = 7, h = 6)
}



####################################################

# HIV PREVALENCE AND CENSUS ELIGIBLE SUSCEPTIBLES #

###################################################

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate)]

# Set to date format
rin[, intdate := as.Date(intdate, format = '%d-%b-%y')]

# keep only first round
rin <- rin[, list(intdate = min(intdate),
                  round = round[intdate == min(intdate)][1],
                  ageyrs = ageyrs[intdate == min(intdate)][1],
                  comm_num = comm_num[intdate == min(intdate)][1]), by = c('study_id', 'sex')]

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

# fit gaussian process
AGEYRS_COVARIATE <- rprev[, sort(unique(AGEYRS))]
rprev <- rprev[order(ROUND, SEX, COMM, AGEYRS)]
rprevfit <- rprev[, {
              fit = rstan::sampling(model,data=list(N = length(AGEYRS),
                                                    x = sort(AGEYRS),
                                                    x_hat = AGEYRS_COVARIATE,
                                                    N_hat = length(AGEYRS_COVARIATE),
                                                    COUNT = COUNT,
                                                    TOTAL_COUNT = TOTAL_COUNT,
                                                    x_in_x_hat = which(AGEYRS %in% AGEYRS_COVARIATE)),
                                    iter=1000,warmup=500,chains=1,
                                    control = list(max_treedepth = 15, adapt_delta = 0.99))
              fit_samples = rstan::extract(fit)
              tmp1 = as.data.table( reshape2::melt(fit_samples[['mu']]))
              setnames(tmp1, 2, 'AGEYRS_INDEX')
              tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('AGEYRS_INDEX')]
              tmp1 = dcast(tmp1, AGEYRS_INDEX ~ q_label, value.var = "q")
              tmp1[, AGEYRS := AGEYRS_COVARIATE[AGEYRS_INDEX]]
              list(M = tmp1[, M], CL = tmp1[, CL], CU = tmp1[, CU], AGEYRS = tmp1[, AGEYRS])
              }, by = c('ROUND', 'SEX', 'COMM')]
rprevfit <- merge(rprevfit, rprev, by = c('ROUND', 'SEX', 'COMM', 'AGEYRS'), all.x = T)

# merge to census eligble and find susceptible
setnames(rprevfit, c('M', 'CL', 'CU'), c('PREVALENCE_PROPORTION', 'PREVALENCE_PROPORTION_CL', 'PREVALENCE_PROPORTION_CU'))
rprevfit[PREVALENCE_PROPORTION < 0 , PREVALENCE_PROPORTION := 0]
rprevfit[, ROUND := substring(ROUND, 3)]
resusc <- merge(re, rprevfit, by = c('COMM', 'AGEYRS', 'SEX', 'ROUND'))
resusc[, SUSCEPTIBLE := ELIGIBLE * (1 - PREVALENCE_PROPORTION)]

# save
write.csv(resusc, file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_prevalence_220411.csv'), row.names = F)
# resusc <- as.data.table(read.csv(file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_prevalence_220411.csv')))

# plot
if(0){
  tmp <- resusc[ROUND == '16']
  tmp[, HIV_NEGATIVE := TOTAL_COUNT - COUNT]
  setnames(tmp, 'COUNT', 'HIV_POSITIVE')
  tmp <- melt.data.table(tmp[, .(SEX, COMM, AGEYRS, HIV_NEGATIVE, HIV_POSITIVE)], id.vars = c('SEX', 'COMM', 'AGEYRS'))
  ggplot(tmp, aes(x = AGEYRS, y = value, fill = variable)) +
    geom_bar(stat = 'identity', position = "stack") +
    labs(y = 'Number of participants in round 16', x = 'Age', fill = '') +
    facet_grid(COMM~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom') +
    scale_y_continuous(expand = c(0,0))
  ggsave(paste0(outdir, '-Participants.png'), w = 7, h = 6)

  tmp <- resusc[COMM == 'inland']
  ggplot(tmp, aes(x = AGEYRS)) +
    geom_line(aes(y = PREVALENCE_PROPORTION)) +
    geom_ribbon(aes(ymin = PREVALENCE_PROPORTION_CL, ymax = PREVALENCE_PROPORTION_CU), alpha = 0.5) +
    geom_point(aes(y = COUNT / TOTAL_COUNT), col = 'darkred') +
    labs(y = 'Prevalence proportion', x = 'Age') +
    facet_grid(ROUND~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom') +
    scale_y_continuous(labels = scales::percent_format())
  ggsave(paste0(outdir, '-PrevalenceProportionGPFit.png'), w = 7, h = 6)

  tmp <- resusc[COMM == 'fishing']
  ggplot(tmp, aes(x = AGEYRS)) +
    geom_line(aes(y = PREVALENCE_PROPORTION)) +
    geom_ribbon(aes(ymin = PREVALENCE_PROPORTION_CL, ymax = PREVALENCE_PROPORTION_CU), alpha = 0.5) +
    geom_point(aes(y = COUNT / TOTAL_COUNT), col = 'darkred') +
    labs(y = 'Prevalence proportion', x = 'Age') +
    facet_grid(ROUND~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom') +
    scale_y_continuous(labels = scales::percent_format())
}

##################

# INCIDENT CASES #

##################

#
# merge to incidence
colnames(incidence) <- toupper(colnames(incidence))
setnames(incidence, 'AGE', 'AGEYRS')
incidence[, COMM := 'inland']
incidence[, SEX := substring(SEX, 1, 1)]
incidence <- incidence[ROUND >= 14]
incidence[, ROUND := as.character(ROUND)]

if(0){
  ggplot(incidence, aes(x = AGEYRS, y = INCIDENCE)) +
    geom_bar(aes(fill = MODEL), stat = 'identity', position = "dodge") +
    geom_errorbar(aes(ymin = LB, ymax = UB),  alpha = 0.5) +
    labs(y = 'Incidence rate per 1 PY in inland community', x = 'Age') +
    facet_grid(ROUND~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-IncidenceEstimate_220524.png'), w = 6, h = 7)
}

# merge to susceptible
dir <- merge(incidence, resusc, by = c('COMM', 'AGEYRS', 'SEX', 'ROUND'))

# find length in years of each round
colnames(df_round) <- toupper(colnames(df_round))
df_round[, ROUND := as.character(ROUND)]
dir <- merge(dir, df_round, by = 'ROUND')
dir[, ROUND_SPANYRS := .year.diff(MAX_SAMPLE_DATE, MIN_SAMPLE_DATE)]

# find incident cases
dir[, INCIDENT_CASES:= SUSCEPTIBLE * ROUND_SPANYRS * INCIDENCE]
dir[, INCIDENT_CASES_UB:= SUSCEPTIBLE * ROUND_SPANYRS * UB]
dir[, INCIDENT_CASES_LB:= SUSCEPTIBLE * ROUND_SPANYRS * LB]

# save
write.csv(dir, file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_incident_cases_220524.csv'), row.names = F)
# dir <- as.data.table(read.csv(file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_incident_cases_220427.csv')))

# plot
if(0){
  ggplot(dir, aes(x = AGEYRS)) +
    geom_bar(aes(y = INCIDENT_CASES, fill = MODEL), stat = 'identity', position = "dodge") +
    geom_errorbar(aes(ymin = INCIDENT_CASES_LB, ymax = INCIDENT_CASES_UB), alpha = 0.5) +
    labs(y = 'Expected number of incident cases \nin inland community', x = 'Age') +
    facet_grid(ROUND~SEX, label = 'label_both') +
    theme_bw() +
    theme(legend.position = 'bottom')
  ggsave(paste0(outdir, '-IncidentCases_0524.png'), w = 7, h = 6)

  tmp <- dir[, list(INCIDENT_CASES = round(sum(INCIDENT_CASES), digits = 1),
                   INCIDENT_CASES_UB = round(sum(INCIDENT_CASES_UB), digits = 1),
                   INCIDENT_CASES_UB = round(sum(INCIDENT_CASES_UB), digits = 1)), by = c('ROUND', 'MODEL')]
  knitr::kable(subset(tmp, select = - c(MODEL)))

}

