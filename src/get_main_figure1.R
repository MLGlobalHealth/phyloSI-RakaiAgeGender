library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
require(lubridate)
library(rstan)
library(gridExtra)
library(lognorm)
library(ggExtra)
library(Hmisc)

# laptop
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir <- '~/git/phyloflows'
  outdir <- '~/Box\ Sync/2021/phyloflows/'

  jobname <- 'new_treatment_cascade'
  stan_model <- 'gp_221201d'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/home/andrea'))
{
  indir <-'~/git/phyloflows'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'test'
  stan_model <- 'gp_221201d'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-stan_model')
  stopifnot(args_line[[7]]=='-jobname')
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  stan_model <- args_line[[6]]
  jobname <- args_line[[8]]
}

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))
outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
if(!dir.exists(dirname(outfile.figures))) dir.create(dirname(outfile.figures))
if(!dir.exists(dirname(outdir.table))) dir.create(dirname(outdir.table))

# indicators -- fixed
only.transmission.after.start.observational.period <- T
only.transmission.before.stop.observational.period <- T
use_number_susceptible_offset <- T
only.one.community <- 'inland'
use_contact_rates_prior <- T

# indicators -- sensitivity analyses
nonparticipants.treated.like.participants <- F
nonparticipants.not.treated <- F
nonparticipants.male.relative.infection <- 1
nonparticipants.female.relative.infection <- 1
remove.pairs.from.rounds <- NULL
use_loess_inc_estimates <- F
pairs_replicates.seed <- NULL
viremic_viral_load_200ml <- F
use_30com_inc_estimates <- F
use_30com_pairs <- F
use_tsi_non_refined <- F

# obtained in script/ 
file.incidence.inland <- file.path(indir, 'data', "Rakai_incpredictions_inland_221107.csv")#central analysis
file.incidence.30com.inland <- file.path(indir, 'data', "Rakai_incpredictions_inland_221119.csv")#sensitivity analysis restricted to community continuoursly surveyed
file.incidence.loess.inland <- file.path(indir, 'data', "Rakai_incpredictions_loess_inland_221116.csv")#sensitivity analysis using estimates obtained with loess regresssion
file.incidence.samples.inland <- file.path(indir, 'data', "Rakai_incpredictions_samples_inland_221107.csv")

# obtained in src/ for analysis
file.path.round.timeline <- file.path(indir, 'data', 'RCCS_round_timeline_220905.RData')
file.eligible.count <- file.path(indir, 'data', 'RCCS_census_eligible_individuals_221116.csv')
file.participation <- file.path(indir, 'data', 'RCCS_participation_221208.csv')
file.prevalence.prop <- file.path(indir, 'fit', 'RCCS_prevalence_estimates_221116.csv')

# obtained in misc/ for analysis
file.pairs <- file.path(indir, 'data', 'pairsdata_toshare_d1_w11_netfrompairs_postponessrem.rds')
file.pairs.nonrefined <- file.path(indir, 'data', 'pairsdata_toshare_d1_w11_netfrompairs_seropairs_sensnoref.rds')

file.treatment.cascade.prop.participants <- file.path(indir, 'fit', "RCCS_treatment_cascade_participants_estimates_221208.csv")
file.treatment.cascade.prop.nonparticipants <- file.path(indir, 'fit', "RCCS_treatment_cascade_nonparticipants_estimates_221208.csv")
file.treatment.cascade.prop.participants.samples <- file.path(indir, 'fit', paste0('RCCS_treatment_cascade_participants_posterior_samples_221208.rds'))
file.treatment.cascade.prop.nonparticipants.samples <- file.path(indir, 'fit', paste0('RCCS_treatment_cascade_nonparticipants_posterior_samples_221208.rds'))

file.treatment.cascade.prop.participants.vl200 <- file.path(indir, 'fit', "RCCS_treatment_cascade_participants_estimates_vl200_221208.csv")
file.treatment.cascade.prop.nonparticipants.vl200 <- file.path(indir, 'fit', "RCCS_treatment_cascade_nonparticipants_estimates_vl200_221208.csv")
file.treatment.cascade.prop.participants.vl200.samples <- file.path(indir, 'fit', paste0('RCCS_treatment_cascade_participants_posterior_samples_vl200_221208.rds')) 
file.treatment.cascade.prop.nonparticipants.vl200.samples <- file.path(indir, 'fit', paste0('RCCS_treatment_cascade_nonparticipants_posterior_samples_vl200_221208.rds')) 

# obtained in misc/ for plots
file.unsuppressed.share <- file.path(indir, 'fit', paste0('RCCS_unsuppressed_share_sex_221208.csv'))
file.unsuppressed_rate_ratio <- file.path(indir, 'fit', paste0('RCCS_unsuppressed_ratio_sex_221208.csv'))
file.prevalence.share <- file.path(indir, 'fit', paste0('RCCS_prevalence_share_sex_221116.csv'))
file.unsuppressed_median_age <-file.path(indir, 'fit', paste0('RCCS_unsuppressed_median_age_221208.csv'))

# sexual partnerships  rates
file.number.sexual.partnerships <- file.path(indir, 'data', paste0('age-age-group-est-cntcts-r15.rds'))
file.sexual.partnerships.rates <- file.path(indir, 'data', paste0('inland_R015_cntcts_rate_1130b.rds'))

# obtained in script/ for plots
file.incidence.samples.inland	<- file.path(indir, 'data', "Rakai_incpredictions_samples_inland_221107.csv")
file.incidence.30com.samples.inland <- file.path(indir, 'data', "Rakai_incpredictions_samples_inland_221119.csv")
file.incidence.loess.samples.inland	<- file.path(indir, 'data', "Rakai_incpredictions_loess_samples_inland_221116.csv")

# stan model
path.to.stan.model <- file.path(indir, 'stan_models', paste0(stan_model, '.stan'))

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'statistics_functions.R'))
source(file.path(indir, 'functions', 'stan_utils.R'))

# load pairs
if(use_tsi_non_refined){
  pairs.all <- read_pairs(file.pairs.nonrefined)
} else{
  pairs.all <- read_pairs(file.pairs)
}

# load round timeline
load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]
# load census eligible ount
eligible_count_smooth <- fread(file.eligible.count)

# load participation (% of census eligible population)
participation <- fread(file.participation)

# load proportion prevalence
proportion_prevalence <- fread(file.prevalence.prop)

# load non-suppressed proportion 
if(viremic_viral_load_200ml){
  treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants.vl200, 
                                              file.treatment.cascade.prop.nonparticipants.vl200)
  if(file.exists(file.treatment.cascade.prop.participants.vl200.samples) & file.exists(file.treatment.cascade.prop.nonparticipants.vl200.samples)){
    treatment_cascade_samples <- read_treatment_cascade_samples(file.treatment.cascade.prop.participants.vl200.samples, 
                                                                file.treatment.cascade.prop.nonparticipants.vl200.samples)
  }
}else{
  treatment_cascade <- read_treatment_cascade(file.treatment.cascade.prop.participants, 
                                              file.treatment.cascade.prop.nonparticipants)
  if(file.exists(file.treatment.cascade.prop.participants.samples) & file.exists(file.treatment.cascade.prop.nonparticipants.samples)){
    treatment_cascade_samples <- read_treatment_cascade_samples(file.treatment.cascade.prop.participants.samples, 
                                                                file.treatment.cascade.prop.nonparticipants.samples)
  }
}

# load incidence estimates 
incidence.inland <- fread(file.incidence.inland)

# for offset
df_estimated_contact_rates <- as.data.table(readRDS(file.sexual.partnerships.rates))#estimated secxual contact rate
if('cntct.rates' %in% colnames(df_estimated_contact_rates))
  setnames(df_estimated_contact_rates, 'cntct.rates', 'cntct.rate')

#for plots
unsuppressed_rate_ratio <- fread(file.unsuppressed_rate_ratio) # sex ratio of unsuppression rate
df_reported_contact <- as.data.table(readRDS(file.number.sexual.partnerships)) # estimated number of sexual contacts
unsuppressed_share <- fread(file.unsuppressed.share) # share of unsuppressed count by sex
infected_share <- fread(file.prevalence.share) # share of infected count by sex
unsuppressed_median_age <-fread(file.unsuppressed_median_age) # median age of unsuppressed

#
# Define start time, end time and cutoff
#

start_first_period_inland <- df_round_inland[round == 'R010', min_sample_date] # "2003-09-26"
stop_first_period_inland <- df_round_inland[round == 'R015', max_sample_date] # "2013-07-05"
start_second_period_inland <-df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"
stop_second_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

stopifnot(start_first_period_inland < stop_first_period_inland)
stopifnot(stop_first_period_inland < start_second_period_inland)
stopifnot(start_second_period_inland < stop_second_period_inland)

df_round <- make.df.round(df_round_inland)

df_period <- make.df.period(start_first_period_inland, stop_first_period_inland, 
                            start_second_period_inland, stop_second_period_inland, 
                            df_round)


#
# Find count eligible susceptible / infected / infected unsuppressed 
# 


# by round
eligible_count_round <- add_susceptible_infected(eligible_count_smooth, proportion_prevalence, participation, 
                                                 nonparticipants.male.relative.infection, 
                                                 nonparticipants.female.relative.infection)
eligible_count_round <- add_infected_unsuppressed(eligible_count_round, treatment_cascade, participation, 
                                                  nonparticipants.treated.like.participants, 
                                                  nonparticipants.not.treated)
eligible_count_round[, table(ROUND, COMM)]


#
# Find incidence cases
#

# by round
incidence_cases_round <- get_incidence_cases_round(incidence.inland, eligible_count_round)
incidence_cases_round[, table(ROUND, COMM)]

# summarise by time period
incidence_cases <- summarise_incidence_cases_period(incidence_cases_round, df_period)
incidence_cases[, table(PERIOD, COMM)]


#
# PREPARE MAPS
#

# prepare age map
df_age <- get.age.map(age_bands_reduced = 4)
df_age_aggregated <- get.age.aggregated.map(c('15-24', '25-34', '35-49'))

# prepare direciton and commuity
df_direction <- get.df.direction()
df_community <- get.df.community()

find_palette_round()

# control linewidth
LWD = .6
ALPHA = .3

make_subplot <- function(DT, M, IL, IU, y_lab='Contribution to incidence cases')
{
    tmp <- copy(DT)
    # icr <- copy(incidence_rates_round.samples)
    tmp <- merge(tmp, df_community, by = 'COMM')
    tmp <- merge(tmp, df_round, by = c('COMM', 'ROUND'))
    tmp[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
    tmp[, SEX_LABEL := 'Women']
    tmp[SEX== 'M', SEX_LABEL := 'Men']

    ggplot(tmp[COMM == 'inland' & round %in% c(10, 12, 14, 16, 18)]) +
        geom_ribbon(aes_string(x = 'AGEYRS', ymin = IL , ymax = IU , fill = 'SEX_LABEL'),  alpha = ALPHA) +
        geom_line(aes_string(x = 'AGEYRS', y = M, col = 'SEX_LABEL'), lwd=LWD) +
        labs(y = y_lab, x = 'Age') +
        facet_grid(.~LABEL_ROUND, scale = 'free_y') +
        theme_bw() +
        scale_color_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
        scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
        theme(legend.position = 'none', 
            strip.background = element_rect(colour="white", fill="white"), 
            legend.title = element_blank(), 
            strip.text = element_text(size = 9.3), 
            axis.title = element_text(size = 12)) + 
        scale_x_continuous(expand = c(0,0), breaks = c(seq(15, 49, 5))) + 
        scale_y_continuous(labels = scales::percent, limits = c(0, NA), expand = expansion(mult = c(0, .05))) + 
        NULL
}

plot_incident_rates_over_time_2 <- function(incidence_cases_round, 
                                            incidence_rates_round.samples,
                                            eligible_count_round,
                                            outdir, outdir.table)
{
    #
    # median age at infection
    ps <- c(0.5, 0.025, 0.975)
    p_labs <- c('M','CL','CU')
    icr <- copy(incidence_rates_round.samples)
    rm(incidence_rates_round.samples)
    icr[COMM == 'fishing', REF.ROUND := 'R015']
    icr[COMM == 'inland', REF.ROUND := 'R010']

    medage <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
    medage[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
    medage[, WEIGHTED_INCIDENCE := INCIDENT_CASES / sum(INCIDENT_CASES), by = c('COMM', 'ROUND', 'SEX', 'iterations')]
    medage <- medage[, list(value = matrixStats::weightedMedian(AGEYRS, WEIGHTED_INCIDENCE )), by = c('iterations', 'COMM', 'ROUND', 'SEX')]
    medage = medage[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('COMM', 'ROUND', 'SEX')]	
    medage = dcast(medage, ... ~ q_label, value.var = "q")


    #
    # incidence rate per person per round

    tmp <- copy(incidence_cases_round)
    tmp <- merge(tmp, df_community, by = 'COMM')
    tmp <- merge(tmp, df_round, by = c('COMM', 'ROUND'))
    tmp[, LABEL_ROUND2 := gsub('(.+)\n.*', '\\1', LABEL_ROUND)]
    tmp[, SEX_LABEL := 'Women']
    tmp[SEX== 'M', SEX_LABEL := 'Men']

    # prepare median age
    max_y_limits <- ifelse(tmp[, max(INCIDENCE)] > 0.021, 2.7, 2.1)
    median_age <-  copy(medage)
    median_age[, SEX_LABEL := 'Women']
    median_age[SEX== 'M', SEX_LABEL := 'Men']
    set.seed(12)
    median_age[ROUND == 'R018' & SEX_LABEL == 'Women', M := M + runif(length(M), 0, 1)]
    median_age[ROUND == 'R018' & SEX_LABEL == 'Men', M := M + runif(length(M), -1, 0)]
    tmp1 <- median_age[COMM == 'inland' & ROUND %in% c('R010', 'R012', 'R014', 'R016', 'R018')]
    tmp1 <- merge(tmp1, df_round, by = c('COMM', 'ROUND'))

    p1 <- ggplot(tmp[COMM == 'inland' & round %in% c(10, 12, 14, 16, 18)]) +
    geom_line(aes(x = AGEYRS, y = INCIDENCE*100, col = SEX_LABEL), lwd=LWD) +
    geom_ribbon(aes(x = AGEYRS, ymin = LB *100, ymax = UB* 100, fill = SEX_LABEL),  alpha = ALPHA) +
    # geom_errorbarh(data = tmp1[SEX_LABEL=='Men' ], aes(y = 0.01, xmin = CL, xmax = CU, col = SEX_LABEL), size =1.5) +
    # geom_errorbarh(data = tmp1[SEX_LABEL=='Men' & ROUND == 'R018'], aes(y = 0.025, xmin = CL, xmax = CU, col = SEX_LABEL), size =1.5) +
    geom_point(data = tmp1[SEX_LABEL=='Men' ], aes(y = 0.08, x = M, fill = SEX_LABEL,  col = SEX_LABEL), shape = 25, size =3,  alpha = 0.5) +
    geom_point(data = tmp1[SEX_LABEL=='Men' ], aes(y = 0.08, x = M, col = SEX_LABEL), shape = 6, size =3) +
    # geom_errorbarh(data = tmp1[SEX_LABEL=='Women' ], aes(y = 0.01, xmin = CL, xmax = CU, col = SEX_LABEL), size =1.5) +
    geom_point(data = tmp1[SEX_LABEL=='Women'], aes(y = 0.08, x = M,  fill = SEX_LABEL,col = SEX_LABEL), shape = 25, size =3,  alpha = 0.5) +
    geom_point(data = tmp1[SEX_LABEL=='Women'], aes(y = 0.08, x = M,  col = SEX_LABEL), shape = 6, size =3) +
    labs(y = 'Incidence rates\nper 100 person-years', x = 'Age') +
    facet_grid(.~LABEL_ROUND, scale = 'free_y') +
    theme_bw() +
    scale_color_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
    theme(legend.position = 'none', 
        strip.background = element_rect(colour="white", fill="white"), 
        legend.title = element_blank(), 
        strip.text = element_text(size = 9.3), 
        axis.title = element_text(size = 12)) + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq(15, 49, 5))) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) + 
    coord_cartesian(ylim= c(0, max_y_limits))
    # ggsave(paste0(outdir, '-data-incidence_rate_round_sex_inland_short.pdf'), w = 8, h = 3.5)

    return(p1)
}

################################
# Load incidence rates samples #
################################

incidence_rates_round.samples <- load_incidence_rates_samples(file.incidence.samples.inland)
ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','IL','IU')

#############################
# CONTRIBUTION TO INCIDENCE #
#############################

# contribution to round incidence
data_1d <- incidence_rates_round.samples[, {
    z <- sum(INCIDENCE.DRAW);
    list(
        SEX=SEX,
        AGEYRS=AGEYRS,
        INCIDENCE_ROUND_TOTAL = z,
        INCIDENCE_CONTRIBUTION = INCIDENCE.DRAW/z
    )} , by=c('ROUND', 'COMM', 'iterations')][, 
    list(
        q = quantile(INCIDENCE_CONTRIBUTION, probs=ps),
        q_label=p_labs),
    by=c('ROUND', 'COMM', 'SEX', 'AGEYRS')] |>
    dcast( ROUND+COMM+SEX+AGEYRS~q_label, value.var='q')

p_1d <- make_subplot(data_1d, M='M', IL='IL', IU='IU', 
    y_lab = "Contribution to incidence rates")
# ggsave('~/Downloads/incidence_rates_contribution.pdf',p_1c, w=8, h=3.5)
# ggsave('~/Downloads/incidence_rates_contribution.png',p_1c, w=8, h=3.5)


###################
# INCIDENCE RATES # 
###################

library(patchwork)
naturemed_reqs()

p_c <- plot_incident_rates_over_time_2(incidence_cases_round, incidence_rates_round.samples, eligible_count_round, outfile.figures, outdir.table)

p_1cd <- (p_c + reqs + theme(axis.title.x=element_blank()))/(p_1d + theme(strip.background = element_blank(), strip.text = element_blank()) + reqs) 

p_1cd
ggsave('~/Downloads/MainFigure1cd.pdf',p_1cd, w=8.2, h=5.3)


#########################################
# EDF5: CONTRIBUTION TO INCIDENCE CASES #
#########################################

# load susceptibles by strata
dsusc <- incidence_cases_round[, .(COMM, ROUND, AGEYRS, SEX, SUSCEPTIBLE)]
data_edf5 <- merge( 
    incidence_rates_round.samples,
    dsusc,
    by=c('COMM', 'ROUND', 'AGEYRS', 'SEX'))

data_edf5 <- data_edf5[, {
    draw <- INCIDENCE.DRAW * SUSCEPTIBLE * ROUND_SPANYRS;
    tot <- sum(draw);
    list(
        SEX=SEX,
        AGEYRS=AGEYRS,
        INCIDENCE_CASES_ROUND_TOTAL=tot,
        INCIDENCE_CASES_CONTRIBUTION=draw/tot
    )
}, by=c('ROUND', 'COMM', 'iterations')][, 
    list(
        q = quantile(INCIDENCE_CASES_CONTRIBUTION, probs=ps),
        q_label=p_labs),
    by=c('ROUND', 'COMM', 'SEX', 'AGEYRS')]|>
    dcast( ROUND+COMM+SEX+AGEYRS~q_label, value.var='q')

p_edf5 <- make_subplot(data_1d, M='M', IL='IL', IU='IU', 
    y_lab = "Contribution to incident cases")

ggsave_nature('~/Downloads/incidence_cases_contribution.pdf',p_edf5+reqs , w=16, h=8)
ggsave_nature('~/Downloads/incidence_cases_contribution.png',p_edf5+reqs , w=16, h=8)
