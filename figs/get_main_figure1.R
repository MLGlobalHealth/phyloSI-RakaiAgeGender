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
gitdir <- here::here()
source(file.path(gitdir, 'config.R'))

# load functions
source(file.path(gitdir.R, 'utils.R'))
source(file.path(gitdir.R.flow, 'summary_functions.R'))
source(file.path(gitdir.R.flow, 'plotting_functions.R'))
source(file.path(gitdir.R.flow, 'statistics_functions.R'))
source(file.path(gitdir.R.flow, 'stan_utils.R'))

# load pairs
pairs.all <- read_pairs(file.pairs)

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
treatment_cascade <- read_treatment_cascade(
    file.treatment.cascade.prop.participants, 
    file.treatment.cascade.prop.nonparticipants)

# get posterior samples 
treatment_cascade_samples <- read_treatment_cascade_samples(
    file.treatment.cascade.prop.participants.samples, 
    file.treatment.cascade.prop.nonparticipants.samples)

# load incidence estimates 
incidence.inland <- fread(file.incidence.inland)

#
# Define start time, end time and cutoff
#

start_first_period_inland <- df_round_inland[round == 'R010', min_sample_date] 
stop_first_period_inland <- df_round_inland[round == 'R015', max_sample_date] 
start_second_period_inland <-df_round_inland[round == 'R016', min_sample_date]
stop_second_period_inland <- df_round_inland[round == 'R018', max_sample_date]
stopifnot(start_first_period_inland < stop_first_period_inland)
stopifnot(stop_first_period_inland < start_second_period_inland)
stopifnot(start_second_period_inland < stop_second_period_inland)

df_round <- make.df.round(df_round_inland)

df_period <- make.df.period(
    start_first_period_inland, 
    stop_first_period_inland, 
    start_second_period_inland, stop_second_period_inland, 
    df_round)


#
# Find count eligible susceptible / infected / infected unsuppressed 
# 


# by round
eligible_count_round <- add_susceptible_infected(
    eligible_count_smooth, 
    proportion_prevalence, 
    participation, 
    nonparticipants.male.relative.infection=1, 
    nonparticipants.female.relative.infection=1)

eligible_count_round <- add_infected_unsuppressed(
    eligible_count_round, 
    treatment_cascade, participation, 
    nonparticipants.treated.like.participants = FALSE, 
    nonparticipants.not.treated = FALSE)



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
        facet_grid(.~LABEL_ROUND, scales = 'free_y') +
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

plot_incident_rates_over_time_2 <- function(
    incidence_cases_round, 
    incidence_rates_round.samples,
    eligible_count_round)
{
    #
    # median age at infection
    ps <- c(M=0.5, CL=0.025, CU=0.975)
    icr <- copy(incidence_rates_round.samples)
    rm(incidence_rates_round.samples)
    icr[COMM == 'fishing', REF.ROUND := 'R015']
    icr[COMM == 'inland', REF.ROUND := 'R010']

    medage <- merge(icr, eligible_count_round, by = c('COMM', 'ROUND', 'AGEYRS', 'SEX'))
    medage[, INCIDENT_CASES := SUSCEPTIBLE * INCIDENCE.DRAW]
    medage[, WEIGHTED_INCIDENCE := INCIDENT_CASES / sum(INCIDENT_CASES), by = c('COMM', 'ROUND', 'SEX', 'iterations')]
    medage <- medage[, list(value = matrixStats::weightedMedian(AGEYRS, WEIGHTED_INCIDENCE )), by = c('iterations', 'COMM', 'ROUND', 'SEX')]
    medage = medage[, list(q= quantile(value, prob=ps, na.rm = T), q_label=names(ps)), by=c('COMM', 'ROUND', 'SEX')]	
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
    facet_grid(.~LABEL_ROUND, scales = 'free_y') +
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
    return(p1)
}

################################
# Load incidence rates samples #
################################

incidence_rates_round.samples <- load_incidence_rates_samples(file.incidence.samples.inland)
ps <- c(`M`=0.5, `IL`=0.025, `IU`=0.975)

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
        q_label=names(ps)),
    by=c('ROUND', 'COMM', 'SEX', 'AGEYRS')] |>
    dcast( ROUND+COMM+SEX+AGEYRS~q_label, value.var='q')

p_1d <- make_subplot(data_1d, M='M', IL='IL', IU='IU', y_lab = "Contribution to incidence rates")
# ggsave('~/Downloads/incidence_rates_contribution.pdf',p_1c, w=8, h=3.5)
# ggsave('~/Downloads/incidence_rates_contribution.png',p_1c, w=8, h=3.5)


###################
# INCIDENCE RATES # 
###################

library(patchwork)
naturemed_reqs()

reqs2 <- theme(
    axis.text = element_text(size=9, family='sans'),
    text=element_text(size=12,family='sans'),
    legend.text=element_text(size=12, family='sans'),
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 12)
)

p_c <- plot_incident_rates_over_time_2(incidence_cases_round, incidence_rates_round.samples, eligible_count_round)
p_1cd <- (p_c + reqs2 + theme(axis.title.x=element_blank()))/(p_1d + reqs2 + theme(strip.background = element_blank(), strip.text = element_blank()) ) 

# 900, 680 = 23.81, 18
ggsave('~/Downloads/MainFigure1cd.pdf',p_1cd, w=23.8 , h=18, units='cm')


#########################################
# EDF3: CONTRIBUTION TO INCIDENCE CASES #
#########################################

# load susceptibles by strata
dsusc <- incidence_cases_round[, .(COMM, ROUND, AGEYRS, SEX, SUSCEPTIBLE)]
data_edf3 <- merge( 
    incidence_rates_round.samples,
    dsusc,
    by=c('COMM', 'ROUND', 'AGEYRS', 'SEX'))

data_edf3 <- data_edf3[, {
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
        q_label=names(ps)),
    by=c('ROUND', 'COMM', 'SEX', 'AGEYRS')]|>
    dcast( ROUND+COMM+SEX+AGEYRS~q_label, value.var='q')

p_edf_contribution <- make_subplot(data_1d, M='M', IL='IL', IU='IU', y_lab = "Contribution to incident cases")

ggsave('~/Downloads/incidence_cases_contribution.pdf',p_edf_contribution+reqs2 , w=20, h=10, units='cm')
ggsave('~/Downloads/incidence_cases_contribution.png',p_edf_contribution+reqs2 , w=20, h=10, units='cm')
