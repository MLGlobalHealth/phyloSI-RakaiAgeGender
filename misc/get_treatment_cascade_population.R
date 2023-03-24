library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "paths.R"))

# TODO: shozen: do you think this would be helpful? 
# library(optparse)
# option_list <- list(
#     make_option(
#         "--outdir",
#         type = "",
#         default = ,
#         help = "",
#         dest= ""
#     ),
# )
# args <- parse_args(OptionParser(option_list = option_list))


# outdir to save figures
if (dir.exists(indir.deepsequence_analyses)) {
  outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')
} else {
  outdir <- 'outputs'
  if (!dir.exists(outdir)) dir.create(outdir);
}

# posterior samples participants
file.exists(file.unsuppressedviralload) |> stopifnot()
file.exists(file.selfreportedart) |> stopifnot()
file.exists(file.unsuppressedviralload.np) |> stopifnot()
file.exists(file.selfreportedart.np) |> stopifnot()
file.exists(file.participation) |> stopifnot()
file.exists(file.spec.sens.art) |> stopifnot()

ps <- c(0.025,0.5,0.975)
qlab <- c('CL','M','CU')


###########################

# COMBINE POSTERIOR SAMPLE

###########################

# load proportion unsuppressed viral load 
uns.p <- as.data.table(readRDS(file.unsuppressedviralload))
uns.np <- as.data.table(readRDS(file.unsuppressedviralload.newly))

# load proportion art
sre.p <- as.data.table(readRDS(file.selfreportedart))
sre.np <- as.data.table(readRDS(file.selfreportedart.newly))

# add label
setnames(uns.p, 'PROP_UNSUPPRESSED_POSTERIOR_SAMPLE', 'PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE')
setnames(uns.np, 'PROP_UNSUPPRESSED_POSTERIOR_SAMPLE', 'PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE')
setnames(sre.p, 'PROP_ART_COVERAGE_POSTERIOR_SAMPLE', 'PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE')
setnames(sre.np, 'PROP_ART_COVERAGE_POSTERIOR_SAMPLE', 'PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE')

# remove crude estimates
set(uns.p, NULL, 'PROP_UNSUPPRESSED_EMPIRICAL', NULL)
set(uns.np, NULL, 'PROP_UNSUPPRESSED_EMPIRICAL', NULL)
set(sre.p, NULL, 'PROP_ART_COVERAGE_EMPIRICAL', NULL)
set(sre.np, NULL, 'PROP_ART_COVERAGE_EMPIRICAL', NULL)

# merge participants (fyi for round < 15, we do not have viral load info)
df.p <- merge(uns.p, sre.p, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'), all.y = T)

# merge non-participants (fyi for round < 15, we do not have viral load info)
df.np <- merge(uns.np, sre.np, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'), all.y = T)

# merge population
df <- merge(df.p, df.np, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'))

# for round 10, set art coverage and unsuppressed as round 11 as questions art related q were introduced from roudn 11
tmp <- df[ROUND == 'R011']
tmp[, ROUND := 'R010']
df <- rbind(tmp, df[ROUND != 'R010'])

# restrict age 
df <- df[AGEYRS > 14 & AGEYRS < 50]

# add participation
participation <- as.data.table(read.csv(file.participation))
participation[, ROUND := paste0('R0', ROUND)]
df <- merge(df, participation, by = c('AGEYRS', 'SEX', 'COMM', 'ROUND'))

# find proportion suppressed
df[, PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]
df[, PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# find suppression rate
df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE / PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE]
df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE / PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# add constraint that suppression rate must at maximum 1
df[SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE > 1, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := 1]
df[SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1]


####################################

# FIND PROPORTION OF ART USE GIVEN INFECTED FOR ROUND 15

####################################

# for round 15 find % on art by using the same suppression rate as round 16
df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE := min(1, PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE / SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE), by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations')]
set(df, NULL, 'SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)

df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
df[ROUND %in% c('R015'), PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE := min(1, PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE / SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE), by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations')]
set(df, NULL, 'SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)

# add constraint that art coverage must be smaller than prop suppression (can occur if prop suppression > art coverage)
df[PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE < PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE, PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]
df[PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE < PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE, PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]


##########################################################################

# FIND PROP SUPPRESSED AMONG HIV POSITIVE PARTICIPANTS FOR ROUNDS WITHOUT VIRAL LOAD MEASUREMENTS

##########################################################################

if(1){
  
  #
  #  USE SENSITIVITY AND SPECIFICITY IN ROUND 15
  #
  
  # load file
  sensitivity_specificity_art <- as.data.table(read.csv(file.spec.sens.art))
  
  # use specificity and sensitivity from round 15 
  spa <- sensitivity_specificity_art[ROUND == 'R015' ]
  
  # select variable of interest
  spa <- spa[, .(SEX, AGEYRS, SPEC_M, SENS_M)]
  
  # merge
  df <- merge(df, spa, by = c('SEX', 'AGEYRS'))
  
  df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE * SPEC_M + (1 - PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE) * (1-SENS_M)]
  df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE * SPEC_M + (1 - PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE) * (1-SENS_M)]
  set(df, NULL, 'SPEC_M', NULL)
  set(df, NULL, 'SENS_M', NULL)

}else{
  
  #
  #  USE PROP SUPPRESSION GIVEN ART UPTAKE IN ROUND 16
  #
  
  # for round <15 find % suppressed by using the same suppression rate as round 16
  df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
  df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
  df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE * SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
  set(df, NULL, 'SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)
  
  df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE[ROUND == 'R016'], by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
  df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
  df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE * SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE, by = c('AGEYRS', 'SEX', 'COMM', 'iterations')]
  set(df, NULL, 'SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016', NULL)
  
}


df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]
df[ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'), PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]

stopifnot(nrow(df[!ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'),SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)
stopifnot(nrow(df[!ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'),SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)


######################################################################################

# ART COVERERAGE and prop suppressed IN POPULATION

######################################################################################

# art coverage in population
stopifnot(nrow(df[PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)
stopifnot(nrow(df[PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)

df[, PROP_ART_COVERAGE_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE * PARTICIPATION + PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE * (1-PARTICIPATION)]

# find suppression in the population
stopifnot(nrow(df[PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)
stopifnot(nrow(df[PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)

df[, PROP_SUPPRESSED_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE * PARTICIPATION + PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE * (1-PARTICIPATION)]


####################################

# FIND PROPORTION OF DIAGNOSED

####################################

# all participants are diagnosed
df[, PROP_DIAGNOSED_PARTICIPANT_POSTERIOR_SAMPLE := 1]

# assume that proportion diagnosed in non -participants = proportion on art
df[, PROP_DIAGNOSED_NONPARTICIPANT_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# find proportion diagnosed in the population
df[, PROP_DIAGNOSED_POSTERIOR_SAMPLE := PROP_DIAGNOSED_PARTICIPANT_POSTERIOR_SAMPLE * PARTICIPATION + PROP_DIAGNOSED_NONPARTICIPANT_POSTERIOR_SAMPLE * (1-PARTICIPATION)]


##########################################################

# FIND PROPORTION OF ON ART AND SUPPRESSED GIVEN DIAGNOSED

##########################################################

# find art coverage given diagnosed in population
df[, PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE := min(1, PROP_ART_COVERAGE_POSTERIOR_SAMPLE / PROP_DIAGNOSED_POSTERIOR_SAMPLE), by = c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations')]
stopifnot(nrow(df[PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE>1]) == 0)

# find suppression given art use in the population (i.e., suppression rate)
df[, PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE / PROP_ART_COVERAGE_POSTERIOR_SAMPLE ]
stopifnot(nrow(df[!ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'),PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE > 1]) == 0)

# find suppression given diagnosed use in the population 
df[, PROP_SUPPRESSION_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE / PROP_DIAGNOSED_POSTERIOR_SAMPLE ]
stopifnot(nrow(df[!ROUND %in% c('R010', 'R011', 'R012', 'R013', 'R014', 'R015S'),PROP_SUPPRESSION_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE > 1]) == 0)


####################################

# SUMMARISE

####################################

tmp <- df[, .(ROUND, SEX, COMM, AGEYRS, iterations, 
              PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE, 
              PROP_SUPPRESSION_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE,
              PROP_SUPPRESSED_POSTERIOR_SAMPLE,
              PROP_ART_COVERAGE_POSTERIOR_SAMPLE,
             PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE, 
              PROP_DIAGNOSED_POSTERIOR_SAMPLE)]
tmp <- melt.data.table(tmp, id.vars= c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'iterations'))
tmp[, variable := gsub('(.+)_POSTERIOR_SAMPLE', '\\1', variable)]
ns = tmp[, list(q= quantile(value, prob=ps, na.rm = T), q_label=qlab), by=c('AGEYRS', 'SEX', 'COMM', 'ROUND', 'variable')]
ns = as.data.table(reshape2::dcast(ns, AGEYRS + SEX + COMM + ROUND + variable ~ q_label, value.var = "q"))

# check all entries are complete
stopifnot(nrow(ns[COMM == 'inland']) == ns[, length(unique(AGEYRS))] * ns[, length(unique(SEX))] * ns[, length(unique(variable))] * ns[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(ns[COMM == 'fishing']) == ns[, length(unique(AGEYRS))] * ns[, length(unique(SEX))] * ns[, length(unique(variable))] * ns[COMM == 'fishing', length(unique(ROUND))])
stopifnot(nrow(df[COMM == 'inland']) == df[, length(unique(AGEYRS))] * df[, length(unique(SEX))] * df[, length(unique(iterations))] * df[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(df[COMM == 'fishing']) == df[, length(unique(AGEYRS))] * df[, length(unique(SEX))] * df[, length(unique(iterations))] * df[COMM == 'fishing', length(unique(ROUND))])


####################################

# SAVE

####################################

# file.name <- file.path(gitdir.fit, paste0('RCCS_treatment_cascade_population_posterior_samples_221208.rds'))
file.name <- file.treatment.cascade 
if(! file.exists(file.name))
{
    cat("Saving file:", file.name, '\n')
    saveRDS(df, file = file.name)
}else{
    cat("File:", file.name, "already exists...\n")
}

file.name <- file.path(gitdir.fit, paste0('RCCS_treatment_cascade_population_estimates_221208.csv'))
if(! file.exists(file.name))
{
    cat("Saving file:", file.name, '\n')
    write.csv(ns, file = file.name, row.names = F)
}else{
    cat("File:", file.name, "already exists...\n")
}

