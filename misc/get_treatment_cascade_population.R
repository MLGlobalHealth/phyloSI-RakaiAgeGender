library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

# outdir to save figures
# outdir to save figures
outdir <- file.path(
  "../phyloSI-RakaiAgeGender-outputs",
  "get_treatment_cascade_participants"
)
if (!dir.exists(outdir)) dir.create(outdir)

# posterior samples participants
file.exists(file.unsuppressedviralload) |> stopifnot()
file.exists(file.selfreportedart) |> stopifnot()
file.exists(file.unsuppressedviralload.newly) |> stopifnot()
file.exists(file.selfreportedart.newly) |> stopifnot()
file.exists(file.participation) |> stopifnot()
file.exists(file.spec.sens.art) |> stopifnot()

ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")

###########################

# COMBINE POSTERIOR SAMPLE

###########################

# load proportion unsuppressed viral load
uns_p <- as.data.table(readRDS(file.unsuppressedviralload))
uns_np <- as.data.table(readRDS(file.unsuppressedviralload.newly))

# load proportion art
sre_p <- as.data.table(readRDS(file.selfreportedart))
sre_np <- as.data.table(readRDS(file.selfreportedart.newly))

# add label
setnames(uns_p, 
         "PROP_UNSUPPRESSED_POSTERIOR_SAMPLE",
         "PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE")
setnames(uns_np,
         "PROP_UNSUPPRESSED_POSTERIOR_SAMPLE",
         "PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE")
setnames(sre_p,
         "PROP_ART_COVERAGE_POSTERIOR_SAMPLE",
         "PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE")
setnames(sre_np,
         "PROP_ART_COVERAGE_POSTERIOR_SAMPLE",
         "PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE")

# remove crude estimates
set(uns_p, NULL, "PROP_UNSUPPRESSED_EMPIRICAL", NULL)
set(uns_np, NULL, "PROP_UNSUPPRESSED_EMPIRICAL", NULL)
set(sre_p, NULL, "PROP_ART_COVERAGE_EMPIRICAL", NULL)
set(sre_np, NULL, "PROP_ART_COVERAGE_EMPIRICAL", NULL)

# merge participants (fyi for round < 15, we do not have viral load info)
groupby_keys <- c("AGEYRS", "SEX", "COMM", "ROUND", "iterations")
df_p <- merge(uns_p, sre_p, by = groupby_keys, all.y = TRUE)

# merge non-participants (fyi for round < 15, we do not have viral load info)
df_np <- merge(uns_np, sre_np, by = groupby_keys, all.y = TRUE)

# merge population
df <- merge(df_p, df_np, by = groupby_keys)

# for round 10, set art coverage and unsuppressed as
# round 11 as questions art related q were introduced from roudn 11
tmp <- df[ROUND == "R011"]
tmp[, ROUND := "R010"]
df <- rbind(tmp, df[ROUND != "R010"])

# restrict age
df <- df[AGEYRS > 14 & AGEYRS < 50]

# add participation
participation <- as.data.table(read.csv(file.participation))
participation[, ROUND := paste0("R0", ROUND)]
df <- merge(df, participation, by = c("AGEYRS", "SEX", "COMM", "ROUND"))

# find proportion suppressed
df[, PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]
df[, PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# find suppression rate
df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE /
                                                       PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE]
df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE /
                                                          PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# add constraint that suppression rate must at maximum 1
df[SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE > 1,
   SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := 1]
df[SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1,
   SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1]


####################################

# FIND PROPORTION OF ART USE GIVEN INFECTED FOR ROUND 15

####################################

# for round 15 find % on art by using the same suppression rate as round 16
df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE[ROUND == "R016"],
    by = c("AGEYRS", "SEX", "COMM", "iterations")]
df[ROUND == "R015", SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016,
    by = c("AGEYRS", "SEX", "COMM", "iterations")]
df[ROUND == "R015", PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE := min(1, PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE / SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE),
    by = c("AGEYRS", "SEX", "COMM", "ROUND", "iterations")]
set(df, NULL, "SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016", NULL)

df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE[ROUND == "R016"],
    by = c("AGEYRS", "SEX", "COMM", "iterations")]
df[ROUND == "R015", SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016,
    by = c("AGEYRS", "SEX", "COMM", "iterations")]
df[ROUND == "R015", PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE := min(1, PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE / SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE),
    by = c("AGEYRS", "SEX", "COMM", "ROUND", "iterations")]
set(df, NULL, "SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016", NULL)

# add constraint that art coverage must be smaller 
# than prop suppression (can occur if prop suppression > art coverage)
df[PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE < PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE,
   PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]
df[PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE < PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE,
   PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]


##########################################################################

# FIND PROP SUPPRESSED AMONG HIV POSITIVE 
# PARTICIPANTS FOR ROUNDS WITHOUT VIRAL LOAD MEASUREMENTS

##########################################################################

if (TRUE) {
  #
  #  USE SENSITIVITY AND SPECIFICITY IN ROUND 15
  #
  # load file
  sensitivity_specificity_art <- as.data.table(read.csv(file.spec.sens.art))
  # use specificity and sensitivity from round 15
  spa <- sensitivity_specificity_art[ROUND == "R015"]

  # select variable of interest
  spa <- spa[, .(SEX, AGEYRS, SPEC_M, SENS_M)]

  # merge
  df <- merge(df, spa, by = c("SEX", "AGEYRS"))
  subset_rounds <- c("R010", "R011", "R012", "R013", "R014", "R015S")
  df[ROUND %in% subset_rounds,
     PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE *
                                                      SPEC_M +
                                                      (1 - PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE) *
                                                      (1 - SENS_M)]
  df[ROUND %in% subset_rounds,
     PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE *
                                                         SPEC_M +
                                                         (1 - PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE) *
                                                         (1 - SENS_M)]
  set(df, NULL, "SPEC_M", NULL)
  set(df, NULL, "SENS_M", NULL)

} else {
  #
  #  USE PROP SUPPRESSION GIVEN ART UPTAKE IN ROUND 16
  #
  # for round <15 find % suppressed by
  # using the same suppression rate as round 16
  subset_rounds <- c("R010", "R011", "R012", "R013", "R014", "R015S")
  groupby_keys <- c("AGEYRS", "SEX", "COMM", "iterations")
  df[, SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE[ROUND == "R016"],
      by = groupby_keys]
  df[ROUND %in% subset_rounds,
     SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016,
      by = groupby_keys]
  df[ROUND %in% subset_rounds,
     PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE *
                                                      SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE,
     by = groupby_keys]
  set(df, NULL, "SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE_R016", NULL)

  df[, SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE[ROUND == "R016"],
      by = groupby_keys]
  df[ROUND %in% subset_rounds,
     SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE := SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016,
     by = groupby_keys]
  df[ROUND %in% subset_rounds,
     PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE *
                                                         SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE,
     by = groupby_keys]
  set(df, NULL, "SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE_R016", NULL)
}

subset_rounds <- c("R010", "R011", "R012", "R013", "R014", "R015S")
groupby_keys <- c("AGEYRS", "SEX", "COMM", "iterations")
df[ROUND %in% subset_rounds,
   PROP_UNSUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE]
df[ROUND %in% subset_rounds,
   PROP_UNSUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE]

stopifnot(
  nrow(df[!ROUND %in% subset_rounds,
          SUPPRESSION_RATE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0
)
stopifnot(
  nrow(df[!ROUND %in% subset_rounds,
          SUPPRESSION_RATE_PARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0
)


########################################

# ART COVERERAGE and prop suppressed IN POPULATION

########################################

# art coverage in population
stopifnot(nrow(df[PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)
stopifnot(nrow(df[PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)

df[, PROP_ART_COVERAGE_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_PARTICIPANTS_POSTERIOR_SAMPLE *
                                           PARTICIPATION +
                                           PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE *
                                           (1 - PARTICIPATION)]

# find suppression in the population
stopifnot(nrow(df[PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)
stopifnot(nrow(df[PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE > 1]) == 0)

df[, PROP_SUPPRESSED_POSTERIOR_SAMPLE := PROP_SUPPRESSED_PARTICIPANTS_POSTERIOR_SAMPLE *
                                         PARTICIPATION +
                                         PROP_SUPPRESSED_NONPARTICIPANTS_POSTERIOR_SAMPLE *
                                         (1 - PARTICIPATION)]


####################################

# FIND PROPORTION OF DIAGNOSED

####################################

# all participants are diagnosed
df[, PROP_DIAGNOSED_PARTICIPANT_POSTERIOR_SAMPLE := 1]

# assume that proportion diagnosed in non -participants = proportion on art
df[, PROP_DIAGNOSED_NONPARTICIPANT_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_NONPARTICIPANTS_POSTERIOR_SAMPLE]

# find proportion diagnosed in the population
df[, PROP_DIAGNOSED_POSTERIOR_SAMPLE := PROP_DIAGNOSED_PARTICIPANT_POSTERIOR_SAMPLE *
                                        PARTICIPATION +
                                        PROP_DIAGNOSED_NONPARTICIPANT_POSTERIOR_SAMPLE *
                                        (1 - PARTICIPATION)]


##########################################################

# FIND PROPORTION OF ON ART AND SUPPRESSED GIVEN DIAGNOSED

##########################################################

# find art coverage given diagnosed in population
df[, PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE := min(1, PROP_ART_COVERAGE_POSTERIOR_SAMPLE / PROP_DIAGNOSED_POSTERIOR_SAMPLE),
    by = c("AGEYRS", "SEX", "COMM", "ROUND", "iterations")]
stopifnot(nrow(df[PROP_ART_COVERAGE_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE > 1]) == 0)

# find suppression given art use in the population (i.e., suppression rate)
df[, PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE /
                                                    PROP_ART_COVERAGE_POSTERIOR_SAMPLE]

subset_rounds <- c("R010", "R011", "R012", "R013", "R014", "R015S")
stopifnot(
  nrow(df[!ROUND %in% subset_rounds,
       PROP_SUPPRESSION_GIVEN_ART_POSTERIOR_SAMPLE > 1]) == 0
)

# find suppression given diagnosed use in the population
df[, PROP_SUPPRESSION_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE /
                                                          PROP_DIAGNOSED_POSTERIOR_SAMPLE]
stopifnot(
  nrow(df[!ROUND %in% subset_rounds,
       PROP_SUPPRESSION_GIVEN_DIAGNOSED_POSTERIOR_SAMPLE > 1]) == 0
)


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
tmp <- melt.data.table(
  data = tmp,
  id.vars = c("AGEYRS", "SEX", "COMM", "ROUND", "iterations")
)
tmp[, variable := gsub("(.+)_POSTERIOR_SAMPLE", "\\1", variable)]
ns <- tmp[, list(q = quantile(value, prob = ps, na.rm = TRUE), q_label = qlab),
          by = c("AGEYRS", "SEX", "COMM", "ROUND", "variable")]
ns <- as.data.table(
  reshape2::dcast(
    data = ns,
    formula = AGEYRS + SEX + COMM + ROUND + variable ~ q_label,
    value.var = "q"
  )
)

# check all entries are complete
stopifnot(
  nrow(ns[COMM == "inland"]) == ns[, uniqueN(AGEYRS) *
                                uniqueN(SEX) *
                                uniqueN(variable)] *
                                ns[COMM == "inland", uniqueN(ROUND)]
)
stopifnot(
  nrow(ns[COMM == "fishing"]) == ns[, uniqueN(AGEYRS) *
                                      uniqueN(SEX) *
                                      uniqueN(variable)] *
                                      ns[COMM == "fishing", uniqueN(ROUND)]
)
stopifnot(
  nrow(df[COMM == "inland"]) == df[, uniqueN(AGEYRS) *
                                     uniqueN(SEX) *
                                     uniqueN(iterations)] *
                                     df[COMM == "inland", uniqueN(ROUND)]
)
stopifnot(
  nrow(df[COMM == "fishing"]) == df[, uniqueN(AGEYRS) *
                                      uniqueN(SEX) *
                                      uniqueN(iterations)] *
                                      df[COMM == "fishing", uniqueN(ROUND)]
)


####################################

# SAVE

####################################

file_name <- file.treatment.cascade
if (! file.exists(file_name) || config$overwrite.existing.files) {
  cat("Saving file:", file_name, "\n")
  saveRDS(df, file = file_name)
} else {
  cat("File:", file_name, "already exists...\n")
}

file_name <- file.treatment.cascade.population
if (! file.exists(file_name) || config$overwrite.existing.files) {
  cat("Saving file:", file_name, "\n")
  write.csv(ns, file = file_name, row.names = FALSE)
} else {
  cat("File:", file_name, "already exists...\n")
}