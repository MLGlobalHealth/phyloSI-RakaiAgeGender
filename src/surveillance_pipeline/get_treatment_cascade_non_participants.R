library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

# outdir directory for stan fit
outdir <- file.path("../phyloSI-RakaiAgeGender-outputs","get_treatment_cascade_non_participants")
if(usr == 'melodiemonod'){
  outdir <- file.path("/Users/melodiemonod/Box Sync/2023//phyloSI-RakaiAgeGender-outputs","get_treatment_cascade_non_participants")
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Specificity and sensitivity ART reporting
ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")

# check files exist
file.exists(c(
  file.unsuppressedviralload.newly ,
  file.selfreportedart.newly,
  file.spec.sens.art))  |> all() |> stopifnot()


###########################

# COMBINE POSTERIOR SAMPLE

###########################

# load proportion unsuppressed viral load
uns <- as.data.table(readRDS(file.unsuppressedviralload.newly))

# load proportion art
sre <- as.data.table(readRDS(file.selfreportedart.newly))

# merge (fyi for round < 15, we do not have viral load info)
join_keys <- c("AGEYRS", "SEX", "COMM", "ROUND", "iterations")
df <- merge(uns, sre, by = join_keys, all.y = TRUE)

# for round 10, set art coverage and unsuppressed as
# round 11 as questions art related q were introduced from roudn 11
tmp <- df[ROUND == "R011"]
tmp[, ROUND := "R010"]
df <- rbind(tmp, df[ROUND != "R010"])

# restrict age
df <- df[AGEYRS > 14 & AGEYRS < 50]

# find proportion suppressed among new registed HIV-positive
df[, PROP_SUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_UNSUPPRESSED_POSTERIOR_SAMPLE]

# find suppression rate  among new registed HIV-positive
df[, SUPPRESSION_RATE_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE /
                                          PROP_ART_COVERAGE_POSTERIOR_SAMPLE]
df[, SUPPRESSION_RATE_EMPIRICAL := (1 - PROP_UNSUPPRESSED_EMPIRICAL) /
                                    PROP_ART_COVERAGE_EMPIRICAL]

# add constraint that suppression rate must at maximum 1
df[SUPPRESSION_RATE_POSTERIOR_SAMPLE > 1,
   SUPPRESSION_RATE_POSTERIOR_SAMPLE := 1]


##############################################################################

# FIND ART COVERERAGE AMONG NEWLY REGISTERED HIV-POSITIVE  fo round 15

##############################################################################

# for round 15 find % on art by using the same suppression rate as round 16
group_keys <- c("AGEYRS", "SEX", "COMM", "iterations")
df[, SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_POSTERIOR_SAMPLE[ROUND == "R016"],
    by = group_keys]
df[, SUPPRESSION_RATE_EMPIRICAL_R016 := SUPPRESSION_RATE_EMPIRICAL[ROUND == "R016"],
    by = group_keys]
df[ROUND == "R015",
   SUPPRESSION_RATE_POSTERIOR_SAMPLE := SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016,
    by = group_keys]
df[ROUND == "R015",
   SUPPRESSION_RATE_EMPIRICAL := SUPPRESSION_RATE_EMPIRICAL_R016,
   by = group_keys]

# find art coverage with the same suppression rate as round 16
df[ROUND == "R015",
   PROP_ART_COVERAGE_POSTERIOR_SAMPLE := min(1, PROP_SUPPRESSED_POSTERIOR_SAMPLE / SUPPRESSION_RATE_POSTERIOR_SAMPLE), # nolint: line_length_linter.
   by = c("AGEYRS", "ROUND", "SEX", "COMM", "iterations")]
set(df, NULL, "SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016", NULL)
set(df, NULL, "SUPPRESSION_RATE_EMPIRICAL_R016", NULL)

# add constraint that art coverage must be smaller
# than prop suppression (can occur if prop suppression > art coverage)
df[PROP_ART_COVERAGE_POSTERIOR_SAMPLE < PROP_SUPPRESSED_POSTERIOR_SAMPLE,
   PROP_ART_COVERAGE_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE]
stopifnot(nrow(df[PROP_ART_COVERAGE_POSTERIOR_SAMPLE > 1]) == 0)

####################################

# FIND DIAGNOSED PROPORTION AMONG NEWLY REGISTERED HIV-POSITIVE

####################################

# assume that proportion diagnosed in non -participants = proportion on art
df[, PROP_DIAGNOSED_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_POSTERIOR_SAMPLE]


##########################################################################

# FIND PROP SUPPRESSED AMONG HIV POSITIVE
# PARTICIPANTS FOR ROUNDS WITHOUT VIRAL LOAD MEASUREMENTS

##########################################################################

if (TRUE) {
  #
  #  USE SENSITIVITY AND SPECIFICITY IN ROUND 15
  #
  # load file
  sensitivity_specificity_art <- fread(file.spec.sens.art)

  # use specificity and sensitivity from round 15
  spa <- sensitivity_specificity_art[ROUND == "R015"]

  # select variable of interest
  spa <- spa[, .(SEX, AGEYRS, SPEC_M, SENS_M)]

  # merge
  df <- merge(df, spa, by = c("SEX", "AGEYRS"))

  df[ROUND %in% c("R010", "R011", "R012", "R013", "R014", "R015S"),
     PROP_SUPPRESSED_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_POSTERIOR_SAMPLE *
                                         SPEC_M + (1 - PROP_ART_COVERAGE_POSTERIOR_SAMPLE) *
                                         (1-SENS_M)]
  set(df, NULL, "SPEC_M", NULL)
  set(df, NULL, "SENS_M", NULL)

} else {
  #
  #  USE PROP SUPPRESSION GIVEN ART UPTAKE IN ROUND 16
  #
  # for round <15 find % suppressed by using
  # the same suppression rate as round 16
  group_keys <- c("AGEYRS", "SEX", "COMM", "iterations")
  df[, SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016 := SUPPRESSION_RATE_POSTERIOR_SAMPLE[ROUND == "R016"],
     by = group_keys]
  df[ROUND %in% c("R010", "R011", "R012", "R013", "R014", "R015S"),
     SUPPRESSION_RATE_POSTERIOR_SAMPLE := SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016,
     by = group_keys]
  df[ROUND %in% c("R010", "R011", "R012", "R013", "R014", "R015S"),
     PROP_SUPPRESSED_POSTERIOR_SAMPLE := PROP_ART_COVERAGE_POSTERIOR_SAMPLE *
                                         SUPPRESSION_RATE_POSTERIOR_SAMPLE,
     by = group_keys]

  set(df, NULL, "SUPPRESSION_RATE_POSTERIOR_SAMPLE_R016", NULL)
}

##########################################################################

# FIND PROP UNSUPPRESSED AND SUPPRESSION RATE

##########################################################################

# find proportion suppressed among HIV positive participants
df[, PROP_UNSUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_POSTERIOR_SAMPLE]

# find prop suppressed given diagnosed among art participants
df[, SUPPRESSION_RATE_POSTERIOR_SAMPLE := PROP_SUPPRESSED_POSTERIOR_SAMPLE /
                                          PROP_ART_COVERAGE_POSTERIOR_SAMPLE]

stopifnot(
  nrow(
    df[!ROUND %in% c("R010", "R011", "R012", "R013", "R014", "R015S"),
        SUPPRESSION_RATE_POSTERIOR_SAMPLE > 1]) == 0
)


####################################

# SUMMARISE

####################################

# melt
df1 <- melt.data.table(
  df,
  id.vars = c("AGEYRS", "SEX", "COMM", "ROUND", "iterations",
              "PROP_UNSUPPRESSED_EMPIRICAL",
              "PROP_ART_COVERAGE_EMPIRICAL",
              "SUPPRESSION_RATE_EMPIRICAL")
)
df1[, variable := gsub("(.+)_POSTERIOR_SAMPLE", "\\1", variable)]

# round the empirical otherwise we get two entries
df1[, PROP_UNSUPPRESSED_EMPIRICAL := round(PROP_UNSUPPRESSED_EMPIRICAL, 7)]
df1[, PROP_ART_COVERAGE_EMPIRICAL := round(PROP_ART_COVERAGE_EMPIRICAL, 7)]
df1[, SUPPRESSION_RATE_EMPIRICAL := round(SUPPRESSION_RATE_EMPIRICAL, 7)]

# summarise
group_keys <- c("AGEYRS", "SEX", "COMM", "ROUND",
                "PROP_UNSUPPRESSED_EMPIRICAL",
                "PROP_ART_COVERAGE_EMPIRICAL",
                "SUPPRESSION_RATE_EMPIRICAL",
                "variable")
ns <- df1[, list(q = quantile(value, prob = ps, na.rm = TRUE),
                q_label = paste0(variable, "_", qlab)),
            by = group_keys]
ns <- as.data.table(
  reshape2::dcast(
    data = ns,
    formula = AGEYRS + SEX + COMM + ROUND +
              PROP_UNSUPPRESSED_EMPIRICAL +
              PROP_ART_COVERAGE_EMPIRICAL +
              SUPPRESSION_RATE_EMPIRICAL ~ q_label,
    value.var = "q"
  )
)

# check all entries are complete
stopifnot(
  nrow(ns[COMM == "inland"]) == ns[, length(unique(AGEYRS))] *
                                ns[, length(unique(SEX))] *
                                ns[COMM == "inland", length(unique(ROUND))]
  )
stopifnot(
  nrow(ns[COMM == "fishing"]) == ns[, length(unique(AGEYRS))] *
                                 ns[, length(unique(SEX))] *
                                 ns[COMM == "fishing", length(unique(ROUND))]
)


####################################

# SAVE

####################################

file_name <- file.treatment.cascade.prop.nonparticipants.samples
if (!file.exists(file_name) || config$overwrite.existing.files) {
  cat("Saving file:", file_name, "\n")
  saveRDS(df, file = file_name)
} else {
  cat("File:", file_name, "already exists...\n")
}

file_name <- file.treatment.cascade.prop.nonparticipants
if (! file.exists(file_name) || config$overwrite.existing.files) {
  cat("Saving file:", file_name, "\n")
  write.csv(ns, file = file_name, row.names = FALSE)
} else {
  cat("File:", file_name, "already exists...\n")
}

