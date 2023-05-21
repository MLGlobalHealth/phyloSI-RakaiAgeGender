library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

# directory to repository
gitdir <- getwd()

# check files exist
file.exists(c(
  file.treatment.cascade ,
  file.prevalence,
  file.eligible.count))  |> all() |> stopifnot()

# load census eligible ount
eligible_count <- fread(file.eligible.count)

# load proportion prevalence
proportion_prevalence <- as.data.table(readRDS(file.prevalence))

# load unsuppressed proportion
treatment_cascade <- as.data.table(readRDS(file.treatment.cascade))


#############################

# FIND INFECTED UNSUPPRESSED

############################

# define round
join_keys <- c("ROUND", "COMM", "AGEYRS", "SEX", "iterations")
df <- merge(proportion_prevalence, treatment_cascade, by = join_keys)
df[, ROUND := gsub("R0(.+)", "\\1", ROUND)]

# merge number of eligible and the prevalence
eligible_count <- eligible_count[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)]
df <- merge(eligible_count, df, by = c("ROUND", "COMM", "AGEYRS", "SEX"))

# find infected
df[, INFECTED := ELIGIBLE * PREVALENCE_POSTERIOR_SAMPLE]

# find infected unsuppressed
df[, PROP_UNSUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_POSTERIOR_SAMPLE]
df[, UNSUPPRESSED := INFECTED * PROP_UNSUPPRESSED_POSTERIOR_SAMPLE]


#####################################################

# UNSUPPRESSED BY 3 AGE GROUP

#####################################################

# aggregated by age groups
df[, INDEX_AGE_GROUP := 3]
df[AGEYRS < 35, INDEX_AGE_GROUP := 2]
df[AGEYRS < 25, INDEX_AGE_GROUP := 1]
df[, AGE_GROUP := c("15-24", "25-34", "35-49")[INDEX_AGE_GROUP]]
df_agg <- df[, list(INFECTED = sum(INFECTED),
                    UNSUPPRESSED = sum(UNSUPPRESSED),
                    ELIGIBLE = sum(ELIGIBLE)),
              by = c("ROUND", "COMM", "AGE_GROUP", "SEX", "iterations")]

# find unsuppressed rate
df_agg[, SUPPRESSED := INFECTED - UNSUPPRESSED]
df_agg[, SUPPRESSION_RATE := SUPPRESSED / INFECTED]
df_agg[, UNSUPPRESSION_RATE := UNSUPPRESSED / INFECTED]

# find unsuppressed rate relative to round 10
df_agg[COMM == "inland", UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == "10"],
       by = c("COMM", "AGE_GROUP", "SEX", "iterations")]
df_agg[COMM == "fishing", UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == "15"],
       by = c("COMM", "AGE_GROUP", "SEX", "iterations")]
df_agg[, UNSUPPRESSION_RATE_REL := UNSUPPRESSION_RATE / UNSUPPRESSION_RATE.REF]

# find ratio of art uptake male to female
df_agg <- dcast(
  data = df_agg,
  formula = ROUND + COMM + AGE_GROUP + iterations ~ SEX,
  value.var = "UNSUPPRESSION_RATE_REL"
)
setnames(
  x = df_agg,
  old = c("M", "F"),
  new = c("UNSUPPRESSION_RATE_REL_M", "UNSUPPRESSION_RATE_REL_F")
)
df_agg[, UNSUPPRESSION_RATE_RATIO := UNSUPPRESSION_RATE_REL_M / UNSUPPRESSION_RATE_REL_F ]

# summarise
ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")
sing_age <- df_agg[, list(q = quantile(UNSUPPRESSION_RATE_RATIO, prob = ps, na.rm = TRUE),
                         q_label = qlab),
                   by = c("ROUND", "COMM", "AGE_GROUP")]
sing_age <- as.data.table(
  reshape2::dcast(sing_age, ... ~ q_label, value.var = "q")
)

# name
setnames(sing_age, qlab, paste0("UNSUPPRESSION_RATE_RATIO_BY_AGE_", qlab))

# plot
ggplot(sing_age, aes(x = ROUND)) +
  geom_point(aes(y = UNSUPPRESSION_RATE_RATIO_BY_AGE_M)) +
  geom_errorbar(aes(ymin = UNSUPPRESSION_RATE_RATIO_BY_AGE_CL,
                    ymax = UNSUPPRESSION_RATE_RATIO_BY_AGE_CU),
                alpha = 0.5) +
  facet_grid(COMM ~ AGE_GROUP) +
  theme_bw()

#####################################################

# FIND SEX RATIO OF INFECTED UNSUPPRESSED ACROSS AGE

#####################################################

df_agg <- df[, list(INFECTED = sum(INFECTED),
                    UNSUPPRESSED = sum(UNSUPPRESSED),
                    ELIGIBLE = sum(ELIGIBLE)),
              by = c("ROUND", "COMM", "SEX", "iterations")]

# find unsuppressed rate
df_agg[, SUPPRESSED := INFECTED - UNSUPPRESSED]
df_agg[, SUPPRESSION_RATE := SUPPRESSED / INFECTED]
df_agg[, UNSUPPRESSION_RATE := UNSUPPRESSED / INFECTED]

# find unsuppressed rate relative to round 10
df_agg[COMM == "inland", UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == "10"],
       by = c("COMM", "SEX", "iterations")]
df_agg[COMM == "fishing", UNSUPPRESSION_RATE.REF := UNSUPPRESSION_RATE[ROUND == "15"],
       by = c("COMM", "SEX", "iterations")]
df_agg[, UNSUPPRESSION_RATE_REL := UNSUPPRESSION_RATE / UNSUPPRESSION_RATE.REF]

# find ratio of art uptake male to female
df_agg <- dcast(
  data = df_agg,
  formula = ROUND + COMM + iterations ~ SEX,
  value.var = "UNSUPPRESSION_RATE_REL"
)
setnames(
  x = df_agg,
  old = c("M", "F"),
  new = c("UNSUPPRESSION_RATE_REL_M", "UNSUPPRESSION_RATE_REL_F")
)
df_agg[, UNSUPPRESSION_RATE_RATIO := UNSUPPRESSION_RATE_REL_M / UNSUPPRESSION_RATE_REL_F]

# summarise
sing <- df_agg[, list(q = quantile(UNSUPPRESSION_RATE_RATIO, prob = ps, na.rm = TRUE),
                      q_label = qlab),
                by = c("ROUND", "COMM")]
sing <- as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

# name
setnames(sing, qlab, paste0("UNSUPPRESSION_RATE_RATIO_RATIO_", qlab))

# plot
ggplot(sing, aes(x = ROUND)) +
  geom_point(aes(y = UNSUPPRESSION_RATE_RATIO_RATIO_M)) +
  geom_errorbar(aes(ymin = UNSUPPRESSION_RATE_RATIO_RATIO_CL,
                ymax = UNSUPPRESSION_RATE_RATIO_RATIO_CU),
                alpha = 0.5) +
  facet_grid(COMM ~ .) +
  theme_bw()

#########################################

# SAVE

#########################################

tmp <- merge(sing_age, sing, by = c("ROUND", "COMM"))
file.name <- file.unsuppressed_rate_ratio
if (! file.exists(file.name) || config$overwrite.existing.files) {
    cat("Saving file:", file.name, "\n")
    write.csv(tmp, file = file.name, row.names = FALSE)
} else {
    cat("File:", file.name, "already exists...\n")
}
