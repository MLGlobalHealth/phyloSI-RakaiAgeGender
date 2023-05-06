library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

outdir <- file.path(
  "../phyloSI-RakaiAgeGender-outputs",
  "get_unsuppressed_share_sex"
)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# files
file.exists(file.treatment.cascade) |> stopifnot()
file.exists(file.prevalence) |> stopifnot()
file.exists(file.eligible.count) |> stopifnot()

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
df <- merge(
  proportion_prevalence,
  treatment_cascade,
  by = c("ROUND", "COMM", "AGEYRS", "SEX", "iterations")
)
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

# FIND SEX SHARE OF INFECTED UNSUPPRESSED ACROSS AGE

#####################################################

# find share of unsuppressed by sex across age
df[, TOTAL_UNSUPPRESSED := sum(UNSUPPRESSED),
    by = c("ROUND", "COMM", "iterations")]
df[, UNSUPPRESSED_SHARE := UNSUPPRESSED / TOTAL_UNSUPPRESSED,
    by = c("ROUND", "COMM", "AGEYRS", "iterations")]

# summarise
ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")
sing_age <- df[, list(q = quantile(UNSUPPRESSED_SHARE, prob = ps, na.rm = TRUE),
                      q_label = qlab),
              by = c("ROUND", "COMM", "SEX", "AGEYRS")]
sing_age <- as.data.table(
  reshape2::dcast(sing_age, ... ~ q_label, value.var = "q")
)

# name
setnames(sing_age, qlab, paste0("UNSUPPRESSED_SHARE_AGE_AND_SEX_", qlab))

#########################################

# FIND SEX SHARE OF INFECTED

#########################################

# find share of unsuppressed by sex across age
tmp <- df[, list(UNSUPPRESSED = sum(UNSUPPRESSED)),
          by = c("ROUND", "COMM", "iterations", "SEX")]
tmp[, UNSUPPRESSED_SHARE := UNSUPPRESSED / sum(UNSUPPRESSED),
    by = c("ROUND", "COMM", "iterations")]

# summarise
sing <- tmp[, list(q = quantile(UNSUPPRESSED_SHARE, prob = ps, na.rm = TRUE),
                   q_label = qlab),
            by = c("ROUND", "COMM", "SEX")]
sing <- as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

# name
setnames(sing, qlab, paste0("UNSUPPRESSED_SHARE_SEX_", qlab))

# merge
sing <- merge(sing_age, sing, by = c("ROUND", "COMM", "SEX"))

#####################################################

# FIND SHARE OF INFECTED UNSUPPRESSED ACROSS AGE BY SEX

#####################################################

# find share of unsuppressed by sex across age
df[, TOTAL_UNSUPPRESSED := sum(UNSUPPRESSED),
    by = c("ROUND", "COMM", "iterations", "SEX")]
df[, UNSUPPRESSED_SHARE := UNSUPPRESSED / TOTAL_UNSUPPRESSED]

# summarise
ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")
sing_age <- df[, list(q = quantile(UNSUPPRESSED_SHARE, prob = ps, na.rm = TRUE),
                      q_label = qlab),
              by = c("ROUND", "COMM", "SEX", "AGEYRS")]
sing_age <- as.data.table(
  reshape2::dcast(sing_age, ... ~ q_label, value.var = "q")
)

# name
setnames(sing_age, qlab, paste0("UNSUPPRESSED_SHARE_AGE_BY_SEX_", qlab))

# merge
sing <- merge(sing_age, sing, by = c("ROUND", "COMM", "SEX", "AGEYRS"))

#########################################

# SAVE

#########################################

tmp <- merge(sing_age, sing, by = c("ROUND", "COMM", "SEX"))
file_name <- file.path(gitdir.fit, "RCCS_unsuppressed_share_sex_221208.csv")
write.csv(tmp, file = file_name, row.names = FALSE)

#########################################

# FIND SEX SHARE OF INFECTED BY AGE GROUP

#########################################

# age groups
df_age_aggregated <- data.table(AGEYRS = df[, sort(unique(AGEYRS))])
df_age_aggregated[, INDEX_AGE_GROUP := 7]
df_age_aggregated[AGEYRS < 45, INDEX_AGE_GROUP := 6]
df_age_aggregated[AGEYRS < 40, INDEX_AGE_GROUP := 5]
df_age_aggregated[AGEYRS < 35, INDEX_AGE_GROUP := 4]
df_age_aggregated[AGEYRS < 30, INDEX_AGE_GROUP := 3]
df_age_aggregated[AGEYRS < 25, INDEX_AGE_GROUP := 2]
df_age_aggregated[AGEYRS < 20, INDEX_AGE_GROUP := 1]
age_grps <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")
df_age_aggregated[, AGE_GROUP := age_grps[INDEX_AGE_GROUP]]

# aggregated by age groups
dfa <- merge(df, df_age_aggregated, by = "AGEYRS")
dfa <- dfa[, list(UNSUPPRESSED = sum(UNSUPPRESSED)),
            by = c("ROUND", "COMM", "AGE_GROUP", "SEX", "iterations")]

# find share of unsuppressed by sex across age
dfa[, UNSUPPRESSED_SHARE := UNSUPPRESSED / sum(UNSUPPRESSED),
    by = c("ROUND", "COMM", "SEX", "iterations")]

# summarise
sing <- dfa[, list(q = quantile(UNSUPPRESSED_SHARE, prob = ps, na.rm = TRUE),
                   q_label = qlab),
            by = c("ROUND", "COMM", "SEX", "AGE_GROUP")]
sing <- as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

# save for round 18 inland
n_digits <- 1
tmp <- sing[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(CL = format(round(CL * 100, n_digits), nsmall = n_digits),
           CU = format(round(CU * 100, n_digits), nsmall = n_digits),
           M = format(round(M * 100, n_digits), nsmall = n_digits))]
tmp[, CL := gsub(" ", "", CL)]
tmp[, CU := gsub(" ", "", CU)]
tmp[, M := gsub(" ", "", M)]
file_name <- file.path(
  outdir,
  "RCCS_shareunsuppressed_age_group_5years_R18_221215.rds"
)
saveRDS(tmp, file = file_name)


#########################################

# FIND SEX SHARE OF INFECTED BY SEX

#########################################


dfa <- df[, list(UNSUPPRESSED = sum(UNSUPPRESSED)),
          by = c("ROUND", "COMM", "SEX", "iterations")]

# find share of unsuppressed by sex across age
dfa[, UNSUPPRESSED_SHARE := UNSUPPRESSED / sum(UNSUPPRESSED),
     by = c("ROUND", "COMM", "iterations")]

# summarise
sing <- dfa[, list(q = quantile(UNSUPPRESSED_SHARE, prob = ps, na.rm = TRUE),
                   q_label = qlab),
            by = c("ROUND", "COMM", "SEX")]
sing <- as.data.table(reshape2::dcast(sing, ... ~ q_label, value.var = "q"))

# save for round 18 inland
n_digits <- 1
tmp <- sing[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(CL = format(round(CL * 100, n_digits), nsmall = n_digits),
           CU = format(round(CU * 100, n_digits), nsmall = n_digits),
           M = format(round(M * 100, n_digits), nsmall = n_digits))]
tmp[, CL := gsub(" ", "", CL)]
tmp[, CU := gsub(" ", "", CU)]
tmp[, M := gsub(" ", "", M)]

file_name <- file.path(outdir, "RCCS_shareunsuppressed_total_R18_221215.rds")
if (! file.exists(file_name) || config$overwrite.existing.files) {
  cat("Saving file:", file_name, "\n")
  saveRDS(tmp, file = file_name)
} else {
  cat("File:", file_name, "already exists...\n")
}
