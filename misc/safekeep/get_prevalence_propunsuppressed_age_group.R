library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(here)

# AB: This only saves things on the HPC,
# but I am not sure the outputs are used in any downstream code

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

file.exists(eligible_count) |> stopifnot()
file.exists(proportion_prevalence) |> stopifnot()
file.exists(treatment_cascade) |> stopifnot()

outdir <- file.path(indir.deepsequence_analyses, "PANGEA2_RCCS", "prevalence_by_gender_loc_age")

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
    by = c("ROUND", "COMM", "AGEYRS", "SEX", "iterations"),
    all.x = TRUE
)
df[, ROUND := gsub("R0(.+)", "\\1", ROUND)]

# merge number of eligible and the prevalence
eligible_count <- eligible_count[, .(ROUND, COMM, AGEYRS, SEX, ELIGIBLE)]
df <- merge(
    eligible_count,
    df,
    by = c("ROUND", "COMM", "AGEYRS", "SEX")
)

# find infected
df[, INFECTED := ELIGIBLE * PREVALENCE_POSTERIOR_SAMPLE]

# find infected unsuppressed
df[, PROP_UNSUPPRESSED_POSTERIOR_SAMPLE := 1 - PROP_SUPPRESSED_POSTERIOR_SAMPLE]
df[, UNSUPPRESSED := INFECTED * PROP_UNSUPPRESSED_POSTERIOR_SAMPLE]


#####################################################

# FIND  PREVALENCE & PROPORTION UNSUPPRESSED BY AGE GROUP
# "15-24", "25-34", "35-49"

#####################################################

# age groups
df_age_aggregated <- data.table(AGEYRS = df[, sort(unique(AGEYRS))])
df_age_aggregated[, INDEX_AGE_GROUP := 3]
df_age_aggregated[AGEYRS < 35, INDEX_AGE_GROUP := 2]
df_age_aggregated[AGEYRS < 25, INDEX_AGE_GROUP := 1]
df_age_aggregated[, AGE_GROUP := c("15-24", "25-34", "35-49")[INDEX_AGE_GROUP]]

# aggregated by age groups
df_1 <- merge(df, df_age_aggregated, by = "AGEYRS")
df_1 <- df_1[, list(
    INFECTED = sum(INFECTED),
    UNSUPPRESSED = sum(UNSUPPRESSED),
    ELIGIBLE = sum(ELIGIBLE)
), by = c("ROUND", "COMM", "AGE_GROUP", "SEX", "iterations")]

# find prevalence by age groups across age
df_1[, PROP_UNSUPPRESSED := UNSUPPRESSED / ELIGIBLE]
df_1[, PREVALENCE := INFECTED / ELIGIBLE]

# summarise unsuppressed
ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")
sing_age <- df_1[, list(q = quantile(na.omit(PROP_UNSUPPRESSED), prob = ps, na.rm = TRUE),
                        q_label = qlab),
                    by = c("ROUND", "COMM", "SEX", "AGE_GROUP")]
sing_age <- as.data.table(
    reshape2::dcast(sing_age, ... ~ q_label, value.var = "q")
)
setnames(sing_age, qlab, paste0("PROP_UNSUPPRESSED_AGE_GROUP_", qlab))

# summarise prevalence
sinf_age <- df_1[, list(q = quantile(PREVALENCE, prob = ps, na.rm = TRUE),
                        q_label = qlab),
                    by = c("ROUND", "COMM", "SEX", "AGE_GROUP")]
sinf_age <- as.data.table(
    reshape2::dcast(sinf_age, ... ~ q_label, value.var = "q")
)
setnames(sinf_age, qlab, paste0("PREVALENCE_AGE_GROUP_", qlab))

# plot
ggplot(sing_age, aes(x = AGE_GROUP, group = SEX)) +
    geom_errorbar(aes(ymin = PROP_UNSUPPRESSED_AGE_GROUP_CL,
                      ymax = PROP_UNSUPPRESSED_AGE_GROUP_CU),
                  alpha = 0.5,
                  width = 0.3,
                  position = position_dodge(width = 0.3)) +
    geom_point(aes(y = PROP_UNSUPPRESSED_AGE_GROUP_M, col = SEX),
              position = position_dodge(width = 0.3)) +
    facet_grid(ROUND ~ COMM) +
    theme_bw()

ggplot(sinf_age, aes(x = AGE_GROUP, group = SEX)) +
    geom_errorbar(aes(ymin = PREVALENCE_AGE_GROUP_CL,
                      ymax = PREVALENCE_AGE_GROUP_CU),
                  alpha = 0.5,
                  width = 0.3,
                  position = position_dodge(width = 0.3)) +
    geom_point(aes(y = PREVALENCE_AGE_GROUP_M, col = SEX),
              position = position_dodge(width = 0.3)) +
    facet_grid(ROUND ~ COMM) +
    theme_bw()

#########################################

# SAVE

#########################################

file_name <- file.unsuppressed.agegroup
if (! file.exists(file_name) || config$overwrite.existing.files) {
    cat("Saving file:", file_name, "\n")
    fwrite(sing_age, file = file_name, row.names = FALSE)

} else {
    cat("File:", file_name, "already exists...\n")
}

file_name <- file.prevalence.agegroup
if (!file.exists(file_name) || config$overwrite.existing.files) {
    cat("Saving file:", file_name, "\n")
    fwrite(sinf_age, file = file_name, row.names = FALSE)

} else {
    cat("File:", file_name, "already exists...\n")
}



#####################################################

# FIND  PREVALENCE &  PROPORTION UNSUPPRESSED BY AGE GROUP
# "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"

#####################################################

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
df_1 <- merge(df, df_age_aggregated, by = "AGEYRS")
df_1 <- df_1[, list(
    INFECTED = sum(INFECTED),
    UNSUPPRESSED = sum(UNSUPPRESSED),
    ELIGIBLE = sum(ELIGIBLE)
), by = c("ROUND", "COMM", "AGE_GROUP", "SEX", "iterations")]

# find prevalence by age groups across age
df_1[, PROP_UNSUPPRESSED := UNSUPPRESSED / INFECTED]
df_1[, PREVALENCE := INFECTED / ELIGIBLE]
df_1[, DIFF_PROP_UNSUPPRESSED := PROP_UNSUPPRESSED[SEX == "M"] - PROP_UNSUPPRESSED[SEX == "F"],
     by = c("ROUND", "COMM", "AGE_GROUP", "iterations")]

# summarise unsuppressed
ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")
sing_age <- df_1[, list(q = quantile(na.omit(PROP_UNSUPPRESSED), prob = ps, na.rm = TRUE),
                        q_label = qlab),
                 by = c("ROUND", "COMM", "SEX", "AGE_GROUP")]
sing_age <- as.data.table(
    reshape2::dcast(sing_age, ... ~ q_label, value.var = "q")
)
setnames(sing_age, qlab, paste0("PROP_UNSUPPRESSED_AGE_GROUP_", qlab))

# summarise prevalence
sinf_age <- df_1[, list(q = quantile(PREVALENCE, prob = ps, na.rm = TRUE),
                        q_label = qlab),
                 by = c("ROUND", "COMM", "SEX", "AGE_GROUP")]
sinf_age <- as.data.table(
    reshape2::dcast(sinf_age, ... ~ q_label, value.var = "q")
)
setnames(sinf_age, qlab, paste0("PREVALENCE_AGE_GROUP_", qlab))

# summarise difference prop unsuppressed
sind_age <- unique(df_1[, .(ROUND, COMM, AGE_GROUP, iterations, DIFF_PROP_UNSUPPRESSED)])
sind_age <- sind_age[, list(q = quantile(DIFF_PROP_UNSUPPRESSED, prob = ps, na.rm = TRUE),
                            q_label = qlab),
                     by = c("ROUND", "COMM", "AGE_GROUP")]
sind_age <- as.data.table(
    reshape2::dcast(sind_age, ... ~ q_label, value.var = "q")
)
setnames(sind_age, qlab, paste0("DIFF_PROP_UNSUPPRESSED_AGE_GROUP_", qlab))


#########################################

# SAVE FOR ROUND 18

#########################################

n_digits <- 1

tmp <- sing_age[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(
    PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CL * 100, n_digits), nsmall = n_digits),
    PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CU * 100, n_digits), nsmall = n_digits),
    PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(PROP_UNSUPPRESSED_AGE_GROUP_M * 100, n_digits), nsmall = n_digits)
)]
tmp[, `:=`(
    PROP_UNSUPPRESSED_AGE_GROUP_CL = gsub(" ", "", PROP_UNSUPPRESSED_AGE_GROUP_CL),
    PROP_UNSUPPRESSED_AGE_GROUP_CU = gsub(" ", "", PROP_UNSUPPRESSED_AGE_GROUP_CU),
    PROP_UNSUPPRESSED_AGE_GROUP_M = gsub(" ", "", PROP_UNSUPPRESSED_AGE_GROUP_M)
)]

file_name <- file.path(outdir, "RCCS_propunsuppressed_age_group_5years_R18_221215.rds")
if (! file.exists(file_name) || config$overwrite.existing.files) {
    cat("Saving file:", file_name, "\n")
    saveRDS(tmp, file = file_name)
} else {
    cat("File:", file_name, "already exists...\n")
}

tmp <- sinf_age[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(
    PREVALENCE_AGE_GROUP_CL = format(round(PREVALENCE_AGE_GROUP_CL * 100, n_digits), nsmall = n_digits),
    PREVALENCE_AGE_GROUP_CU = format(round(PREVALENCE_AGE_GROUP_CU * 100, n_digits), nsmall = n_digits),
    PREVALENCE_AGE_GROUP_M = format(round(PREVALENCE_AGE_GROUP_M * 100, n_digits), nsmall = n_digits)
)]
tmp[, `:=`(
    PREVALENCE_AGE_GROUP_CL = gsub(" ", "", PREVALENCE_AGE_GROUP_CL),
    PREVALENCE_AGE_GROUP_CU = gsub(" ", "", PREVALENCE_AGE_GROUP_CU),
    PREVALENCE_AGE_GROUP_M = gsub(" ", "", PREVALENCE_AGE_GROUP_M)
)]

file_name <- file.path(
    outdir,
    "RCCS_prevalence_age_group_5years_R18_221215.rds"
)
if (!file.exists(file_name) || config$overwrite.existing.files) {
    cat("Saving file:", file_name, "\n")
    saveRDS(tmp, file = file_name)
} else {
    cat("File:", file_name, "already exists...\n")
}


tmp <- sind_age[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL * 100, n_digits), nsmall = n_digits),
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU * 100, n_digits), nsmall = n_digits),
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M * 100, n_digits), nsmall = n_digits)
)]

tmp[, `:=`(
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL = gsub(" ", "", DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL),
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU = gsub(" ", "", DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU),
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M = gsub(" ", "", DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M)
)]

file_name <- file.path(
    outdir,
    "RCCS_diffpropunsuppressed_age_group_5years_R18_221215.rds"
)
if (! file.exists(file_name) || config$overwrite.existing.files) {
    cat("Saving file:", file_name, "\n")
    saveRDS(tmp, file = file_name)
} else {
    cat("File:", file_name, "already exists...\n")
}



#####################################################

# FIND PREVALENCE & PROPORTION UNSUPPRESSED ACROSS AGE

#####################################################

# aggregated across ages
df_1 <- df[, list(
    INFECTED = sum(INFECTED),
    UNSUPPRESSED = sum(UNSUPPRESSED),
    ELIGIBLE = sum(ELIGIBLE)
), by = c("ROUND", "COMM", "SEX", "iterations")]

# find prevalence and unsuppressed
df_1[, PROP_UNSUPPRESSED := UNSUPPRESSED / INFECTED]
df_1[, PREVALENCE := INFECTED / ELIGIBLE]
df_1[, DIFF_PROP_UNSUPPRESSED := PROP_UNSUPPRESSED[SEX == "M"] - PROP_UNSUPPRESSED[SEX == "F"],
     by = c("ROUND", "COMM", "iterations")]

# summarise unsuppressed
ps <- c(0.025, 0.5, 0.975)
qlab <- c("CL", "M", "CU")
sing_t <- df_1[, list(q = quantile(na.omit(PROP_UNSUPPRESSED), prob = ps, na.rm = TRUE),
                      q_label = qlab),
                by = c("ROUND", "COMM", "SEX")]
sing_t <- as.data.table(reshape2::dcast(sing_t, ... ~ q_label, value.var = "q"))
setnames(sing_t, qlab, paste0("PROP_UNSUPPRESSED_AGE_GROUP_", qlab))

# summarise prevalence
sinf_t <- df_1[, list(q = quantile(PREVALENCE, prob = ps, na.rm = TRUE),
                      q_label = qlab),
                by = c("ROUND", "COMM", "SEX")]
sinf_t <- as.data.table(reshape2::dcast(sinf_t, ... ~ q_label, value.var = "q"))
setnames(sinf_t, qlab, paste0("PREVALENCE_AGE_GROUP_", qlab))

# summarise difference prop unsuppressed
sind_t <- unique(df_1[, .(ROUND, COMM, iterations, DIFF_PROP_UNSUPPRESSED)])
sind_t <- sind_t[, list(q = quantile(DIFF_PROP_UNSUPPRESSED, prob = ps, na.rm = TRUE),
                        q_label = qlab),
                 by = c("ROUND", "COMM")]
sind_t <- as.data.table(reshape2::dcast(sind_t, ... ~ q_label, value.var = "q"))
setnames(sind_t, qlab, paste0("DIFF_PROP_UNSUPPRESSED_AGE_GROUP_", qlab))


#########################################

# SAVE FOR ROUND 18

#########################################

n_digits <- 1

tmp <- sing_t[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(
    PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CL * 100, n_digits), nsmall = n_digits),
    PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(PROP_UNSUPPRESSED_AGE_GROUP_CU * 100, n_digits), nsmall = n_digits),
    PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(PROP_UNSUPPRESSED_AGE_GROUP_M * 100, n_digits), nsmall = n_digits)
)]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_CL := gsub(" ", "", PROP_UNSUPPRESSED_AGE_GROUP_CL)]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_CU := gsub(" ", "", PROP_UNSUPPRESSED_AGE_GROUP_CU)]
tmp[, PROP_UNSUPPRESSED_AGE_GROUP_M := gsub(" ", "", PROP_UNSUPPRESSED_AGE_GROUP_M)]
file_name <- file.path(outdir, "RCCS_propunsuppressed_total_R18_221215.rds")
if (!file.exists(file_name) | config$overwrite.existing.files) {
    cat("saving file:", file_name, "\n")
    saveRDS(tmp, file = file_name)
} else {
    cat("File:", file_name, "already exists...\n")
}

tmp <- sinf_t[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(
    PREVALENCE_AGE_GROUP_CL = format(round(PREVALENCE_AGE_GROUP_CL * 100, n_digits), nsmall = n_digits),
    PREVALENCE_AGE_GROUP_CU = format(round(PREVALENCE_AGE_GROUP_CU * 100, n_digits), nsmall = n_digits),
    PREVALENCE_AGE_GROUP_M = format(round(PREVALENCE_AGE_GROUP_M * 100, n_digits), nsmall = n_digits)
)]
tmp[, PREVALENCE_AGE_GROUP_CL := gsub(" ", "", PREVALENCE_AGE_GROUP_CL)]
tmp[, PREVALENCE_AGE_GROUP_CU := gsub(" ", "", PREVALENCE_AGE_GROUP_CU)]
tmp[, PREVALENCE_AGE_GROUP_M := gsub(" ", "", PREVALENCE_AGE_GROUP_M)]
file_name <- file.path(outdir, "RCCS_prevalence_total_R18_221215.rds")
if (!file.exists(file_name) | config$overwrite.existing.files) {
    cat("Saving file:", file_name, "\n")
    saveRDS(tmp, file = file_name)
} else {
    cat("File:", file_name, "already exists...\n")
}

tmp <- sind_t[ROUND == "18" & COMM == "inland"]
tmp[, `:=`(
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL * 100, n_digits), nsmall = n_digits),
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU * 100, n_digits), nsmall = n_digits),
    DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M = format(round(DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M * 100, n_digits), nsmall = n_digits)
)]

tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL := gsub(" ", "", DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CL)]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU := gsub(" ", "", DIFF_PROP_UNSUPPRESSED_AGE_GROUP_CU)]
tmp[, DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M := gsub(" ", "", DIFF_PROP_UNSUPPRESSED_AGE_GROUP_M)]
file_name <- file.path(outdir, "RCCS_diffpropunsuppressed_total_R18_221215.rds")
if (!file.exists(file_name) | config$overwrite.existing.files) {
    cat("Saving file:", file_name, "\n")
    saveRDS(tmp, file = file_name)
} else {
    cat("File:", file_name, "already exists...\n")
}
