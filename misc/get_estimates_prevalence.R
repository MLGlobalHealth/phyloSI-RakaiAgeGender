library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library(haven)
library(here)
library(yaml)

# directory to repository
gitdir <- here()

# Load paths
source(file.path(gitdir, "config.R"))

# Output directory to save stan fits and figures
outdir <- file.path(
  "../phyloSI-RakaiAgeGender-outputs",
  "get_estimates_prevalence"
)
if (!dir.exists(outdir)) dir.create(outdir)

# Path to stan model
path_stan <- file.path(gitdir.misc, "stan_models", "binomial_gp.stan")

# Read Stan configurations
model_config <- read_yaml(file.path(gitdir.misc, "stan_models", "config.yml"))

# Load count of participants by hiv status
rprev <- fread(path.count.hivpositive)

# Load nature med requirements
source(file.path(gitdir.functions, "plotting_functions.R"))
naturemed_reqs()

#################################

# PLOT #

#################################

if (TRUE) {
  tmp <- copy(rprev)
  tmp[, Negative := TOTAL_COUNT - COUNT]
  setnames(tmp, "COUNT", "Positive")
  tmp <- melt.data.table(
    data = tmp,
    id.vars = c("ROUND", "COMM", "SEX", "AGEYRS", "TOTAL_COUNT")
  )
  tmp <- tmp[!(ROUND == "R015S" & COMM == "inland")]
  tmp[, ROUND := gsub("R0(.+)", "\\1", ROUND)]
  tmp[, ROUND_LABEL := paste0("Round ", ROUND)]
  tmp[, SEX_LABEL := "Female"]
  tmp[SEX == "M", SEX_LABEL := "Male"]
  tmp[, COMM_LABEL := "Fishing\n communities"]
  tmp[COMM == "inland", COMM_LABEL := "Inland\n communities"]

  # plot
  tmp_subset <- tmp[!ROUND %in% c("06", "07", "08", "09")]
  p <- ggplot(tmp_subset, aes(x = AGEYRS, y = value)) +
    geom_bar(aes(fill = variable), stat = "identity") +
    labs(x = "Age", y = "Count participants", fill = "HIV status") +
    facet_grid(ROUND_LABEL ~ COMM_LABEL + SEX_LABEL) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.background = element_rect(colour = "white", fill = "white"),
          strip.text = element_text(size = rel(1))) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  p

  save_path <- file.path(outdir, "count_participants_by_gender_loc_age.png")
  ggsave(p, file = save_path, w = 8, h = 10)
}



########################

# FIND EMPIRICAL PREVENCE #

########################

setnames(
  rprev,
  old = c("COMM", "SEX", "AGEYRS"),
  new = c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL")
)
rprev[, LOC := as.integer(LOC_LABEL == "fishing")]
rprev[, SEX := as.integer(SEX_LABEL == "M")]
rprev[, AGE := AGE_LABEL - 14L]
rprev[, ROW_ID := seq_len(nrow(rprev))]

# find empirical proportions
rprev[, EMPIRICAL_PREVALENCE := COUNT / TOTAL_COUNT,
      by = c("ROUND", "LOC", "SEX", "AGE")] # prevalence


########################

# FIND SMOOTH PREVENCE #

########################

# find smooth proportion
rounds <- c("R010", "R011", "R012", "R013",
            "R014", "R015", "R016", "R017", "R018")
for (round in rounds) {
  dtbl <- copy(rprev[ROUND == round])

  stopifnot(length(round) == 1)
  cat("Fitting stan model for round ", round, "\n")

  # account for unobserved entries
  tmp <- data.table(
    expand.grid(
      LOC = c(0, 1),
      SEX = c(0, 1),
      AGE_LABEL = rprev[, sort(unique(AGE_LABEL))]
    )
  )
  dtbl <- merge(dtbl, tmp, by = c("LOC", "SEX", "AGE_LABEL"), all.y = TRUE)
  dtbl[is.na(COUNT), COUNT := 0]
  dtbl[is.na(TOTAL_COUNT), TOTAL_COUNT := 0]
  dtbl <- dtbl[order(SEX, LOC, AGE_LABEL)]

  # predicts age
  x_predict <- seq(rprev[, min(AGE_LABEL)], rprev[, max(AGE_LABEL) + 1], 0.5)

  # make stan data
  stan_data <- list(
    x_predict = x_predict,
    y_observed_00 = dtbl[SEX == 0 & LOC == 0, COUNT],
    y_observed_10 = dtbl[SEX == 1 & LOC == 0, COUNT],
    y_observed_01 = dtbl[SEX == 0 & LOC == 1, COUNT],
    y_observed_11 = dtbl[SEX == 1 & LOC == 1, COUNT],
    total_observed_00 = dtbl[SEX == 0 & LOC == 0, TOTAL_COUNT],
    total_observed_10 = dtbl[SEX == 1 & LOC == 0, TOTAL_COUNT],
    total_observed_01 = dtbl[SEX == 0 & LOC == 1, TOTAL_COUNT],
    total_observed_11 = dtbl[SEX == 1 & LOC == 1, TOTAL_COUNT],
    alpha_hyper_par_00 = 2,
    alpha_hyper_par_10 = 2,
    alpha_hyper_par_01 = 2,
    alpha_hyper_par_11 = 2
  )
  stan_data$N_predict <- length(stan_data$x_predict)
  stan_data$observed_idx <- which(stan_data$x_predict %% 1 == 0.5)
  stan_data$N_observed <- length(stan_data$observed_idx)
  stan_data$rho_hyper_par_00 <- diff(range(stan_data$x_predict)) / 3
  stan_data$rho_hyper_par_10 <- diff(range(stan_data$x_predict)) / 3
  stan_data$rho_hyper_par_01 <- diff(range(stan_data$x_predict)) / 3
  stan_data$rho_hyper_par_11 <- diff(range(stan_data$x_predict)) / 3

  # load stan model
  stan_model <- stan_model(path_stan, model_name = "gp_all")

  # run and save model
  fit <- sampling(
    stan_model,
    data = stan_data,
    iter = model_config$iter,
    warmup = model_config$warmup,
    chains = model_config$chains,
    cores = model_config$cores,
    control = list(
      max_treedepth = model_config$control$max_treedepth,
      adapt_delta = model_config$control$adapt_delta
    )
  )

  file_name <- paste0("hivprevalence_gp_stanfit_round",
                     gsub("R0", "", round),
                     "_221116.rds")
  file_name <- file.path(outdir, file_name)
  saveRDS(fit, file = file_name)
}

# load results
rounds <- c(10:15, 16:18)
nsinf <- vector(mode = "list", length = length(rounds))
nspred <- vector(mode = "list", length = length(rounds))
nsinf_samples <- vector(mode = "list", length = length(rounds))
convergence_list <- vector(mode = "list", length = length(rounds))
for (i in seq_along(rounds)) {
  round <- rounds[i]
  dtbl <- copy(rprev[ROUND == paste0("R0", round)])

  # age to predict
  x_predict <- seq(rprev[, min(AGE_LABEL)], rprev[, max(AGE_LABEL) + 1], 0.5)

  # load samples
  filename <- paste0("hivprevalence_gp_stanfit_round", round, "_221116.rds")
  fit <- readRDS(file.path(outdir, filename))
  re <- rstan::extract(fit)

  #
  # Find rhat and neff
  #

  sum_fit <- summary(fit)
  neff <- na.omit(sum_fit$summary[, 9])
  rhat <- na.omit(sum_fit$summary[, 10])
  convergence <- data.table(ROUND = round, neff = neff, rhat = rhat)


  #
  #	summarise estimated prevalence by sex and age
  #

  # extract estimates
  ps <- c(0.025, 0.5, 0.975)
  qlab <- c("CL", "M", "CU")
  tmp <- as.data.table(reshape2::melt(re$p_predict_00))
  tmp[, `:=`(SEX = 0, LOC = 0)]
  tmp1 <- as.data.table(reshape2::melt(re$p_predict_10))
  tmp1[, `:=`(SEX = 1, LOC = 0)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(re$p_predict_01))
  tmp1[, `:=`(SEX = 0, LOC = 1)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(re$p_predict_11))
  tmp1[, `:=`(SEX = 1, LOC = 1)]
  tmp <- rbind(tmp, tmp1)

  tmp[, AGE_LABEL := x_predict[Var2]]
  set(tmp, NULL, "Var2", NULL)

  # summarise
  nsinf_by_age <- tmp[, list(q = quantile(value, prob = ps, na.rm = TRUE),
                             q_label = qlab),
                      by = c("SEX", "LOC", "AGE_LABEL")]
  nsinf_by_age <- as.data.table(
    reshape2::dcast(nsinf_by_age, ... ~ q_label, value.var = "q")
  )

  #
  #	summarise predicted prevalence by sex and age
  #

  # merge to total count
  tmp <- merge(
    tmp,
    dtbl[, .(SEX, LOC, AGE_LABEL, TOTAL_COUNT)],
    by = c("SEX", "LOC", "AGE_LABEL"),
    all.x = TRUE
  )

  # predict count and predict prevalence
  tmp[!is.na(TOTAL_COUNT), COUNT_PREDICT := rbinom(1, TOTAL_COUNT, value),
      by = c("SEX", "LOC", "AGE_LABEL", "iterations")]
  tmp[, PREVALENCE_PREDICT := COUNT_PREDICT / TOTAL_COUNT]

  # summarise
  nspred_by_age <- tmp[, list(q = quantile(PREVALENCE_PREDICT, prob = ps, na.rm = TRUE), # nolint: line_length_linter.
                              q_label = qlab),
                      by = c("SEX", "LOC", "AGE_LABEL")]
  nspred_by_age <- as.data.table(
    reshape2::dcast(nspred_by_age, ... ~ q_label, value.var = "q")
  )

  #
  # POSTPROCESING
  #

  # merge to data
  subset_vars <- c("SEX", "SEX_LABEL", "LOC",
                   "LOC_LABEL", "AGE_LABEL", "EMPIRICAL_PREVALENCE")
  join_keys <- c("SEX", "LOC", "AGE_LABEL")
  dtbl_subset <- subset(dtbl, select = subset_vars)
  nsinf_by_age <- merge(dtbl_subset, nsinf_by_age, by = join_keys)
  nsinf_samples_by_age <- merge(dtbl_subset, tmp, by = join_keys)
  nspred_by_age <- merge(dtbl_subset, nspred_by_age, by = join_keys)

  # load change of var name
  set(nsinf_by_age, NULL, "SEX", NULL)
  set(nsinf_by_age, NULL, "LOC", NULL)
  setnames(
    nsinf_by_age,
    old = c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "M", "CL", "CU"),
    new = c("COMM", "SEX", "AGEYRS",
            "PREVALENCE_M",
            "PREVALENCE_CL",
            "PREVALENCE_CU")
  )
  nsinf_by_age[, ROUND := paste0("R0", round)]

  # load change of var name
  set(nspred_by_age, NULL, "SEX", NULL)
  set(nspred_by_age, NULL, "LOC", NULL)
  setnames(
    nspred_by_age,
    old = c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "M", "CL", "CU"),
    new = c("COMM", "SEX", "AGEYRS",
            "PREVALENCE_M",
            "PREVALENCE_CL",
            "PREVALENCE_CU")
  )
  nspred_by_age[, ROUND := paste0("R0", round)]

  # load change of var name
  set(nsinf_samples.by.age, NULL, "SEX", NULL)
  set(nsinf_samples.by.age, NULL, "LOC", NULL)
  set(nsinf_samples.by.age, NULL, "COUNT_PREDICT", NULL)
  set(nsinf_samples.by.age, NULL, "TOTAL_COUNT", NULL)
  set(nsinf_samples.by.age, NULL, "PREVALENCE_PREDICT", NULL)
  setnames(
    nsinf_samples.by.age,
    old = c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "value"),
    new = c("COMM", "SEX", "AGEYRS", "PREVALENCE_POSTERIOR_SAMPLE"))
  nsinf_samples.by.age[, ROUND := paste0("R0", round)]

  # keep
  nsinf[[i]] <- nsinf_by_age
  nspred[[i]] <- nspred_by_age
  nsinf_samples[[i]] <- nsinf_samples.by.age
  convergence_list[[i]] <- convergence
}
nsinf <- do.call("rbind", nsinf)
nspred <- do.call("rbind", nspred)
nsinf_samples <- do.call("rbind", nsinf_samples)
convergence <- do.call("rbind", convergence_list)

# check that all entries are complete
stopifnot(
  nrow(nsinf[COMM == "inland"]) == nsinf[, length(unique(AGEYRS))] *
                                   nsinf[, length(unique(SEX))] *
                                   nsinf[COMM == "inland",
                                         length(unique(ROUND))]
)
stopifnot(
  nrow(nsinf[COMM == "fishing"]) == nsinf[, length(unique(AGEYRS))] *
                                    nsinf[, length(unique(SEX))] *
                                    nsinf[COMM == "fishing",
                                          length(unique(ROUND))]
)


##############################

# PLOT PREVALENCE #

##############################

# PREDICTED PREVALENCE
tmp <- copy(nspred)
tmp <- tmp[!(ROUND == "R015S" & COMM == "inland")]
tmp[, ROUND := gsub("R0(.+)", "\\1", ROUND)]
tmp[, ROUND_LABEL := paste0("Round ", ROUND)]
tmp[, SEX_LABEL := "Women"]
tmp[SEX == "M", SEX_LABEL := "Men"]
tmp[, COMM_LABEL := "Fishing\n communities"]
tmp[COMM == "inland", COMM_LABEL := "Inland\n communities"]

ggplot(tmp[COMM == "inland"], aes(x = AGEYRS)) +
  geom_ribbon(aes(ymin = PREVALENCE_CL, ymax = PREVALENCE_CU, fill = SEX_LABEL),
              alpha = 0.5) +
  geom_line(aes(y = PREVALENCE_M, col = SEX_LABEL, linetype = "Prediction")) +
  geom_point(aes(y = EMPIRICAL_PREVALENCE, col = SEX_LABEL, shape = "Data"),
             alpha = 0.7) +
  facet_grid(ROUND_LABEL ~ .) +
  theme_bw() +
  scale_color_manual(values = c("Men" = "royalblue3", "Women" = "deeppink")) +
  scale_fill_manual(values = c("Men" = "lightblue3", "Women" = "lightpink1")) +
  theme(legend.position = "bottom",
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = rel(1)),
        panel.spacing.y = unit(0.7, "lines")) +
  labs(x = "Age", y = "HIV prevalence in RCCS participants",
       col = "", fill = "",
       shape = "", linetype = "") +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, .5), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2),
         color = guide_legend(order = 3), fill = guide_legend(order = 3))

ggsave(
  file = file.path(outdir, "smooth_predicted_prevalence_221116.png"),
  w = 5,
  h = 10
)


# ESTIMATED PREVALENCE
tmp <- copy(nsinf)
tmp <- tmp[!(ROUND == "R015S" & COMM == "inland")]
tmp[, ROUND := gsub("R0(.+)", "\\1", ROUND)]
tmp[, ROUND_LABEL := paste0("Round ", ROUND)]
tmp[, SEX_LABEL := "Women"]
tmp[SEX == "M", SEX_LABEL := "Men"]
tmp[, COMM_LABEL := "Fishing\n communities"]
tmp[COMM == "inland", COMM_LABEL := "Inland\n communities"]

ggplot(tmp[COMM == "inland"], aes(x = AGEYRS)) +
  geom_ribbon(aes(ymin = PREVALENCE_CL, ymax = PREVALENCE_CU, fill = SEX_LABEL),
                  alpha = 0.7) +
  geom_line(aes(y = PREVALENCE_M, col = SEX_LABEL, linetype = "Fit")) +
  geom_point(aes(y = EMPIRICAL_PREVALENCE, col = SEX_LABEL, shape = "Data"),
             alpha = 0.7) +
  facet_wrap(~ROUND_LABEL) +
  theme_bw() +
  scale_color_manual(values = c("Men" = "royalblue3", "Women" = "deeppink")) +
  scale_fill_manual(values = c("Men" = "lightblue3", "Women" = "lightpink1")) +
  theme(legend.position = "bottom",
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = rel(1)),
        panel.spacing.y = unit(0.7, "lines")) +
  labs(x = "Age", y = "HIV prevalence in RCCS participants",
       col = "", fill = "",
       shape = "", linetype = "") +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, .5),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2),
         color = guide_legend(order = 3), fill = guide_legend(order = 3)) + reqs

ggsave(
  file = file.path(outdir, "smooth_estimated_prevalence_221116.pdf"),
  w = 7,
  h = 7
)

###########################

# STATISTICS FOR PAPER #

###########################

# get proportion of predicted prevalence inside credible interval
stats <- list()
tmp <- nspred[COMM == "inland" & !is.na(EMPIRICAL_PREVALENCE)]
tmp[, within.CI := data.table::between(EMPIRICAL_PREVALENCE,
                                       PREVALENCE_CL,
                                       PREVALENCE_CU)]
stats[["within.CI"]] <- tmp[, paste0(round(mean(within.CI) * 100, 2))]

# get lowest rhat and lowest neff
stats[["min_neff"]] <- convergence[, round(min(neff))]
stats[["max_rhat"]] <- convergence[, round(max(rhat), 4)]


#########

# SAVE #

#########

file_name <- file.prevalence.prop
if (!file.exists(file_name) || config$overwrite.existing.files) {
  cat("Saving file:", file_name, "\n")
  write.csv(nsinf, file = file_name, row.names = FALSE)
} else {
    cat("File:", file_name, "already exists...\n")
}

file_name <- file.prevalence
if (!file.exists(file_name) || config$overwrite.existing.files) {
  cat("Saving file:", file_name, "\n")
  saveRDS(nsinf_samples, file = file_name)
} else {
  cat("File:", file_name, "already exists...\n")
}

file_name <- file.path(
  outdir,
  "RCCS_prevalence_model_fit_convergence_221116.RDS"
)
saveRDS(stats, file = file_name)