library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library(haven)
library(here)
library(yaml)

# directory of the repository
gitdir <- here()

# load file paths
source(file.path(gitdir, "config.R"))

# outdir directory for stan fit
outdir <- file.path("../phyloSI-RakaiAgeGender-outputs","get_estimates_art_coverage_non_participants")
if(usr == 'melodiemonod'){
  outdir <- file.path("/Users/melodiemonod/Box Sync/2023//phyloSI-RakaiAgeGender-outputs","get_estimates_art_coverage_non_participants")
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

file.exists(c(
  path_stan_binomialgp ,
  path_binomialgp_model_config ,
  path.newly.registered.art ))  |> all() |> stopifnot()

# find count of newly registered participants who reported art use
rart <- fread(path.newly.registered.art)

# config stan
binomialgp_model_config <- read_yaml(path_binomialgp_model_config)

#################################

# PLOT  #

#################################

# plot
if (1) {
  tmp <- copy(rart)
  tmp[, `Do not use` := TOTAL_COUNT - COUNT]
  setnames(tmp, "COUNT", "Use")
  tmp <- melt.data.table(
    tmp,
    id.vars = c("ROUND", "COMM", "SEX", "AGEYRS", "TOTAL_COUNT")
  )
  tmp[, variable := factor(variable, levels = c("Use", "Do not use"))]
  tmp <- tmp[!(ROUND == "R015S" & COMM == "inland")]
  tmp[, ROUND := gsub("R0(.+)", "\\1", ROUND)]
  tmp[, ROUND_LABEL := paste0("Round ", ROUND)]
  tmp[, SEX_LABEL := "Women"]
  tmp[SEX == "M", SEX_LABEL := "Men"]
  tmp[, COMM_LABEL := "Fishing\n communities"]
  tmp[COMM == "inland", COMM_LABEL := "Inland\n communities"]
  tmp <- tmp[AGEYRS > 14 & AGEYRS < 50]

  # plot
  tmp_subset <- tmp[!ROUND %in% c("06", "07", "08", "09", "10") & COMM == "inland"]
  p <- ggplot(tmp_subset, aes(x = AGEYRS, y = value)) +
    geom_bar(aes(fill = variable), stat = "identity") +
    labs(x = "Age", y = "Count newly registered HIV-positive participants", fill = "") +
    facet_grid(ROUND_LABEL ~ SEX_LABEL) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.background = element_rect(colour = "white", fil = "white"),
          strip.text = element_text(size = rel(1))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = c("#90B77D", "#425F57"),
                      labels = c("Reported ART use", "Did not report ART use"))

  save_path <- file.path(
     outdir,
    "count_selfreportedart_by_gender_loc_age_newlyregistered_230703.pdf"
  )
  ggsave(p, file = save_path, w = 7, h = 9)
}


########################

# FIND CRUDE PROPORTION #

########################

setnames(
  rart,
  c("COMM", "SEX", "AGEYRS"),
  c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL")
)
rart[, LOC := as.integer(LOC_LABEL == "fishing")]
rart[, SEX := as.integer(SEX_LABEL == "M")]
rart[, AGE := AGE_LABEL - 14L]
rart[, ROW_ID := seq_len(nrow(rart))]

# find empirical proportions
rart[, PROP_ART_COVERAGE_EMPIRICAL := COUNT / TOTAL_COUNT,
      by = c("ROUND", "LOC", "SEX", "AGE")]


########################

# FIND SMOOTH PROPORTION #

########################

# find smooth proportion
rounds <- c("R010", "R011", "R012", "R013", "R014",
            "R015", "R016", "R017", "R018")
for (r in rounds) {
  dtbl <- copy(rart[ROUND == r])
  dtbl <- dtbl[order(SEX, LOC, AGE_LABEL)]

  stopifnot(length(r) == 1)
  cat("Fitting stan model for round ", r, "\n")

  # account for unobserved entries
  tmp <- data.table(
    expand.grid(LOC = c(0, 1),
    SEX = c(0, 1),
    AGE_LABEL = rart[, sort(unique(AGE_LABEL))])
  )
  dtbl <- merge(tmp, dtbl, by = c("LOC", "SEX", "AGE_LABEL"), all.x = TRUE)
  dtbl[is.na(COUNT), COUNT := 0]
  dtbl[is.na(TOTAL_COUNT), TOTAL_COUNT := 0]
  dtbl <- dtbl[order(SEX, LOC, AGE_LABEL)]

  # predicts age
  x_predict <- seq(rart[, min(AGE_LABEL)], rart[, max(AGE_LABEL) + 1], 0.5)

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
  stan_model <- stan_model(path_stan_binomialgp, model_name = "gp_all")

  # run and save model
  fit <- sampling(
    stan_model,
    data = stan_data,
    iter = binomialgp_model_config$iter,
    warmup = binomialgp_model_config$warmup,
    chains = binomialgp_model_config$chains,
    cores = binomialgp_model_config$cores,
    control = list(
      max_treedepth = binomialgp_model_config$control$max_treedepth,
      adapt_delta = binomialgp_model_config$control$adapt_delta
    )
  )
  filename <- paste0("art_gp_stanfit_round",
                     gsub("R0", "", r),
                     "_newlyregistered_230703.rds")
  saveRDS(fit, file = file.path(outdir, filename))
}

# load results
rounds <- c(10:15, "16", "17", "18")
nsinf <- vector(mode = "list", length = length(rounds))
nsinf_samples <- vector(mode = "list", length = length(rounds))
nspred <- vector(mode = "list", length = length(rounds))
convergence_list <- vector(mode = "list", length = length(rounds))

for (i in seq_along(rounds)) {
  round <- rounds[i]
  dtbl <- copy(rart[ROUND == paste0("R0", round)])

  # account for unobserved age entries but not loc
  tmp <- data.table(
    expand.grid(
      LOC = dtbl[, unique(LOC)],
      SEX = c(0, 1),
      AGE_LABEL = rart[, sort(unique(AGE_LABEL))],
      ROUND = dtbl[, unique(ROUND)])
  )

  rart_unique <- unique(rart[LOC %in% dtbl[, unique(LOC)], .(LOC, SEX, SEX_LABEL, LOC_LABEL)]) # nolint: line_length_linter.
  tmp <- merge(tmp, rart_unique, by = c("LOC", "SEX"), all.x = TRUE)

  join_keys <- c("LOC", "SEX", "AGE_LABEL", "SEX_LABEL", "LOC_LABEL", "ROUND")
  dtbl <- merge(tmp, dtbl, by = join_keys, all.x = TRUE)
  dtbl[is.na(COUNT), COUNT := 0]
  dtbl[is.na(TOTAL_COUNT), TOTAL_COUNT := 0]
  dtbl <- dtbl[order(SEX, LOC, AGE_LABEL)]

  # age to predict
  x_predict <- seq(rart[, min(AGE_LABEL)], rart[, max(AGE_LABEL) + 1], 0.5)

  # load samples
  filename <- paste0(
    "art_gp_stanfit_round",
    round,
    "_newlyregistered_230703.rds"
  )

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
  #	summarise estimated art use by sex and age
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

  nsinf_by_age <- tmp[, list(q = quantile(value, prob = ps, na.rm = TRUE),
                            q_label = qlab),
                      by = c("SEX", "LOC", "AGE_LABEL")]
  nsinf_by_age <- as.data.table(
    reshape2::dcast(nsinf_by_age, ... ~ q_label, value.var = "q")
  )

  #
  #	summarise predicted art use by sex and age
  #

  # merge to total count
  tmp <- merge(
    tmp, dtbl[, .(SEX, LOC, AGE_LABEL, TOTAL_COUNT)],
    by = c("SEX", "LOC", "AGE_LABEL"),
    all.x = TRUE
  )

  # predict count and predict art use
  tmp[!is.na(TOTAL_COUNT), COUNT_PREDICT := rbinom(1, TOTAL_COUNT, value),
      by = c("SEX", "LOC", "AGE_LABEL", "iterations")]
  tmp[, ART_USE_PREDICT := COUNT_PREDICT / TOTAL_COUNT]

  # summarise
  nspred_by_age <- tmp[, list(q = quantile(ART_USE_PREDICT, prob = ps, na.rm = TRUE), # nolint: line_length_linter.
                             q_label = qlab),
                      by = c("SEX", "LOC", "AGE_LABEL")]
  nspred_by_age <- as.data.table(
    reshape2::dcast(nspred_by_age, ... ~ q_label, value.var = "q")
  )

  # sub-sample the last 9500 iterations
  it <- data.table(iterations = tmp[, sort(unique(iterations))])
  it[, iterations_rev := max(iterations):1]
  tmp <- merge(it, tmp, by = "iterations")
  tmp <- tmp[iterations_rev %in% 1:9500]
  tmp[, iterations := iterations - min(iterations) + 1]
  set(tmp, NULL, "iterations_rev", NULL)

  #
  # POSTPROCESING
  #

  # Merge to data
  var_names <- c("SEX", "SEX_LABEL", "LOC", "LOC_LABEL",
                 "AGE_LABEL", "PROP_ART_COVERAGE_EMPIRICAL")
  dtbl_subset <- subset(dtbl, select = var_names)
  join_keys <- c("SEX", "LOC", "AGE_LABEL")

  nsinf_by_age <- merge(dtbl_subset, nsinf_by_age, by = join_keys)
  nsinf_samples_by_age <- merge(dtbl_subset, tmp, by = join_keys)
  nspred_by_age <- merge(dtbl_subset, nspred_by_age, by = join_keys)

  # load change of var name
  set(nsinf_by_age, NULL, "SEX", NULL)
  set(nsinf_by_age, NULL, "LOC", NULL)
  setnames(
    nsinf_by_age,
    c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "M", "CL", "CU"),
    c("COMM", "SEX", "AGEYRS",
      "PROP_ART_COVERAGE_M",
      "PROP_ART_COVERAGE_CL",
      "PROP_ART_COVERAGE_CU")
  )
  nsinf_by_age[, ROUND := paste0("R0", round)]

  # load change of var name
  set(nspred_by_age, NULL, "SEX", NULL)
  set(nspred_by_age, NULL, "LOC", NULL)

  old_names <- c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "M", "CL", "CU")
  new_names <- c("COMM", "SEX", "AGEYRS",
                 "PROP_ART_COVERAGE_M",
                 "PROP_ART_COVERAGE_CL",
                 "PROP_ART_COVERAGE_CU")
  setnames(nspred_by_age, old_names, new_names)

  nspred_by_age[, ROUND := paste0("R0", round)]

  # load change of var name
  set(nsinf_samples_by_age, NULL, "SEX", NULL)
  set(nsinf_samples_by_age, NULL, "LOC", NULL)
  set(nsinf_samples_by_age, NULL, "TOTAL_COUNT", NULL)
  set(nsinf_samples_by_age, NULL, "COUNT_PREDICT", NULL)
  set(nsinf_samples_by_age, NULL, "ART_USE_PREDICT", NULL)
  set(nsinf_samples_by_age, NULL, "Var2", NULL)

  old_names <- c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "value")
  new_names <- c("COMM", "SEX", "AGEYRS", "PROP_ART_COVERAGE_POSTERIOR_SAMPLE")
  setnames(nsinf_samples_by_age, old_names, new_names)

  nsinf_samples_by_age[, ROUND := paste0("R0", round)]

  # keep
  nsinf[[i]] <- nsinf_by_age
  nspred[[i]] <- nspred_by_age
  nsinf_samples[[i]] <- nsinf_samples_by_age
  convergence_list[[i]] <- convergence
}

nsinf <- do.call("rbind", nsinf)
nsinf_samples <- do.call("rbind", nsinf_samples)
nspred <- do.call("rbind", nspred)
convergence <- do.call("rbind", convergence_list)

# check all entries are complete
stopifnot(
  nrow(nsinf[COMM == "inland"]) == nsinf[, length(unique(AGEYRS))] *
                                   nsinf[, length(unique(SEX))] *
                                   nsinf[COMM == "inland", length(unique(ROUND))] # nolint: line_length_linter.
)
stopifnot(
  nrow(nsinf[COMM == "fishing"]) == nsinf[, length(unique(AGEYRS))] *
                                    nsinf[, length(unique(SEX))] *
                                    nsinf[COMM == "fishing", length(unique(ROUND))] # nolint: line_length_linter.
)


###########################

# STATISTICS FOR PAPER #

###########################

# get proportion of predicted art use inside credible interval
stats <- list()
tmp <- nspred[!is.na(PROP_ART_COVERAGE_EMPIRICAL)]
tmp[, within.CI := data.table::between(PROP_ART_COVERAGE_EMPIRICAL,
                                       PROP_ART_COVERAGE_CL,
                                       PROP_ART_COVERAGE_CU)]
stats[["within.CI"]] <- tmp[, paste0(round(mean(within.CI) * 100, 2))]

# get lowest rhat and lowest neff
stats[["min_neff"]] <- convergence[, round(min(neff))]
stats[["max_rhat"]] <- convergence[, round(max(rhat), 4)]


#########

# SAVE #

#########

# samples
file.name <- file.selfreportedart.newly
if (!file.exists(file.name) || config$overwrite.existing.files) {
  saveRDS(nsinf_samples, file = file.name)
}

# stats
file.name <- file.path(outdir, "RCCS_art_model_fit_newlyregistered_230703.RDS")
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(stats, file = file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}
cat("\n Done \n")


