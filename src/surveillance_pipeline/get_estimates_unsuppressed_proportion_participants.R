library(data.table)
library(ggplot2)
library(Hmisc)
library(rstan)
library(here)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

# outdir directory for stan fit
outdir <- file.path("../phyloSI-RakaiAgeGender-outputs","get_estimates_unsuppressed_proportion_participants")
if(usr == 'melodiemonod'){
  outdir <- file.path("/Users/melodiemonod/Box Sync/2023//phyloSI-RakaiAgeGender-outputs","get_estimates_unsuppressed_proportion_participants")
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

file.exists(c(
  path_stan_binomialgp ,
  path.count.unsupp ))  |> all() |> stopifnot()

# Load count of participants with unsuppressed viral loads
vla <- fread(path.count.unsupp)

##########################################

# PLOT #

##########################################

if (1) {
  tmp <- vla[, .(ROUND, LOC_LABEL, SEX_LABEL, AGE_LABEL, HIV_N, VLNS_N)]
  tmp[, `Non viremic` := HIV_N - VLNS_N]
  setnames(tmp, "VLNS_N", "Viremic")
  tmp <- melt.data.table(
    tmp,
    id.vars = c("ROUND", "LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "HIV_N")
  )
  tmp[, variable := factor(variable, levels = c("Non viremic", "Viremic"))]
  setnames(tmp, "LOC_LABEL", "COMM")
  setnames(tmp, "SEX_LABEL", "SEX")
  tmp[, ROUND := as.character(ROUND)]
  tmp[ROUND == "15.5", ROUND := "15S"]
  tmp <- tmp[ROUND != "15S"]
  tmp[, ROUND_LABEL := paste0("Round ", ROUND)]
  tmp <- tmp[!(ROUND == "15S")]
  tmp[, SEX_LABEL := "Women"]
  tmp[SEX== "M", SEX_LABEL := "Men"]
  tmp[, COMM_LABEL := "Fishing\n communities"]
  tmp[COMM == "inland", COMM_LABEL := "Inland\n communities"]
  tmp <- tmp[AGE_LABEL > 14 & AGE_LABEL < 50]

  # plot
  p <- ggplot(tmp[COMM == "inland"], aes(x = AGE_LABEL, y = value)) +
    geom_bar(aes(fill = variable), stat = "identity") +
    labs(x = "Age", y = "Count HIV-positive participants", fill = "") +
    facet_grid(ROUND_LABEL ~ SEX_LABEL) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.background = element_rect(colour = "white", fill = "white"),
          strip.text = element_text(size = rel(1))) +
    scale_fill_manual(values = c("#9F73AB", "#432C7A"),
                      labels = c("Viral load <= 1,000 copies/mL", "Viral load > 1,000 copies/mL")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(p,
         file = file.path(outdir,
                          "count_unsuppressed_by_gender_loc_age_221104.pdf"),
         w = 7, h = 5.2)
}


##########################################

# FIND UNSUPPRESSED PROPORTION SMOOTH ESTIMATE #

##########################################

# find smooth proportion
for (round in 15:18) {
  dtbl <- copy(vla[ROUND == round])
  stopifnot(length(round) == 1)
  cat("Fitting stan model for round ", round, "\n")

  # predicts age
  x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL) + 1], 0.5)

  # make stan data
  stan_data <- list(
    x_predict = x_predict,
    y_observed_00 = dtbl[SEX == 0 & LOC == 0, HIV_N - VLNS_N],
    y_observed_10 = dtbl[SEX == 1 & LOC == 0, HIV_N - VLNS_N],
    y_observed_01 = dtbl[SEX == 0 & LOC == 1, HIV_N - VLNS_N],
    y_observed_11 = dtbl[SEX == 1 & LOC == 1, HIV_N - VLNS_N],
    total_observed_00 = dtbl[SEX == 0 & LOC == 0, HIV_N],
    total_observed_10 = dtbl[SEX == 1 & LOC == 0, HIV_N],
    total_observed_01 = dtbl[SEX == 0 & LOC == 1, HIV_N],
    total_observed_11 = dtbl[SEX == 1 & LOC == 1, HIV_N],
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
    iter = 10e3,
    warmup = 5e2,
    chains = 1,
    control = list(max_treedepth = 15, adapt_delta = 0.999)
  )
  filename <- paste0("220729f_notsuppAmongInfected_gp_stan_round",
                     round,
                     "_vl_1000.rds")
  filename <- file.path(outdir, filename)
  if (!file.exists(filename) || config$overwrite.existing.files) {
      saveRDS(fit, file = filename)
  }
}

###############

# LOAD RESULTS

###############

rounds <- 15:18
nsinf <- vector(mode = "list", length = length(rounds))
nsinf_samples <-  vector(mode = "list", length = length(rounds))
nspred <- vector(mode = "list", length = length(rounds))
convergence_list <- vector(mode = "list", length = length(rounds))
for (i in seq_along(rounds)) {

  round <- rounds[i]
  dtbl <- copy(vla[ROUND == round])

  # age to predict
  x_predict <- seq(vla[, min(AGE_LABEL)], vla[, max(AGE_LABEL) + 1], 0.5)

  # load samples
  filename <- paste0("220729f_notsuppAmongInfected_gp_stan_round",
                     round,
                     "_vl_1000.rds")
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
  #	summarise estimated unsuppressed by sex and age
  #

  # extract estimates
  ps <- c(0.025, 0.5, 0.975)
  qlab <- c("CL", "M", "CU")
  tmp <- as.data.table(reshape2::melt(1 - re$p_predict_00))
  tmp[, `:=`(SEX = 0, LOC = 0)]
  tmp1 <- as.data.table(reshape2::melt(1 - re$p_predict_10))
  tmp1[, `:=`(SEX = 1, LOC = 0)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(1 - re$p_predict_01))
  tmp1[, `:=`(SEX = 0, LOC = 1)]
  tmp <- rbind(tmp, tmp1)
  tmp1 <- as.data.table(reshape2::melt(1 - re$p_predict_11))
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
  #	summarise predicted unsuppressed by sex and age
  #

  # merge to total count
  tmp <- merge(tmp, dtbl[, .(SEX, LOC, AGE_LABEL, HIV_N)],
               by = c("SEX", "LOC", "AGE_LABEL"),
               all.x = TRUE)

  # predict count and predict unsuppressed
  tmp[!is.na(HIV_N), COUNT_PREDICT := rbinom(1, HIV_N, value),
      by = c("SEX", "LOC", "AGE_LABEL", "iterations")]
  tmp[, UNSUPPRESSED_PREDICT := COUNT_PREDICT / HIV_N]

  # summarise
  nspred_by_age <- tmp[, list(q = quantile(UNSUPPRESSED_PREDICT,
                                           prob = ps,
                                           na.rm = TRUE),
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


  #
  # POSTPROCESING
  #

  # merge to data
  nsinf_by_age <- merge(
    subset(dtbl, select = c(SEX,
                            SEX_LABEL,
                            LOC,
                            LOC_LABEL,
                            AGE_LABEL,
                            EMPIRICAL_VLNS_IN_HIV)),
    nsinf_by_age,
    by = c("SEX", "LOC", "AGE_LABEL")
  )
  nsinf_samples_by_age <- merge(
    subset(dtbl, select = c(SEX,
                            SEX_LABEL,
                            LOC,
                            LOC_LABEL,
                            AGE_LABEL,
                            EMPIRICAL_VLNS_IN_HIV)),
    tmp,
    by = c("SEX", "LOC", "AGE_LABEL")
  )
  nspred_by_age <- merge(
    subset(dtbl, select = c(SEX,
                            SEX_LABEL,
                            LOC,
                            LOC_LABEL,
                            AGE_LABEL,
                            EMPIRICAL_VLNS_IN_HIV)),
    nspred_by_age,
    by = c("SEX", "LOC", "AGE_LABEL")
  )

  # change of var name
  set(nsinf_by_age, NULL, "SEX", NULL)
  set(nsinf_by_age, NULL, "LOC", NULL)
  setnames(
    nsinf_by_age,
    c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "EMPIRICAL_VLNS_IN_HIV",
      "M", "CL", "CU"),
    c("COMM", "SEX", "AGEYRS", "PROP_UNSUPPRESSED_EMPIRICAL",
      "PROP_UNSUPPRESSED_M", "PROP_UNSUPPRESSED_CL", "PROP_UNSUPPRESSED_CU")
  )
  nsinf_by_age[, ROUND := paste0("R0", round)]

  # load change of var name
  set(nspred_by_age, NULL, "SEX", NULL)
  set(nspred_by_age, NULL, "LOC", NULL)
  setnames(nspred_by_age,
    c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "M", "CL", "CU"),
    c("COMM", "SEX", "AGEYRS",
      "PROP_UNSUPPRESSED_M", "PROP_UNSUPPRESSED_CL", "PROP_UNSUPPRESSED_CU")
  )
  nspred_by_age[, ROUND := paste0("R0", round)]

  # load change of var name
  set(nsinf_samples_by_age, NULL, "SEX", NULL)
  set(nsinf_samples_by_age, NULL, "LOC", NULL)
  set(nsinf_samples_by_age, NULL, "iterations_rev", NULL)
  set(nsinf_samples_by_age, NULL, "HIV_N", NULL)
  set(nsinf_samples_by_age, NULL, "Var2", NULL)
  set(nsinf_samples_by_age, NULL, "UNSUPPRESSED_PREDICT", NULL)
  set(nsinf_samples_by_age, NULL, "COUNT_PREDICT", NULL)
  setnames(
    nsinf_samples_by_age,
    c("LOC_LABEL", "SEX_LABEL", "AGE_LABEL", "EMPIRICAL_VLNS_IN_HIV", "value"),
    c("COMM", "SEX", "AGEYRS",
      "PROP_UNSUPPRESSED_EMPIRICAL", "PROP_UNSUPPRESSED_POSTERIOR_SAMPLE")
  )
  nsinf_samples_by_age[, ROUND := paste0("R0", round)]

  # keep
  nsinf[[i]] <- nsinf_by_age
  nsinf_samples[[i]] <- nsinf_samples_by_age
  nspred[[i]] <- nspred_by_age
  convergence_list[[i]] <- convergence
}
nsinf <- do.call("rbind", nsinf)
nsinf_samples <- do.call("rbind", nsinf_samples)
nspred <- do.call("rbind", nspred)
convergence <- do.call("rbind", convergence_list)


###########################

# STATISTICS FOR PAPER #

###########################

# get proportion of predicted art use inside credible interval
stats <- list()
tmp <- nspred[COMM == "inland" & !is.na(EMPIRICAL_VLNS_IN_HIV)]

tmp[, within.CI := between(EMPIRICAL_VLNS_IN_HIV,
                           PROP_UNSUPPRESSED_CL,
                           PROP_UNSUPPRESSED_CU)]
stats[["within.CI"]] <- tmp[, paste0(round(mean(within.CI) * 100, 2))]

# get lowest rhat and lowest neff
stats[["min_neff"]] <- convergence[, round(min(neff))]
stats[["max_rhat"]] <- convergence[, round(max(rhat), 4)]


#########

# SAVE #

#########

file.name <- file.unsuppressedviralload
if (!file.exists(file.name) || config$overwrite.existing.files) {
  saveRDS(nsinf_samples, file = file.name)
}


file.name <- file.path(outdir, "RCCS_nonsuppressed_proportion_model_fit_221101.RDS")
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  saveRDS(stats, file = file.name)
}else{
  cat("\n Output file", file.name, "already exists\n")
}
cat("\n Done \n")


