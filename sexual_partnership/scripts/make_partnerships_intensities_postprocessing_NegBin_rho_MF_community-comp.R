# Postprocessing for simulation age gender ----
# Preamble ----
# This script aims to analyse the outputs from the stan model (assess mixing and convergence) and do postprocessing for the outputs

# Load the required packages
require(data.table)
require(ggplot2)
require(cmdstanr)
require(bayesplot)
require(posterior)
require(tidyverse)
#
#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                        help = "Print extra output [default]"),
  optparse::make_option("--seed", type = "integer", default = 18L,
                        help = "Random number seed [default %default]",
                        dest = "seed"),
  optparse::make_option("--pkg_dir", type = "character", default = NA_character_,
                        help = "Absolute file path to package directory, used as long we don t build an R package [default]",
                        dest = "prj.dir"),
  optparse::make_option("--out_dir_base", type = "character", default = NA_character_,
                        help = "Absolute file path to base directory where all output is stored [default]",
                        dest = "out.dir")
)
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))


# Define input arguments that can be changed by users----
# args <- list()
# args$seed <- 18L

# define prj.dir and out.dir
tmp <- Sys.info()
if (tmp["user"] == "yc2819" & grepl("hpc.ic.ac.uk",tmp["nodename"])) # outdir yu
{
  if (is.na(args$out.dir))
  {
    args$out.dir <- "/rds/general/user/yc2819/home/github/phyloflows/sexual_partnership/results"
  }
  if (is.na(args$prj.dir))
  {
    args$prj.dir <- "/rds/general/user/yc2819/home/github/phyloflows/sexual_partnership"
  }
}
# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
  args$out.dir <- here::here()
}

# Source functions
source( file.path(args$prj.dir,'functions','plot_functions.R') )

# read simulated contact intensities, input args, simulated data to model from workspace
args$in.dir <- args$out.dir
infiles <- list.files(args$in.dir, pattern = 'area', full.names = TRUE)
outfile.base <- gsub('_stan_fit.rds','',infiles[grepl('stan_fit.rds',infiles)])

cat('\nFound files\n', paste0(infiles, collapse = '\n'))

# Load workspace----
tmp <- infiles[grepl('image.RData',infiles)]
cat('\n Reading workspace from file: ', tmp , '\n')
load( tmp, tmp_env <- new.env())
tmp_list <- as.list.environment(tmp_env)
reported.partnerships <- tmp_list$reported.partnerships

cat("\n Loading input args...")
tmp <- tmp_list$args
rm(tmp_list)
for (x in names(tmp))
{
  if (!(x %in% names(args)))
  {
    args[[x]] <- tmp[[x]]
  }
  # stopifnot(!x %in% names(args))
}

args$if_save_fig <- TRUE
cat("\n Input args are \n ")
str(args)

# Load the Stan output----
tmp <- infiles[grepl('stan_fit.rds',infiles)]
cat("\n Loading fitted data from ", tmp, "\n")
model_fit <- readRDS( file = tmp )

# Need seed for posterior predictions
set.seed(args$seed)

# Print the model used
cat("\nModel is used in stan model: ", args$stan_model, "\n")

# HMC samplers analysis ----
if ( grepl("BSGP", args$stan_model))
{
  su.target.vars <- c('log_random_effect_baseline', 'log_random_effect_community', 'log_random_effect_community_sd', 'overdispersion','gp_rho_age1','gp_rho_age2','gp_alpha','z', 'f_mf')
  model_id <- paste0('NegBin_BSGP_knots-', args$n_knots, '-comm-', args$area)
}
if ( grepl("HSGP", args$stan_model))
{
  su.target.vars <- c('log_random_effect_baseline', 'log_random_effect_community', 'log_random_effect_community_sd', 'overdispersion','gp_rho_age1','gp_rho_age2','gp_alpha','z', 'f_mf')
  model_id <- paste0('NegBin_HSGP_c-',args$hsgp_boundary_inflation*100,"_m-",args$hsgp_m, '-comm-', args$area)
}

# Rhat Neff----
su <- as.data.table(posterior::summarise_draws(
  model_fit$draws(
    variables = su.target.vars,
    inc_warmup = FALSE)
))
cat("\nThe minimal nESS is: ",su[,min(ess_bulk)], "\nThe maximal rhat is: ",su[,max(rhat)])
tmp <- file.path(paste0(outfile.base,"_convergence.csv"))
cat("\nWrite convergences to file ", tmp)
write.csv(su, file = tmp, row.names = TRUE)

# Diagnostics----
diagnostics <- posterior::summarise_draws(posterior::as_draws_df(
  model_fit$sampler_diagnostics()), ~quantile(.x, probs = c(0.01, 0.5, 0.99)
  ))
diagnostics <- rbind(diagnostics, c('min_ess_bulk', su[,min(ess_bulk)],su[,min(ess_bulk)],su[,min(ess_bulk)]))
diagnostics <- rbind(diagnostics, c('max_rhat', su[,max(rhat)],su[,max(rhat)],su[,max(rhat)]))

tmp <- file.path(paste0(outfile.base,"_diagnostics.csv"))
cat("\nWrite diagnostics to file ", tmp)
write.csv(diagnostics, file = tmp, row.names = TRUE)

# Worst 5 traces----
# traces must be done WITH warmup (so always use $sample with save_warmup = TRUE)
tmp <- su[order(ess_bulk),][1:5,][,variable]
bayesplot::color_scheme_set("mix-blue-red")
p <- bayesplot::mcmc_trace(
  model_fit$draws(variables = tmp, inc_warmup = TRUE),
  pars = tmp,
  n_warmup = model_fit$metadata()$iter_warmup, # TODO please note
  facet_args = list(ncol = 1, strip.position = "left")
) +
  theme(legend.position = 'bottom')
tmp <- file.path(paste0(outfile.base,"_worst5traces.pdf"))
cat("\nWrite worst 5 traces to file ", tmp)
ggsave(file = tmp, p, width = 10, height = 12)
tmp <- file.path(paste0(outfile.base,"_worst5traces.png"))
cat("\nWrite worst 5 traces to file ", tmp)
ggsave(file = tmp, p, width = 10, height = 12)

# Draw the df baseline and community-specific baseline----
p <- model_fit$draws() %>%
  as_draws_df() %>%
  as_tibble() %>%
  select(starts_with("log_random"))  %>%
  setNames(c('log_baseline', paste0('log_baseline_comm_', unique(reported.partnerships$part.comm_num)), 'log_basline_comm_sd')) %>%
  mcmc_areas()
tmp <- file.path(paste0(outfile.base,"_posteriordist.pdf"))
cat("\nSave posterior distribution to file ", tmp)
ggsave(file = tmp, p, width = 4, height = 1 * (1 + length(unique(reported.partnerships$part.comm_num))), limitsize = FALSE)
tmp <- file.path(paste0(outfile.base,"_posteriordist.png"))
cat("\nSave posterior distribution to file ", tmp)
ggsave(file = tmp, p, width = 4, height = 1 * (1 + length(unique(reported.partnerships$part.comm_num))), limitsize = FALSE)

# Computing time----
compute.time <- model_fit$metadata()$time
tmp <- file.path(paste0(outfile.base,"_computing-time.csv"))
cat("\nWrite computing time to file ", tmp)
write.csv(compute.time, file = tmp, row.names = TRUE)

# Extract Monte Carlo samples for postprocessing----
cat("\nExtract Monte Carlo samples for postprocessing ...\n")
su.target.vars <- c('log_random_effect_baseline','log_random_effect_community', 'f_mf','overdispersion')
pd <- model_fit$draws(variables = su.target.vars, inc_warmup = FALSE)

select.chains <- seq_along(dimnames(pd)[['chain']])
iters <- seq(from = 1, to = model_fit$metadata()[['iter_sampling']], length.out = ceiling(1e4 / length(select.chains)))
# iters <- seq_len(model_fit$metadata()[['iter_sampling']])
po <- list()
tmp <- pd[,,which(grepl('f_mf',dimnames(pd)[[3]]))]
po$f_mf <- unname(apply(tmp[iters,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('log_random_effect_baseline',dimnames(pd)[[3]]))]
po$log_random_effect_baseline <- unname(apply(tmp[iters,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('log_random_effect_community',dimnames(pd)[[3]]))]
po$log_random_effect_community <- unname(apply(tmp[iters,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('overdispersion',dimnames(pd)[[3]]))]
po$overdispersion <- as.vector(unname(apply(tmp[iters,select.chains,], 3, rbind)))

pd <- NULL
gc()

pds.quantiles <- c(.025,.25,.5,.75,.975)
pds.quantilelabels <- c('CL','IL','M','IU','CU')
# pds.MF <- 1 # from Stan model file
# Extract age and gender specific contact intensities----
cat("\nExtract age and gender specific contact intensities ...\n")
pd <- as.data.table(reshape2::melt(po$f_mf))
setnames(pd, 1:3, c('iterations','mf_mat_idx','value'))
tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                    log_random_effect_baseline = po$log_random_effect_baseline
                    )
setnames(tmp2, 1:2, c('iterations','log_random_effect_baseline'))
pd <- merge(pd, tmp2, by = 'iterations')
pd[, 'DUMMY' := 1]
tmp <- data.table( DUMMY = 1L, c = seq_len(ncol(po$log_random_effect_community)))
pd <- merge(tmp, pd, by = c('DUMMY'), allow.cartesian = TRUE)
set(pd, NULL, 'DUMMY', NULL)

tmp2 <- as.data.table(reshape2::melt(po$log_random_effect_community))
setnames(tmp2, 1:3, c('iterations','c','log_random_effect_community'))

pd_save <- merge(pd, tmp2, by = c('iterations', 'c'))

tmp2 <- subset(reported.partnerships, part.sex == 'M', select = c(mf_mat_idx,c, pop)) #, gamma))
pd <- merge(pd_save, tmp2, by = c('mf_mat_idx', 'c'))
set(pd, NULL, 'value_m', pd[, exp(log_random_effect_baseline + log_random_effect_community + value + log(pop))])
set(pd, NULL, 'value_m_all', pd[, exp(log_random_effect_baseline + log_random_effect_community + value + log(pop) + log(gamma))])
tmp <- pd[,
          list(sd = sd(value_m),
               value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
               stat = pds.quantilelabels),
          by = c('mf_mat_idx', 'c')
          ]
tmp <- dcast.data.table(tmp, mf_mat_idx + c + sd~stat, value.var = 'value_m')
tmp2 <- subset(reported.partnerships, part.sex == 'M', select = c(part.sex, part.age, cont.age, mf_mat_idx, capped_cntcts, capped_cntcts_in, part, pop, rho, part.comm_num, c))
pds <- merge(tmp, tmp2, by = c('mf_mat_idx', 'c'))
tmp3 <- pd[,
          list(sd = sd(value_m_all),
               value_m_all = quantile(value_m_all, p = pds.quantiles, na.rm = TRUE),
               stat = pds.quantilelabels),
          by = c('mf_mat_idx', 'c')
          ]
tmp3 <- dcast.data.table(tmp3, mf_mat_idx +c + sd~stat, value.var = 'value_m_all')
pds.r <- merge(tmp3, tmp2, by = c('mf_mat_idx', 'c'))

# and only now we repeat for fm
tmp2 <- subset(reported.partnerships, part.sex == 'F', select = c(mf_mat_idx, pop, c, gamma))
pd <- merge(pd_save, tmp2, by = c('mf_mat_idx', 'c'))

set(pd, NULL, 'value_m', pd[, exp(log_random_effect_baseline + log_random_effect_community + value + log(pop))])

set(pd, NULL, 'value_m_all', pd[, exp(log_random_effect_baseline + value + log(pop) + log(gamma))])
tmp <- pd[,
          list(sd = sd(value_m),
               value_m = quantile(value_m, p = pds.quantiles, na.rm = TRUE),
               stat = pds.quantilelabels),
          by = c('mf_mat_idx', 'c')
          ]
tmp <- dcast.data.table(tmp, mf_mat_idx + c + sd~stat, value.var = 'value_m')
tmp2 <- subset(reported.partnerships, part.sex == 'F', select = c(part.sex, part.age, cont.age, c, part.comm_num, mf_mat_idx, capped_cntcts, capped_cntcts_in, part, pop, rho))
tmp <- merge(tmp, tmp2, by = c('mf_mat_idx', 'c'))
tmp3 <- pd[,
          list(sd = sd(value_m_all),
               value_m_all = quantile(value_m_all, p = pds.quantiles, na.rm = TRUE),
               stat = pds.quantilelabels),
          by = c('mf_mat_idx', 'c')
          ]
tmp3 <- dcast.data.table(tmp3, mf_mat_idx + sd~stat, value.var = 'value_m_all')
tmp3 <- merge(tmp3, tmp2, by = c('mf_mat_idx', 'c'))
pds <- rbind(pds, tmp, use.names = TRUE, fill = TRUE)
pds.r <- rbind(pds.r, tmp3, use.names = TRUE, fill = TRUE)

pds[, cntct_intensity_rho_scaled_empirical := capped_cntcts_in/part/rho]
pds[, cntct_intensity_empirical := capped_cntcts_in/part]
pds.r[, cntct_intensity_rho_scaled_empirical := capped_cntcts_in/part/rho]
pds.r[, cntct_intensity_empirical := capped_cntcts_in/part]
pds.r[, cntct_intensity_rho_scaled_empirical := capped_cntcts/part/rho]
pds.r[, cntct_intensity_empirical := capped_cntcts/part]

tmp <- file.path(paste0(outfile.base, "_contact-intensities-within-comm.rds"))
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds, file = tmp)
tmp <- file.path(paste0(outfile.base, "_contact-intensities-all.rds"))
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds.r, file = tmp)

# Visualize age and gender specific contact intensities----
cat("\nVisualize age and gender specific contact intensities ...\n")
pds <- pds[part.age %in% 15:49]
# pds.r <- pds.r[part.age %in% 15:49]

cols <- c("empirical" = "black","rho_scaled_empirical" = "purple")
shapes <- c("empirical" = "*","rho_scaled_empirical" = "+")

# Plot age posterior distribution----
if (args$if_save_fig)
{
  cat("\nPlot age posterior distribution ...")
  # we just visualize the contact intensities for male-female
  pm <- plot_age_contact_intensity_MF(pds[c == 1])
  pm <- pm + ggtitle(paste0('Age dist using ', model_id, ' community ', unique(pds[c == 1]$part.comm_num), " - MF" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-within-comm-MF.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-within-comm-MF.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}

if (args$if_save_fig)
{
  cat("\nPlot age posterior distribution with empirical contact intensities ...")
  # we just visualize the contact intensities for male-female
  pm <- plot_age_contact_intensity_MF_empirical(pds[c == 1])
  pm <- pm + ggtitle(paste0('Age dist with empirical contact intensities using ', model_id,' community ', unique(pds[c == 1]$part.comm_num), " - MF" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-MF-within-comm-with-empirical.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-MF-within-comm-with-empirical.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}



if (args$if_save_fig)
{
  cat("\nPlot age posterior distribution ...")
  # we just visualize the contact intensities for female-male
  pm <- plot_age_contact_intensity_FM(pds[c == 1])
  pm <- pm + ggtitle(paste0('Age dist using ', model_id,' community ', unique(pds[c == 1]$part.comm_num), " - FM" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-within-comm-FM.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-within-comm-FM.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}

if (args$if_save_fig)
{
  cat("\nPlot age posterior distribution with empirical contact intensities ...")
  # we just visualize the contact intensities for female-male
  pm <- plot_age_contact_intensity_FM_empirical(pds[c == 1])
  pm <- pm + ggtitle(paste0('Age dist with empirical contact intensities using ', model_id,' community ', unique(pds[c == 1]$part.comm_num), " - FM" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-FM-within-comm-with-empirical.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-FM-within-comm-with-empirical.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}
# Plot the selected age contact intensities for comparison ----
if ( args$if_save_fig)
{
  tmp <- pds[part.age %in% c(10, 20, 30, 40)]
  tmp$part.age <- as.character(tmp$part.age)
  tmp[, part.age := paste0("participants: ", part.age)]
  tmp[, id := model_id]
  tmp[, age_gender := paste0(part.age, '-', part.sex)]

  pm <- ggplot(tmp, aes(x = cont.age, col = age_gender)) +
    geom_step(aes(y = M)) +
    geom_point(aes(y = cntct_intensity_empirical), col = "black", shape = "*", size = 2) +
    geom_point(aes(y = cntct_intensity_rho_scaled_empirical), col = "purple", shape = "+", size = 2) +

    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = age_gender), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    facet_grid(id + paste0('Community ID: ', part.comm_num) ~age_gender, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") +
    scale_colour_manual(name = "Empirical values",values = cols) +
    scale_shape_manual(name = "Empirical values",values = shapes) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab("Age of contacted individuals") +
    ylab("Posterior mean contact intensities") +
    ggtitle(paste0('Selected age contact intensities dist within community with empirical intensities using ', model_id))

  pm <- pm + theme(axis.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   strip.text = element_text(size = 16)
                   # axis.text.x = element_text(angle = 45)
  )

  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-within-comm-comparison.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6 * length(unique(reported.partnerships$c)), limitsize = FALSE)
  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-within-comm-comparison.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6 * length(unique(reported.partnerships$c)), limitsize = FALSE)
}

# Selected age dist comparison between first 4 communities----
if ( args$if_save_fig)
{
  tmp <- pds[part.age %in% c(10, 20, 30, 40)]
  tmp$part.age <- as.character(tmp$part.age)
  tmp[, part.age := paste0("participants: ", part.age)]
  tmp[, id := model_id]
  tmp[, age_gender := paste0(part.age, '-', part.sex)]
  tmp <- tmp[c %in% 1:4]
  tmp$part.comm_num <- as.character(tmp$part.comm_num)
  pm <- ggplot(tmp, aes(x = cont.age, col = part.comm_num)) +
    geom_step(aes(y = M)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.comm_num), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    facet_grid(id ~age_gender, scales = "free_y") +
    theme_bw() +
    scale_colour_brewer(palette = "Set1") +
    guides(scales = "free_y") +
    guides(color = guide_legend("Community ID")) +
    theme(legend.position = "bottom") +
    xlab("Age of contacted individuals") +
    ylab("Posterior mean contact intensities by community") +
    ggtitle(paste0('Selected age contact intensities dist within community with empirical intensities using ', model_id))

  pm <- pm + theme(axis.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   strip.text = element_text(size = 16)
                   # axis.text.x = element_text(angle = 45)
  )

  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-within-comm-comparison-comm.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6)
  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-within-comm-comparison-comm.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6)
}
# Plot contact intensities patterns----
if (args$if_save_fig)
{
  cat("\nPlot contact intensities patterns ...")
  pds[, M := ifelse(M == 0, NaN, M)]
  pm <- plot_contact_intensity_matrix(pds[c == 1]) + labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "contact intensity" )
  pm <- pm + ggtitle(paste0('Posterior mean contacts intensities using ', model_id, ' community ', unique(pds[c == 1]$part.comm_num) ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-within-comm-patterns.pdf"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-within-comm-patterns.png"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
}

if (args$if_save_fig & 0)
{
  cat("\nPlot contact intensities patterns ...")
  pm <- plot_contact_intensity_sd_matrix(pds[c == 1]) + labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "standard deviation of the contact intensity" )
  pm <- pm + ggtitle(paste0('Standard deviation of the posterior contacts intensities using ', model_id, ' community ', unique(pds[c == 1]$part.comm_num) ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-sd-within-comm-patterns.pdf"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-sd-within-comm-patterns.png"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
}

# show the comparison with the empirical intensities
tmp.emp <- subset(pds, select = c('part.sex', 'part.age', 'cont.age', 'cntct_intensity_empirical', 'part.comm_num', 'c'))
setnames(tmp.emp, 'cntct_intensity_empirical', 'M')
set(tmp.emp, NULL, "label", "empirical value")
set(pds, NULL, "label", model_id)
dm <- rbind(pds, tmp.emp, fill = TRUE)
if (args$if_save_fig)
{
  pm <- plot_contact_intensity_matrix_empirical(dm[c == 1]) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "contact intensity" )

  pm <- pm +  ggtitle(paste0('Contact intensities patterns comparison using ', model_id, ' community ', unique(dm[c == 1]$part.comm_num)  ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns-within-comm-with-truth.pdf"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns-within-comm-with-truth.png"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
}

# show the comparison with the rho scaled empirical intensities
tmp.emp <- subset(pds, select = c('part.sex', 'part.age', 'cont.age', 'cntct_intensity_rho_scaled_empirical', 'c', 'part.comm_num'))
setnames(tmp.emp, 'cntct_intensity_rho_scaled_empirical', 'M')
set(tmp.emp, NULL, "label", "scaled empirical value")
set(pds, NULL, "label", model_id)
dm <- rbind(pds, tmp.emp, fill = TRUE)
if (args$if_save_fig)
{
  pm <- plot_contact_intensity_matrix_empirical(dm[c == 1]) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "contact intensity" )

  pm <- pm +  ggtitle(paste0('Contact intensities patterns comparison with rho scaled empirical intensities using ', model_id, ' community ', unique(dm[c == 1]$part.comm_num)  ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns-within-comm-with-rho-truth.pdf"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns-within-comm-with-rho-truth.png"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
}

# Plot intensities within and without commu----
if (args$if_save_fig & 0)
{
  cat("\nPlot age posterior distribution ...")
  # we just visualize the contact intensities for all communities for male-female
  pm <- plot_age_contact_intensity_MF(pds.r)
  pm <- pm + ggtitle(paste0('Age dist for all communities using ', model_id, " - MF" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-MF.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-MF.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}

if (args$if_save_fig & 0)
{
  cat("\nPlot age posterior distribution with empirical contact rates ...")
  # we just visualize the contact rates for male-female
  pm <- plot_age_contact_intensity_MF_empirical(pds.r)
  pm <- pm + ggtitle(paste0('Age dist with empirical contact intensities for all communities using ', model_id, " - MF" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-MF-with-empirical.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-MF-with-empirical.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}



if (args$if_save_fig & 0)
{
  cat("\nPlot age posterior distribution ...")
  # we just visualize the contact intensities for female-male
  pm <- plot_age_contact_intensity_FM(pds.r)
  pm <- pm + ggtitle(paste0('Age dist for all communities using ', model_id, " - FM" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-FM.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-FM.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}

if (args$if_save_fig & 0)
{
  # Participants and the contacted individuals aged 0 - 17 are missing in the BICs survey

  cat("\nPlot age posterior distribution with empirical contact rates ...")
  # we just visualize the contact intensities for female-male
  pm <- plot_age_contact_intensity_FM_empirical(pds.r)
  pm <- pm + ggtitle(paste0('Age dist with empirical contact intensities for all communities using ', model_id, " - FM" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-FM-with-empirical.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-all-comm-FM-with-empirical.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}
# Plot contact patterns all communities----
if (args$if_save_fig & 0)
{
  cat("\nPlot contact rate patterns ...")
  pds.r[, M := ifelse(M == 0, NaN, M)]
  pm <- plot_contact_intensity_matrix(pds.r) + labs(x = "Age of participant", y = "Age of contact", fill = "contact rate" )
  pm <- pm + theme(legend.text = element_text(angle = 45, vjust = 0.5, hjust = 1))
  pm <- pm + ggtitle(paste0('Posterior mean contacts intensities for all communities using ', model_id ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-all-comm-patterns.pdf"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-all-comm-patterns.png"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
}

if (0 & args$if_save_fig)
{
  cat("\nPlot the standard deviation of the contact rate patterns ...")
  pm <- plot_contact_rate_sd_matrix(pds.r) + labs(x = "Age of participant", y = "Age of contact", fill = "sd of the contact rate" )
  pm <- pm + ggtitle(paste0('Standard deviation of the posterior contacts rates using ', model_id )) +
    theme(legend.text = element_text(angle = 45, vjust = 0.5, hjust = 1))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-all-comm-patterns.pdf"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-all-comm-patterns.png"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
}


if (args$if_save_fig & 0)
{
  tmp.emp <- subset(pds.r, select = c('part.sex', 'part.age', 'cont.age', 'cntct_intensity_empirical'))
  setnames(tmp.emp, 'cntct_intensity_empirical', 'M')
  set(tmp.emp, NULL, "label", "empirical value")
  set(pds.r, NULL, "label", model_id)
  dm <- rbind(pds.r, tmp.emp, fill = TRUE)

  pm <- plot_contact_intensity_matrix_empirical(dm)
  pm <- pm +  ggtitle(paste0('Contact intensities patterns comparison for all communities using ', model_id  ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-all-comm-patterns-with-truth.pdf"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-all-comm-patterns-with-truth.png"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
}


# Compare the intensities within and all community ----
# Plot the selected age contact intensities for comparison ----
if ( args$if_save_fig & 0 )
{
  pds$label <- 'within-community'
  pds.r$label <- 'all'
  tmp <- rbind(pds, pds.r)

  tmp <- tmp[part.age %in% c(10, 20, 30, 40)]
  tmp$part.age <- as.character(tmp$part.age)
  tmp[, part.age := paste0("participants: ", part.age)]
  tmp[, id := model_id]
  tmp[, age_gender := paste0(part.age, '-', part.sex)]

  pm <- ggplot(tmp, aes(x = cont.age, col = label)) +
    geom_step(aes(y = M)) +
    # geom_point(aes(y = cntct_intensity_empirical), col = "black", shape = "*", size = 2) +
    # geom_point(aes(y = cntct_intensity_rho_scaled_empirical), col = "purple", shape = "+", size = 2) +
    #
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = label), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    facet_grid(id ~age_gender, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") +
    # scale_colour_manual(name = "Empirical values",values = cols) +
    # scale_shape_manual(name = "Empirical values",values = shapes) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlab("Age of contacted individuals") +
    ylab("Posterior mean contact intensities") +
    ggtitle(paste0('Seletcted age contact intensities dist comparison within or all communities using ', model_id  ))

  pm <- pm + theme(axis.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   strip.text = element_text(size = 16)
                   # axis.text.x = element_text(angle = 45)
  )

  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-comm-comparison.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6)
  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-comm-comparison.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6)
}


# Extract age and gender specific marginal contact intensities----
#   here we will process the Monte Carlo output in batches to avoid handling tables with >100m rows
cat("\nExtract age and gender specific marginal contact intensities ...\n")
# dcm <- subset(reported.partnerships, part.sex != cont.sex)
# setkey(dcm, part.sex, part.age)
tmp <- unique(subset(reported.partnerships, part.age <= 49, select = c(part.sex, part.age, c)))
setkey(tmp, part.sex, part.age)
tmp[, m_index := seq_len(length(part.age)), by = c('c')]
dcm <- merge(reported.partnerships[part.age <= 49], tmp, by = c('part.sex', 'part.age', 'c'))
pds.midx <- unique(tmp$m_index)

pds <- vector('list',length(pds.midx))
for (i in pds.midx)
{
  cat("\nProcessing ",i," of ",length(pds.midx))
  gender_i <- subset(dcm, m_index == i )[, part.sex][1]
  mf_mat_idx <- unique(subset(dcm, m_index == i & cont.sex != gender_i)[, mf_mat_idx])

  if (gender_i == 'F')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline
    )
    setnames(tmp2, 2, c('log_random_effect_baseline'))
    pd <- merge(pd, tmp2, by = 'iterations')
    pd[, 'DUMMY' := 1]
    tmp <- data.table( DUMMY = 1L, c = seq_len(ncol(po$log_random_effect_community)))
    pd <- merge(tmp, pd, by = c('DUMMY'), allow.cartesian = TRUE)
    set(pd, NULL, 'DUMMY', NULL)

    tmp2 <- as.data.table(reshape2::melt(po$log_random_effect_community))
    setnames(tmp2, 1:3, c('iterations','c','log_random_effect_community'))
    pd <- merge(pd, tmp2, by = c('iterations', 'c'))
    # pd[, log_contact_rate := log_random_effect_baseline + log_random_effect_community + f_mf]
    tmp <- subset(dcm, m_index == i, select = c(m_index, part.sex, cont.age, part.age, mf_mat_idx, pop, c))
    pd <- merge(pd, tmp, by = c('mf_mat_idx', 'c'))
    pd[, log_contact_intensity_within := log_random_effect_baseline + log_random_effect_community + f_mf + log(pop) ]
    # pd[, log_contact_intensity_all := log_random_effect_baseline + f_mf + log(pop) + log(gamma)]

    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','log_random_effect_community','f_mf','pop'), NULL)
    pds[[i]] <- copy(pd)

    pds[[i]] <- pds[[i]][, .(log_contact_intensity_within = log(sum(exp(log_contact_intensity_within), na.rm = TRUE)))#,
                             # log_contact_intensity_all = log(sum(exp(log_contact_intensity_all))))
                         , by = c('iterations','part.sex','part.age', 'c')]
  }

  if (gender_i == 'M')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline
    )
    setnames(tmp2, 2, c('log_random_effect_baseline'))
    pd <- merge(pd, tmp2, by = 'iterations')
    pd[, 'DUMMY' := 1]
    tmp <- data.table( DUMMY = 1L, c = seq_len(ncol(po$log_random_effect_community)))
    pd <- merge(tmp, pd, by = c('DUMMY'), allow.cartesian = TRUE)
    set(pd, NULL, 'DUMMY', NULL)

    tmp2 <- as.data.table(reshape2::melt(po$log_random_effect_community))
    setnames(tmp2, 1:3, c('iterations','c','log_random_effect_community'))
    pd <- merge(pd, tmp2, by = c('iterations', 'c'))
    # pd[, log_contact_rate := log_random_effect_baseline + log_random_effect_community + f_mf]

    tmp <- subset(dcm, m_index == i, select = c(m_index, part.sex, cont.age, part.age, mf_mat_idx, pop, c))
    pd <- merge(pd, tmp, by = c('mf_mat_idx', 'c'))
    pd[, log_contact_intensity_within := log_random_effect_baseline + log_random_effect_community + f_mf + log(pop) ]
    # pd[, log_contact_intensity_all := log_random_effect_baseline + f_mf + log(pop) + log(gamma)]
    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','f_mf','pop', 'log_random_effect_community'), NULL)
    pds[[i]] <- copy(pd)

    # marginalise with log sum exp
    pds[[i]] <- pds[[i]][, .(log_contact_intensity_within = log(sum(exp(log_contact_intensity_within), na.rm = TRUE))), #,
                             # log_contact_intensity_all = log(sum(exp(log_contact_intensity_all))))
                         by = c('iterations','part.sex','part.age', 'c')]
  }
}
pds <- data.table::rbindlist( pds )

# now summarise
# pds.all <- pds[,
#            list(log_contact_intensity_all = quantile(log_contact_intensity_all, p = pds.quantiles, na.rm = TRUE),
#                 stat = pds.quantilelabels),
#            by = c('part.sex','part.age')
# ]
# pds.all[, contact_intensity := exp(log_contact_intensity_all)]
# pds.all <- dcast.data.table(pds.all, part.sex+part.age~stat, value.var = 'contact_intensity')

pds <- pds[,
               list(log_contact_intensity_within = quantile(log_contact_intensity_within, p = pds.quantiles, na.rm = TRUE),
                    stat = pds.quantilelabels),
               by = c('part.sex','part.age', 'c')
]
pds[, contact_intensity := exp(log_contact_intensity_within)]
pds <- dcast.data.table(pds, part.sex+part.age+c~stat, value.var = 'contact_intensity')

tmp <- dcm[part > 0,
           list(cntct_intensity_empirical = sum(capped_cntcts_in / part, na.rm = TRUE),
                cntct_intensity_rho_scaled_empirical = sum(capped_cntcts_in / part/rho, na.rm = TRUE)),
           by = c('part.sex','part.age', 'c')
]

# tmp2 <- dcm[part > 0,
#            list(cntct_intensity_empirical = sum(capped_cntcts / part),
#                 cntct_intensity_rho_scaled_empirical = sum(capped_cntcts / part/rho)),
#            by = c('part.sex','part.age')
# ]
# pds.all <- merge(pds.all, tmp, by = c('part.sex','part.age'), all.x = TRUE)
# tmp <- file.path(paste0(outfile.base,"_marginal-contact-intensities-all.rds"))
# cat("\nWrite marginal contact intensities to file ", tmp)
# saveRDS(pds.all, file = tmp)

pds <- merge(pds, tmp, by = c('part.sex','part.age', 'c'), all.x = TRUE)
tmp <- file.path(paste0(outfile.base,"_marginal-contact-intensities-with.rds"))
cat("\nWrite marginal contact intensities to file ", tmp)
saveRDS(pds, file = tmp)
# Visualize age and gender specific contact intensities----
cat("\nVisualize marginal contact intensities ...\n")
pds <- pds[part.age %in% 15:49]
pds <- merge(pds, unique(subset(reported.partnerships, select = c('c', 'part.comm_num'))), by = c('c'))
# pds.all <- pds.all[part.age %in% 15:49]

# Plot the marginal contact intensities within community ----
if (args$if_save_fig)
{
  pm <- ggplot(pds, aes(x = part.age)) +
    geom_step(aes(y = M, col = part.sex)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.sex), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = cntct_intensity_rho_scaled_empirical), col = "purple", linetype = "dashed") +

    geom_step(aes(y = cntct_intensity_empirical), col = "black", linetype = "dashed") +
    facet_grid(paste0("community ID: ", part.comm_num) ~ part.sex) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age of contacting individual', y = 'Marginal contact intensity') +
    theme_bw() +
    scale_colour_manual(name = "Empirical values",values = cols) +
    scale_shape_manual(name = "Empirical values",values = shapes) +
    ggtitle(paste0('Posterior age marginal distribution within community using ', model_id  )) +
    theme(legend.position = "bottom")#, strip.background = element_rect(fill = "transparent"))
  pm <- pm + scale_color_manual(name = 'empirical type',
                                # breaks = c('empirical value', 'rho scaled empirical value'),
                                values = c('empirical value' = 'black', 'rho scaled empirical value' = 'purple'))
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities-within.pdf"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6 * length(unique(reported.partnerships$c)), limitsize = FALSE)
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities-within.png"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6 * length(unique(reported.partnerships$c)), limitsize = FALSE)
}

# Plot the marginal contact intensities within community ----
if (args$if_save_fig & 0)
{
  pm <- ggplot(pds.all, aes(x = part.age)) +
    geom_step(aes(y = M, col = part.sex)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.sex), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = cntct_intensity_rho_scaled_empirical), col = "purple", linetype = "dashed") +

    geom_step(aes(y = cntct_intensity_empirical), col = "black", linetype = "dashed") +
    facet_grid(~part.sex) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age of contacting individual', y = 'Marginal contact intensity') +
    theme_bw() +
    scale_colour_manual(name = "Empirical values",values = cols) +
    scale_shape_manual(name = "Empirical values",values = shapes) +
    ggtitle(paste0('Posterior age marginal distribution for all communities using ', model_id  )) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "transparent"))
  pm <- pm + scale_color_manual(name = 'empirical type',
                                # breaks = c('empirical value', 'rho scaled empirical value'),
                                values = c('empirical value' = 'black', 'rho scaled empirical value' = 'purple'))
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities-all.pdf"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 10)
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities-all.png"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 10)
}

# Plot the marginal contact intensities community comparison ----
if (args$if_save_fig & 0)
{
  pds$label <- 'within-comm'
  pds.all$label <- 'all'
  tmp <- rbind(pds, pds.all)
  pm <- ggplot(tmp, aes(x = part.age)) +
    geom_step(aes(y = M, col = label)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = label), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = cntct_intensity_rho_scaled_empirical), col = "purple", linetype = "dashed") +

    geom_step(aes(y = cntct_intensity_empirical), col = "black", linetype = "dashed") +
    facet_grid(~part.sex) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age of contacting individual', y = 'Marginal contact intensity') +
    theme_bw() +
    scale_colour_manual(name = "Empirical values",values = cols) +
    scale_shape_manual(name = "Empirical values",values = shapes) +
    ggtitle(paste0('Posterior age marginal distribution within and all community comparison using ', model_id  )) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "transparent"))
  pm <- pm + scale_color_manual(name = 'empirical type',
                                # breaks = c('empirical value', 'rho scaled empirical value'),
                                values = c('empirical value' = 'black', 'rho scaled empirical value' = 'purple'))
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities-comparison.pdf"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 10)
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities-comparison.png"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 10)
}

# Posterior predictive check on age and gender specific contact intensities----
#   here we will process the Monte Carlo output in batches to avoid handling tables with >100m rows
cat("\nPosterior predictive check on age and gender specific contact intensities ...\n")
dcm <- subset(reported.partnerships, part > 0 & part.age <= 49)
setkey(dcm, part.sex, part.age)
tmp <- unique(subset(dcm, select = c(part.sex, part.age, c)))
tmp[, m_index := seq_len(length(part.age)), by = 'c']
dcm <- merge(dcm, tmp, by = c('part.sex', 'part.age', 'c'))
pds.midx <- unique(tmp$m_index)

pds <- vector('list',length(pds.midx))
pds_marg <- vector('list',length(pds.midx))
for (i in pds.midx)
{
  #i <- 45
  cat("\nProcessing ",i," of ",length(pds.midx))
  gender_i <- subset(dcm, m_index == i )[, part.sex][1]
  mf_mat_idx <- unique(subset(dcm, m_index == i & cont.sex != gender_i)[, mf_mat_idx])

  if (gender_i == 'F')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline
    )
    setnames(tmp2, 2, c('log_random_effect_baseline'))
    pd <- merge(pd, tmp2, by = 'iterations')
    pd[, 'DUMMY' := 1]
    tmp <- data.table( DUMMY = 1L, c = seq_len(ncol(po$log_random_effect_community)))
    pd <- merge(tmp, pd, by = c('DUMMY'), allow.cartesian = TRUE)
    set(pd, NULL, 'DUMMY', NULL)

    tmp2 <- as.data.table(reshape2::melt(po$log_random_effect_community))
    setnames(tmp2, 1:3, c('iterations','c','log_random_effect_community'))
    pd <- merge(pd, tmp2, by = c('iterations', 'c'))
    tmp <- subset(dcm, m_index == i, select = c(m_index, part.sex, part.age, cont.age, mf_mat_idx, part, pop, rho, c))
    pd <- merge(pd, tmp, by = c('mf_mat_idx', 'c'))
    pd[, log_exp_contacts := log_random_effect_baseline + log_random_effect_community + f_mf + log(part) + log(pop) + log(rho)]
    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','log_random_effect_community','f_mf','part', 'pop', 'rho'), NULL)
    pds[[i]] <- copy(pd)
    pds_marg[[i]] <- pds[[i]][, .(log_exp_contacts = log(sum(exp(log_exp_contacts), na.rm = TRUE))), by = c('iterations','part.sex','part.age', 'c')]

  }
  if (gender_i == 'M')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline)
    setnames(tmp2, 2, c('log_random_effect_baseline'))
    pd <- merge(pd, tmp2, by = 'iterations')
    pd[, 'DUMMY' := 1]
    tmp <- data.table( DUMMY = 1L, c = seq_len(ncol(po$log_random_effect_community)))
    pd <- merge(tmp, pd, by = c('DUMMY'), allow.cartesian = TRUE)
    set(pd, NULL, 'DUMMY', NULL)
    tmp2 <- as.data.table(reshape2::melt(po$log_random_effect_community))
    setnames(tmp2, 1:3, c('iterations','c','log_random_effect_community'))
    pd <- merge(pd, tmp2, by = c('iterations', 'c'))

    tmp <- subset(dcm, m_index == i, select = c(m_index, part.sex, part.age, cont.age, mf_mat_idx, part, pop, rho, c))
    pd <- merge(pd, tmp, by = c('mf_mat_idx', 'c'))
    pd[, log_exp_contacts := log_random_effect_baseline + log_random_effect_community + f_mf + log(part) + log(pop) + log(rho) ]
    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','f_mf','part','pop','rho','log_random_effect_community'), NULL)
    pds[[i]] <- copy(pd)
    pds_marg[[i]] <- pds[[i]][, .(log_exp_contacts = log(sum(exp(log_exp_contacts), na.rm = TRUE))), by = c('iterations','part.sex','part.age', 'c')]


  }
  tmp2 <- data.table( iterations = seq_len(length(po$overdispersion)),
                      overdispersion = po$overdispersion
  )
  pds[[i]] <- merge(pds[[i]], tmp2, by = 'iterations')

  # conversion Stan <-> R is 1-prob == 1/(beta+1) and size == alpha
  pds[[i]][, beta := 1/overdispersion]
  pds[[i]][, cntcts_ppr := rnbinom( nrow(pds[[i]]),
                                    size = pds[[i]][, exp(log_exp_contacts) * beta ],
                                    prob = pds[[i]][, beta/(beta + 1)]
  )
  ]


  pds_marg[[i]] <- merge(pds_marg[[i]], tmp2, by = 'iterations')

  # conversion Stan <-> R is 1-prob == 1/(beta+1) and size == alpha
  pds_marg[[i]][, beta := 1/overdispersion]


  pds_marg[[i]][, cntcts_ppr := rnbinom( nrow(pds_marg[[i]]),
                                    size = pds_marg[[i]][, exp(log_exp_contacts) * beta ],
                                    prob = pds_marg[[i]][, beta/(beta + 1)]
  )
  ]
  # summarise

  pds[[i]] <- pds[[i]][,
                       list(cntcts_ppr = quantile(cntcts_ppr, p = pds.quantiles, na.rm = TRUE),
                            stat = pds.quantilelabels),
                       by = c('part.sex', 'part.age', 'cont.age', 'c')
  ]
}
pds <- data.table::rbindlist( pds )
pds <- dcast.data.table(pds, part.sex+part.age+cont.age+c~stat, value.var = 'cntcts_ppr')
tmp <- subset(dcm, select = c('part.sex','part.age','cont.age','capped_cntcts','capped_cntcts_in', 'c', 'part.comm_num'))
pds <- merge(pds, tmp, by = c('part.sex','part.age','cont.age','c'))

tmp <- file.path(paste0(outfile.base,"_posterior-predictions.rds"))
cat("\nWrite posterior predictions to file ", tmp)
saveRDS(pds, file = tmp)

pds[, if_out := (capped_cntcts_in < CL | CU < capped_cntcts_in)]
tmp <- mean(pds[, as.integer(if_out)], na.rm = TRUE)
tmp <- paste0(round(tmp*100, 2), '%')
cat("\nProportion of data points outside the 95% posterior prediction interval = ", tmp)

tmp.m <- mean(pds[part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
tmp.m <- paste0(round(tmp.m*100, 2), '%')
cat("\nProportion of male reported data points outside the 95% posterior prediction interval = ", tmp.m)

# Plot the posterior predictive check----
# A simple plot showing the gender age maps (x-axis age, y-axis age, panel for gender combination)
# and as value show "red" if data outside 95% prediction interval, "grey" if inside, and "transparent" if no data
cat("\n Plot the posterior predictive check contacts patterns whether the contact data is outside 95% prediction interval ...")
if (args$if_save_fig)
{
  pm <- vector('list',length(unique(reported.partnerships$c)))
  names(pm) <- paste0('Community ID: ', unique(reported.partnerships$part.comm_num))
  for (i in seq_along(unique(reported.partnerships$c)))
  {
    tmp <- mean(pds[c == i, as.integer(if_out)], na.rm = TRUE)
    tmp <- paste0(round(tmp*100, 2), '%')

    tmp.m <- mean(pds[c == i & part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
    tmp.m <- paste0(round(tmp.m*100, 2), '%')

    pm[[i]] <- ggplot(pds[c == i]) +
    geom_tile(aes(x = part.age, y = cont.age, fill = if_out)) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "If the true contact data outside 95% P.I." ) +
    coord_equal() +
    facet_grid( ~ paste0("Participants' sex: ", part.sex)) +
    scale_fill_manual(values = c("grey", "red"), na.value = 'transparent') +
    theme_bw() +
    theme(legend.position = 'bottom') +
      ggtitle(paste0('Community ', unique(reported.partnerships[c == i]$part.comm_num), ': ', tmp, ' of data points outside the 95% posterior prediction interval',
                     '\n', tmp.m, ' of Male reported data outside the 95% posterior prediction interval'))

    tmp <- file.path(paste0(outfile.base,'_posterior-predictive-check-contact-patterns-', unique(reported.partnerships[c == i]$part.comm_num), '.pdf'))
    cat("\nSave posterior predictive check contact intensities to file ", tmp)
    ggsave(file = tmp, pm[[i]], width = 12, height = 6 )
    tmp <- file.path(paste0(outfile.base,'_posterior-predictive-check-contact-patterns-', unique(reported.partnerships[c == i]$part.comm_num), '.png'))
    cat("\nSave posterior predictive check contact intensities to file ", tmp)
    ggsave(file = tmp, pm[[i]], width = 12, height = 6 )

  }
  p <- ggpubr::ggarrange(plotlist = pm,
                         nrow = length(unique(reported.partnerships$c)),
                         widths = 12,
                         heights = 6 * length(unique(reported.partnerships$c)),
                         common.legend = TRUE,
                         legend = "none"
                         )
  p <- p %>%
    gridExtra::grid.arrange(ggpubr::get_legend(pm[[1]]), ncol = 1,  widths = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns.pdf"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, p, width = 12, height = 6 * length(unique(reported.partnerships$c)))
  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns.png"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, p, width = 12, height = 6 * length(unique(reported.partnerships$c)))

}
cat("\n Plot the  comparison between the actual reported contacts and the estimated reported contact for each community ...")
if (args$if_save_fig)
{
  tmp <- subset(pds, select = c('part.age', 'part.sex', 'cont.age', 'capped_cntcts_in', 'c'))
  setnames(tmp, 'capped_cntcts_in', 'M')
  set(tmp, NULL, "label", "reported contacts")
  set(pds, NULL, "label", 'estimated contacts')
  dm <- rbind(pds, tmp, fill = TRUE)
  pds[, diff := M - capped_cntcts_in]
  set(pds, NULL, "label", 'estimated - reported')


  pm <- vector('list',length(unique(reported.partnerships$c)))
  names(pm) <- paste0('Community ID: ', unique(reported.partnerships$part.comm_num))
  pmav <- vector('list',length(unique(reported.partnerships$c)))
  names(pmav) <- paste0('Community ID: ', unique(reported.partnerships$part.comm_num))

  for (i in seq_along(unique(reported.partnerships$c)))
  {

    pm[[i]] <- ggplot(dm[c == i]) +
      geom_tile(aes(x = part.age, y = cont.age, fill = M)) +
      labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "The estimated or actual contacts" ) +
      coord_equal() +
      facet_grid(label ~ paste0("Participants' sex: ", part.sex)) +
      scale_fill_continuous(type = "viridis", na.value = 'transparent') +
      theme_bw() +
      theme(legend.position = 'bottom') +
      ggtitle(paste0('Comparison between the actual reported contacts and the estimated contacts within the community ', unique(reported.partnerships[c == i]$part.comm_num)))

    tmp <- file.path(paste0(outfile.base,'_posterior-predictive-check-contact-patterns-comparison-', unique(reported.partnerships[c == i]$part.comm_num), '.pdf'))
    cat("\nSave posterior predictive check contact intensities to file ", tmp)
    ggsave(file = tmp, pm[[i]], width = 12, height = 12 )
    tmp <- file.path(paste0(outfile.base,'_posterior-predictive-check-contact-patterns-comparison-', unique(reported.partnerships[c == i]$part.comm_num), '.png'))
    cat("\nSave posterior predictive check contact intensities to file ", tmp)
    ggsave(file = tmp, pm[[i]], width = 12, height = 12 )

    pmav[[i]] <- ggplot(pds[c == i]) +
      geom_tile(aes(x = part.age, y = cont.age, fill = diff)) +
      labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "Estimated contacts - reported contacts" ) +
      coord_equal() +
      facet_grid( ~ paste0("Participants' sex: ", part.sex)) +
      scale_fill_continuous(type = "viridis", na.value = 'transparent') +
      theme_bw() +
      theme(legend.position = 'bottom') +
      ggtitle(paste0('Comparison between the actual reported contacts and the estimated contacts within the community ', unique(reported.partnerships[c == i]$part.comm_num)))

    tmp <- file.path(paste0(outfile.base,'_posterior-predictive-check-contact-patterns-diff-', unique(reported.partnerships[c == i]$part.comm_num), '.pdf'))
    cat("\nSave posterior predictive check contact intensities to file ", tmp)
    ggsave(file = tmp, pmav[[i]], width = 12, height = 6 )
    tmp <- file.path(paste0(outfile.base,'_posterior-predictive-check-contact-patterns-diff-', unique(reported.partnerships[c == i]$part.comm_num), '.png'))
    cat("\nSave posterior predictive check contact intensities to file ", tmp)
    ggsave(file = tmp, pmav[[i]], width = 12, height = 6 )


  }
  p <- ggpubr::ggarrange(plotlist = pm,
                         nrow = length(unique(reported.partnerships$c)),
                         widths = 12,
                         heights = 6 * length(unique(reported.partnerships$c)),
                         common.legend = TRUE,
                         legend = "none"
  )
  p <- p %>%
    gridExtra::grid.arrange(ggpubr::get_legend(pm[[1]]), ncol = 1,  widths = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns.pdf"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, p, width = 12, height = 6 * length(unique(reported.partnerships$c)))
  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns.png"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, p, width = 12, height = 6 * length(unique(reported.partnerships$c)))

}
pds <- pds[, list(capped_cntcts_in = mean(capped_cntcts_in, na.rm = TRUE),
                            CL = mean(CL, na.rm = TRUE),
                            CU = mean(CU, na.rm = TRUE),
                            M = mean(M, na.rm = TRUE)),
                     by = c('part.age', 'part.sex', 'cont.age')]
pds[, if_out := (capped_cntcts_in < CL | CU < capped_cntcts_in)]
tmp <- mean(pds[, as.integer(if_out)], na.rm = TRUE)
tmp <- paste0(round(tmp*100, 2), '%')
cat("\nProportion of data points outside the 95% posterior prediction interval = ", tmp)

tmp.m <- mean(pds[part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
tmp.m <- paste0(round(tmp.m*100, 2), '%')
cat("\nProportion of male reported data points outside the 95% posterior prediction interval = ", tmp.m)

if (args$if_save_fig)
{
  pm <- ggplot(pds) +
    geom_tile(aes(x = part.age, y = cont.age, fill = if_out)) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "If the true contact data outside 95% P.I." ) +
    coord_equal() +
    facet_grid( ~ paste0("Participants' sex: ", part.sex)) +
    scale_fill_manual(values = c("red", "grey"), na.value = 'transparent') +
    theme_bw() +
    theme(legend.position = 'bottom') +
    ggtitle(paste0('Posterior predictive check contact patterns using ',model_id, '\nwith ', tmp, ' of data points outside the 95% posterior prediction interval',
                   '\nwith ', tmp.m, ' of Male reported data points outside the 95% posterior prediction interval'))

  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns-mean.pdf"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6 )
  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns-mean.png"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6)
}

# Plot for the comparison/discrepancy between the estimated reports and the actual reports
if (args$if_save_fig)
{
  tmp <- subset(pds, select = c('part.age', 'part.sex', 'cont.age', 'capped_cntcts_in'))
  setnames(tmp, 'capped_cntcts_in', 'M')
  set(tmp, NULL, "label", "reported contacts within the community")
  set(pds, NULL, "label", 'estimated reported contacts within the community')
  dm <- rbind(pds, tmp, fill = TRUE)

  pm <- ggplot(dm) +
    geom_tile(aes(x = part.age, y = cont.age, fill = if_out)) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "If the true contact data outside 95% P.I." ) +
    coord_equal() +
    facet_grid( ~ paste0("Participants' sex: ", part.sex)) +
    scale_fill_manual(values = c("red", "grey"), na.value = 'transparent') +
    theme_bw() +
    theme(legend.position = 'bottom') +
    ggtitle(paste0('Posterior predictive check contact patterns using ',model_id, '\nwith ', tmp, ' of data points outside the 95% posterior prediction interval',
                   '\nwith ', tmp.m, ' of Male reported data points outside the 95% posterior prediction interval'))

  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns-mean.pdf"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6 )
  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns-mean.png"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6)
}




# PPC for marginal
pds_marg <- data.table::rbindlist( pds_marg )
pds_marg <- pds_marg[,
                     list(cntcts_ppr = quantile(cntcts_ppr, p = pds.quantiles),
                          stat = pds.quantilelabels),
                     by = c('part.sex', 'part.age','c')
]
pds_marg <- dcast.data.table(pds_marg, part.sex+part.age+c~stat, value.var = 'cntcts_ppr')

tmp <- dcm[part > 0,
           list(y_marginal = sum(capped_cntcts_in, na.rm = TRUE),
                m_marginal = sum(capped_cntcts_in/part, na.rm = TRUE),
                part = unique(part)),
           by = c('part.sex','part.age','c')
]
pds_marg <- merge(pds_marg, tmp, by = c('part.sex','part.age','c'), all.x = TRUE)
set(pds_marg, NULL, 'CL.m', pds_marg[, CL/part])
set(pds_marg, NULL, 'IL.m', pds_marg[, IL/part])
set(pds_marg, NULL, 'M.m', pds_marg[, M/part])
set(pds_marg, NULL, 'IU.m', pds_marg[, IU/part])
set(pds_marg, NULL, 'CU.m', pds_marg[, CU/part])

tmp <- file.path(paste0(outfile.base,"_posterior-predictions-marginal.rds"))
cat("\nWrite marginal posterior predictions to file ", tmp)
saveRDS(pds_marg, file = tmp)

pds_marg[, if_out := (y_marginal < CL | CU < y_marginal)]
tmp <- mean(pds_marg[, as.integer(if_out)], na.rm = TRUE)
tmp <- paste0(round(tmp*100, 2), '%')
cat("\nProportion of data points outside the 95% posterior prediction interval = ", tmp)

tmp.m <- mean(pds_marg[part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
tmp.m <- paste0(round(tmp.m*100, 2), '%')
cat("\nProportion of male reported data outside the 95% posterior prediction interval = ", tmp.m)


cat("\n Plot the posterior predictive check marginal age contacts whether the contact data is outside 95% prediction interval ...")
if (args$if_save_fig)
{
  pm <- vector('list',length(unique(reported.partnerships$c)))
  names(pm) <- paste0('Community ID: ', unique(reported.partnerships$part.comm_num))
  for (i in seq_along(unique(reported.partnerships$c)))
  {
    tmp <- mean(pds_marg[c == i, as.integer(if_out)], na.rm = TRUE)
    tmp <- paste0(round(tmp*100, 2), '%')
    tmp.m <- mean(pds_marg[c == i & part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
    tmp.m <- paste0(round(tmp.m*100, 2), '%')

    pm[[i]] <- ggplot(pds_marg[c == i], aes(x = part.age)) +
      geom_step(aes(y = M, col = part.sex)) +
      pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.sex), linetype = "dotted",alpha = 0.2,
                                  colour = "transparent", show.legend = FALSE) +

      geom_step(aes(y = y_marginal), col = "black", linetype = "dashed") +
      facet_grid(paste0('Community ID:' part.comm_num) ~part.sex) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x = 'Age of contacting individual', y = 'Marginal predictied contacts') +
      theme_bw() +
      ggtitle(paste0(tmp, ' of data points outside the 95% posterior prediction interval',
              '\n', tmp.m, ' of Male reported data outside the 95% posterior prediction interval')) +
      theme(legend.position = "bottom")

  }
  p <- ggpubr::ggarrange(plotlist = pm,
                         nrow = length(unique(reported.partnerships$c)),
                         common.legend = TRUE,
                         legend = "none") %>%
    gridExtra::grid.arrange(ggpubr::get_legend(pm[[1]]), nrow = length(pm))
  # p <- annotate_figure(p,
  #                      #bottom = text_grob("Date",  hjust = 0.5, size = 25, vjust = -6),
  #                      left = text_grob("Change from the baseline", size = 25, rot = 90, vjust = 1))

  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact.pdf"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, p, width = 12, height = 6 * length(pm), limitsize = FALSE)
  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact.png"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, p, width = 12, height = 6 * length(pm), limitsize = FALSE)
}

pds_marg <- pds_marg[, list(y_marginal = mean(y_marginal),
                CL = mean(CL),
                CU = mean(CU),
                M = mean(M)),
         by = c('part.age', 'part.sex')]
pds_marg[, if_out := (y_marginal < CL | CU < y_marginal)]
tmp <- mean(pds_marg[, as.integer(if_out)], na.rm = TRUE)
tmp <- paste0(round(tmp*100, 2), '%')
cat("\nProportion of data points outside the 95% posterior prediction interval = ", tmp)

tmp.m <- mean(pds_marg[part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
tmp.m <- paste0(round(tmp.m*100, 2), '%')
cat("\nProportion of male reported data outside the 95% posterior prediction interval = ", tmp.m)

cat("\n Plot the posterior predictive check marginal age contacts whether the contact data is outside 95% prediction interval ...")
if (args$if_save_fig)
{
  pm <- ggplot(pds_marg, aes(x = part.age)) +
    geom_step(aes(y = M.m, col = part.sex)) +
    pammtools:::geom_stepribbon(aes(ymin = CL.m, ymax = CU.m, fill = part.sex), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +

    geom_step(aes(y = m_marginal), col = "black", linetype = "dashed") +
    facet_grid(part.com_num~part.sex) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age of contacting individual', y = 'Marginal predicted contact intensities') +
    theme_bw() +
    ggtitle(paste0('Posterior predictive check age marginal contact intensities using ',model_id, '\nwith ', tmp, ' of empirical contact intensities outside the 95% posterior prediction interval',
                   '\nwith ', tmp.m, ' of male reported empirical contact intensities outside the 95% posterior prediction interval')) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "transparent"))
  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact-intensities-mean.pdf"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6 * length(unique(reported.partnerships$c)), limitsize = FALSE)
  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact-intensities-mean.png"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 6 * length(unique(reported.partnerships$c)), limitsize = FALSE)
}

cat("\nDone\n")
