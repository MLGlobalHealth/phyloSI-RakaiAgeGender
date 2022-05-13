# Postprocessing for simulation age gender ----
# Preamble ----
# This script aims to analyse the outputs from the stan model (assess mixing and convergence) and do postprocessing for the outputs

# Load the required packages
require(data.table)
require(ggplot2)
require(cmdstanr)
require(bayesplot)
require(posterior)
#
#
# Define input arguments that can be changed by users
#
# option_list <- list(
#   optparse::make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
#                         help = "Print extra output [default]"),
#   optparse::make_option("--seed", type = "integer", default = 18L,
#                         help = "Random number seed [default %default]",
#                         dest = "seed"),
#   optparse::make_option("--pkg_dir", type = "character", default = NA_character_,
#                         help = "Absolute file path to package directory, used as long we don t build an R package [default]",
#                         dest = "prj.dir"),
#   optparse::make_option("--out_dir_base", type = "character", default = NA_character_,
#                         help = "Absolute file path to base directory where all output is stored [default]",
#                         dest = "out.dir")
# )
# args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))


# Define input arguments that can be changed by users----
args <- list()
args$seed <- 18L

# define prj.dir and out.dir
tmp <- Sys.info()
if (tmp["user"] == "yc2819" & grepl("hpc.ic.ac.uk",tmp["nodename"])) # outdir yu
{

    args$out.dir <- "/rds/general/user/yc2819/home/github/phyloflows/sexual_partnership/results"

    args$prj.dir <- "/rds/general/user/yc2819/home/github/phyloflows/sexual_partnership"

}

# Source functions
source( file.path(args$prj.dir,'functions','plot_functions.R') )

# read simulated contact intensities, input args, simulated data to model from workspace
args$in.dir <- args$out.dir
infiles <- list.files(args$in.dir, pattern = 'comm', full.names = TRUE)
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
  su.target.vars <- c('log_random_effect_baseline','gp_rho_age1','gp_rho_age2','gp_alpha','z', 'f_mf')
  model_id <- paste0('BSGP_knots-', args$n_knots, '-comm-', args$community)
}
if ( grepl("HSGP", args$stan_model))
{
  su.target.vars <- c('log_random_effect_baseline','gp_rho_age1','gp_rho_age2','gp_sigma','z', 'f_mf')
  model_id <- paste0('HSGP_c-',args$hsgp_boundary_inflation*100,"_m-",args$hsgp_m, '-comm-', args$community)
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

# Computing time----
compute.time <- model_fit$metadata()$time
tmp <- file.path(paste0(outfile.base,"_computing-time.csv"))
cat("\nWrite computing time to file ", tmp)
write.csv(compute.time, file = tmp, row.names = TRUE)

# Extract Monte Carlo samples for postprocessing----
cat("\nExtract Monte Carlo samples for postprocessing ...\n")
su.target.vars <- c('log_random_effect_baseline','f_mf')
pd <- model_fit$draws(variables = su.target.vars, inc_warmup = FALSE)

select.chains <- seq_along(dimnames(pd)[['chain']])
iters <- seq(from = 1, to = model_fit$metadata()[['iter_sampling']], length.out = ceiling(1e4 / length(select.chains)))
# iters <- seq_len(model_fit$metadata()[['iter_sampling']])
po <- list()
tmp <- pd[,,which(grepl('f_mf',dimnames(pd)[[3]]))]
po$f_mf <- unname(apply(tmp[iters,select.chains,], 3, rbind))
tmp <- pd[,,which(grepl('log_random_effect_baseline',dimnames(pd)[[3]]))]
po$log_random_effect_baseline <- unname(apply(tmp[iters,select.chains,], 3, rbind))

pd <- NULL
gc()

pds.quantiles <- c(.025,.25,.5,.75,.975)
pds.quantilelabels <- c('CL','IL','M','IU','CU')
pds.MF <- 1 # from Stan model file

# Extract age and gender specific contact intensities----
# note log(T) are just constants and exp is monotonic,
#   so in this special case we can summarise the quantiles of "log_random_effect_baseline + f" already now
#   to avoid large mem costs
cat("\nExtract age and gender specific contact intensities ...\n")
pd <- as.data.table(reshape2::melt(po$f_mf))
setnames(pd, 1:3, c('iterations','mf_mat_idx','value'))
tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                    log_random_effect_baseline = po$log_random_effect_baseline[,pds.MF]
                    )
pd <- merge(pd, tmp2, by = 'iterations')
tmp2 <- subset(reported.partnerships, part.sex == 'M' & cont.sex == 'F', select = c(mf_mat_idx, T))
pd <- merge(pd, tmp2, by = 'mf_mat_idx')
set(pd, NULL, 'value', pd[, exp(log_random_effect_baseline + value + log(T))])
tmp <- pd[,
          list(sd = sd(value),
               value = quantile(value, p = pds.quantiles),
               stat = pds.quantilelabels),
          by = c('mf_mat_idx')
          ]
tmp <- dcast.data.table(tmp, mf_mat_idx + sd~stat, value.var = 'value')
tmp2 <- subset(reported.partnerships, part.sex == 'M' & cont.sex == 'F', select = c(part.sex, cont.sex, part.age, cont.age, mf_mat_idx, y, N, T))
pds <- merge(tmp, tmp2, by = 'mf_mat_idx')

# and only now we repeat for fm
pd <- as.data.table(reshape2::melt(po$f_mf))
setnames(pd, 1:3, c('iterations','mf_mat_idx','value'))
tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                    log_random_effect_baseline = po$log_random_effect_baseline[,pds.MF]
                    )
pd <- merge(pd, tmp2, by = 'iterations')
tmp2 <- subset(reported.partnerships, part.sex == 'F' & cont.sex == 'M', select = c(mf_mat_idx, T))
pd <- merge(pd, tmp2, by = 'mf_mat_idx')

set(pd, NULL, 'value', pd[, exp(log_random_effect_baseline + value + log(T))])
tmp <- pd[,
          list(sd = sd(value),
               value = quantile(value, p = pds.quantiles),
               stat = pds.quantilelabels),
          by = c('mf_mat_idx')
          ]
tmp <- dcast.data.table(tmp, mf_mat_idx + sd~stat, value.var = 'value')
tmp2 <- subset(reported.partnerships, part.sex == 'F' & cont.sex == 'M', select = c(part.sex, cont.sex, part.age, cont.age, mf_mat_idx, y, N, T))
tmp <- merge(tmp, tmp2, by = 'mf_mat_idx')
pds <- rbind(pds, tmp, use.names = TRUE, fill = TRUE)

pds[, cntct_intensity_empirical := y/N]
tmp <- file.path(paste0(outfile.base, "_contact-intensities.rds"))
cat("\nWrite contact intensities to file ", tmp)
saveRDS(pds, file = tmp)

# Visualize age and gender specific contact intensities----
cat("\nVisualize age and gender specific contact intensities ...\n")
# Plot age posterior distribution----
if (args$if_save_fig)
{
  cat("\nPlot age posterior distribution ...")
  # we just visualize the contact intensities for male-female
  pm <- plot_age_contact_intensity_MF(pds)
  pm <- pm + ggtitle(paste0('Age dist using ', model_id, " - MF" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-MF.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-MF.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}

if (args$if_save_fig)
{
  # Participants and the contacted individuals aged 0 - 17 are missing in the BICs survey

  cat("\nPlot age posterior distribution with empirical contact intensities ...")
  # we just visualize the contact intensities for male-female
  pm <- plot_age_contact_intensity_MF_empirical(pds)
  pm <- pm + ggtitle(paste0('Age dist with emprical contact intensities using ', model_id, " - MF" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-MF-with-empirical.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-MF-with-empirical.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}



if (args$if_save_fig)
{
  cat("\nPlot age posterior distribution ...")
  # we just visualize the contact intensities for female-male
  pm <- plot_age_contact_intensity_FM(pds)
  pm <- pm + ggtitle(paste0('Age dist using ', model_id, " - FM" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-FM.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-FM.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}

if (args$if_save_fig)
{
  # Participants and the contacted individuals aged 0 - 17 are missing in the BICs survey

  cat("\nPlot age posterior distribution with empirical contact intensities ...")
  # we just visualize the contact intensities for female-male
  pm <- plot_age_contact_intensity_FM_empirical(pds)
  pm <- pm + ggtitle(paste0('Age dist with emprical contact intensities using ', model_id, " - FM" ))

  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-FM-with-empirical.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_age-posterior-dist-contact-intensities-FM-with-empirical.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)

}
# Plot the selected age contact intensites for comparison ----
# make this figure just for one dataset
if ( args$if_save_fig)
{
  tmp <- pds[part.age %in% c(10, 20, 30, 40)]
  tmp$part.age <- as.character(tmp$part.age)
  tmp[, part.age := paste0("participants: ", part.age)]
  tmp[, id := model_id]
  tmp[, age_gender := paste0(part.age, '-', part.sex)]

  pm <- ggplot(tmp, aes(x = cont.age, col = age_gender)) +
    geom_step(aes(y = M)) +
    geom_step(aes(y = cntct_intensity_empirical), col = "black", linetype = "dashed") +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = age_gender), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    facet_grid(id ~age_gender, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") +
    theme(legend.position = "none", legend.title = element_blank()) +
    xlab("Age of contacted individuals") +
    ylab("Posterior mean contact intensities")
  pm <- pm + theme(axis.text = element_text(size = 14),
                   axis.title = element_text(size = 14),
                   strip.text = element_text(size = 16)
                   # axis.text.x = element_text(angle = 45)
  )

  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-comparison.pdf"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6)
  tmp <- file.path(paste0(outfile.base,"_selected-age-posterior-contact-intensities-comparison.png"))
  cat("\nSave age posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 22, height = 6)
}

# Plot contact intensities patterns----
if (args$if_save_fig)
{
  cat("\nPlot contact intensities patterns ...")
  pds[, M := ifelse(M == 0, NaN, M)]
  pm <- plot_contact_intensity_matrix(pds) + labs(x = "Age of participant", y = "Age of contact", fill = "contact intensity" )
  pm <- pm + ggtitle(paste0('Posterior mean contacts intensities using ', model_id ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns.pdf"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns.png"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
}

if (args$if_save_fig)
{
  cat("\nPlot contact intensities patterns ...")
  pm <- plot_contact_intensity_sd_matrix(pds) + labs(x = "Age of participant", y = "Age of contact", fill = "standard deviation of the contact intensity" )
  pm <- pm + ggtitle(paste0('Standard deviation of the posterior contacts intensities using ', model_id )) +
    theme(legend.text = element_text(angle = 45, vjust = 0.5, hjust = 1))

  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-sd-patterns.pdf"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-sd-patterns.png"))
  cat("\nSave posterior contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 12)
}

tmp.emp <- subset(pds, select = c('part.sex', 'cont.sex', 'part.age', 'cont.age', 'cntct_intensity_empirical'))
setnames(tmp.emp, 'cntct_intensity_empirical', 'M')
set(tmp.emp, NULL, "label", "empirical value")
set(pds, NULL, "label", model_id)
dm <- rbind(pds, tmp.emp, fill = TRUE)
if (args$if_save_fig)
{
  pm <- plot_contact_intensity_matrix_emprical(dm)
  pm <- pm +  ggtitle(paste0('Contact intensities patterns comparison using ', model_id  ))
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns-with-truth.pdf"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
  tmp <- file.path(paste0(outfile.base,"_posterior-contact-intensities-patterns-with-truth.png"))
  cat("\nSave posterior contact intensities with ground truth to file ", tmp)
  ggsave(file = tmp, pm, width = 12, height = 8)
}


# Extract age and gender specific marginal contact intensities----
#   here we will process the Monte Carlo output in batches to avoid handling tables with >100m rows
cat("\nExtract age and gender specific marginal contact intensities ...\n")
dcm <- subset(reported.partnerships, part.sex != cont.sex)
setkey(dcm, part.sex, part.age)
tmp <- unique(subset(dcm, select = c(part.sex, part.age)))
tmp[, m_index := seq_len(nrow(tmp))]
dcm <- merge(dcm, tmp, by = c('part.sex', 'part.age'))
pds.midx <- tmp$m_index

pds <- vector('list',length(pds.midx))
for (i in pds.midx)
{
  cat("\nProcessing ",i," of ",length(pds.midx))
  gender_i <- subset(dcm, m_index == i )[, part.sex][1]
  mf_mat_idx <- subset(dcm, m_index == i & cont.sex != gender_i)[, mf_mat_idx]

  if (gender_i == 'F')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline[,pds.MF]
    )
    pd <- merge(pd, tmp2, by = 'iterations')
    pd[, log_contact_rate := log_random_effect_baseline + f_mf]
    tmp <- subset(dcm, m_index == i, select = c(m_index, cont.sex, part.sex, cont.age, part.age, mf_mat_idx, T))
    pd <- merge(pd, tmp, by = 'mf_mat_idx')
    pd[, log_contact_intensity := log_random_effect_baseline + f_mf + log(T)]
    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','f_mf','log_contact_rate','T'), NULL)
    pds[[i]] <- copy(pd)

    pds[[i]] <- pds[[i]][, .(log_contact_intensity = log(sum(exp(log_contact_intensity)))), by = c('iterations','part.sex','part.age')]
  }
  if (gender_i == 'M')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline[,pds.MF]
    )
    pd <- merge(pd, tmp2, by = 'iterations')
    pd[, log_contact_rate := log_random_effect_baseline + f_mf]
    tmp <- subset(dcm, m_index == i, select = c(m_index, cont.sex, part.sex, cont.age, part.age, mf_mat_idx, T))
    pd <- merge(pd, tmp, by = 'mf_mat_idx')
    pd[, log_contact_intensity := log_random_effect_baseline + f_mf + log(T)]
    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','f_mf','log_contact_rate','T'), NULL)
    pds[[i]] <- copy(pd)

    # marginalise with log sum exp
    pds[[i]] <- pds[[i]][, .(log_contact_intensity = log(sum(exp(log_contact_intensity)))), by = c('iterations','part.sex','part.age')]
  }
}
pds <- data.table::rbindlist( pds )

# now summarise
pds <- pds[,
           list(log_contact_intensity = quantile(log_contact_intensity, p = pds.quantiles),
                stat = pds.quantilelabels),
           by = c('part.sex','part.age')
]
pds[, contact_intensity := exp(log_contact_intensity)]
pds <- dcast.data.table(pds, part.sex+part.age~stat, value.var = 'contact_intensity')

tmp <- dcm[N > 0,
           list(cntct_intensity_empirical = sum(y / N)),
           by = c('part.sex','part.age')
]
pds <- merge(pds, tmp, by = c('part.sex','part.age'), all.x = TRUE)
tmp <- file.path(paste0(outfile.base,"_marginal-contact-intensities.rds"))
cat("\nWrite marginal contact intensities to file ", tmp)
saveRDS(pds, file = tmp)

# Visualize age and gender specific contact intensities----
cat("\nVisualize marginal contact intensities ...\n")
# Plot the marginal contact intensities ----
if (args$if_save_fig)
{
  pm <- ggplot(pds, aes(x = part.age)) +
    geom_step(aes(y = M, col = part.sex)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.sex), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +

    # TODO: not consistent with the empirical value
    geom_step(aes(y = cntct_intensity_empirical), col = "black", linetype = "dashed") +
    facet_grid(~part.sex) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age of contacting individual', y = 'Contact intensity') +
    theme_bw() +
    ggtitle(paste0('Posterior age marginal distribution using ', model_id  )) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "transparent"))
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities.pdf"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
  tmp <- file.path(paste0(outfile.base,"_age-marginal-contact-intensities.png"))
  cat("\nSave mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
}

# Posterior predictive check on age and gender specific contact intensities----
#   here we will process the Monte Carlo output in batches to avoid handling tables with >100m rows
cat("\nPosterior predictive check on age and gender specific contact intensities ...\n")
dcm <- subset(reported.partnerships, N > 0 & part.sex != cont.sex)
setkey(dcm, part.sex, part.age)
tmp <- unique(subset(dcm, select = c(part.sex, part.age)))
tmp[, m_index := seq_len(nrow(tmp))]
dcm <- merge(dcm, tmp, by = c('part.sex', 'part.age'))
pds.midx <- tmp$m_index

pds <- vector('list',length(pds.midx))
pds_marg <- vector('list',length(pds.midx))

for (i in pds.midx)
{
  #i <- 45
  cat("\nProcessing ",i," of ",length(pds.midx))
  gender_i <- subset(dcm, m_index == i )[, part.sex][1]
  mf_mat_idx <- subset(dcm, m_index == i & cont.sex != gender_i)[, mf_mat_idx]

  if (gender_i == 'F')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline[,pds.MF]
    )
    pd <- merge(pd, tmp2, by = 'iterations')
    tmp <- subset(dcm, m_index == i, select = c(m_index, part.sex, cont.sex, part.age, cont.age, mf_mat_idx, T, N))
    pd <- merge(pd, tmp, by = 'mf_mat_idx')
    pd[, log_exp_contacts := log_random_effect_baseline + f_mf + log(N) + log(T)]
    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','f_mf','T','N'), NULL)
    pds[[i]] <- copy(pd)
    pds_marg[[i]] <- pds[[i]][, .(log_exp_contacts = log(sum(exp(log_exp_contacts)))), by = c('iterations','part.sex','part.age')]

  }
  if (gender_i == 'M')
  {
    pd <- as.data.table(reshape2::melt(po$f_mf[,mf_mat_idx]))
    setnames(pd, 1:3, c('iterations','col_idx','f_mf'))
    tmp2 <- data.table(col_idx = seq_along(mf_mat_idx), mf_mat_idx = mf_mat_idx)
    pd <- merge(pd, tmp2, by = 'col_idx')
    tmp2 <- data.table( iterations = seq_len(nrow(po$log_random_effect_baseline)),
                        log_random_effect_baseline = po$log_random_effect_baseline[,pds.MF])
    pd <- merge(pd, tmp2, by = 'iterations')
    tmp <- subset(dcm, m_index == i, select = c(m_index, part.sex, cont.sex, part.age, cont.age, mf_mat_idx, T, N))
    pd <- merge(pd, tmp, by = 'mf_mat_idx')
    pd[, log_exp_contacts := log_random_effect_baseline + f_mf + log(N) + log(T)]
    set(pd, NULL, c('m_index','col_idx','mf_mat_idx','log_random_effect_baseline','f_mf','T','N'), NULL)
    pds[[i]] <- copy(pd)
    pds_marg[[i]] <- pds[[i]][, .(log_exp_contacts = log(sum(exp(log_exp_contacts)))), by = c('iterations','part.sex','part.age')]


  }
  pds[[i]][, cntcts_ppr := rpois( nrow(pds[[i]]),
                                  lambda = pds[[i]][, exp(log_exp_contacts) ]
  )
  ]
  pds_marg[[i]][, cntcts_ppr := rpois( nrow(pds_marg[[i]]),
                                  lambda = pds_marg[[i]][, exp(log_exp_contacts) ]
  )
  ]
  # summarise
  pds[[i]] <- pds[[i]][,
                       list(cntcts_ppr = quantile(cntcts_ppr, p = pds.quantiles),
                            stat = pds.quantilelabels),
                       by = c('part.sex', 'cont.sex', 'part.age', 'cont.age')
  ]

}
pds <- data.table::rbindlist( pds )
pds <- dcast.data.table(pds, part.sex+cont.sex+part.age+cont.age~stat, value.var = 'cntcts_ppr')
tmp <- subset(dcm, select = c('part.sex','cont.sex','part.age','cont.age','y'))
pds <- merge(pds, tmp, by = c('part.sex','cont.sex','part.age','cont.age'))

tmp <- file.path(paste0(outfile.base,"_posterior-predictions.rds"))
cat("\nWrite posterior predictions to file ", tmp)
saveRDS(pds, file = tmp)

pds[, if_out := (y < CL | CU < y)]
tmp <- mean(pds[, as.integer(if_out)], na.rm = TRUE)
tmp <- paste0(round(tmp*100, 2), '%')
cat("\nProportion of data points outside the 95% posterior prediction interval = ", tmp)

tmp.m <- mean(pds[part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
tmp.m <- paste0(round(tmp.m*100, 2), '%')
cat("\nProportion of male reported data points outside the 95% posterior prediction interval = ", tmp.m)

# Plot the posterior predictive check----
# I think just a simple plot showing the gender age maps (x-axis age, y-axis age, panel for gender combination)
# and as value show "red" if data outside 95% prediction interval, "green" if inside, and "transparent" if no data
cat("\n Plot the posterior predictive check contacts patterns whether the contact data is outside 95% prediction interval ...")
if (args$if_save_fig)
{
  set(pds, NULL,"gender_label", pds[, paste0(part.sex,' (Participant)')])
  set(pds, NULL,"contact_gender_label", pds[, paste0(cont.sex,' (Contacted)')])

  pm <- ggplot(pds) +
    geom_tile(aes(x = part.age, y = cont.age, fill = if_out)) +
    labs(x = "Age of contacting individuals", y = "Age of contacted individuals", fill = "If the true contact data outside 95% P.I." ) +
    coord_equal() +
    facet_grid(contact_gender_label ~ gender_label) +
    scale_fill_manual(values = c("red", "green"), na.value = 'transparent') +
    theme_bw() +
    theme(legend.position = 'bottom') +
    ggtitle(paste0('Posterior predictive check contact patterns using ',model_id, '\nwith ', tmp, ' of data points outside the 95% posterior prediction interval',
                   '\nwith ', tmp.m, ' of Male reported data points outside the 95% posterior prediction interval'))


  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns.pdf"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
  tmp <- file.path(paste0(outfile.base,"_posterior-predictive-check-contact-patterns.png"))
  cat("\nSave posterior predictive check contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
}


# PPC for marginal
pds_marg <- data.table::rbindlist( pds_marg )
pds_marg <- pds_marg[,
                     list(cntcts_ppr = quantile(cntcts_ppr, p = pds.quantiles),
                          stat = pds.quantilelabels),
                     by = c('part.sex', 'part.age')
]
pds_marg <- dcast.data.table(pds_marg, part.sex+part.age~stat, value.var = 'cntcts_ppr')

tmp <- dcm[N > 0,
           list(y_marginal = sum(y),
                m_marginal = sum(y/N)),
           by = c('part.sex','part.age', 'N')
]
pds_marg <- merge(pds_marg, tmp, by = c('part.sex','part.age'), all.x = TRUE)
set(pds_marg, NULL, 'CL.m', pds_marg[, CL/N])
set(pds_marg, NULL, 'IL.m', pds_marg[, IL/N])
set(pds_marg, NULL, 'M.m', pds_marg[, M/N])
set(pds_marg, NULL, 'IU.m', pds_marg[, IU/N])
set(pds_marg, NULL, 'CU.m', pds_marg[, CU/N])

tmp <- file.path(paste0(outfile.base,"_posterior-predictions-marginal.rds"))
cat("\nWrite marginal posterior predictions to file ", tmp)
saveRDS(pds, file = tmp)

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
    geom_step(aes(y = M, col = part.sex)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.sex), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +

    geom_step(aes(y = y_marginal), col = "black", linetype = "dashed") +
    facet_grid(~part.sex) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age of contacting individual', y = 'Marginal predictied contacts') +
    theme_bw() +
    ggtitle(paste0('Posterior predictive check age marginal contact using ',model_id, '\nwith ', tmp, ' of data points outside the 95% posterior prediction interval',
                   '\nwith ', tmp.m, ' of Male reported data outside the 95% posterior prediction interval')) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "transparent"))
  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact.pdf"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact.png"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
}

pds_marg[, if_out := (m_marginal < CL.m | CU.m < m_marginal)]
tmp <- mean(pds_marg[, as.integer(if_out)], na.rm = TRUE)
tmp <- paste0(round(tmp*100, 2), '%')
cat("\nProportion of empirical contact intensities outside the 95% posterior prediction interval = ", tmp)
tmp.m <- mean(pds_marg[part.sex == 'M', as.integer(if_out)], na.rm = TRUE)
tmp.m <- paste0(round(tmp.m*100, 2), '%')
cat("\nProportion of Male reported empirical contact intensities outside the 95% posterior prediction interval = ", tmp.m)


cat("\n Plot the posterior predictive check marginal age contacts whether the contact data is outside 95% prediction interval ...")
if (args$if_save_fig)
{
  pm <- ggplot(pds_marg, aes(x = part.age)) +
    geom_step(aes(y = M.m, col = part.sex)) +
    pammtools:::geom_stepribbon(aes(ymin = CL.m, ymax = CU.m, fill = part.sex), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +

    geom_step(aes(y = m_marginal), col = "black", linetype = "dashed") +
    facet_grid(~part.sex) +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Age of contacting individual', y = 'Marginal predicted contact intensities') +
    theme_bw() +
    ggtitle(paste0('Posterior predictive check age marginal contact intensities using ',model_id, '\nwith ', tmp, ' of empirical contact intensities outside the 95% posterior prediction interval',
                   '\nwith ', tmp.m, ' of male reported empirical contact intensities outside the 95% posterior prediction interval')) +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "transparent"))
  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact-intensities.pdf"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
  tmp <- file.path(paste0(outfile.base,"_ppc-age-marginal-contact-intensities.png"))
  cat("\nSave posterior predictive mariginal contact intensities to file ", tmp)
  ggsave(file = tmp, pm, width = 6, height = 5)
}

cat("\nDone\n")
