library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(reshape2)

#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                        help = "Print extra output [default]"),
  optparse::make_option("--seed", type = "integer", default = 18L,
                        help = "Random number seed [default %default]",
                        dest = "seed"),
  optparse::make_option("--data_set", type = "character", default = "RCCS_reported_partnership_220505.csv",
                        help = "The data set to be used [default \"%default\"]",
                        dest = "data_set"),
  optparse::make_option("--community", type = "character", default = "inland",
                        help = "The participant community[default \"%default\"]",
                        dest = "community"),
  optparse::make_option("--spline_degree", type = "integer", default = 3L,
                        help = "The degree of the B splines[default \"%default\"]",
                        dest = "spline_degree"),
  optparse::make_option("--n_knots", type = "integer", default = 30,
                        help = "The number of knots [default \"%default\"]",
                        dest = "n_knots"),
  optparse::make_option("--pkg_dir", type = "character", default = NA_character_,
                        help = "Absolute file path to package directory, used as long we don t build an R package [default]",
                        dest = "prj.dir"),
  optparse::make_option("--out_dir_base", type = "character", default = NA_character_,
                        help = "Absolute file path to base directory where all output is stored [default]",
                        dest = "out.dir")
)
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
# #
# # Define input arguments that can be changed by users----
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

args$indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
# args$community <- 'inland'
# args$n_knots <- 30
# args$spline_degree <- 3
# args$data_set = 'RCCS_reported_partnership_220505.csv'

args$infile.reported.partnerships <- file.path(args$indir.deepsequencedata, 'RCCS_R15_R18', args$data_set)
args$stan_model <- file.path('stan_models','NegBin_BSGP_rho_MF.stan')
args$outfile.prefix <- paste0('NegBin_BSGP-knots-', args$n_knots, '-comm-', args$community)

# Create directories if needed----
ifelse(!dir.exists(args$out.dir), dir.create(args$out.dir), FALSE)

cat("\n Input args are \n ")
str(args)

# Source functions
source( file.path(args$prj.dir,'functions','GP-functions.R') )

# load data----
cat("\n Loading the partnership data at individual level...")
reported.partnerships <- as.data.table(read.csv(file.path(args$infile.reported.partnerships)))

# part.age: participant age
# part.sex: participant sex

# part.T: Census eligible individuals in the age/sex/round/comm of the participant
# Z number of sexual intercourse that the participant have had in the part year
# parter_age_X: Reported age of the partner in the Xth sexual intercourse that occurred in the past year.

# part.round: participant round
# part.comm: participant community = inland, fishing

#
reported.partnerships <- as.data.table(read.csv(file.path(args$infile.reported.partnerships)))
reported.partnerships <- subset(reported.partnerships, part.comm == args$community)

# 16327

# Clean the data----
# remove those reported level not a specific number
# noted that there will be '>3' character in the Z col
reported.partnerships <- reported.partnerships[-which(grepl('>', reported.partnerships$Z)),]

# Obtain the participant size (part) and total reported partnerships (total_cntcts)----
ds <- subset(reported.partnerships, select = c('part.age', 'part.sex', 'part.round', 'Z'))
ds <- ds[,
         list(part = length(Z),
              total_cntcts = sum(as.numeric(Z))),
         by = c('part.age', 'part.sex', 'part.round')]
# check if there is missing participant groups
# and keep all categories with zero participants, setting these to 0
tmp <- as.data.table(expand.grid(part.age = 15:49, part.sex = c('M','F')))
ds <- merge(tmp, ds, by = c('part.age', 'part.sex'), all.x = TRUE)
set(ds, ds[, which(is.na(part))], "part", 0L)

# process the contacted detailed informations----
# and keep all combinations, setting these to NA
dc <- subset(reported.partnerships, select = c('pt_id', 'partner_age_1', 'partner_age_2', 'partner_age_3', 'partner_age_4'))
tmp <- as.data.table(reshape2::melt(dc, id = c('pt_id')))
setnames(tmp, 'value', 'cont.age')
set(tmp, NULL, 'variable', NULL)
tmp <- tmp[-which(is.na(tmp$cont.age)),]
# nrow(tmp)
# 15443
dc <- subset(reported.partnerships, select = c('part.age', 'part.sex', 'part.round', 'part.comm','pt_id'))
dc <- merge(dc, tmp, by = 'pt_id')
dc <- dc[, list(capped_cntcts = length(pt_id)),
         by = c('part.age', 'part.sex', 'cont.age', 'part.round', 'part.comm')]

# check if all combinations are included
# Noted that here the part.round only contains R015
tmp <- as.data.table(expand.grid(part.sex = c("M", "F"), part.age = 15:49, cont.age = 15:49))
tmp[, part.round := unique(dc$part.round)]
tmp[, part.comm := unique(dc$part.comm)]

dc <- merge(tmp, dc, by = c('part.age', 'part.sex', 'cont.age', 'part.round', 'part.comm'), all.x = TRUE)
dc[, cont.sex := ifelse(part.sex == 'M', 'F', 'M')]

# combine the total contacts
dc <- merge(dc, ds, by = c('part.age', 'part.sex', 'part.round'))
# if the part is not NA, the capped contacts will be
set(dc, dc[, which(is.na(capped_cntcts) & part > 0L)], "capped_cntcts", 0L)

# Get the population size(pop) ----
pop <- unique(subset(reported.partnerships, select = c('part.age', 'part.sex', 'part.round', 'part.T')))
setnames(pop, c('part.age', 'part.sex', 'part.T'), c('cont.age', 'cont.sex', 'pop'))
reported.partnerships <- merge(dc, pop, by = c('cont.age', 'cont.sex', 'part.round'), all.x = TRUE)

# Generate Stan data----
cat("\n Generating stan data ...")
setkey(reported.partnerships, cont.sex, part.sex, cont.age, part.age)
reported.partnerships[, index := 1:nrow(reported.partnerships)]

# apply rho to account for the missing reports
tmp <- reported.partnerships[,
                             list(rho = sum(capped_cntcts, na.rm = TRUE) / unique(total_cntcts)),
                             by = c("part.sex", "part.age", "part.round")
]
reported.partnerships <- merge(reported.partnerships, tmp, by = c("part.sex","part.age","part.round"), all.x = TRUE)

# We only interesting in heterosexuals partnerships
# make obs_2_MF_matrix index
tmp <- subset(reported.partnerships,
              cont.sex == 'F' & part.sex == 'M',
              select = c(cont.sex, part.sex, cont.age, part.age)
)
setkey(tmp, cont.sex, part.sex, cont.age, part.age)
tmp[, mf_mat_idx := seq_len(nrow(tmp))]
tmp2 <- copy(tmp)
setnames(tmp2, c('cont.age','part.age','cont.sex','part.sex'), c('part.age','cont.age','part.sex','cont.sex'))
tmp <- rbind(tmp, tmp2)
reported.partnerships <- merge(reported.partnerships, tmp, by = c('cont.sex','part.sex','cont.age','part.age'), all.x = TRUE)

stan_data <- list()
stan_data$A <- length(unique(reported.partnerships$part.age))
stan_data$age1 <- unique(reported.partnerships$part.age)
stan_data$age2 <- unique(reported.partnerships$cont.age)

# contact data
stan_data$Nmf <- nrow(reported.partnerships[part.sex == 'M' & cont.sex == 'F' & part > 0L,])
stan_data$ymf <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & part > 0L, capped_cntcts]
stan_data$ymf_rowmajor_matrix_index <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & part > 0L, mf_mat_idx]

# offsets
stan_data$log_pop_mf <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & part > 0L, log(pop)]
stan_data$log_participants_mf <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & part > 0L, log(part)]
stan_data$log_rho_mf <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & part > 0L, log(rho)]

# B-splines projected GP parameters
tmp <- stan_data$age1[seq(1, stan_data$A, length.out = args$n_knots)]
stan_data$M_age1 <- length(tmp) + args$spline_degree - 1
stan_data$basis_age1 <- bsplines(stan_data$age1, tmp, args$spline_degree)
stan_data$idx_basis_age1 <- seq_len(stan_data$M_age1)
tmp <- stan_data$age2[seq(1, stan_data$A, length.out = args$n_knots)]
stan_data$M_age2 <- length(tmp) + args$spline_degree - 1
stan_data$basis_age2 <- bsplines(stan_data$age2, tmp, args$spline_degree)
stan_data$idx_basis_age2 <- seq_len(stan_data$M_age2)

tmp <- paste0(file.path(args$out.dir, args$outfile.prefix), "_image.RData")
cat("\n Save image prior to running Stan to file ", tmp , "\n")
save.image( file = tmp )

# Run Stan model---
cat("\n Run Stan model ...")
tmp <- file.path(args$prj.dir, args$stan_model)
model <- cmdstanr::cmdstan_model(stan_file = tmp)

if (0) # for debugging/testing
{
  model_fit <- model$sample(
    data = stan_data,
    seed = args$seed,
    chains = 1,
    parallel_chains = 1,
    refresh = 1e2,
    iter_warmup = 1e3,
    iter_sampling = 1e4,
    max_treedepth = 13,
    save_warmup = TRUE
    )
}

model_fit <- model$sample(
  data = stan_data,
  seed = args$seed,
  chains = 2,
  parallel_chains = 2,
  refresh = 1e2,
  iter_warmup = 1e3,
  iter_sampling = 1e4,
  max_treedepth = 13,
  save_warmup = TRUE
)

tmp <- paste0(file.path(args$out.dir, args$outfile.prefix), "_stan_fit.rds")
cat("\n Save fitted data to file ", tmp , "\n")
model_fit$save_object(file = tmp)

cat("\nDone\n.")
