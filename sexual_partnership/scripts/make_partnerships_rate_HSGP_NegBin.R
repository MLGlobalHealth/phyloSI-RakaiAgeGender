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
  optparse::make_option("--hsgp_boundary_inflation", type = "double", default = 1.2,
                        help = "The boundary inflation of the HSGP prior in any dimension [default \"%default\"]",
                        dest = "hsgp_boundary_inflation"),
  optparse::make_option("--hsgp_m", type = "integer", default = 20,
                        help = "The number of the HSGP basis functions in any dimension [default \"%default\"]",
                        dest = "hsgp_m"),
  optparse::make_option("--pkg_dir", type = "character", default = NA_character_,
                        help = "Absolute file path to package directory, used as long we don t build an R package [default]",
                        dest = "prj.dir"),
  optparse::make_option("--out_dir_base", type = "character", default = NA_character_,
                        help = "Absolute file path to base directory where all output is stored [default]",
                        dest = "out.dir")
)
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
#
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
# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
  args$out.dir <- here::here()
}

args$indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
# args$community <- 'fishing'
# args$hsgp_boundary_inflation <- 1.2
# args$hsgp_m <- 20
args$infile.reported.partnerships <- file.path(args$indir.deepsequencedata, 'RCCS_R15_R18', args$data_set)
args$stan_model <- file.path('stan_models','NegBin_HSGP.stan')
args$outfile.prefix <- paste0("NegBin_HSGP_c-",args$hsgp_boundary_inflation*100,"_m-",args$hsgp_m,'-comm-', args$community)

# Create directories if needed----
ifelse(!dir.exists(args$out.dir), dir.create(args$out.dir), FALSE)

cat("\n Input args are \n ")
str(args)

# Source functions
source( file.path(args$prj.dir,'functions','GP-functions.R') )

# load data----
reported.partnerships <- as.data.table(read.csv(file.path(args$infile.reported.partnerships)))

# cont.age: contacted age
# part.age: participant age
# part.sex: participant sex
# cont.sex: contacted sex
# y: number of partnership reported by participant
# N: Number of participants
# T: population count of contacted
# U = N x T

# round: participant round
# comm: participant community = inland, fishing

# Generate Stan data----
cat("\n Generating stan data ...")
reported.partnerships <- reported.partnerships[part.comm == args$community]
setkey(reported.partnerships, cont.sex, part.sex, cont.age, part.age)
reported.partnerships[, index := 1:nrow(reported.partnerships)]

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
stan_data$Nmf <- nrow(reported.partnerships[part.sex == 'M' & cont.sex == 'F' & N > 0L,])
stan_data$Nfm <- nrow(reported.partnerships[part.sex == 'F' & cont.sex == 'M' & N > 0L,])
stan_data$ymf <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & N > 0L, y]
stan_data$yfm <- reported.partnerships[part.sex == 'F' & cont.sex == 'M' & N > 0L, y]
stan_data$ymf_rowmajor_matrix_index <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & N > 0L, mf_mat_idx]
stan_data$yfm_rowmajor_matrix_index <- reported.partnerships[part.sex == 'F' & cont.sex == 'M' & N > 0L, mf_mat_idx]

# offsets
stan_data$log_pop_mf <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & N > 0L, log(T)]
stan_data$log_pop_fm <- reported.partnerships[part.sex == 'F' & cont.sex == 'M' & N > 0L, log(T)]
stan_data$log_participants_mf <- reported.partnerships[part.sex == 'M' & cont.sex == 'F' & N > 0L, log(N)]
stan_data$log_participants_fm <- reported.partnerships[part.sex == 'F' & cont.sex == 'M' & N > 0L, log(N)]

# HSGP input variables
stan_data$c_age1 <- args$hsgp_boundary_inflation
stan_data$c_age2 <- args$hsgp_boundary_inflation
stan_data$M_age1 <- args$hsgp_m
stan_data$M_age2 <- args$hsgp_m

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
  chains = 4,
  parallel_chains = 4,
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
