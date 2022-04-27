library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(reshape2)

# Define input arguments that can be changed by users----
args <- list()
args$seed <- 18L

# define prj.dir and out.dir
tmp <- Sys.info()
if (tmp["user"] == "yc2819" & grepl("hpc.ic.ac.uk",tmp["nodename"])) # outdir yu
{
  args$out.dir <- "/rds/general/user/yc2819/home/github/phyloflows"
  args$prj.dir <- "/rds/general/user/yc2819/home/github/phyloflows"
}else
{
  args$prj.dir <- here::here()
  args$out.dir <- here::here()
}
args$indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'

args$infile.reported.partnerships <- file.path(args$indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_reported_partnership_220426.csv')
args$stan_model <- file.path('sexual_partnership','stan_models','Bspline_gp.stan')
args$outfile.prefix <- 'BSGP-inland'

# Source functions
source( file.path(args$prj.dir,'sexual_partnership','functions','GP-functions.R') )

# load data----
reported.partnerships <- as.data.table(read.csv(file.path(args$prj.dir, args$infile.reported.partnerships)))

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
reported.partnerships <- reported.partnerships[comm == 'inland']
setkey(reported.partnerships, cont.sex, part.sex, cont.age, part.age)
reported.partnerships[, index := 1:nrow(reported.partnerships)]

# We only interesting in heterosexuals partnerships
# make obs_2_MF_matrix index
tmp <- subset(reported.partnerships,
              cont.sex == 'Female' & part.sex == 'Male',
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
stan_data$Nmf <- nrow(reported.partnerships[part.sex == 'Male' & cont.sex == 'Female' & N > 0L,])
stan_data$Nfm <- nrow(reported.partnerships[part.sex == 'Female' & cont.sex == 'Male' & N > 0L,])
stan_data$ymf <- reported.partnerships[part.sex == 'Male' & cont.sex == 'Female' & N > 0L, cntcts_data]
stan_data$yfm <- reported.partnerships[part.sex == 'Female' & cont.sex == 'Male' & N > 0L, cntcts_data]
stan_data$ymf_rowmajor_matrix_index <- reported.partnerships[part.sex == 'Male' & cont.sex == 'Female' & N > 0L, mf_mat_idx]
stan_data$yfm_rowmajor_matrix_index <- reported.partnerships[part.sex == 'Female' & cont.sex == 'Male' & N > 0L, mf_mat_idx]

# offsets
stan_data$log_offset_mf <- reported.partnerships[part.sex == 'Male' & cont.sex == 'Female' & N > 0L, log(U)]
stan_data$log_offset_fm <- reported.partnerships[part.sex == 'Female' & cont.sex == 'Male' & N > 0L, log(U)]

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
    iter_warmup = 5e2,
    iter_sampling = 1e3,
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
