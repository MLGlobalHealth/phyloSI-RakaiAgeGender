library(rstan)
library(data.table)	

jobname <- 'diagonalprior'
stan_model <- 'gp_220317'
DEBUG <- F

indir <- "/rds/general/user/mm3218/home/git/phyloflows"
outdir <- paste0("/rds/general/user/mm3218/home/projects/2021/phyloflows/", stan_model, '-', jobname)

if(0)
{
  indir <- '~/git/phyloflows'
  outdir <- file.path('~/Box\ Sync/2021/phyloflows/', paste0(stan_model, '-', jobname))
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-stan_model')
  stopifnot(args_line[[7]]=='-jobname')
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  stan_model <- args_line[[6]]
  jobname <- args_line[[8]]
}

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))

# load stan data
path.to.stan.data <- paste0(outfile, "-stanin_",jobname,".RData")
load(path.to.stan.data)

# make stan model
path.to.stan.model <- file.path(indir, 'stan_models', paste0(stan_model, '.stan'))
model = rstan::stan_model(path.to.stan.model)

# run stan model
# stan_init <- list()
# stan_init$rho_gp1 = rep(1.25, stan_data$N_group)
# stan_init$rho_gp2 = rep(1.25, stan_data$N_group)
# stan_init$alpha = 0
# stan_init$nu = rep(0, stan_data$N_group)

if(DEBUG){
  fit <- sampling(model, data = stan_data, iter = 10, warmup = 5, chains=1, thin=1)
}else{
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- sampling(model, data = stan_data,
                  iter = 3000, warmup = 500, chains=4, thin=1, seed = 5,
                  verbose = FALSE, control = list(adapt_delta = 0.999,max_treedepth=15))
}

# sum = summary(fit)
# sum$summary[which(sum$summary[,9] < 100),]

file = paste0(outfile, "-stanout_", jobname, ".rds")
cat("Save file ", file)
saveRDS(fit,file = file)
     
