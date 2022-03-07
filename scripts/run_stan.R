library(rstan)
library(data.table)	

jobname <- '2014_IpriorGP'
.stan_model <- 'gp_220108b'
DEBUG <- F
JOBID = 12
lab <- paste0("OnlyHTX_TRUE_threshold_0.5", '_jobname_', jobname)

.indir <- "/rds/general/user/mm3218/home/git/phyloflows"
datadir <- "/rds/general/user/mm3218/home/projects/2021/phyloflows"

if(0)
{
  .indir <- '~/git/phyloflows'
  datadir <- '~/Box\ Sync/2021/phyloflows/'
}
if(Sys.info()[['user']] == 'andrea')
{
  .indir <- '~/git/phyloflows'
  datadir <- '~/Documents/Box/2021/phyloflows'
}

.outdir <- file.path(datadir, lab)
datadir <- file.path(datadir, lab)

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-datadir')
  stopifnot(args_line[[5]]=='-outdir')
  stopifnot(args_line[[7]]=='-stan_model')
  stopifnot(args_line[[9]]=='-JOBID')
  stopifnot(args_line[[11]]=='-lab')
  .indir <- args_line[[2]]
  datadir <- args_line[[4]]
  .outdir <- args_line[[6]]
  .stan_model <- args_line[[8]]
  .JOBID <- args_line[[10]]
  lab <- args_line[[12]]
}

path.to.stan.data <- file.path(datadir, paste0("stanin_",lab,".RData"))
load(path.to.stan.data)

indir <- .indir; outdir <- .outdir; 
path.to.stan.model <- file.path(indir, 'stan_models', paste0(.stan_model, '.stan'))

# run stan model
model = rstan::stan_model(path.to.stan.model)

stan_init <- list()
stan_init$rho_gp1 = rep(1.25, stan_data$N_group)
stan_init$rho_gp2 = rep(1.25, stan_data$N_group)
stan_init$alpha = 0
stan_init$nu = rep(0, stan_data$N_group)

if(DEBUG){
  fit <- sampling(model, data = stan_data, iter = 10, warmup = 5, chains=1, thin=1, init=rep(list(stan_init),1))
}else{
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- sampling(model, data = stan_data,
                  iter = 3000, warmup = 500, chains=4, thin=1, seed = 5,
                  verbose = FALSE, control = list(adapt_delta = 0.999,max_treedepth=15),
                  init=rep(list(stan_init), 4))
}

sum = summary(fit)
sum$summary[which(sum$summary[,9] < 100),]

.JOBID <- 'TEST'
file = file.path(outdir, paste0(.stan_model,'-', .JOBID), paste0(.stan_model,'-', .JOBID, '_', lab, '.rds'))
cat("Save file ", file)
saveRDS(fit,file = file)
     
