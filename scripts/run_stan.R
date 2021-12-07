library(rstan)
library(data.table)	

.indir <- '/rds/general/user/mm3218/home/git/phyloflows'
.outdir <- '/rds/general/user/mm3218/home/projects/2021/phyloflows'

if(0){
  .indir <- '~/git/phyloflows'
  .outdir <- '~/Box\ Sync/2021/RCCS/outputs'
}

lab <- "MRC_FALSE_OnlyHTX_TRUE_threshold_0.5"
stan_model <- 'gp_211207'
DEBUG <- F

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-stan_model')
  stopifnot(args_line[[7]]=='-JOBID')
  stopifnot(args_line[[9]]=='-lab')
  .indir <- args_line[[2]]
  .outdir <- args_line[[4]]
  stan_model <- args_line[[6]]
  JOBID <- as.numeric(args_line[[8]])
  lab <- as.numeric(args_line[[10]])
}

path.to.stan.data <- file.path(.outdir, lab, paste0("stanin_",lab,".RData"))
load(path.to.stan.data)

indir <- .indir; outdir <- .outdir; outdir.lab <- file.path(outdir, lab)
path.to.stan.model <- file.path(indir, 'stan_models', paste0(stan_model, '.stan'))

# run stan model
model = rstan::stan_model(path.to.stan.model)

if(DEBUG){
  fit <- sampling(model, data = stan_data, iter = 100, warmup = 50, chains=1, thin=1)
}else{
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  fit <- sampling(model, data = stan_data,
                  iter = 3000, warmup = 500, chains=4, thin=1, seed = 5,
                  verbose = FALSE, control = list(adapt_delta = 0.999,max_treedepth=15))
}

sum = summary(fit)
sum$summary[which(sum$summary[,9] < 100),]

file = file.path(outdir.lab, paste0(stan_model,'-', JOBID, '_', lab, '.rds'))
saveRDS(fit,file = file)
     