library(data.table)
library(ggplot2)
library(loo)

indir <- '~/git/phyloflows'
# indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
# indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
outdir.all <- "/rds/general/user/mm3218/home/projects/2021/phyloflows/"

# models
jobname <- c('central')
stan_model <- c( 'gp_221201d', 'gp_221208c', 'gp_221208d', 'gp_221208e', 'gp_221208f', 'gp_221208a', 'gp_221208b')
# stan_model <- c( 'gp_221201d', 'gp_221208c', 'gp_221208d', 'gp_221208e', 'gp_221208a', 'gp_221208b')

outdir.here <- file.path(outdir.all, paste0(stan_model, '-', jobname))
N <- length(stan_model)
# data_in <- file.path(outdir.here, paste0(stan_model, '-', jobname, '-stanin_', jobname, '.RData'))

# labels
label_models <- c('Central model', paste0('Alternative model: ', stan_model[2:N]))


# statistics_prediction.RDS  
# WAIC.rds 
# LOO.rds   

prettyround <- function(x,n) format(round(x, n), nsmall = n)

#
# load LOO
#

paths <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-LOO.rds'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  output[[i]] <- readRDS(paths[i])
}
loo_compare(output)


#
# load prediction 
#

paths <- file.path(outdir.here, 'tables', paste0(stan_model, '-', jobname, '-statistics_prediction.RDS'))
output <- vector(mode = 'list', length = N)
for(i in 1:N){
  output[[i]] <- as.data.table(readRDS(paths[i]))
  output[[i]][, type := label_models[i]]
}
output <- rbindlist(output)
output[, type := factor(type, levels= label_models)]
output <- output[order(type)]
output[, pairs_MAE := format(round(pairs_MAE, 4), nsmall = 4)]
output[, incidence_MAE := format(round(incidence_MAE, 5), nsmall = 5)]

# save
saveRDS(output, file = file.path(outdir.here[1], 'tables', paste0(stan_model[1], '-', jobname[1], '-sensitivity-model_fit.rds')))

