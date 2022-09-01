library(rstan)

indir.repository <- '~/git/phyloflows'

# load solver
file.stan.model <- file.path(indir.repository, 'misc', 'stan_models', 'gp_prior_tune.stan')
model = rstan::stan_model(file.stan.model)

# choose lower and upper bound for lenghscale
stan_data <- list(lower_b = 1, upper_b = (50-15)/2)

# fit
fit <- sampling(model, data=stan_data, iter=1, warmup=0, chains=1,
                seed=5838298, algorithm="Fixed_param")

# find solutions
sol <- extract(fit)
print(exp(sol$y_sol))#3.172883 8.696633
a = exp(sol$y_sol[1]); b = exp(sol$y_sol[2])

# plot
library(invgamma)
library(extraDistr)
x = seq(0, 40, 0.05)

plot(x, invgamma::dinvgamma(x,a, b), type = 'l')
lines(x, extraDistr::dtnorm(x, 0, diff(c(15, 49))/3, a = 0), col = 'red')
lines(x, invgamma::dinvgamma(x,2, 2), col = 'green')


