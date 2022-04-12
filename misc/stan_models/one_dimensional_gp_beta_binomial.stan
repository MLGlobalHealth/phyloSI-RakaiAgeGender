functions{
  vector gp(real[] x, 
            real alpha, real rho,
            real delta0, 
            vector eta) {
    int N = size(x);
    matrix[N, N] K = cov_exp_quad(x, alpha, rho) + diag_matrix(rep_vector(delta0, N));
    matrix[N, N] L_K = cholesky_decompose(K);
    return L_K * eta;
  }
}

data {
  int<lower=1> N;
  real x[N];
  int COUNT[N];
  int TOTAL_COUNT[N];
  
  int<lower=1> N_hat;
  real x_hat[N_hat];
  int x_in_x_hat[N];
}

transformed data {
  real delta0 = 1e-9;
}

parameters {
  real<lower=0> rho;
  real<lower=0> sigma;
  real<lower=0> nu;
  vector[N_hat] eta;
}

transformed parameters {
   vector[N_hat] f = gp(x_hat, sigma, rho, delta0, eta);
   vector<lower=0, upper=1>[N_hat] mu = inv_logit(f);
   vector<lower=0>[N_hat] alpha = mu * nu;
   vector<lower=0>[N_hat] beta = (1 - mu) * nu;
}

model {
  rho ~ inv_gamma(5, 5);
  nu ~ exponential(1);
  sigma ~ cauchy(0,1);
  eta ~ std_normal();

  COUNT ~ beta_binomial(TOTAL_COUNT, alpha[x_in_x_hat], beta[x_in_x_hat]);
}

generated quantities {
  real theta_hat[N_hat] = beta_rng(alpha, beta);
}
