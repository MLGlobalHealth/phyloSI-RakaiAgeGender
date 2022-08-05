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
  real X[N];
  int COUNT_00[N];
  int COUNT_01[N];
  int COUNT_10[N];
  int COUNT_11[N];
  int TOTAL_COUNT_00;
  int TOTAL_COUNT_01;
  int TOTAL_COUNT_10;
  int TOTAL_COUNT_11;
  real rho_hyper_par_1;
  real rho_hyper_par_2;
}

transformed data {
  real delta0 = 1e-9;
}

parameters {
  real alpha_00;
  real alpha_01;
  real alpha_10;
  real alpha_11;
  real<lower=0> rho_00;
  real<lower=0> rho_01;
  real<lower=0> rho_10;
  real<lower=0> rho_11;
  real<lower=0> sigma_00;
  real<lower=0> sigma_01;
  real<lower=0> sigma_10;
  real<lower=0> sigma_11;
  vector[N] eta_00;
  vector[N] eta_01;
  vector[N] eta_10;
  vector[N] eta_11;
}

transformed parameters {
  vector[N] f_00 = gp(X, sigma_00, rho_00, delta0, eta_00);
  vector[N] f_01 = gp(X, sigma_01, rho_01, delta0, eta_01);
  vector[N] f_10 = gp(X, sigma_10, rho_10, delta0, eta_10);
  vector[N] f_11 = gp(X, sigma_11, rho_11, delta0, eta_11);
  vector[N] log_lambda_00 = alpha_00 + f_00;
  vector[N] log_lambda_01 = alpha_01 + f_01;
  vector[N] log_lambda_10 = alpha_10 + f_10;
  vector[N] log_lambda_11 = alpha_11 + f_11;
  vector[N] lambda_00 = exp(log_lambda_00);
  vector[N] lambda_01 = exp(log_lambda_01);
  vector[N] lambda_10 = exp(log_lambda_10);
  vector[N] lambda_11 = exp(log_lambda_11);
}

model {
  alpha_00 ~ normal(0, 100);
  alpha_01 ~ normal(0, 100);
  alpha_10 ~ normal(0, 100);
  alpha_11 ~ normal(0, 100);
  
  rho_00 ~ inv_gamma(rho_hyper_par_1, rho_hyper_par_2);
  rho_01 ~ inv_gamma(rho_hyper_par_1, rho_hyper_par_2);
  rho_10 ~ inv_gamma(rho_hyper_par_1, rho_hyper_par_2);
  rho_11 ~ inv_gamma(rho_hyper_par_1, rho_hyper_par_2);
  
  sigma_00 ~ normal(0, 10);
  sigma_01 ~ normal(0, 10);
  sigma_10 ~ normal(0, 10);
  sigma_11 ~ normal(0, 10);
  
  eta_00 ~ std_normal();
  eta_01 ~ std_normal();
  eta_10 ~ std_normal();
  eta_11 ~ std_normal();

  COUNT_00 ~ poisson_log(log_lambda_00);
  COUNT_01 ~ poisson_log(log_lambda_01);
  COUNT_10 ~ poisson_log(log_lambda_10);
  COUNT_11 ~ poisson_log(log_lambda_11);

  TOTAL_COUNT_00 ~ poisson(sum(lambda_00));
  TOTAL_COUNT_01 ~ poisson(sum(lambda_01));
  TOTAL_COUNT_10 ~ poisson(sum(lambda_10));
  TOTAL_COUNT_11 ~ poisson(sum(lambda_11));
}

