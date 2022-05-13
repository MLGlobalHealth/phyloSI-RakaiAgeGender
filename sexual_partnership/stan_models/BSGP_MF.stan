functions
{
  matrix kron_mvprod(matrix A, matrix B, matrix V)
  {
    return transpose(A*transpose(B*V));
  }

  matrix gp(int N_rows, int N_columns, array[] real rows_idx, array[] real columns_index,
            real delta0,
            real alpha_gp,
            real rho_gp1, real rho_gp2,
            matrix z1)
  {
    matrix[N_rows,N_columns] GP;

    matrix[N_rows, N_rows] K1;
    matrix[N_rows, N_rows] L_K1;

    matrix[N_columns, N_columns] K2;
    matrix[N_columns, N_columns] L_K2;

    // compute the Exponentiated quadratic covariance function:
    K1 = gp_exp_quad_cov(rows_idx, alpha_gp, rho_gp1) + diag_matrix(rep_vector(delta0, N_rows));
    K2 = gp_exp_quad_cov(columns_index, alpha_gp, rho_gp2) + diag_matrix(rep_vector(delta0, N_columns));
    // Cholesky Decomposition K1 = L_K1 L_K1^T
    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);

    GP = kron_mvprod(L_K2, L_K1, z1);

    return(GP);
  }
}

data
{
  int<lower=1> A; // number of age inputs
  array[A] real age1; // age of contacting individuals
  array[A] real age2; // age of contacted individuals

  int<lower=0> Nmf; // number of observed male-male contacts
  // int<lower=0> Nfm;

  array[Nmf] int<lower=0> ymf; // the reported male-male contacts
  // array[Nfm] int<lower=0> yfm;

  array[Nmf] int<lower=1> ymf_rowmajor_matrix_index; // the row-major index of the observations in the full contact matrix
  // array[Nfm] int<lower=1> yfm_rowmajor_matrix_index;

  vector[Nmf] log_participants_mf; // the offset terms
  // vector[Nfm] log_participants_fm;
  vector[Nmf] log_pop_mf;
  // vector[Nfm] log_pop_fm;

  // B spline basis function
  int M_age1; // number of B-Splines basis functions for age of participants
  int M_age2; // number of B-Splines basis functions for age of contacted individuals
  matrix[M_age1, A] basis_age1; // B-splines basis functions on the age of participants (rows)
  matrix[M_age2, A] basis_age2; // B-splines basis functions on the age of contacted individuals (columns)
  array[M_age1] real idx_basis_age1; // index of the B-splines basis functions rows
  array[M_age2] real idx_basis_age2; // index of the B-splines basis functions columns
}

transformed data
{
  // int<lower=0> N = Nmf + Nfm; // all observations
  int<lower=1> A_squared = A*A;
  real gp_delta = 1e-9; // GP nugget
  int MF = 1; // index of gender combinations in parameter vectors
  // int FM = 2;
  // array[N] int y;
  // // append data
  // y = append_array(ymf, yfm);
}

parameters
{
  vector[2] log_random_effect_baseline;

  real<lower=0> gp_rho_age1; // GP hyperparameters
  real<lower=0> gp_rho_age2;
  real<lower=0> gp_alpha;

  matrix[M_age1, M_age2] z; // GP variables
}

transformed parameters
{
  matrix[M_age1, M_age2] f_beta_mf; // GPs on spline coefficients

  vector[A_squared] f_mf;// B-splines GP approx

  // GP prior for male-female and female-male contact rates
  f_beta_mf = gp(
    M_age1,
    M_age2,
    idx_basis_age1,
    idx_basis_age2,
    gp_delta,
    gp_alpha,
    gp_rho_age1,
    gp_rho_age2,
    z
    );
  f_mf = to_vector(transpose( (basis_age1') * f_beta_mf * basis_age2 ));
}

model
{
  // GP priors
  target += inv_gamma_lpdf(gp_rho_age1 | 5, 5);
  target += inv_gamma_lpdf(gp_rho_age2 | 5, 5);
  target += cauchy_lpdf(gp_alpha | 0, 1);
  target += std_normal_lpdf( to_vector(z) );

  // baseline
  target += normal_lpdf(log_random_effect_baseline | 0, 10);

  // model in local scope
  {
    vector[Nmf] log_mu = log_random_effect_baseline[MF] + f_mf[ymf_rowmajor_matrix_index] + log_pop_mf + log_participants_mf;
    target += poisson_lpmf( ymf | exp(log_mu));
  }
}

generated quantities
{
  array[Nmf] real log_lik;

  {
    //local scope

    //pointwise log likelihood
    vector[Nmf] log_mu = log_random_effect_baseline[MF] + f_mf[ymf_rowmajor_matrix_index] + log_pop_mf + log_participants_mf;
    for(i in 1:Nmf)
    {
      log_lik[i] = poisson_lpmf( ymf[i] | exp(log_mu[i]) );
    }
  }
}
