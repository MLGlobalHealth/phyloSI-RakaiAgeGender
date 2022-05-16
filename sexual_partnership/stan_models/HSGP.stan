functions
{
  vector diagSPD_EQ(real alpha, real rho, real L, int M)
  {
    return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }

  vector diagSPD_Matern32(real alpha, real rho, real L, int M)
  {
    return 2*alpha * (sqrt(3)/rho)^1.5 * inv((sqrt(3)/rho)^2 + ((pi()/2/L) * linspaced_vector(M, 1, M))^2);
  }

  matrix PHI(int N, int M, real L, vector x)
  {
    return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
  }

  matrix kron_mvprod(matrix A, matrix B, matrix V)
  {
    return (B*V) * transpose(A);
  }
}

data
{
  int<lower=1> A; // number of age inputs
  vector[A] age1; // age of contacting individuals
  vector[A] age2; // age of contacted individuals

 int<lower=0> Nmf; // number of observed male-male contacts
  int<lower=0> Nfm;

  array[Nmf] int<lower=0> ymf; // the reported male-male contacts
  array[Nfm] int<lower=0> yfm;

  array[Nmf] int<lower=1> ymf_rowmajor_matrix_index; // the row-major index of the observations in the full contact matrix
  array[Nfm] int<lower=1> yfm_rowmajor_matrix_index;

  vector[Nmf] log_participants_mf; // the offset terms
  vector[Nfm] log_participants_fm;
  vector[Nmf] log_pop_mf;
  vector[Nfm] log_pop_fm;

  // HSGP arguments
  real<lower=0> c_age1; // factor c to determine the boundary value L for age of participants (age 1)
  int<lower=1> M_age1; // number of basis functions
  real<lower=0> c_age2; // factor c to determine the boundary value L for age of contacted individuals (age 2)
  int<lower=1> M_age2; // number of basis functions
}

transformed data
{
  int<lower=0> N = Nmf + Nfm; // all observations
  int<lower=1> A_squared = A*A;
  real gp_delta = 1e-9; // GP nugget
  int MF = 1; // index of gender combinations in parameter vectors
  int FM = 2;
  array[N] int y;
 
  vector[A] age1_std; // standardised age of contacting individuals
  vector[A] age2_std; // standardised age of contacted individuals
  real L_age1;
  real L_age2;
  matrix[A, M_age1] PHI_age1;
  matrix[A, M_age2] PHI_age2;

 // append data
  y = append_array(ymf, yfm);

  // standardise age inputs
  age1_std = (age1 - mean(age1))/sd(age1);
  age2_std = (age2 - mean(age2))/sd(age2);

  // HSGP basis functions in both age dimensions, will be different if age1 different from age2
  L_age1 = c_age1*max(age1_std);
  L_age2 = c_age2*max(age2_std);
  PHI_age1 = PHI(A, M_age1, L_age1, age1_std);
  PHI_age2 = PHI(A, M_age2, L_age2, age2_std);
}

parameters
{
  vector[2] log_random_effect_baseline;
  real<lower=0> overdispersion;

 real<lower=0> gp_rho_age1; // GP hyperparameters
  real<lower=0> gp_rho_age2;
  real<lower=0> gp_alpha;

  matrix[M_age1, M_age2] z; // HSGP basis function coefficients for mm, ff, mf (which is symmetric to fm)
}

transformed parameters
{
  vector[A_squared] f_mf;// B-splines GP approx

  {
    // local scope
    vector[M_age1] hsgp_sqrt_spd_eq_mf_age1;// sqare root of spectral densities for mf in each dimension age1 age2
    vector[M_age2] hsgp_sqrt_spd_eq_mf_age2;

    // square root of spectral densities
    hsgp_sqrt_spd_eq_mf_age1 = diagSPD_EQ( gp_alpha, gp_rho_age1, L_age1, M_age1);
    hsgp_sqrt_spd_eq_mf_age2 = diagSPD_EQ( gp_alpha, gp_rho_age2, L_age2, M_age2);

    // Kronecker HSGP approximations for mf

    f_mf = to_vector(
      kron_mvprod(
        diag_post_multiply(PHI_age1, hsgp_sqrt_spd_eq_mf_age1),
        diag_post_multiply(PHI_age2, hsgp_sqrt_spd_eq_mf_age2),
        z
        )
      );
  }
}

model
{
  // GP priors
  target += inv_gamma_lpdf(gp_rho_age1 | 5, 5);
  target += inv_gamma_lpdf(gp_rho_age2 | 5, 5);
  target += cauchy_lpdf(gp_alpha | 0, 1);
  target += std_normal_lpdf( to_vector(z) );

  // overdispersion
  target += exponential_lpdf(overdispersion | 1);

  // baseline
  target += normal_lpdf(log_random_effect_baseline | 0, 10);

  // model in local scope
  {
    vector[N] log_mu =
      append_row(log_random_effect_baseline[MF] + f_mf[ymf_rowmajor_matrix_index] + log_pop_mf + log_participants_mf,
      log_random_effect_baseline[FM] + f_mf[yfm_rowmajor_matrix_index] + log_pop_fm + log_participants_fm
      );
    target += poisson_lpmf( y | exp(log_mu) );
  }
}

generated quantities
{
  array[N] real y_pred;
  {
    //local scope

    //pointwise log likelihood
    vector[N] log_mu =
      append_row(log_random_effect_baseline[MF] + f_mf[ymf_rowmajor_matrix_index] + log_pop_mf + log_participants_mf,
      log_random_effect_baseline[FM] + f_mf[yfm_rowmajor_matrix_index] + log_pop_fm + log_participants_fm
      );
    for(i in 1:N)
    {
      y_pred[i] = poisson_lpmf( y[i] | exp(log_mu[i]));
    }
  }
}
