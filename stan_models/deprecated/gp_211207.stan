functions {
  matrix kronecker_prod(matrix A, matrix B) {
    matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
    int m;
    int n;
    int p;
    int q;
    m = rows(A);
    n = cols(A);
    p = rows(B);
    q = cols(B);
    for (i in 1:m) {
      for (j in 1:n) {
        int row_start;
        int row_end;
        int col_start;
        int col_end;
        row_start = (i - 1) * p + 1;
        row_end = (i - 1) * p + p;
        col_start = (j - 1) * q + 1;
        col_end = (j - 1) * q + q;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }

  matrix kron_mvprod(matrix A, matrix B, matrix V) 
  {
      return B*V*(A');
  }
    
  matrix gp_3D(int N_1D, int N_2D, int N_3D,
               real[] idx_1D, real[] idx_2D, real[] idx_3D,
               real delta0,
               real alpha_gp,
               real rho_gp1, real rho_gp2, real rho_gp3,
               matrix z1)
  {
    
    matrix[N_1D,N_2D*N_3D] GP;
    
    matrix[N_1D, N_1D] K1;
    matrix[N_1D, N_1D] L_K1;
    
    matrix[N_2D, N_2D] K2;
    matrix[N_2D, N_2D] L_K2;
    
    matrix[N_3D, N_3D] K3;
    matrix[N_3D, N_3D] L_K3;
    
    K1 = cov_exp_quad(idx_1D, pow(alpha_gp, 1.0/3), rho_gp1) + diag_matrix(rep_vector(delta0, N_1D));
    K2 = cov_exp_quad(idx_2D, pow(alpha_gp, 1.0/3), rho_gp2) + diag_matrix(rep_vector(delta0, N_2D));
    K3 = cov_exp_quad(idx_3D, pow(alpha_gp, 1.0/3), rho_gp3) + diag_matrix(rep_vector(delta0, N_3D));

    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);
    L_K3 = cholesky_decompose(K3);
    
    GP = kron_mvprod(kronecker_prod(L_K3, L_K2), L_K1, z1);

    return(GP);
  }
}

data {
  int<lower=1> N_group; // number of directions
  int<lower=1> N_per_group; // age-age-time entries
  int<lower=1> N_non_missing_per_group; // age-age-time entries with observation
  int y[N_non_missing_per_group, N_group]; // count of transmissions for each age-age entry
  int<lower=0,upper=1> is_mf[N_group]; // is the direction mf
	int idx_obs[N_non_missing_per_group, N_group]; // coordinates 
	
	int A; // number of ages
	int T; // number of times
	
	//splines
  int num_basis_1D;
  int num_basis_2D;
  int num_basis_3D;
  matrix[num_basis_1D, A] BASIS_1D; 
  matrix[num_basis_2D, A] BASIS_2D; 
  matrix[num_basis_3D, T] BASIS_3D; 
  
  // GP
  real IDX_BASIS_1D[num_basis_1D];
  real IDX_BASIS_2D[num_basis_2D];
  real IDX_BASIS_3D[num_basis_3D];
}

transformed data
{   
  real delta0 = 1e-9;  
  int num_basis_2D_times_3D = num_basis_2D*num_basis_3D;
}

parameters {
  real<lower=0> rho_gp1[N_group];
  real<lower=0> rho_gp2[N_group];
  real<lower=0> rho_gp3[N_group];
  real<lower=0> alpha_gp[N_group];
  real nu;
  real mu;
  matrix[num_basis_1D,num_basis_2D_times_3D] z1[N_group];
}

transformed parameters {
  vector[N_per_group] log_lambda[N_group];  
  vector[N_per_group] f[N_group];
  matrix[num_basis_1D,num_basis_2D_times_3D] beta[N_group]; 
  
  for(i in 1:N_group){
    beta[i] = gp_3D(num_basis_1D, num_basis_2D, num_basis_3D,
                    IDX_BASIS_1D, IDX_BASIS_2D, IDX_BASIS_3D,
                    delta0,
                    alpha_gp[i], rho_gp1[i], rho_gp2[i], rho_gp3[i], 
                    z1[i]);

    f[i] = to_vector((BASIS_1D') * beta[i] * kronecker_prod(BASIS_3D, BASIS_2D));

    log_lambda[i] = mu + f[i];
    
    if(is_mf[i])
      log_lambda[i] += nu;
  }
  
}

model {
  alpha_gp ~ cauchy(0,1);
  rho_gp1 ~ inv_gamma(5, 5);
  rho_gp2 ~ inv_gamma(5, 5);
  rho_gp3 ~ inv_gamma(5, 5);
  mu ~ normal(0, 5);
  nu ~ normal(0, 5);
  
  for(i in 1:num_basis_1D){
    for(j in 1:num_basis_2D_times_3D){
      z1[:,i,j] ~ normal(0,1);
    }
  } 
      
  for (i in 1:N_group){
    y[:,i] ~ poisson_log(log_lambda[i][idx_obs[:,i]]);
  }
}

generated quantities{
  // int y_predict[N_per_group,N_group];
  real log_lik[N_non_missing_per_group,N_group];
  for(i in 1:N_group){
    for(j in 1:N_non_missing_per_group){
      // y_predict[j,i] = poisson_log_rng(log_lambda[i][j]);
      log_lik[j,i] = poisson_log_lpmf(y[j,i]| log_lambda[i][idx_obs[j,i]]);
    }
  }
}

