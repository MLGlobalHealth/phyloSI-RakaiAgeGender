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
  int<lower=1> N_per_group_reduced; // age-age-time entries with observation
  int y[N_per_group_reduced, N_group]; // count of transmissions for each age-age entry
  int<lower=0,upper=1> is_mf[N_group]; // is the direction mf
	
	// ages and time array
	int A; // number of ages
	int T; // number of times
	
	// reduced map
	int max_per_group;
	int<lower=-1,upper=N_per_group> MAP_TO_REDUCED[N_per_group_reduced, max_per_group];
	
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
  vector[N_per_group_reduced] log_lambda_reduced[N_group];  
  
  for(i in 1:N_group){
    vector[N_per_group] f;
    matrix[num_basis_1D,num_basis_2D_times_3D] beta;  
    vector[N_per_group] lambda;  
    vector[N_per_group_reduced] lambda_reduced;  
    
    beta = gp_3D(num_basis_1D, num_basis_2D, num_basis_3D,
                  IDX_BASIS_1D, IDX_BASIS_2D, IDX_BASIS_3D,
                  delta0,
                  alpha_gp[i], rho_gp1[i], rho_gp2[i], rho_gp3[i], 
                  z1[i]);
                  
    f = to_vector((BASIS_1D') * beta * kronecker_prod(BASIS_3D, BASIS_2D)); 
    
    log_lambda[i] = mu + f;
    
    if(is_mf[i])
      log_lambda[i] += nu;
    
    lambda = exp(log_lambda[i]);
    for(j in 1:N_per_group_reduced){
      lambda_reduced[j] = 0;
      for(k in MAP_TO_REDUCED[j,:]){
        if(k != -1){
          lambda_reduced[j] += lambda[k];
        }
      }
    }
    

    log_lambda_reduced[i] = log(lambda_reduced);
  }
  
  
}

model {
  alpha_gp ~ cauchy(0,1);
  rho_gp1 ~ inv_gamma(5, 5);
  rho_gp2 ~ inv_gamma(5, 5);
  rho_gp3 ~ inv_gamma(5, 5);
  mu ~ normal(0, 1);
  nu ~ normal(0, 1);
  
  for(i in 1:num_basis_1D){
    for(j in 1:num_basis_2D_times_3D){
      z1[:,i,j] ~ normal(0,1);
    }
  } 
      
  for (i in 1:N_group){
    y[:,i] ~ poisson_log(log_lambda_reduced[i]);
  }
}

generated quantities{
  real log_lik[N_per_group_reduced,N_group];
  for(i in 1:N_group){
    for(j in 1:N_per_group_reduced){
      log_lik[j,i] = poisson_log_lpmf(y[j,i]| log_lambda_reduced[i][j]);
    }
  }
}

