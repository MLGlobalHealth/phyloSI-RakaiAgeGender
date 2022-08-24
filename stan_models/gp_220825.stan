functions {
    matrix kron_mvprod(matrix A, matrix B, matrix V) 
    {
        return transpose(A*transpose(B*V));
    }
    
  matrix gp(int N_rows, int N_columns, real[] rows_idx, real[] columns_index,
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
    
    K1 = cov_exp_quad(rows_idx, sqrt(alpha_gp), rho_gp1) + diag_matrix(rep_vector(delta0, N_rows));
    K2 = cov_exp_quad(columns_index, sqrt(alpha_gp), rho_gp2) + diag_matrix(rep_vector(delta0, N_columns));

    L_K1 = cholesky_decompose(K1);
    L_K2 = cholesky_decompose(K2);
    
    GP = kron_mvprod(L_K2, L_K1, z1);

    return(GP);
  }
  
  vector gp_1D(int N, real[] idx, 
            real delta0,
            real alpha_gp, 
            real rho_gp, 
            vector z)
  {
    
    vector[N] GP;
    
    matrix[N, N] K = cov_exp_quad(idx, sqrt(alpha_gp), rho_gp) + diag_matrix(rep_vector(delta0, N));
    matrix[N, N] L_K = cholesky_decompose(K);
    
    GP = L_K * z;

    return(GP);
  }
}

data {
  int<lower=1> N_DIRECTION; // number of directions
  int<lower=1> N_COMMUNITY; // number of community
  int<lower=1> N_PERIOD; // number of period
  int<lower=1> N_ROUND; // number of period
  int<lower=1> N_PER_GROUP; // number of age-age entries
  int<lower=1> N_AGE; // N_AGE*N_AGE = N_PER_GROUP
  int y[N_PER_GROUP, N_DIRECTION, N_COMMUNITY, N_PERIOD]; // count of transmissions for each age-age entry
  int z[N_AGE, N_DIRECTION, N_COMMUNITY, N_ROUND]; // count of transmissions for each age-age entry
  vector[N_PER_GROUP] log_offset[N_DIRECTION, N_COMMUNITY, N_ROUND]; 
  vector[N_PER_GROUP] log_prop_sampling[N_DIRECTION, N_COMMUNITY, N_PERIOD]; 
  int sampling_index_y[N_PER_GROUP, N_DIRECTION, N_COMMUNITY, N_PERIOD]; 
  int n_sampling_index_y[N_DIRECTION, N_COMMUNITY, N_PERIOD];
  int map_age_source[N_PER_GROUP];
  int map_age_recipient[N_PER_GROUP];
  int map_round_period[N_ROUND];
  
	//splines
	int number_rows; // = N_AGE
	int number_columns;  // = N_AGE
  int num_basis_rows;
  int num_basis_columns;
  matrix[num_basis_rows, number_rows] BASIS_ROWS; 
  matrix[num_basis_columns, number_columns] BASIS_COLUMNS; 
  
  // GP
  real IDX_BASIS_ROWS[num_basis_rows];
  real IDX_BASIS_COLUMNS[num_basis_columns];
  
}

transformed data
{   
  real delta0 = 1e-9;  
  matrix[N_AGE, N_PER_GROUP] matrix_map_age_recipient = rep_matrix(0.0, N_AGE, N_PER_GROUP);
  int N_ROUND_PER_PERIOD[N_PERIOD] = rep_array(0, N_PERIOD);
  for(i in 1:N_PER_GROUP){
    matrix_map_age_recipient[map_age_recipient[i], i] = 1;
  }
  for(t in 1:N_ROUND){
    N_ROUND_PER_PERIOD[map_round_period[t]] += 1;
  }
  print("N_ROUND_PER_PERIOD", N_ROUND_PER_PERIOD);
}


parameters {
  real log_beta_baseline;
  
  real log_beta_baseline_contrast_community;
  real log_beta_baseline_contrast_direction;
  
  real log_beta_baseline_contrast_round[N_ROUND - 1, N_DIRECTION,N_COMMUNITY];
  real<lower=0> sigma_beta_baseline_contrast_round[N_COMMUNITY];
  
  real<lower=0> rho_gp_community[N_DIRECTION];
  real<lower=0> alpha_gp_community[N_DIRECTION];
  vector[num_basis_rows] z_community[N_DIRECTION];
  
  real<lower=0> rho_gp_period[N_DIRECTION,N_COMMUNITY];
  real<lower=0> alpha_gp_period[N_DIRECTION,N_COMMUNITY];
  vector[num_basis_rows] z_period[N_DIRECTION,N_COMMUNITY];
  
  real<lower=0> rho_gp1[N_DIRECTION];
  real<lower=0> rho_gp2[N_DIRECTION];
  real<lower=0> alpha_gp[N_DIRECTION];
  matrix[num_basis_rows,num_basis_columns] z1[N_DIRECTION];
}

transformed parameters {
  vector[N_PER_GROUP] lambda[N_DIRECTION, N_COMMUNITY, N_PERIOD] = rep_array(rep_vector(0.0, N_PER_GROUP), N_DIRECTION, N_COMMUNITY, N_PERIOD);  
  vector[N_PER_GROUP] log_lambda[N_DIRECTION, N_COMMUNITY, N_PERIOD];  
  vector[N_PER_GROUP] log_lambda_latent[N_DIRECTION, N_COMMUNITY, N_ROUND]; 
  vector[N_PER_GROUP] lambda_latent[N_DIRECTION, N_COMMUNITY, N_ROUND]; 
  vector[N_AGE] log_lambda_latent_recipient[N_DIRECTION, N_COMMUNITY, N_ROUND]; 
  vector[N_PER_GROUP] log_beta[N_DIRECTION, N_COMMUNITY, N_ROUND];  
  vector[N_PER_GROUP] log_beta_direction_contrast[N_DIRECTION];
  vector[N_PER_GROUP] log_beta_community_contrast[N_DIRECTION];
  vector[N_PER_GROUP] log_beta_period_contrast[N_DIRECTION,N_COMMUNITY];
  matrix[N_ROUND - 1, N_PER_GROUP] log_beta_round_contrast[N_DIRECTION,N_COMMUNITY];
  matrix[num_basis_rows,num_basis_columns] low_rank_gp_direction[N_DIRECTION]; 
  
  // start with baseline
  log_beta = rep_array(rep_vector(log_beta_baseline, N_PER_GROUP), N_DIRECTION, N_COMMUNITY, N_ROUND);
  
  for(i in 1:N_DIRECTION){
    
    // find direction contrast
    low_rank_gp_direction[i] = gp(num_basis_rows, num_basis_columns, IDX_BASIS_ROWS, IDX_BASIS_COLUMNS, delta0,
              alpha_gp[i], rho_gp1[i], rho_gp2[i], z1[i]);
    log_beta_direction_contrast[i] = to_vector(((BASIS_ROWS') * low_rank_gp_direction[i] * BASIS_COLUMNS)');
    if(i == 2){
      log_beta_direction_contrast[i] = log_beta_direction_contrast[i] + rep_vector(log_beta_baseline_contrast_direction, N_PER_GROUP);
    }
    
     // find community contrast
     log_beta_community_contrast[i] = rep_vector(log_beta_baseline_contrast_community, N_PER_GROUP) + 
         (BASIS_ROWS' * gp_1D(num_basis_rows, IDX_BASIS_ROWS, delta0, alpha_gp_community[i], rho_gp_community[i], z_community[i]))[map_age_source];;

    
    for(j in 1:N_COMMUNITY){
      
          // find period contrast
    log_beta_period_contrast[i,j] = (BASIS_ROWS' * gp_1D(num_basis_rows, IDX_BASIS_ROWS, delta0, alpha_gp_period[i,j], rho_gp_period[i,j], z_period[i,j]))[map_age_source];
    log_beta_round_contrast[i,j] = rep_matrix(to_vector(log_beta_baseline_contrast_round[:,i,j]), N_PER_GROUP) 
                                    + append_row(rep_matrix(rep_row_vector(0.0, N_PER_GROUP), (N_ROUND_PER_PERIOD[1] - 1) ), 
                                                 rep_matrix(to_row_vector(log_beta_period_contrast[i,j]), N_ROUND_PER_PERIOD[N_PERIOD])) ;

      for(k in 1:N_ROUND){
        
        // add direction contrast
        log_beta[i,j,k] += log_beta_direction_contrast[i];
        
        // add community contrast
        if(j == 1){
          log_beta[i,j,k] += log_beta_community_contrast[i]; 
        }
        
        // add round contrast
        if(k > 1){
          log_beta[i,j,k] += to_vector(log_beta_round_contrast[i,j][(k-1),:]); 
        }
        
        // add offset
        log_lambda_latent[i,j,k] = log_beta[i,j,k] + log_offset[i,j,k];
        lambda_latent[i,j,k] = exp(log_lambda_latent[i,j,k]);
        log_lambda_latent_recipient[i,j,k] = log(matrix_map_age_recipient * lambda_latent[i,j,k]);
        
        // aggregate by period and add sampling probability
        lambda[i,j,map_round_period[k]] += exp(log_lambda_latent[i,j,k] + log_prop_sampling[i,j,map_round_period[k]]);
      }
      

    }

  }
  
  log_lambda = log(lambda);
  
}

model {
  log_beta_baseline ~ normal(0, 10);
  log_beta_baseline_contrast_community ~ normal(0, 10);
  log_beta_baseline_contrast_direction ~ normal(0, 10);
  
  sigma_beta_baseline_contrast_round~ cauchy(0,1);
  
  alpha_gp ~ cauchy(0,1);
  rho_gp1 ~ inv_gamma(2, 2);
  rho_gp2 ~ inv_gamma(2, 2);
  
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
      z1[:,i,j] ~ normal(0,1);
    }
  }
  
  for (i in 1:N_DIRECTION){
    
    alpha_gp_period[i] ~ cauchy(0,1);
    rho_gp_period[i] ~ inv_gamma(2, 2);
    
    alpha_gp_community[i] ~ cauchy(0,1);
    rho_gp_community[i] ~ inv_gamma(2, 2);
    z_community[i] ~ normal(0,1);
    
    log_beta_baseline_contrast_round[1,i,:] ~ normal(0, 10);
    for(k in 2:(N_ROUND - 1)){
      log_beta_baseline_contrast_round[k,i,:] ~ normal(log_beta_baseline_contrast_round[k - 1,i,:], sigma_beta_baseline_contrast_round);
    }
  
    for (j in 1:N_COMMUNITY){
        z_period[i,j] ~ normal(0,1);
          
       for (k in 1:N_ROUND){
        z[:,i,j,k] ~ poisson_log(log_lambda_latent_recipient[i,j,k]);
      }
      
      for (p in 1:N_PERIOD){
        y[sampling_index_y[1:n_sampling_index_y[i,j,p],i,j,p],i,j,p] ~ poisson_log(log_lambda[i,j,p][sampling_index_y[1:n_sampling_index_y[i,j,p],i,j,p]]);
      }
    }
  }
}

generated quantities{
  int y_predict[N_PER_GROUP, N_DIRECTION, N_COMMUNITY, N_PERIOD];
  int z_predict[N_PER_GROUP, N_DIRECTION, N_COMMUNITY, N_ROUND];
  real log_lik[N_PER_GROUP, N_DIRECTION, N_COMMUNITY, N_ROUND];

  for(i in 1:N_DIRECTION){
     for(j in 1:N_COMMUNITY){
       for(n in 1:N_PER_GROUP){
         for(k in 1:N_ROUND){
           z_predict[n,i,j,k] = poisson_log_rng(log_lambda_latent[i,j,k][n]);
           log_lik[n,i,j,k] = poisson_log_lpmf( z[map_age_source[n],i,j,k] | log_lambda_latent_recipient[i,j,k][map_age_source[n]] ) / N_AGE;
            if(sampling_index_y[n,i,j,map_round_period[k]] != -1){
              log_lik[n,i,j,k] += poisson_log_lpmf( y[n,i,j,map_round_period[k]] | log_lambda[i,j,map_round_period[k]][n] ) / N_ROUND_PER_PERIOD[map_round_period[k]];
            }
         }
         for(p in 1:N_PERIOD){
           y_predict[n,i,j,p] = poisson_log_rng(log_lambda[i,j,p][n]);
         }
       }
     }
  }
}





