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
  int<lower=1> N_PERIOD; // number of period
  int<lower=1> N_ROUND; // number of round inland
  int<lower=1> N_PER_GROUP; // number of age-age entries
  int<lower=1> N_AGE; // N_AGE*N_AGE = N_PER_GROUP
  int y[N_PER_GROUP, N_DIRECTION, N_PERIOD]; // count of transmissions for each age-age entry
  real<lower=0> ir[N_AGE, N_DIRECTION, N_ROUND]; // incidence rates 
  real ir_lognorm_mean[N_AGE, N_DIRECTION, N_ROUND]; // incidence rates 
  real<lower=0> ir_lognorm_sd[N_AGE, N_DIRECTION, N_ROUND]; // incidence rates 
  vector[N_PER_GROUP] log_offset[N_DIRECTION, N_ROUND]; // offset including prop susceptible, number of unsuppressed
  vector[N_PER_GROUP] log_offset_time[N_DIRECTION, N_ROUND];  // offset time
  vector[N_PER_GROUP] log_offset_susceptible[N_DIRECTION, N_ROUND];  // log number of susceptible
  vector[N_PER_GROUP] log_prop_sampling[N_DIRECTION, N_PERIOD];  // offset probability of sampling
  int sampling_index_y[N_PER_GROUP, N_DIRECTION, N_PERIOD]; 
  int n_sampling_index_y[N_DIRECTION, N_PERIOD];
  int map_age_source[N_PER_GROUP];
  int map_age_recipient[N_PER_GROUP];
  int map_round_period[N_ROUND];
  int N_ROUND_PER_PERIOD[N_PERIOD];
  int N_OBS;
    
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
  real log_ir[N_AGE, N_DIRECTION, N_ROUND] = log(ir); 
  matrix[N_AGE, N_PER_GROUP] matrix_map_age_recipient = rep_matrix(0.0, N_AGE, N_PER_GROUP);
  for(i in 1:N_PER_GROUP){
    matrix_map_age_recipient[map_age_recipient[i], i] = 1;
  }
}


parameters {
  real log_beta_baseline;
  real log_beta_baseline_contrast_direction;
  real log_beta_baseline_contrast_round[N_ROUND - 1];
  real log_beta_baseline_contrast_period;
  
  real<lower=0> rho_gp_round_inland[N_ROUND - 1,N_DIRECTION];
  real<lower=0> alpha_gp_round_inland[N_ROUND - 1,N_DIRECTION];
  vector[num_basis_rows] z_round_inland[N_ROUND - 1, N_DIRECTION];
  
  real<lower=0> rho_gp1[N_DIRECTION];
  real<lower=0> rho_gp2[N_DIRECTION];
  real<lower=0> alpha_gp[N_DIRECTION];
  matrix[num_basis_rows,num_basis_columns] z1[N_DIRECTION];
}

transformed parameters {
  vector[N_PER_GROUP] lambda[N_DIRECTION, N_PERIOD] = rep_array(rep_vector(0.0, N_PER_GROUP), N_DIRECTION, N_PERIOD);  
  vector[N_PER_GROUP] log_lambda[N_DIRECTION, N_PERIOD];  
  vector[N_PER_GROUP] log_lambda_latent[N_DIRECTION, N_ROUND]; 
  vector[N_PER_GROUP] lambda_latent[N_DIRECTION, N_ROUND]; 
  vector[N_PER_GROUP] log_lambda_latent_peryear[N_DIRECTION, N_ROUND]; 
  vector[N_PER_GROUP] lambda_latent_peryear[N_DIRECTION, N_ROUND]; 
  vector[N_AGE] log_lambda_latent_recipient[N_DIRECTION, N_ROUND]; 
  vector[N_AGE] log_lambda_latent_peryear_recipient[N_DIRECTION, N_ROUND]; 
  vector[N_AGE] lambda_latent_peryear_recipient[N_DIRECTION, N_ROUND]; 
  vector[N_PER_GROUP] log_beta[N_DIRECTION, N_ROUND];  
  vector[N_PER_GROUP] log_beta_direction_contrast[N_DIRECTION];
  vector[N_PER_GROUP] log_beta_period_contrast[N_DIRECTION];
  matrix[N_ROUND - 1, N_PER_GROUP] log_beta_round_contrast[N_DIRECTION];
  matrix[num_basis_rows,num_basis_columns] low_rank_gp_direction[N_DIRECTION]; 
  vector[N_PER_GROUP] log_beta_baseline_contrast_round_inland[N_ROUND - 1, N_DIRECTION] = rep_array(rep_vector(0.0, N_PER_GROUP), N_ROUND - 1, N_DIRECTION);

  
  // start with baseline
  log_beta = rep_array(rep_vector(log_beta_baseline, N_PER_GROUP), N_DIRECTION, N_ROUND);
  
  for(i in 1:N_DIRECTION){
    
  // find direction contrast
  low_rank_gp_direction[i] = gp(num_basis_rows, num_basis_columns, IDX_BASIS_ROWS, IDX_BASIS_COLUMNS, delta0,
            alpha_gp[i], rho_gp1[i], rho_gp2[i], z1[i]);
  log_beta_direction_contrast[i] = to_vector(((BASIS_ROWS') * low_rank_gp_direction[i] * BASIS_COLUMNS)');
  if(i == 2){
    log_beta_direction_contrast[i] +=  rep_vector(log_beta_baseline_contrast_direction, N_PER_GROUP);
  }
    
  // find period contrast
  log_beta_period_contrast[i] = rep_vector(log_beta_baseline_contrast_period, N_PER_GROUP);
  
  // add period contrast to round contrast
  log_beta_round_contrast[i] = append_row(rep_matrix(rep_row_vector(0.0, N_PER_GROUP), N_ROUND_PER_PERIOD[1] - 1 ), 
                                              rep_matrix(to_row_vector(log_beta_period_contrast[i]), N_ROUND_PER_PERIOD[N_PERIOD]));

    for(k in 1:N_ROUND){
      
      // find round contrast
      if(k > 1 && k <= N_ROUND){
        log_beta_baseline_contrast_round_inland[k-1,i] = (BASIS_ROWS' * gp_1D(num_basis_rows, IDX_BASIS_ROWS, delta0, alpha_gp_round_inland[k-1,i], rho_gp_round_inland[k-1,i], z_round_inland[k-1,i]))[map_age_recipient];
        log_beta_round_contrast[i][k-1,:] += to_row_vector(log_beta_baseline_contrast_round_inland[k-1,i]) ;
        log_beta_round_contrast[i][k-1,:] += log_beta_baseline_contrast_round[k-1];
      }

      // add direction contrast
      log_beta[i,k] += log_beta_direction_contrast[i];
      
      // add round contrast
      if(k > 1){
        log_beta[i,k] += to_vector(log_beta_round_contrast[i][(k-1),:]); 
      }
        
      // lambda aggregated by age of the recipient per year per susceptible
      log_lambda_latent_peryear[i,k] = log_beta[i,k] + log_offset[i,k] - log_offset_susceptible[i,k];
      lambda_latent_peryear[i,k] = exp(log_lambda_latent_peryear[i,k]);
      log_lambda_latent_peryear_recipient[i,k] = log(matrix_map_age_recipient * lambda_latent_peryear[i,k]) ;
      lambda_latent_peryear_recipient[i,k]= exp(log_lambda_latent_peryear_recipient[i,k]);
      
      // lambda aggregated by age of the recipient
      log_lambda_latent[i,k] = log_beta[i,k] + log_offset[i,k] + log_offset_time[i,k];
      lambda_latent[i,k] = exp(log_lambda_latent[i,k]);
      log_lambda_latent_recipient[i,k] = log(matrix_map_age_recipient * lambda_latent[i,k]);
  
      // lambda aggregated by period and add sampling probability
      lambda[i,map_round_period[k]] += exp(log_lambda_latent[i,k] + log_prop_sampling[i,map_round_period[k]]);
      
      
    }
  }
  
  log_lambda = log(lambda);
}

model {
  
  //
  // PRIORS
  //
  
  // basline 
  log_beta_baseline ~ normal(0, 5);
  
  // direction, round and period contrasts
  log_beta_baseline_contrast_direction ~ normal(0, 1);
  log_beta_baseline_contrast_round ~ normal(0, 1);
  log_beta_baseline_contrast_period ~ normal(0, 1);
  
  // hyperparameters baseline surface on the age of the source and recipient 
  alpha_gp ~ cauchy(0,1);
  rho_gp1 ~ inv_gamma(2, 2);
  rho_gp2 ~ inv_gamma(2, 2);
  
  // baseline surface on the age of the source and recipient, standardised
  for(i in 1:num_basis_rows){
    for(j in 1:num_basis_columns){
      z1[:,i,j] ~ normal(0,1);
    }
  }

  for (i in 1:N_DIRECTION){

    for (k in 1:N_ROUND){
      
      if(k > 1){
        // hyperparameters round contrast over the age of the recipient
        rho_gp_round_inland[k-1,i] ~ inv_gamma(2, 2);
        alpha_gp_round_inland[k-1,i] ~ cauchy(0,1);
    
        // round contrast on the age of the recipient, standardised
        z_round_inland[k-1, i] ~ normal(0,1);
      }
    }
  }
        
  
  //
  // LIKELIHOOD
  //
  
  for (i in 1:N_DIRECTION){
    
    // likelihood on the incidence rate
    for (k in 1:N_ROUND){

      lambda_latent_peryear_recipient[i,k] ~ lognormal(ir_lognorm_mean[:,i,k], ir_lognorm_sd[:,i,k]);
    }
    
    // likelihood on the phylo pairs
    for (p in 1:N_PERIOD){
      if(n_sampling_index_y[i,p] > 0){
        y[sampling_index_y[1:n_sampling_index_y[i,p],i,p],i,p] ~ poisson_log(log_lambda[i,p][sampling_index_y[1:n_sampling_index_y[i,p],i,p]]);
      }
    }
  }
}

generated quantities{
  int y_predict[N_PER_GROUP, N_DIRECTION, N_PERIOD];
  int z_predict[N_PER_GROUP, N_DIRECTION, N_ROUND];
  real ir_predict[N_AGE, N_DIRECTION, N_ROUND];
  real log_lik[N_OBS];

  {
    int index = 1;
    for(i in 1:N_DIRECTION){
        
      // save the predicted incidence rate in another format
      for(n in 1:N_AGE){
        for(k in 1:N_ROUND){
            ir_predict[n,i,k] = lambda_latent_peryear_recipient[i,k][n];
        }
      }
      
        
      for(n in 1:N_PER_GROUP){
        
        for(k in 1:N_ROUND){
          
          // predict total transmissions
          z_predict[n,i,k] = poisson_log_rng(log_lambda_latent[i,k][n]);
          
        }
        
        for(p in 1:N_PERIOD){
          
          // save log likelihood values on the phylo pairs
          if(sampling_index_y[n,i,p] != -1){
            log_lik[index] = poisson_log_lpmf( y[n,i,p] | log_lambda[i,p][n] ) ;
            index = index + 1;
          }
          
          // predict detected transmissions
          y_predict[n,i,p] = poisson_log_rng(log_lambda[i,p][n]);
                    
        }
        
      }

    }
  }
}





