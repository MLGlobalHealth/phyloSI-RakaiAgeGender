
data{	
	int<lower=1> N_predict;
  	real x_predict[N_predict];
  	int<lower=1> N_observed;
  	int<lower=1, upper=N_predict> observed_idx[N_observed];
  	int y_observed_00[N_observed];
	int y_observed_10[N_observed];
	int y_observed_01[N_observed];	
	int y_observed_11[N_observed];	
	int total_observed_00[N_observed];
	int total_observed_10[N_observed];
	int total_observed_01[N_observed];
	int total_observed_11[N_observed];
  	real<lower=0> rho_hyper_par_00;
	real<lower=0> rho_hyper_par_10;
	real<lower=0> rho_hyper_par_01;
	real<lower=0> rho_hyper_par_11;  	
  	real<lower=0> alpha_hyper_par_00;
	real<lower=0> alpha_hyper_par_10;
	real<lower=0> alpha_hyper_par_01;
	real<lower=0> alpha_hyper_par_11;
}

parameters {
	real<lower=0> rho_00;
	real<lower=0> rho_10;
	real<lower=0> rho_01;
	real<lower=0> rho_11;
	real<lower=0> alpha_00;
	real<lower=0> alpha_10;
	real<lower=0> alpha_01;
	real<lower=0> alpha_11;
	real sex0_loc0;
	real sex1_loc0;
	real sex0_loc1;
	real sex1_loc1;  	
  	vector[N_predict] f_tilde_00;
	vector[N_predict] f_tilde_10;
	vector[N_predict] f_tilde_01;
	vector[N_predict] f_tilde_11;
}

transformed parameters {
  	matrix[N_predict, N_predict] L_cov;
  	vector[N_predict] logit_p_predict_00;
	vector[N_predict] logit_p_predict_10;
	vector[N_predict] logit_p_predict_01;
	vector[N_predict] logit_p_predict_11;
	// GP for 00 and 01 (women)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_00, rho_00) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_00 = sex0_loc0 + L_cov * f_tilde_00;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_01, rho_01) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_01 = sex0_loc1 + L_cov * f_tilde_01;
	// GP for 10 and 10 (men)
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_10, rho_10) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_10 = sex1_loc0 + L_cov * f_tilde_10;
	L_cov = cholesky_decompose(cov_exp_quad(x_predict, alpha_11, rho_11) + diag_matrix(rep_vector(1e-10, N_predict)));
	logit_p_predict_11 = sex1_loc1 + L_cov * f_tilde_11;
}

model {
	rho_00 ~ normal(0, rho_hyper_par_00);
	rho_10 ~ normal(0, rho_hyper_par_10);
	rho_01 ~ normal(0, rho_hyper_par_01);
	rho_11 ~ normal(0, rho_hyper_par_11);  	
  	alpha_00 ~ normal(0, alpha_hyper_par_00);
	alpha_10 ~ normal(0, alpha_hyper_par_10);
	alpha_01 ~ normal(0, alpha_hyper_par_01);
	alpha_11 ~ normal(0, alpha_hyper_par_11);
  	sex0_loc0 ~ normal( 0 , 10 );
	sex0_loc1 ~ normal( 0 , 10 );
	sex1_loc0 ~ normal( 0 , 10 );
	sex1_loc1 ~ normal( 0 , 10 );
  	f_tilde_00 ~ normal(0, 1);
	f_tilde_01 ~ normal(0, 1);
	f_tilde_10 ~ normal(0, 1);
	f_tilde_11 ~ normal(0, 1);
  	y_observed_00 ~ binomial_logit(total_observed_00, logit_p_predict_00[observed_idx] );
	y_observed_01 ~ binomial_logit(total_observed_01, logit_p_predict_01[observed_idx] );
	y_observed_10 ~ binomial_logit(total_observed_10, logit_p_predict_10[observed_idx] );
	y_observed_11 ~ binomial_logit(total_observed_11, logit_p_predict_11[observed_idx] );
}

generated quantities {
  	vector[N_predict] p_predict_00;
	vector[N_predict] p_predict_01;
	vector[N_predict] p_predict_10;
	vector[N_predict] p_predict_11;
	p_predict_00 = inv_logit(logit_p_predict_00);
	p_predict_01 = inv_logit(logit_p_predict_01);  
	p_predict_10 = inv_logit(logit_p_predict_10);
	p_predict_11 = inv_logit(logit_p_predict_11);
}			
