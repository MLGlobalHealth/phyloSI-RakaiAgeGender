functions {
  vector tail_delta(vector y, vector theta, real[] x_r, int[] x_i) {
    vector[2] deltas;
    deltas[1] = inv_gamma_cdf(theta[1], exp(y[1]), exp(y[2])) - 0.01;
    deltas[2] = 1 - inv_gamma_cdf(theta[2], exp(y[1]), exp(y[2])) - 0.01;
    return deltas;
  }
}

data{
  real lower_b;
  real upper_b;
}

transformed data {
  vector[2] y_guess = [log(10), log(20)]';
  vector[2] theta = [lower_b, upper_b]';
  vector[2] y;
  real x_r[0];
  int x_i[0];

  y = algebra_solver(tail_delta, y_guess, theta, x_r, x_i);
}

generated quantities {
  vector[2] y_sol = y;
}
