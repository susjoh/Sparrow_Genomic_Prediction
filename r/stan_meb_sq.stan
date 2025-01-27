data {
  int<lower=0> N;                    // Number of observations
  int<lower=0> M;                    // Number of exact observations
  vector[M] X_exact;                 // Predictor observed perfectly
  vector[N-M] X_obs;                 // Predictor observed with Berkson error
  vector[N] Y;                       // Response variable
  real<lower=0> pc_par;
}

parameters {
  real alpha;                        // Intercept
  real beta_X;                       // Coefficient for X
  real beta_X2;                      // Coefficient for X^2
  real<lower=0> sigma;               // Residual sd
  real<lower=0> sigma_error;         // Berkson error sd
  vector[N-M] X_raw;                 // Std.norm noise in X
}

transformed parameters {

  vector[N-M] X_lat = X_obs + sigma_error * X_raw;  // Non-centered parameterization

}

model {
  // Priors
  alpha ~ normal(0, 1e1);
  beta_X ~ normal(0, 1e1);
  beta_X2 ~ normal(0, 1e1);
  sigma ~ exponential(pc_par); // PC-prior
  sigma_error ~ exponential(pc_par); // PC-prior

  // Non-centered prior for latent X to improve convergence
  X_raw ~ std_normal();

  // Likelihood for Y: Using the latent true X values
  Y ~ normal(alpha + beta_X * append_row(X_exact, X_lat) + beta_X2 * square(append_row(X_exact, X_lat)), sigma);
}

// generated quantities {
//   vector[N] Y_pred;
//   vector[N] X_squared = square(X); // Precompute square(X) once
//
//   for (n in 1:N) {
//     Y_pred[n] = normal_rng(alpha + beta_X * X[n] + beta_X2 * X_squared[n], sigma);
//   }
// }
