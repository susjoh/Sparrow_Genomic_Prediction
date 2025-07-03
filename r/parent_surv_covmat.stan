data {
  int<lower=0> N;                        // Number of observations
  vector[N] age_q1;                     // Age
  vector[N] age_q2;                     // Age
  vector[N] f;                           // Inbreeding
  int<lower=0> N_ye;                     // Number of levels in year random effect
  int<lower=0,upper=N_ye> ye_idx[N];
  int<lower=0> N_ll;                     // Number of levels in last locality random effect
  int<lower=0,upper=N_ll> ll_idx[N];
  int<lower=0> N_id;                     // Number of levels in identity random effect
  int<lower=0,upper=N_id> id_idx[N];
  int<lower=0> N_par;                     // Number of levels in identity random effect
  int<lower=0,upper=N_par> par_idx[N];
  vector[N_par] bv_mean;                  // Posterior means of breeding values
  matrix[N_par, N_par] bv_covmat_chol;     // Chol decomp of covariance for bvs
  real bv_mean_std;                  // Constant used to standardize the vector of breeding values
  real bv_sd_std;                    // Constant used to standardize the vector of breeding values
  int<lower=0,upper=1> Y[N];             // Response variable (yearly survival)
  real<lower=0> exp_rate;                // Rate in exponential priors
  // Parameters for the priors on coefficients:
  real alpha_prior_mean;
  real<lower=0> beta_prior_sd;
}

transformed data {
  // Center/scale predictors
  real mean_age_q1 = mean(age_q1);
  real mean_age_q2 = mean(age_q2);
  real sd_age_q1 = sd(age_q1);
  real sd_age_q2 = sd(age_q2);
  vector[N] age_q1_std = (age_q1 - mean_age_q1) / sd_age_q1;
  vector[N] age_q2_std = (age_q2 - mean_age_q2) / sd_age_q2;
  real age_q1_const = mean_age_q1 / sd_age_q1; // Constant used in gen. quant. block
  real age_q2_const = mean_age_q2 / sd_age_q2; // Constant used in gen. quant. block
  real mean_f = mean(f);
  real sd_f = sd(f);
  vector[N] f_std = (f - mean_f) / sd_f;
  real f_const = mean_f / sd_f; // Constant used in gen. quant. block

  // Constants need for gen. quant. block
  real bv_std_const1 = bv_mean_std / bv_sd_std;
  real bv_std_const2 = square(bv_std_const1);
  real bv_std_const3 = square(bv_sd_std);
  real bv_std_const4 = 2 * bv_mean_std / bv_std_const3;

  // Matrix containing fixed predictors
  matrix[N, 4] fixed_pred_mat = append_col(append_col(append_col(rep_vector(1, N), age_q1_std), age_q2_std), f_std);
}

parameters {
  vector[N_ye] z_ye;              // Std.normal noise for year random effect
  vector[N_ll] z_ll;              // Std.normal noise for last locality random effect
  vector[N_id] z_id;              // Std.normal noise for identity random effect
  // vector[N] z_res;                // Std.normal noise for residual random effect
  vector[N_par] z_bv;              // Std.norm noise in bv
  real z_alpha;
  real z_beta_bv;
  real z_beta_bv2;
  real z_beta_age_q1;
  real z_beta_age_q2;
  real z_beta_f;
  real<lower=0,upper=1> sigma_ye_raw;
  real<lower=0,upper=1> sigma_ll_raw;
  real<lower=0,upper=1> sigma_id_raw;
  // real<lower=0,upper=1> sigma_res_raw;
}

transformed parameters {
  real<lower=0> sigma_ye = -log(sigma_ye_raw) / exp_rate;
  real<lower=0> sigma_ll = -log(sigma_ll_raw) / exp_rate;
  real<lower=0> sigma_id = -log(sigma_id_raw) / exp_rate;
  // real<lower=0> sigma_res = -log(sigma_res_raw) / exp_rate;
  vector[N_ye] ye = z_ye * sigma_ye;                // Levels in year random effect
  vector[N_ll] ll = z_ll * sigma_ll;                // Levels in last locality random effect
  vector[N_id] id = z_id * sigma_id;                // Levels in identity random effect
  // vector[N] res = z_res * sigma_res;             // Levels in residual random effect

  // Full bv vector (non-centered parameterization berkson errored GP results)
  real alpha_std = alpha_prior_mean + beta_prior_sd * z_alpha;
  real beta_bv_std = beta_prior_sd * z_beta_bv;
  real beta_bv2_std = beta_prior_sd * z_beta_bv2 / sqrt(2);
  real beta_age_q1_std = beta_prior_sd * z_beta_age_q1;
  real beta_age_q2_std = beta_prior_sd * z_beta_age_q2;
  real beta_f_std = beta_prior_sd * z_beta_f;

  vector[N_par] bv_lat = bv_mean + bv_covmat_chol * z_bv;

  // Auxilliary variables
  // Vector of regression coefficients
  vector[6] beta_vec = [alpha_std, beta_age_q1_std, beta_age_q2_std, beta_f_std, beta_bv_std, beta_bv2_std]';
  // Random intercepts
  vector[N] rand_inter;
  // Predictor matrix
  matrix[N, 6] x_mat;
  vector[N] bv_lat_full;
  for (i in 1:N) {
    bv_lat_full[i] = (bv_lat[par_idx[i]] - bv_mean_std) / bv_sd_std;
  }
  x_mat = append_col(append_col(fixed_pred_mat, bv_lat_full), square(bv_lat_full));
  for (i in 1:N) {
    rand_inter[i] = ye[ye_idx[i]] + ll[ll_idx[i]] + id[id_idx[i]]; // + res[i];
  }
}

model {
  // Priors
  sigma_ye_raw ~ uniform(0, 1);
  sigma_ll_raw ~ uniform(0, 1);
  sigma_id_raw ~ uniform(0, 1);
  // sigma_res_raw ~ uniform(0, 1);

  // Non-centered parameterizations
  z_ye ~ std_normal();
  z_ll ~ std_normal();
  z_id ~ std_normal();
  // z_res ~ std_normal();
  z_bv ~ std_normal();

  z_alpha ~ std_normal();
  z_beta_bv ~ std_normal();
  z_beta_bv2 ~ std_normal();
  z_beta_age_q1 ~ std_normal();
  z_beta_age_q2 ~ std_normal();
  z_beta_f ~ std_normal();

  // Likelihood
  Y ~ bernoulli_logit_glm(x_mat, rand_inter, beta_vec);
}

generated quantities {
  real alpha = alpha_std
  - beta_bv_std * bv_std_const1
  + beta_bv2_std * bv_std_const2
  - beta_age_q1_std * age_q1_const
  - beta_age_q2_std * age_q2_const
  - beta_f_std * f_const;

  real beta_bv = beta_bv_std / bv_sd_std - beta_bv2_std * bv_std_const4;
  real beta_bv2 = beta_bv2_std / bv_std_const3;
  real beta_age_q1 = beta_age_q1_std / sd_age_q1;
  real beta_age_q2 = beta_age_q2_std / sd_age_q2;
  real beta_f = beta_f_std / sd_f;

  // Posterior predictions
  int y_rep[N]; // Posterior predictive samples
  y_rep = bernoulli_logit_glm_rng(x_mat, rand_inter, beta_vec);
}
