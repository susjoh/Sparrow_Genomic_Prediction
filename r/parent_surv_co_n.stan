data {
  int<lower=0> N;                        // Number of observations
  vector[N] age_q1;                     // Age
  vector[N] age_q2;                     // Age
  vector[N] f;                           // Inbreeding
  vector[N] co_n;                           // number cc measurements
  int<lower=0> N_ye;                     // Number of levels in year random effect
  int<lower=0,upper=N_ye> ye_idx[N];
  int<lower=0> N_ll;                     // Number of levels in last locality random effect
  int<lower=0,upper=N_ll> ll_idx[N];
  int<lower=0> N_id;                     // Number of levels in identity random effect
  int<lower=0,upper=N_id> id_idx[N];
  int<lower=0> N_par;                     // Number of levels in parent random effect
  int<lower=0,upper=N_par> par_idx[N];
  int<lower=0,upper=N_par> id_to_par_idx[N_id];
  vector[N_par] bv_mean;                  // Posterior means of breeding values
  matrix[N_par, N_par] bv_covmat_chol;     // Chol decomp of covariance for bvs
  real bv_mean_std;                  // Constant used to standardize the vector of breeding values
  real bv_sd_std;                    // Constant used to standardize the vector of breeding values
  int<lower=0,upper=1> Y[N];             // Response variable (yearly survival)
  real<lower=0> exp_rate_ars;                // Rate in exponential priors
  // Parameters for the priors on coefficients:
  real alpha_prior_mean_ars;
  real<lower=0> beta_prior_sd_ars;
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
  real mean_co_n = mean(co_n);
  real sd_co_n = sd(co_n);
  vector[N] co_n_std = (co_n - mean_co_n) / sd_co_n;
  real co_n_const = mean_co_n / sd_co_n; // Constant used in gen. quant. block

  // Constants need for gen. quant. block
  real bv_std_const1 = bv_mean_std / bv_sd_std;
  real bv_std_const2 = square(bv_std_const1);
  real bv_std_const3 = square(bv_sd_std);
  real bv_std_const4 = 2 * bv_mean_std / bv_std_const3;

  // Matrix containing fixed predictors
  matrix[N, 5] fixed_pred_mat = append_col(append_col(append_col(append_col(rep_vector(1, N), age_q1_std), age_q2_std), f_std), co_n_std);
}

parameters {
  vector[N_ye] z_ye;              // Std.normal noise for year random effect
  vector[N_ll] z_ll;              // Std.normal noise for last locality random effect
  vector[N_id] z_id;              // Std.normal noise for identity random effect
  vector[N_par] z_par;            // Std.normal noise for parent random effect
  vector[N_par] z_bv;             // Std.norm noise in bv
  real z_alpha;
  real z_beta_bv;
  real z_beta_bv2;
  real z_beta_age_q1;
  real z_beta_age_q2;
  real z_beta_f;
  real z_beta_co_n;
  real<lower=0,upper=1> sigma_ye_raw;
  real<lower=0,upper=1> sigma_ll_raw;
  real<lower=0,upper=1> sigma_id_raw;
  real<lower=0,upper=1> sigma_par_raw;
}

transformed parameters {
  real<lower=0> sigma_ye = -log(sigma_ye_raw) / exp_rate_ars;
  real<lower=0> sigma_ll = -log(sigma_ll_raw) / exp_rate_ars;
  real<lower=0> sigma_id = -log(sigma_id_raw) / exp_rate_ars;
  real<lower=0> sigma_par = -log(sigma_par_raw) / exp_rate_ars;
  vector[N_ye] ye = z_ye * sigma_ye;                // Levels in year random effect
  vector[N_ll] ll = z_ll * sigma_ll;                // Levels in last locality random effect
  vector[N_id] id_offset = z_id * sigma_id;                // Levels in identity random effect
  vector[N_par] par = z_par * sigma_par;             // Levels in parent random effect
  vector[N_id] id;


  // Full bv vector (non-centered parameterization berkson errored GP results)
  real alpha_std = alpha_prior_mean_ars + beta_prior_sd_ars * z_alpha;
  real beta_bv_std = beta_prior_sd_ars * z_beta_bv;
  real beta_bv2_std = beta_prior_sd_ars * z_beta_bv2 / sqrt(2);
  real beta_age_q1_std = beta_prior_sd_ars * z_beta_age_q1;
  real beta_age_q2_std = beta_prior_sd_ars * z_beta_age_q2;
  real beta_f_std = beta_prior_sd_ars * z_beta_f;
  real beta_co_n_std = beta_prior_sd_ars * z_beta_co_n;

  vector[N_par] bv_lat = bv_mean + bv_covmat_chol * z_bv;

  // Auxilliary variables
  // Vector of regression coefficients
  vector[7] beta_vec = [alpha_std, beta_age_q1_std, beta_age_q2_std, beta_f_std, beta_co_n_std, beta_bv_std, beta_bv2_std]';
  // Random intercepts
  vector[N] rand_inter;
  // Predictor matrix
  matrix[N, 7] x_mat;
  vector[N] bv_lat_full;
  for (i in 1:N) {
    bv_lat_full[i] = (bv_lat[par_idx[i]] - bv_mean_std) / bv_sd_std;
  }
  x_mat = append_col(append_col(fixed_pred_mat, bv_lat_full), square(bv_lat_full));
  // nested effect
  for (i in 1:N_id) {
    id[i] = par[id_to_par_idx[i]] + id_offset[i];
  }
  for (i in 1:N) {
    rand_inter[i] = ye[ye_idx[i]] + ll[ll_idx[i]] + id[id_idx[i]];
  }
}

model {
  // Priors
  sigma_ye_raw ~ uniform(0, 1);
  sigma_ll_raw ~ uniform(0, 1);
  sigma_id_raw ~ uniform(0, 1);
  sigma_par_raw ~ uniform(0, 1);

  // Non-centered parameterizations
  z_ye ~ std_normal();
  z_ll ~ std_normal();
  z_id ~ std_normal();
  z_par ~ std_normal();
  z_bv ~ std_normal();

  z_alpha ~ std_normal();
  z_beta_bv ~ std_normal();
  z_beta_bv2 ~ std_normal();
  z_beta_age_q1 ~ std_normal();
  z_beta_age_q2 ~ std_normal();
  z_beta_f ~ std_normal();
  z_beta_co_n ~ std_normal();

  // Likelihood
  Y ~ bernoulli_logit_glm(x_mat, rand_inter, beta_vec);
}

generated quantities {
  real alpha = alpha_std
  - beta_bv_std * bv_std_const1
  + beta_bv2_std * bv_std_const2
  - beta_age_q1_std * age_q1_const
  - beta_age_q2_std * age_q2_const
  - beta_f_std * f_const
  - beta_co_n_std * co_n_const;

  real beta_bv = beta_bv_std / bv_sd_std - beta_bv2_std * bv_std_const4;
  real beta_bv2 = beta_bv2_std / bv_std_const3;
  real beta_age_q1 = beta_age_q1_std / sd_age_q1;
  real beta_age_q2 = beta_age_q2_std / sd_age_q2;
  real beta_f = beta_f_std / sd_f;
  real beta_co_n = beta_co_n_std / sd_co_n;

  // Posterior predictions
  int y_rep[N]; // Posterior predictive samples
  y_rep = bernoulli_logit_glm_rng(x_mat, rand_inter, beta_vec);
}
