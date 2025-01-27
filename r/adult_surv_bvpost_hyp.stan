data {
  int<lower=0> N;                        // Number of observations
  vector[N] age;                         // Age
  vector[N] f;                           // Inbreeding
  int<lower=0> N_ye;                     // Number of levels in year random effect
  int<lower=0,upper=N_ye> ye_idx[N];
  int<lower=0> N_ll;                     // Number of levels in last locality random effect
  int<lower=0,upper=N_ll> ll_idx[N];
  int<lower=0> N_id;                     // Number of levels in identity random effect
  int<lower=0,upper=N_id> id_idx[N];
  vector[N_id] bv_mean;                  // Posterior means of breeding values
  vector<lower=0>[N_id] bv_sd;           // Posterior sd.s of breeding values
  int<lower=0,upper=1> Y[N];             // Response variable (yearly survival)
  real<lower=0> exp_rate;                // Rate in exponential priors
  // Parameters for the priors on coefficients:
  real<lower=0> hyperprior_shape;
  real<lower=0> hyperprior_rate;
  real<lower=0> quad_scale;
  real alpha_prior_mean;
}

transformed data {
  // Center/scale predictors
  real mean_age = mean(age);
  real sd_age = sd(age);
  real sd_age2 = square(sd_age); // Constant used in gen. quant. block
  real age_const = mean_age / sd_age; // Constant used in gen. quant. block
  real age_const2 = mean_age / square(sd_age); // Constant used in gen. quant. block
  real age2_const = square(age_const); // Constant used in gen. quant. block
  vector[N] age_std = (age - mean_age) / sd_age;
  real mean_f = mean(f);
  real sd_f = sd(f);
  vector[N] f_std = (f - mean_f) / sd_f;
  real f_const = mean_f / sd_f; // Constant used in gen. quant. block

  // Matrix containing fixed predictors
  matrix[N, 4] fixed_pred_mat = append_col(append_col(append_col(rep_vector(1, N),
                                                                 age_std),
                                                      square(age_std)),
                                           f_std);
}

parameters {
  vector[N_ye] z_ye;              // Std.normal noise for year random effect
  vector[N_ll] z_ll;              // Std.normal noise for last locality random effect
  vector[N_id] z_id;              // Std.normal noise for identity random effect
  vector[N] z_res;                // Std.normal noise for residual random effect
  vector[N_id] z_bv;              // Std.norm noise in bv
  real z_alpha;
  real z_beta_bv;
  real z_beta_bv2;
  real z_beta_age;
  real z_beta_age2;
  real z_beta_f;
  real<lower=0,upper=1> sigma_ye_raw;
  real<lower=0,upper=1> sigma_ll_raw;
  real<lower=0,upper=1> sigma_id_raw;
  real<lower=0,upper=1> sigma_res_raw;
  // real<lower=0> sigma_res;     // Residual sd
  // real<lower=0> sigma_ye;         // Hatch year sd
  // real<lower=0> sigma_ll;         // Last locality sd
  // real<lower=0> sigma_id;         // Individual sd
  // real<lower=0> sigma_res;     // Residual sd
  real<lower=0> beta_prior_sd;    // Prior sd for regression coeffs
}

transformed parameters {
  real<lower=0> sigma_ye = -log(sigma_ye_raw) / exp_rate;
  real<lower=0> sigma_ll = -log(sigma_ll_raw) / exp_rate;
  real<lower=0> sigma_id = -log(sigma_id_raw) / exp_rate;
  real<lower=0> sigma_res = -log(sigma_res_raw) / exp_rate;
  vector[N_ye] ye = z_ye * sigma_ye;                // Levels in year random effect
  vector[N_ll] ll = z_ll * sigma_ll;                // Levels in last locality random effect
  vector[N_id] id = z_id * sigma_id;                // Levels in identity random effect
  vector[N] res = z_res * sigma_res;             // Levels in residual random effect

  // Full bv vector (non-centered parameterization berkson errored GP results)
  real alpha_std = alpha_prior_mean + beta_prior_sd * z_alpha;
  real beta_bv_std = beta_prior_sd * z_beta_bv;
  real beta_bv2_std = (beta_prior_sd / quad_scale) * z_beta_bv2;
  real beta_age_std = beta_prior_sd * z_beta_age;
  real beta_age2_std = (beta_prior_sd / quad_scale) * z_beta_age2;
  real beta_f_std = beta_prior_sd * z_beta_f;

  vector[N_id] bv_lat = bv_mean + bv_sd .* z_bv;
  vector[N] bv_lat_full;
  real bv_lat_full_mean;
  real bv_lat_full_sd;
  for (i in 1:N) {
    bv_lat_full[i] = bv_lat[id_idx[i]];
  }
  bv_lat_full_mean = mean(bv_lat_full);
  bv_lat_full_sd = sd(bv_lat_full);

}

model {
  // Auxilliary variables
  // Vector of regression coefficients
  vector[6] beta_vec = [alpha_std, beta_age_std, beta_age2_std, beta_f_std, beta_bv_std, beta_bv2_std]';
  // Random intercepts
  vector[N] rand_inter;
  // Standardized breeding values
  vector[N] bv_std = (bv_lat_full - bv_lat_full_mean) / bv_lat_full_sd;
  // Predictor matrix
  matrix[N, 6] x_mat = append_col(append_col(fixed_pred_mat, bv_std), square(bv_std));

  for (i in 1:N) {
    rand_inter[i] = ye[ye_idx[i]] + ll[ll_idx[i]] + id[id_idx[i]]; // + res[i];
  }

  // Priors
  beta_prior_sd ~ gamma(hyperprior_shape, hyperprior_rate); // hyperprior
  // sigma_ye ~ exponential(exp_rate);          // PC-prior for hatch year sd
  // sigma_ll ~ exponential(exp_rate);          // PC-prior for last locality sd
  // sigma_id ~ exponential(exp_rate);          // PC-prior for individual sd
  sigma_ye_raw ~ uniform(0, 1);
  sigma_ll_raw ~ uniform(0, 1);
  sigma_id_raw ~ uniform(0, 1);
  sigma_res_raw ~ uniform(0, 1);

  // Non-centered parameterizations
  z_alpha ~ std_normal();
  z_ye ~ std_normal();
  z_ll ~ std_normal();
  z_id ~ std_normal();
  z_res ~ std_normal();
  z_bv ~ std_normal();
  z_beta_bv ~ std_normal();
  z_beta_bv2 ~ std_normal();
  z_beta_age ~ std_normal();
  z_beta_age2 ~ std_normal();
  z_beta_f ~ std_normal();

  // Likelihood
  Y ~ bernoulli_logit_glm(x_mat, rand_inter, beta_vec);
}

generated quantities {
  real alpha = alpha_std
  - beta_bv_std * bv_lat_full_mean / bv_lat_full_sd
  + beta_bv2_std * square(bv_lat_full_mean / bv_lat_full_sd)
  - beta_age_std * age_const
  + beta_age2_std * age2_const
  - beta_f_std * f_const;

  real beta_bv = beta_bv_std / bv_lat_full_sd - 2 * beta_bv2_std * bv_lat_full_mean / square(bv_lat_full_sd);
  real beta_bv2 = beta_bv2_std / square(bv_lat_full_sd);
  real beta_age = beta_age_std / sd_age - 2 * beta_age2_std * age_const2;
  real beta_age2 = beta_age2_std / sd_age2;
  real beta_f = beta_f_std / sd_f;
}
