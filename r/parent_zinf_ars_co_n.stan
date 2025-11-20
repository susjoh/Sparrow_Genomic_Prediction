data {
  int<lower=0> N;                    // Number of observations
  vector[N] age_q1;                  // Age (orth. poly of degree 1)
  vector[N] age_q2;                  // Age (orth. poly of degree 2)
  vector[N] f;                       // Inbreeding
  vector[N] co_n;                    // numer of cc measurements
  int<lower=0> N_ye;                 // Number of levels in year random effect
  int<lower=0,upper=N_ye> ye_idx[N];
  int<lower=0> N_ll;                 // Number of levels in last locality random effect
  int<lower=0,upper=N_ll> ll_idx[N];
  int<lower=0> N_id;                 // Number of levels in identity random effect
  int<lower=0,upper=N_id> id_idx[N];
  int<lower=0> N_par;                 // Number of levels in identity random effect
  int<lower=0,upper=N_par> par_idx[N];
  vector[N_par] bv_mean;              // Posterior means of breeding values
  matrix[N_par, N_par] bv_covmat_chol;     // Chol decomp of covariance for bvs
  real bv_mean_std;                  // Constant used to standardize the vector of breeding values
  real bv_sd_std;                    // Constant used to standardize the vector of breeding values
  int<lower=0> Y[N];                 // Response variable (yearly number of offspring)
  // Rate in exponential priors
  real<lower=0> exp_rate_ars;
  real<lower=0> exp_rate_zi;
  // real<lower=0> phi_inv_rate;
  // Parameters for the priors on coefficients:
  real alpha_prior_mean_ars;
  real<lower=0> beta_prior_sd_ars;
  real alpha_zi_prior_mean;
  real<lower=0> beta_zi_prior_sd;
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
  // Std.normal noise for random effect levels
  vector[N_ye] z_ye;
  vector[N_ll] z_ll;
  vector[N_id] z_id;
  vector[N_par] z_par;
  // vector[N_ye] z_ye_zi;
  // vector[N_ll] z_ll_zi;
  // vector[N_id] z_id_zi;
  // vector[N_par] z_par_zi;
  // Std.norm noise in breeding value
  vector[N_par] z_bv;
  // Uniform noise for
  // real<lower=0,upper=1> phi_inv_raw;
  // Uniform noise for random effects standard deviations:
  real<lower=0,upper=1> sigma_ye_raw;
  real<lower=0,upper=1> sigma_ll_raw;
  real<lower=0,upper=1> sigma_id_raw;
  real<lower=0,upper=1> sigma_par_raw;
  // real<lower=0,upper=1> sigma_ye_zi_raw;
  // real<lower=0,upper=1> sigma_ll_zi_raw;
  // real<lower=0,upper=1> sigma_id_zi_raw;
  // real<lower=0,upper=1> sigma_par_zi_raw;
  // real<lower=0,upper=1> sigma_res_zi_raw;
  // Std.normal noise for count component coefficients:
  real z_alpha;
  real z_beta_bv;
  real z_beta_bv2;
  real z_beta_age_q1;
  real z_beta_age_q2;
  real z_beta_f;
  real z_beta_co_n;
  // Std.normal noise for zero-inflation component coefficients:
  real z_alpha_zi;
  // real z_beta_zi_bv;
  // real z_beta_zi_bv2;
  // real z_beta_zi_age_q1;
  // real z_beta_zi_age_q2;
  // real z_beta_zi_f;
}

transformed parameters {
  //  Transform random effect standard deviations from uniform to exponential
  real<lower=0> sigma_ye = -log(sigma_ye_raw) / exp_rate_ars;
  real<lower=0> sigma_ll = -log(sigma_ll_raw) / exp_rate_ars;
  real<lower=0> sigma_id = -log(sigma_id_raw) / exp_rate_ars;
  real<lower=0> sigma_par = -log(sigma_par_raw) / exp_rate_ars;
  // real<lower=0> sigma_ye_zi = -log(sigma_ye_zi_raw) / exp_rate_ars;
  // real<lower=0> sigma_ll_zi = -log(sigma_ll_zi_raw) / exp_rate_ars;
  // real<lower=0> sigma_id_zi = -log(sigma_id_zi_raw) / exp_rate_ars;
  // real<lower=0> sigma_par_zi = -log(sigma_par_zi_raw) / exp_rate_ars;
  // Levels in random effects:
  vector[N_ye] ye = z_ye * sigma_ye;
  vector[N_ll] ll = z_ll * sigma_ll;
  vector[N_id] id = z_id * sigma_id;
  vector[N_par] par = z_par * sigma_par;
  // vector[N_ye] ye_zi = z_ye_zi * sigma_ye_zi;
  // vector[N_ll] ll_zi = z_ll_zi * sigma_ll_zi;
  // vector[N_id] id_zi = z_id_zi * sigma_id_zi;
  // vector[N_par] par_zi = z_par_zi * sigma_par_zi;
  // 1 / (overdispersion parameter)
  // real<lower=0> phi_inv = -log(phi_inv_raw) / phi_inv_rate;
  // Overdispersion parameter
  // real<lower=0> phi = 1 / phi_inv;

  // Count component coefficients
  real alpha_std = alpha_prior_mean_ars + beta_prior_sd_ars * z_alpha;
  real beta_bv_std = beta_prior_sd_ars * z_beta_bv;
  real beta_bv2_std = beta_prior_sd_ars * z_beta_bv2 / sqrt(2);
  real beta_age_q1_std = beta_prior_sd_ars * z_beta_age_q1;
  real beta_age_q2_std = beta_prior_sd_ars * z_beta_age_q2;
  real beta_f_std = beta_prior_sd_ars * z_beta_f;
  real beta_co_n_std = beta_prior_sd_ars * z_beta_co_n;
  // Standardized zero-inflation regression coefficients
  real alpha_zi_std = alpha_zi_prior_mean + beta_zi_prior_sd * z_alpha_zi;
  // real beta_zi_bv_std = beta_zi_prior_sd * z_beta_zi_bv;
  // real beta_zi_bv2_std = beta_zi_prior_sd * z_beta_zi_bv2 / sqrt(2);
  // real beta_zi_age_q1_std = beta_zi_prior_sd * z_beta_zi_age_q1;
  // real beta_zi_age_q2_std = beta_zi_prior_sd * z_beta_zi_age_q2;
  // real beta_zi_f_std = beta_zi_prior_sd * z_beta_zi_f;

  // Full bv vector (non-centered parameterization berkson errored GP results)
  vector[N_par] bv_lat = bv_mean + bv_covmat_chol * z_bv;

  // Auxilliary variables
  // Vector of regression coefficients
  vector[7] beta_vec = [alpha_std, beta_age_q1_std, beta_age_q2_std, beta_f_std, beta_co_n_std, beta_bv_std, beta_bv2_std]';
  // Random intercepts
  vector[N] rand_inter;
  // Predictor matrix
  matrix[N, 7] x_mat;
  // Zero-inflated linear predictor
  vector[N] logit_theta;
  vector[N] theta;

  vector[N] bv_lat_full;
  for (i in 1:N) {
    bv_lat_full[i] = (bv_lat[par_idx[i]] - bv_mean_std) / bv_sd_std;
  }

  x_mat = append_col(append_col(fixed_pred_mat, bv_lat_full), square(bv_lat_full));

  for (i in 1:N) {
    rand_inter[i] = ye[ye_idx[i]] + ll[ll_idx[i]] + id[id_idx[i]] + par[par_idx[i]];
  }

  // Zero-inflation linear predictor
  for (i in 1:N) {
    logit_theta[i] = alpha_zi_std;
    // + beta_zi_bv_std * bv_lat_full[i]
    // + beta_zi_bv2_std * square(bv_lat_full[i])
    // + beta_zi_age_q1_std * age_q1_std[i]
    // + beta_zi_age_q2_std * age_q2_std[i]
    // + beta_zi_f_std * f_std[i]
    // + ye_zi[ye_idx[i]] + ll_zi[ll_idx[i]] + id_zi[id_idx[i]]+ par_zi[par_idx[i]];
  }
  theta = inv_logit(logit_theta);
}

model {

  // Priors
  // Non-centered parameterizations to improve convergence
  z_bv ~ std_normal();

  z_ye ~ std_normal();
  z_ll ~ std_normal();
  z_id ~ std_normal();
  z_par ~ std_normal();

  // z_ye_zi ~ std_normal();
  // z_ll_zi ~ std_normal();
  // z_id_zi ~ std_normal();
  // z_par_zi ~ std_normal();

  // phi_inv_raw ~ uniform(0, 1);
  sigma_ye_raw ~ uniform(0, 1);
  sigma_ll_raw ~ uniform(0, 1);
  sigma_id_raw ~ uniform(0, 1);
  sigma_par_raw ~ uniform(0, 1);
  // sigma_ye_zi_raw ~ uniform(0, 1);
  // sigma_ll_zi_raw ~ uniform(0, 1);
  // sigma_id_zi_raw ~ uniform(0, 1);
  // sigma_par_zi_raw ~ uniform(0, 1);

  z_alpha ~ std_normal();
  z_beta_bv ~ std_normal();
  z_beta_bv2 ~ std_normal();
  z_beta_age_q1 ~ std_normal();
  z_beta_age_q2 ~ std_normal();
  z_beta_f ~ std_normal();
  z_beta_co_n ~ std_normal();

  // z_alpha_zi ~ std_normal();
  // z_beta_zi_bv ~ std_normal();
  // z_beta_zi_bv2 ~ std_normal();
  // z_beta_zi_age_q1 ~ std_normal();
  // z_beta_zi_age_q2 ~ std_normal();
  // z_beta_zi_f ~ std_normal();

  // Likelihood
  for (i in 1:N) {
    if (Y[i] == 0) {
      target += log_sum_exp(log(theta[i]),
      log1m(theta[i])
      + poisson_log_glm_lpmf(Y[i] | to_matrix(x_mat[i]), rand_inter[i], beta_vec));
    } else {
      target += log1m(theta[i])
      + poisson_log_glm_lpmf(Y[i] | to_matrix(x_mat[i]), rand_inter[i], beta_vec);
    }
  }
}

generated quantities {
  // Back-transform regression parameters to original scale
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

  real alpha_zi = alpha_zi_std;
  // - beta_zi_bv_std * bv_std_const1
  // + beta_zi_bv2_std * bv_std_const2
  // - beta_zi_age_q1_std * age_q1_const
  // - beta_zi_age_q2_std * age_q2_const
  // - beta_zi_f_std * f_const;

  // real beta_zi_bv = beta_zi_bv_std / bv_sd_std - beta_zi_bv2_std * bv_std_const4;
  // real beta_zi_bv2 = beta_zi_bv2_std / bv_std_const3;
  // real beta_zi_age_q1 = beta_zi_age_q1_std / sd_age_q1;
  // real beta_zi_age_q2 = beta_zi_age_q2_std / sd_age_q2;
  // real beta_zi_f = beta_zi_f_std / sd_f;


  // Posterior predictions
  int y_rep[N]; // Posterior predictive samples
  int<lower=0,upper=1> zinf[N];
  vector[N] count_linpred = x_mat * beta_vec + rand_inter; // Linear predictor
  for (i in 1:N) {
    zinf[i] = bernoulli_rng(theta[i]);
    if (zinf[i] == 1) {
      y_rep[i] = 0;
    } else {
      y_rep[i] = poisson_log_rng(count_linpred[i]);
    }
  }
}
