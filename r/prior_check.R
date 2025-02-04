ars <- function(i, data, alpha_mean, beta_sd, sigma_rate, phi_inv_rate) {

  alpha_std  <- rnorm(1, mean = alpha_mean, sd = beta_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age2_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_f_std <- rnorm(1, mean = 0, sd = beta_sd)

  phi_inv <- rexp(1, rate = phi_inv_rate)

  sigma_ye <- rexp(1, rate = sigma_rate)
  sigma_ll <- rexp(1, rate = sigma_rate)
  sigma_id <- rexp(1, rate = sigma_rate)

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)

  phi <- 1 / phi_inv

  age_q <- poly(data$age, degree = 2)

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = (bv_lat[data$id_idx] - data$bv_mean_std) / data$bv_sd_std

  eta1 <- alpha_std +
    beta_bv_std * bv_lat_full +
    beta_bv2_std * bv_lat_full^2 +
    beta_age_std * scale(age_q[, 1]) +
    beta_age2_std * scale(age_q[, 2]) +
    beta_f_std * scale(data$f)
  eta2 <- ye_full +
    ll_full +
    id_full

  eta <- eta1 + eta2

  mu <- exp(eta)
  y <- rnbinom(n = data$N, mu = mu, size = phi)

  c(summary(y),
    sd_y = sd(y),
    prop = var(eta1) / var(eta),
    alpha_std = alpha_std,
    beta_bv_std = beta_bv_std,
    beta_bv2_std = beta_bv2_std,
    beta_age_std = beta_age_std,
    beta_age2_std = beta_age2_std,
    beta_f_std = beta_f_std,
    phi_inv = phi_inv,
    phi = phi,
    sigma_ye = sigma_ye,
    sigma_ll = sigma_ll,
    sigma_id = sigma_id
  )
}

data <- tar_read(stan_data_adult_f)

res <- lapply(1:2e3,
              ars,
              data = data,
              alpha_mean = data$sum_recruit_log_mean,
              sigma_rate = 10,
              beta_sd = 0.05, #sqrt(2 / 10^2),
              phi_inv_rate = mean(data$sum_recruit)^2 /
                (var(data$sum_recruit) - mean(data$sum_recruit))) %>%
  do.call(rbind, .);summary(res)
plot(density(res[, "Mean"]))
plot(density(res[, "Max."]))

surv <- function(i, data, alpha_mean, beta_sd, sigma_rate) {

  alpha_std  <- rnorm(1, mean = alpha_mean, sd = beta_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age2_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_f_std <- rnorm(1, mean = 0, sd = beta_sd)

  sigma_ye <- rexp(1, rate = sigma_rate)
  sigma_ll <- rexp(1, rate = sigma_rate)
  sigma_id <- rexp(1, rate = sigma_rate)
  sigma_res <- rexp(1, rate = sigma_rate)

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)
  res <-rnorm(n = data$N, 0, sigma_res)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = (bv_lat[data$id_idx] - data$bv_mean_std) / data$bv_sd_std

  age_q <- poly(data$age, degree = 2)

  eta1 <- alpha_std +
    beta_bv_std * bv_lat_full +
    beta_bv2_std * bv_lat_full^2 +
    beta_age_std * scale(age_q[, 1]) +
    beta_age2_std * scale(age_q[, 2]) +
    beta_f_std * scale(data$f)
  eta2 <- ye_full +
    ll_full +
    id_full +
    res

  eta <- eta1 + eta2

  p <- 1 / (1 + exp(-eta))
  y <- rbinom(n = data$N, p = p, size = 1)

  c(summary(y),
    # summary(p),
    sd_y = sd(y),
    prop = var(eta1) / var(eta))
}

data <- tar_read(stan_data_adult_f)

res <- lapply(1:2e3,
              surv,
              data = data,
              alpha_mean = data$survival_logit_mean,
              sigma_rate = 40,
              beta_sd = sqrt(2 / 40^2)
              ) %>%
  do.call(rbind, .);summary(res)
