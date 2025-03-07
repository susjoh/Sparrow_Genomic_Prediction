ars <- function(i,
                data,
                alpha_mean,
                beta_sd,
                sigma_rate,
                alpha_zi_mean,
                beta_zi_sd,
                sigma_rate_zi# ,phi_inv_rate,
                ) {

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)
  bv_lat_full = (bv_lat[data$id_idx] - data$bv_mean_std) / data$bv_sd_std

  # Linear predictor in count component
  alpha_std  <- rnorm(1, mean = alpha_mean, sd = beta_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age2_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_f_std <- rnorm(1, mean = 0, sd = beta_sd)

  sigma_ye <- rexp(1, rate = sigma_rate)
  sigma_ll <- rexp(1, rate = sigma_rate)
  sigma_id <- rexp(1, rate = sigma_rate)

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]

  eta1 <- alpha_std +
    beta_bv_std * bv_lat_full +
    beta_bv2_std * bv_lat_full^2 +
    beta_age_std * scale(data$age_q1) +
    beta_age2_std * scale(data$age_q2) +
    beta_f_std * scale(data$f)
  eta2 <- ye_full +
    ll_full +
    id_full

  eta <- eta1 + eta2

  lambda <- exp(eta)

  # Zero-inflation linear predictor
  alpha_zi_std  <- rnorm(1, mean = alpha_zi_mean, sd = beta_zi_sd)
  beta_zi_bv_std <- rnorm(1, mean = 0, sd = beta_zi_sd)
  beta_zi_bv2_std <- rnorm(1, mean = 0, sd = beta_zi_sd)
  beta_zi_age_std <- rnorm(1, mean = 0, sd = beta_zi_sd)
  beta_zi_age2_std <- rnorm(1, mean = 0, sd = beta_zi_sd)
  beta_zi_f_std <- rnorm(1, mean = 0, sd = beta_zi_sd)

  sigma_ye_zi <- rexp(1, rate = sigma_rate_zi)
  sigma_ll_zi <- rexp(1, rate = sigma_rate_zi)
  sigma_id_zi <- rexp(1, rate = sigma_rate_zi)

  ye_zi <- rnorm(n = data$N_ye, 0, sigma_ye_zi)
  ll_zi <- rnorm(n = data$N_ll, 0, sigma_ll_zi)
  id_zi <- rnorm(n = data$N_id, 0, sigma_id_zi)

  ye_zi_full = ye_zi[data$ye_idx]
  ll_zi_full = ll_zi[data$ll_idx]
  id_zi_full = id_zi[data$id_idx]

  eta1_zi <- alpha_zi_std +
    beta_zi_bv_std * bv_lat_full +
    beta_zi_bv2_std * bv_lat_full^2 +
    beta_zi_age_std * scale(data$age_q1) +
    beta_zi_age2_std * scale(data$age_q2) +
    beta_zi_f_std * scale(data$f)
  eta2_zi <- ye_zi_full +
    ll_zi_full +
    id_zi_full

  eta_zi <- eta1_zi + eta2_zi

  prob_zi <- 1 / (1 + exp(-eta_zi))

  # overdisperion
  # phi_inv <- rexp(1, rate = phi_inv_rate)
  # phi <- 1 / phi_inv
  zinf <- rbinom(n = data$N, prob = prob_zi, size = 1)
  y <- zinf * rpois(n = data$N, lambda = lambda)#, size = phi)

  # c(summary(y),
  #   sd_y = sd(y),
  #   prop = var(eta1) / var(eta),
  #   alpha_std = alpha_std,
  #   beta_bv_std = beta_bv_std,
  #   beta_bv2_std = beta_bv2_std,
  #   beta_age_std = beta_age_std,
  #   beta_age2_std = beta_age2_std,
  #   beta_f_std = beta_f_std,
  #   # phi_inv = phi_inv,
  #   # phi = phi,
  #   sigma_ye = sigma_ye,
  #   sigma_ll = sigma_ll,
  #   sigma_id = sigma_id
  # )
  y
}

data <- tar_read(stan_data_adult_ss_f)

res <- lapply(X = 1:2e3,
              FUN = ars,
              data = data,
              alpha_mean = log(mean(c(data$sum_recruit[data$sum_recruit != 0],
                                      data$sum_recruit))),
              beta_sd = 0.1,
              sigma_rate = 1 / 0.1,
              alpha_zi_mean = log(mean(data$sum_recruit == 0) /
                                    (1 - mean(data$sum_recruit == 0))),
              beta_zi_sd = 0.5,
              sigma_rate_zi = 1 / 0.5
              # phi_inv_rate = mean(data$sum_recruit)^2 /
                # (var(data$sum_recruit) - mean(data$sum_recruit))
              ) %>%
  do.call(rbind, .)
summary(res)
plot(density(res[, "Mean"]))
plot(density(res[, "Max."]))

surv <- function(i, data, alpha_mean, beta_sd, sigma_rate) {

  alpha_std  <- rnorm(1, mean = alpha_mean, sd = beta_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd / sqrt(2))
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
  res <- rnorm(n = data$N, 0, sigma_res)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = (bv_lat[data$id_idx] - data$bv_mean_std) / data$bv_sd_std

  eta1 <- alpha_std +
    beta_bv_std * bv_lat_full +
    beta_bv2_std * bv_lat_full^2 +
    beta_age_std * as.numeric(scale(data$age_q1)) +
    beta_age2_std * as.numeric(scale(data$age_q2)) +
    beta_f_std * as.numeric(scale(data$f))
  eta2 <- ye_full +
    ll_full +
    id_full +
    res

  eta <- eta1 + eta2

  p <- 1 / (1 + exp(-eta))
  y <- rbinom(n = data$N, p = p, size = 1)

  c(# summary(y),
    # summary(p),
    summary(eta),
    1.96 * sd(eta),
    # sd_y = sd(y),
    prop = var(eta1) / var(eta))
}

data <- tar_read(stan_data_adult_f)

res <- lapply(1:2e3,
              surv,
              data = data,
              alpha_mean = data$survival_logit_mean,
              sigma_rate = sqrt(1 / 0.5^2),
              beta_sd = 0.5) %>%
  do.call(rbind, .);summary(res)
