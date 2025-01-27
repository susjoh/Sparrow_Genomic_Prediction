# data <- tar_read(stan_data_adult_ars_f)
data <- tar_read(stan_data_adult_surv_f)

pc_rate <- function(U, a = 0.05) {
  -log(a) / U
}

ars <- function(i, data, beta_sd, sigma_rate, sigma_err, phi_inv_rate) {
  alpha_std  <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd / 4)
  beta_age_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age2_std <- rnorm(1, mean = 0, sd = beta_sd / 4)
  beta_f_std <- rnorm(1, mean = 0, sd = beta_sd)

  phi_inv <- rexp(1, rate = phi_inv_rate)

  sigma_ye <- rexp(1, rate = sigma_rate)
  sigma_ll <- rexp(1, rate = sigma_rate)
  sigma_id <- rexp(1, rate = sigma_rate)

  sigma_error <- rgamma(n = 1, shape = data$gam_shape, rate = data$gam_rate)

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)

  bv_raw <-  rnorm(n = data$N - data$M, 0, 1)

  bv_lat <- data$bv_obs + sigma_error * bv_raw
  bv <- c(data$bv_exact, bv_lat)
  phi <- 1 / phi_inv

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]

  eta <- alpha_std +
    beta_bv_std * scale(bv) +
    beta_bv2_std * scale(bv)^2 +
    beta_age_std * scale(data$age) +
    beta_age2_std * scale(data$age)^2 +
    beta_f_std * scale(data$f) +
    ye_full +
    ll_full +
    id_full

  mu <- exp(eta)
  # print(summary(mu))
  y <- rnbinom(n = data$N, mu = mu, size = phi)

  c(summary(y),
    var_y = var(y),
    sd_y = sd(y),
    var_eta = var(eta),
    sd_eta = sd(eta),
    # alpha_std = alpha_std,
    # beta_bv_std = beta_bv_std,
    # beta_bv2_std = beta_bv2_std,
    # beta_age_std = beta_age_std,
    # beta_age2_std = beta_age2_std,
    # beta_f_std = beta_f_std,
    phi_inv = phi_inv,
    phi = phi,
    sigma_error = sigma_error,
    sigma_ye = sigma_ye,
    sigma_ll = sigma_ll,
    sigma_id = sigma_id)
}

summary(log(data$Y[data$Y > 0]))

pc_rate(U = 0.3) # -> rate 10 exps

var(log(data$Y[data$Y > 0])) # -> 0.05 betas

res <- lapply(1:1e3, ars, data = data,
              beta_sd = 5e-2,
              sigma_rate = pc_rate(U = 5e-2 * 6),
              phi_inv_rate = 10) %>%
  do.call(rbind, .) %>%
  round(., 2)

head(res, 20)

plot(density(res[, "Mean"]))
plot(density(res[, "Median"]))
plot(density(res[, "Max."]))

colMeans(res)

surv <- function(i, data, beta_sd, sigma_rate) {
  alpha_std  <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd / 4)
  beta_age_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age2_std <- rnorm(1, mean = 0, sd = beta_sd / 4)
  beta_f_std <- rnorm(1, mean = 0, sd = beta_sd)

  rate <- 1.29 / (1.29 - sqrt(3.5))^2
  shape <- rate * 1.29

  sigma_error <- rgamma(n = 1, shape = shape, rate = rate)

  sigma_ye <- rexp(1, rate = sigma_rate)
  sigma_ll <- rexp(1, rate = sigma_rate)
  sigma_id <- rexp(1, rate = sigma_rate)
  # sigma_res <- rexp(1, rate = sigma_rate)

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)
  # res <- rnorm(n = data$N, 0, sigma_res)

  bv_raw <-  rnorm(n = data$N - data$M, 0, 1)

  bv_lat <- data$bv_obs + sigma_error * bv_raw
  bv <- c(data$bv_exact, bv_lat)

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]

  eta <- alpha_std +
    beta_bv_std * scale(bv) +
    beta_bv2_std * scale(bv)^2 +
    beta_age_std * scale(data$age) +
    beta_age2_std * scale(data$age)^2 +
    beta_f_std * scale(data$f) +
    ye_full +
    ll_full +
    id_full #+
  # res

  p <- 1 / (1 + exp(-eta))
  # print(summary(p))
  y <- rbinom(n = data$N, p = p, size = 1)

  c(summary(y),
    var_y = var(y),
    sd_y = sd(y),
    var_eta = var(eta),
    sd_eta = sd(eta),
    alpha_std = alpha_std,
    beta_bv_std = beta_bv_std,
    beta_bv2_std = beta_bv2_std,
    beta_age_std = beta_age_std,
    beta_age2_std = beta_age2_std,
    beta_f_std = beta_f_std,
    sigma_error = sigma_error,
    sigma_ye = sigma_ye,
    sigma_ll = sigma_ll,
    sigma_id = sigma_id
    # sigma_res = sigma_res
  )
}

p_hat <- mean(data$Y)

logit_p_hat <- log(p_hat / (1 - p_hat))

var_p_hat <- 1 / (p_hat * (1 - p_hat)) #-> 0.8 betas (divide by 5 since there are 5)

# pc_rate(U = 1) # -> rate 10 exps
#
# var(log(data$Y[data$Y > 0])) # -> 0.05 betas
#
res <- lapply(1:1e3, surv, data = data,
              beta_sd = 0.8,
              sigma_rate = pc_rate(U = 0.8 * 6)) %>%
  do.call(rbind, .) %>%
  round(., 2)

# plot(density(res[, "Mean"]))
# plot(density(res[, "Median"]))
#
# plot(density(res[, ""]))
#
# colMeans(res)

ars_bvpost <- function(i, data, alpha_sd, alpha_mean, r, s, sigma_rate, phi_inv_rate) {

  beta_sd <- rgamma(1, shape = s, rate = r)
  alpha_std  <- rnorm(1, mean = alpha_mean, sd = alpha_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd / sqrt(2))
  beta_age_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age2_std <- rnorm(1, mean = 0, sd = beta_sd / sqrt(2))
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

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = bv_lat[data$id_idx]

  eta1 <- alpha_std +
    beta_bv_std * scale(bv_lat_full) +
    beta_bv2_std * scale(bv_lat_full)^2 +
    beta_age_std * scale(data$age) +
    beta_age2_std * scale(data$age)^2 +
    beta_f_std * scale(data$f)
  eta2 <- ye_full +
    ll_full +
    id_full

  eta <- eta1 + eta2

  mu <- exp(eta)
  # print(summary(mu))
  y <- rnbinom(n = data$N, mu = mu, size = phi)

  c(summary(y),
    sd_y = sd(y),
    mean_eta = mean(eta),
    mean_eta1 = mean(eta2),
    mean_eta2 = mean(eta1),
    sd_eta = sd(eta),
    sd_eta1 = sd(eta1),
    sd_eta2 = sd(eta2),
    prop = var(eta1) / var(eta),
    # alpha_std = alpha_std,
    # # beta_bv_std = beta_bv_std,
    # beta_bv2_std = beta_bv2_std,
    # beta_age_std = beta_age_std,
    # beta_age2_std = beta_age2_std,
    # beta_f_std = beta_f_std,
    phi_inv = phi_inv
    # phi = phi,
    # sigma_ye = sigma_ye,
    # sigma_ll = sigma_ll,
    # sigma_id = sigma_id
  )
}

data <- tar_read(stan_data_adult_ars_bvpost_f)

bv_lat_full <- rnorm(n = data$N_id,
                     mean = data$bv_mean,
                     sd = data$bv_sd)[data$id_idx]

a <- 2.5 # -> About 1/3 explained by quadratic effects
sqrt((scale(bv_lat_full)^4 / a^2 + scale(data$age)^4 / a^2) /
       (scale(bv_lat_full)^2 + scale(bv_lat_full)^4 / a^2 + scale(data$age)^2 + scale(data$age)^4 / a^2 + scale(data$f)^2)) %>%
  summary()


num <- sqrt(1 + scale(bv_lat_full)^2 + scale(bv_lat_full)^4 / 2^2 + scale(data$age)^2 + scale(data$age)^4 / 2^2 + scale(data$f)^2) %>%
  mean()

sr <- 25

res <- lapply(1:2e3, ars_bvpost, data = data,
              s = 2, r = 40, alpha_sd = 1, alpha_mean = log(mean(data$Y)),
              # Ensure the fixed effects explain about 2 times more sd than the random effects (given a =2)
              sigma_rate = pc_rate(U = 1e-1),
              # PC-prior on phi_inv with mean equal to observed overdispersion
              phi_inv_rate = mean(data$Y)^2 / (var(data$Y) - mean(data$Y))) %>%
  do.call(rbind, .)
summary(res)


surv_bvpost <- function(i, data, beta_sd, sigma_rate) {
  alpha_std  <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_bv2_std <- rnorm(1, mean = 0, sd = beta_sd / 2)
  beta_age_std <- rnorm(1, mean = 0, sd = beta_sd)
  beta_age2_std <- rnorm(1, mean = 0, sd = beta_sd / 2)
  beta_f_std <- rnorm(1, mean = 0, sd = beta_sd)

  sigma_ye <- rexp(1, rate = sigma_rate)
  sigma_ll <- rexp(1, rate = sigma_rate)
  sigma_id <- rexp(1, rate = sigma_rate)

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)


  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = bv_lat[data$id_idx]

  eta1 <- alpha_std +
    beta_bv_std * scale(bv_lat_full) +
    beta_bv2_std * scale(bv_lat_full)^2 +
    beta_age_std * scale(data$age) +
    beta_age2_std * scale(data$age)^2 +
    beta_f_std * scale(data$f)
  eta2 <- ye_full +
    ll_full +
    id_full

  eta <- eta1 + eta2

  p <- 1 / (1 + exp(-eta))
  # print(summary(p))
  y <- rbinom(n = data$N, p = p, size = 1)

  c(summary(y),
    sd_y = sd(y),
    mean_eta = mean(eta),
    mean_eta1 = mean(eta2),
    mean_eta2 = mean(eta1),
    sd_eta = sd(eta),
    sd_eta1 = sd(eta1),
    sd_eta2 = sd(eta2)
    # alpha_std = alpha_std,
    # # beta_bv_std = beta_bv_std,
    # beta_bv2_std = beta_bv2_std,
    # beta_age_std = beta_age_std,
    # beta_age2_std = beta_age2_std,
    # beta_f_std = beta_f_std,
    # sigma_ye = sigma_ye,
    # sigma_ll = sigma_ll,
    # sigma_id = sigma_id
  )
}

data <- tar_read(stan_data_adult_surv_bvpost_f)

sr <- 80

res <- lapply(1:2e3, surv_bvpost, data = data,
              # Ensure the fixed effects explain about 2 times more sd than the random effects (given a =2)
              sigma_rate = sr, # pc_rate(U = 1e-1),
              beta_sd = (2 * sqrt(3)) / ((sr) * num)) %>%
  do.call(rbind, .)

summary(res)


ars_bvpost_sim <- function(i, data) {

  alpha_std <- alpha <- -0.077466411
  beta_bv_std <- beta_bv <-  0.1
  beta_bv2_std <- beta_bv2 <- -0.5
  beta_age_std <- beta_age <- 0.182861276
  beta_age2_std <- beta_age2 <- -0.044908613
  beta_f_std <- beta_f <- -0.111374375

  sigma_ye <- 0.367186607
  sigma_ll <- 0.411729306
  sigma_id <- 0.599423570

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)

  phi <- 1.670098105

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = bv_lat[data$id_idx]

  eta1 <-
    # alpha_std +
    # beta_bv_std * scale(bv_lat_full) +
    # beta_bv2_std * scale(bv_lat_full)^2 +
    # beta_age_std * scale(data$age) +
    # beta_age2_std * scale(data$age)^2 +
    # beta_f_std * scale(data$f)
    alpha +
    beta_bv * bv_lat_full +
    beta_bv2 * (bv_lat_full)^2 +
    beta_age * data$age +
    beta_age2 * (data$age)^2 +
    beta_f * data$f
  eta2 <- ye_full +
    ll_full +
    id_full

  eta <- eta1 + eta2

  mu <- exp(eta)
  rnbinom(n = data$N, mu = mu, size = phi)
}

surv_bvpost_sim <- function(i, data) {

  alpha_std <- alpha <- -0.077466411
  beta_bv_std <- beta_bv <-  0.1
  beta_bv2_std <- beta_bv2 <- -0.5
  beta_age_std <- beta_age <- 0.182861276
  beta_age2_std <- beta_age2 <- -0.044908613
  beta_f_std <- beta_f <- -0.111374375

  sigma_ye <- 0.367186607
  sigma_ll <- 0.411729306
  sigma_id <- 0.599423570

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = bv_lat[data$id_idx]

  eta1 <-
    # alpha_std +
    # beta_bv_std * scale(bv_lat_full) +
    # beta_bv2_std * scale(bv_lat_full)^2 +
    # beta_age_std * scale(data$age) +
    # beta_age2_std * scale(data$age)^2 +
    # beta_f_std * scale(data$f)
    alpha +
    beta_bv * bv_lat_full +
    beta_bv2 * (bv_lat_full)^2 +
    beta_age * data$age +
    beta_age2 * (data$age)^2 +
    beta_f * data$f
  eta2 <- ye_full +
    ll_full +
    id_full

  eta <- eta1 + eta2

  p <- 1 / (1 + exp(-eta))
  # print(summary(p))
  rbinom(n = data$N, p = p, size = 1)
}
