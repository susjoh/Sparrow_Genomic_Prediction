library(rstan)
# library(INLA)
library(ggplot2)
library(qs)
library(magrittr)

n = 1e4
beta = 6
beta2 <- 4
prec.u = 1 / 0.4^2
prec.y = 1 / 0.3^2
prior.beta = c(0, 0.0001)
prior.prec.u = c(10, 9/prec.u)
prior.prec.y = c(10, 9/prec.y)
## heteroscedastic scaling
s = c(rep(1e8, 200), rep(1, 800))

# formula1 = y ~ f(mix, model="meb", scale=s, values=mix,
#                  hyper = list(
#                    beta = list(
#                      prior = "gaussian",
#                      param = prior.beta,
#                      fixed = FALSE
#                    ),
#                    prec.u = list(
#                      prior = "pc.prec",
#                      param = c(1, 0.05),
#                      initial = log(prec.u),
#                      fixed = TRUE
#                    )
#                  )
# ) + f(mix2, model="meb", scale=s, values=mix2,
#       hyper = list(
#         beta = list(
#           prior = "gaussian",
#           param = prior.beta,
#           fixed = FALSE
#         ),
#         prec.u = list(
#           prior = "pc.prec",
#           param = c(1, 0.05),
#           initial = log(prec.u),
#           fixed = FALSE)))
#
# formula2 = y ~ f(mix, model="meb", scale=s, values=mix,
#                  hyper = list(
#                    beta = list(
#                      prior = "gaussian",
#                      param = prior.beta,
#                      fixed = FALSE
#                    ),
#                    prec.u = list(
#                      prior = "pc.prec",
#                      param = c(1, 0.05),
#                      initial = log(prec.u),
#                      fixed = TRUE)))
#
# formula3 <- y ~ mix + mix2

set.seed(1)
# Simulate data
w = rnorm(n)
u <- rnorm(n, sd = 1/sqrt(prec.u))
## true but unobserved covariate
x <- w + u
m <- 2661
mix = c(x[1:m], w[(m + 1):n])
age <- sample(1:5, size = n, replace = TRUE) %>% scale()
y = 1 + beta*x + beta2*x^2 + rnorm(n, sd = 1/sqrt(prec.y))
f = rnorm(n, mean = 1, sd = 0.04)
n_ye <- 22
ye <- rnorm(n_ye, sd = 0.3)
ye_idx <- sample(size = n, x = n_ye, replace = TRUE)
n_ll <- 12
ll <- rnorm(n_ll, sd = 0.7)
ll_idx <- sample(size = n, x = n_ll, replace = TRUE)
n_id <- 5268
id <- rnorm(n_id, sd = 0.2)
id_idx <- c(1:n_id, sample(size = n-n_id, x = n_id, replace = TRUE))
mu <- exp(0.1 + 0.2*x + 0.1*x^2 + 0.1 * age + 0.15 * age^2 + f * 0.5 +
            ye[ye_idx] +
            ll[ll_idx] +
            id[id_idx])

y_neg_bin <- rnbinom(n = n, mu = mu, size = 1.3) # size is phi when mu is set

# data = data.frame(y,
#                   mix,
#                   mix2 = mix^2,
#                   s,
#                   id = 1:n,
#                   w,
#                   u,
#                   x,
#                   x2 = x^2,
#                   x_sqrt = sqrt(x),
#                   age,
#                   age2 = age^2,
#                   mu,
#                   f,
#                   y_neg_bin)

# r1 = inla(formula1, data = data,
#           family = "gaussian",
#           control.family = list(
#             hyper = list(
#               prec = list(param = prior.prec.y,
#                           fixed = FALSE))))
#
# r3 = inla(formula3, data = data,
#           family = "gaussian",
#           control.family = list(
#             hyper = list(
#               prec = list(param = prior.prec.y,
#                           fixed = FALSE))))
#
# data.frame(r1b1 = r1$summary.hyperpar["MEB beta for mix", "0.025quant"] < beta &
#              beta < r1$summary.hyperpar["MEB beta for mix", "0.975quant"],
#            r1b2 = r1$summary.hyperpar["MEB beta for mix2", "0.025quant"] < beta2 &
#              beta < r1$summary.hyperpar["MEB beta for mix2", "0.975quant"],
#            r3b1 = r3$summary.fixed["mix", "0.025quant"] < beta &
#              beta < r3$summary.fixed["mix", "0.975quant"],
#            r3b2 = r3$summary.fixed["mix2", "0.025quant"] < beta2 &
#              beta < r3$summary.fixed["mix2", "0.975quant"])

options(mc.cores = 8)
rstan_options(auto_write = TRUE)

# stan_dat <- list(N = n,
#                  M = m,
#                  X_exact = x[1:m],
#                  X_obs = w[(m + 1):n],
#                  Y = y,
#                  pc_par = -log(0.05) / (sd(y)))
#
# fit <- stan(file = "r/stan_meb_sq.stan",
#             data = stan_dat,
#             iter = 1e7,
#             warmup = 1e4,
#             chains = 8,
#             thin = 100,
#             verbose = TRUE)
#
# qsave(fit, file = "stan_result")

stan_dat_nb <- list(N = n,
                    M = m,
                    N_ll = n_ll,
                    N_ye = n_ye,
                    N_id = n_id,
                    ye_idx = ye_idx,
                    ll_idx = ll_idx,
                    id_idx = id_idx,
                    X_exact = x[1:m],
                    X_obs = w[(m + 1):n],
                    Y = y_neg_bin,
                    age = age[,1],
                    f = f)

fit_nb <- stan(file = "r/adult_ars.stan",
               data = stan_dat_nb,
               iter = 1.1e3,
               warmup = 1e2,
               chains = 8,
               cores = 8,
               control = list(stepsize = 0.95),
               model_name = "asdf",
               thin = 10)

qsave(fit_nb, file = "stan_res_nb_random")
