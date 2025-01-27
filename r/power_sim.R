library(magrittr)
library(dplyr)
library(ggplot2)

classical <- function(i,
                      cor_x,
                      x_var = 1,
                      n_perf = 1354, # number phenogone females
                      n_noise = 10000,
                      effect_size = 0.01,
                      intercept = 0) {
  # Generate true covariate
  x <- rnorm(n_perf + n_noise, sd = sqrt(x_var))

  # Generate response
  y <- intercept +  x * effect_size + rnorm(n_perf + n_noise, sd = 1)

  df <- data.frame(x, y)

  # Model using all data
  mod <- lm(y ~ x, data = df)

  # Model using just data perfectly observed
  mod_perf <- lm(y ~ x, data = df[1:n_perf, ])

  # Covariate observed with noise
  df$x_noise <- x +
    rnorm(n_perf + n_noise, sd = sqrt(x_var * (1 / cor_x^2 - 1)))
  # cor(x_noise, x) # check correctness of corr

  # Model using just noisy obs
  mod_noise <- lm(y ~ x_noise, data = df[-(1:n_perf), ])

  df$x_mixed <- df$x_noise
  df$x_mixed[1:n_perf] <- df$x[1:n_perf]

  mod_mixed <- lm(y ~ x_mixed, data = df)

  # Corrected estimate
  data.frame(
    i = i,
    cor_x = cor_x,
    effect_size = effect_size,
    coef_perf = coef(mod_perf)[2],
    coef_mixed = coef(mod_mixed)[2],
    coef_mixed_corrected = coef(mod_mixed)[2] / cor_x^2,
    coef_noise = coef(mod_noise)[2],
    coef_noise_corrected = coef(mod_noise)[2] / cor_x^2,
    sign_perf = summary(mod_perf)$coefficients["x", "Pr(>|t|)"],
    sign_mixed = summary(mod_mixed)$coefficients["x_mixed", "Pr(>|t|)"],
    sign_noise = summary(mod_noise)$coefficients["x_noise", "Pr(>|t|)"]#,
    # sign_perf = summary(mod_perf)$coefficients["x", "Pr(>|t|)"]
  )
}

berkson <- function(i,
                    cor_x,
                    x_var = 1,
                    n_perf = 1354, # number phenogone females
                    n_noise = 10000,
                    effect_size = 0.3,
                    intercept = 0) {
  # Generate measurements
  w <- rnorm(n_perf + n_noise, sd = sqrt(x_var / (1 / cor_x^2 - 1)))
  # True covariate
  x <- w + rnorm(n_perf + n_noise, sd = sqrt(x_var))

  # cor(x, w) # OK

  # Generate response
  y <- intercept +  x * effect_size + rnorm(n_perf + n_noise, sd = 1)

  df <- data.frame(x, y, w)

  # Model using true covariate
  mod <- lm(y ~ x, data = df)

  # Model using just data perfectly observed
  mod_perf <- lm(y ~ x, data = df[1:n_perf, ])

  # Model using just observation
  mod_obs <- lm(y ~ w, data = df)

  df$x_mixed <- df$w
  df$x_mixed[1:n_perf] <- df$x[1:n_perf]

  # Model using just observation
  mod_mixed <- lm(y ~ x_mixed, data = df)

  # Corrected estimate
  data.frame(
    i = i,
    cor_x = cor_x,
    effect_size = effect_size,
    coef_perf = coef(mod_perf)[2],
    coef_obs = coef(mod_obs)[2],
    coef_mixed = coef(mod_mixed)[2],
    sign_perf = summary(mod_perf)$coefficients["x", "Pr(>|t|)"],
    sign_obs = summary(mod_obs)$coefficients["w", "Pr(>|t|)"],
    sign_mixed = summary(mod_mixed)$coefficients["x_mixed", "Pr(>|t|)"]
  )
}

combs <- expand.grid(i = 1:100,
                     cor_x = seq(0.1, 0.6, length = 6),
                     effect_size = c(0.01, 0.025, 0.05, 0.1,  0.3))

res_class <- mapply(FUN = classical,
                    i = combs$i,
                    cor_x = combs$cor_x,
                    effect_size = combs$effect_size,
                    SIMPLIFY = FALSE) %>%
  do.call(what = rbind, .)

res_berk <- mapply(FUN = berkson,
                    i = combs$i,
                    cor_x = combs$cor_x,
                    effect_size = combs$effect_size,
                    SIMPLIFY = FALSE) %>%
  do.call(what = rbind, .)

res_avg_class <- res_class %>%
  group_by(cor_x, effect_size) %>%
  summarise(coef_perf_mean = mean(coef_perf),
            coef_mixed_mean = mean(coef_mixed),
            coef_noise_mean = mean(coef_noise),
            coef_noise_corrected_mean = mean(coef_noise_corrected),
            coef_perf_sd = sd(coef_perf),
            coef_mixed_sd = sd(coef_mixed),
            coef_noise_sd = sd(coef_noise),
            coef_noise_corrected_sd = sd(coef_noise_corrected),
            sign_prop_mixed = mean(sign_mixed < 0.05),
            sign_prop_noise = mean(sign_noise < 0.05),
            sign_prop_perf = mean(sign_perf < 0.05),
            n = n())

res_avg_berk <- res_berk %>%
  group_by(cor_x, effect_size) %>%
  summarise(coef_perf_mean = mean(coef_perf),
            coef_mixed_mean = mean(coef_mixed),
            coef_obs_mean = mean(coef_obs),
            coef_perf_sd = sd(coef_perf),
            coef_mixed_sd = sd(coef_mixed),
            coef_obs_sd = sd(coef_obs),
            sign_prop_mixed = mean(sign_mixed < 0.05),
            sign_prop_obs = mean(sign_obs < 0.05),
            sign_prop_perf = mean(sign_perf < 0.05),
            n = n())

summary(res_class)

ggplot(data = res_class[res_class$effect_size == 0.1, ],
       aes(group = cor_x, x = cor_x, y = coef_mixed)) +
  geom_boxplot()

ggplot(data = res_class[res_class$effect_size == 0.1, ],
       aes(group = cor_x, x = cor_x, y = coef_mixed_corrected)) +
  geom_boxplot()

ggplot(data = res_class[res_class$effect_size == 0.1, ],
       aes(group = cor_x, x = cor_x, y = coef_noise_corrected)) +
  geom_boxplot()

ggplot(data = res_class[res_class$effect_size == 0.1, ],
       aes( y = coef_perf)) +
  geom_boxplot()

ggplot(data = res_class, aes(x = cor_x, group = cor_x, y = sign_mixed)) +
  geom_boxplot() +
  facet_wrap(vars(effect_size)) +
  geom_smooth(method = "lm")

ggplot(data = res_avg_class, aes(x = cor_x, y = sign_prop_mixed)) +
  geom_point() +
  facet_wrap(vars(effect_size)) +
  geom_smooth(method = "lm")

ggplot(data = res_avg_class, aes(x = cor_x, y = sign_prop_perf)) +
  geom_point() +
  facet_wrap(vars(effect_size)) +
  geom_smooth(method = "lm")

ggplot(data = res_avg_class, aes()) +
  facet_grid(cor_x ~ effect_size) +
  geom_text(aes(label = sign_prop_noise, x = 1, y = 1))

ggplot(data = res_avg_class, aes()) +
  facet_grid(cor_x ~ effect_size) +
  geom_text(aes(label = sign_prop_mixed, x = 1, y = 1))

ggplot(data = res_avg_class, aes()) +
  facet_grid(cor_x ~ effect_size) +
  geom_text(aes(label = sign_prop_perf, x = 1, y = 1))

# Improvement
ggplot(data = res_avg_class) +
  facet_grid(rows = vars(cor_x),
             cols = vars(effect_size),
             switch = "both") +
  geom_rect(aes(xmin = 0,
                ymin = 0,
                xmax = 1,
                ymax = 1,
                fill = sign_prop_mixed - sign_prop_perf)) +
  coord_equal() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    limits = c(-1, 1),
    name  = "") +
  scale_y_continuous(c(0, 1), expand = c(0, 0), name = "Accuracy in x") +
  scale_x_continuous(c(0, 1), expand = c(0, 0), name = "Effect size") +
  geom_text(aes(label = round(sign_prop_mixed - sign_prop_perf, 2),
                x = 0.5,
                y = 0.5),
            size = 8) +
  theme(panel.border = element_rect(fill = NA),
        axis.text = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.ticks = element_blank(),
        strip.text.y = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.title = element_text(size = 18)) +
  labs(title = "Improvement in effect discovery rate",
       subtitle = "when including noisy covariate measurements")
