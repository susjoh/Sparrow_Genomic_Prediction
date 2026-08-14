library(ggplot2)
ms <- 5e2
n_repl <- 1e4

foo <- function(i, max_samp) {
  samp <- sample(20, size = max_samp, replace = TRUE)
  samp_1 <- (samp == 1)
  cs_s1 <- cumsum(samp_1)
  prop_1 <- cs_s1 / seq_len(max_samp)
  prop_1
}

res <- sapply(seq_len(n_repl), foo, max_samp = ms)

res_df <- as.data.frame(res)
res_df$mean <- apply(res, 1, mean)
res_df$q_025 <- apply(res, 1, function(x) quantile(x, probs = 0.025))
res_df$q_975 <- apply(res, 1, function(x) quantile(x, probs = 0.975))

res_df$q_025[250]
res_df$q_975[250]

res_df_l <- tidyr::pivot_longer(res_df, cols = starts_with("V"))
res_df_l$x <- rep(seq_len(ms), each = n_repl)

ggplot(dat = res_df_l) +
  geom_line(aes(x = x, y = value, group = factor(name)), alpha = 0.01) +
  geom_line(aes(x = x, y = mean), color = "red") +
  geom_line(aes(x = x, y = q_025), color = "blue") +
  geom_line(aes(x = x, y = q_975), color = "blue")
