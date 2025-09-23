sca <- 2
library(targets)
library(ggplot2)
library(ggpubr)
source("r/crossover_gp_inla_func.R")

ars_bv_pred_plot_poster_f=
  plot_lines_posterior(df = tar_read(ars_covmat_bv_preds_and_marg_f)$df_pred,
                       xlab = "",
                       ylab = "",
                       title = "",
                       lw = sca, leg.x = 0.75,
                       bs = 11 * sca)

ars_bv_pred_plot_poster_m=
  plot_lines_posterior(df = tar_read(ars_covmat_bv_preds_and_marg_m)$df_pred,
                       xlab = "",
                       ylab = "",
                       title = "",
                       lw = sca, leg.x = 0.75,
                       bs = 11 * sca)

surv_bv_pred_plot_poster_f=
  plot_lines_posterior(df = tar_read(surv_covmat_bv_preds_and_marg_f)$df_pred,
                       xlab = "",
                       ylab = "",
                       title = "",
                       lw = sca, leg.x = 0.75,
                       bs = 11 * sca)

surv_bv_pred_plot_poster_m=
  plot_lines_posterior(df = tar_read(surv_covmat_bv_preds_and_marg_m)$df_pred,
                       xlab = "",
                       ylab = "",
                       title = "",
                       lw = sca, leg.x = 0.75,
                       bs = 11 * sca)

nest_bv_pred_plot_poster_f=
  plot_lines_posterior(df = tar_read(nest_bv_preds_and_marg_f)$df_pred,
                       xlab = "",
                       ylab = "",
                       title = "",
                       lw = sca, leg.x = 0.75,
                       bs = 11 * sca)

nest_bv_pred_plot_poster_m=
  plot_lines_posterior(df = tar_read(nest_bv_preds_and_marg_m)$df_pred,
                       xlab = "",
                       ylab = "",
                       title = "",
                       lw = sca, leg.x = 0.75,
                       bs = 11 * sca)

plt <- ggarrange(ars_bv_pred_plot_poster_f ,
          ars_bv_pred_plot_poster_m,
          surv_bv_pred_plot_poster_f,
          surv_bv_pred_plot_poster_m,
          nest_bv_pred_plot_poster_f,
          nest_bv_pred_plot_poster_m,
          common.legend = FALSE,
          ncol = 2,
          nrow = 3)

ggsave_path(paste0("fasdf.pdf"),
            plot = plt,
            width = 14,
            height = 21,
            device = "pdf")

ggsave_path(paste0("fasdf.png"),
            plot = plt,
            width = 14,
            height = 21,
            device = "png")
