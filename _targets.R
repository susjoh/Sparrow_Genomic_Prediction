# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(crew)
library(crew.cluster)
library(tibble)

controller_local <- crew_controller_local(name = "my_local_controller",
                                          workers = 6)

controller_slurm <- crew_controller_slurm(
  name = "my_slurm_controller",
  workers = 50,
  seconds_idle = 120,
  tasks_max = 8,
  options_cluster = crew_options_slurm(
    script_lines = paste("#SBATCH --account=share-nv-bio \n module load",
                         "R/4.4.2-gfbf-2024a R-bundle-CRAN/2024.11-foss-2024a",
                         "CMake/3.29.3-GCCcore-13.3.0"),
    log_output = "Jobs/%A.log",
    memory_gigabytes_per_cpu = 6,
    cpus_per_task = 16,
    time_minutes = 60 * 24 * 20, # minutes * hours * days
    partition = "CPUQ",
    verbose = TRUE)
)

# Set target options:
tar_option_set(
  # Packages that your targets need for their tasks.
  packages = c("tibble",
               "data.table",
               "stats",
               "parallel",
               "ggplot2",
               "GGally",
               "dplyr",
               "qs2",
               "magrittr",
               "tools",
               "INLA",
               "rstan",
               "bayesplot",
               "tidyr",
               "mvtnorm",
               "ggpubr",
               "MCMCglmm",
               "HDInterval",
               "rstan",
               "bayesplot"),
  format = "qs", # Optionally set the default storage format. qs is fast.
  error = "continue", # produce result even if the target errored
  resources = tar_resources(
    crew = tar_resources_crew(controller = "my_slurm_controller")),
  memory = "transient",
  garbage_collection = TRUE,
  iteration = "list",
  retrieval = "worker",
  storage = "worker",
  controller = controller_slurm
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source(files = "r/crossover_gp_inla_func.R")

# Filepaths that might need to be changed

plink_path <- "PLINK/plink_linux" # path to plink program

values_fitmod <- tibble(
  mod = c("ars_adult",
          "surv_adult",
          "ars_parent",
          "surv_parent",
          "nest"),
  gp_data_func = rlang::syms(c("make_adult_gp_data",
                               "make_adult_gp_data",
                               "make_parent_gp_data",
                               "make_parent_gp_data",
                               "make_nestling_gp_data")),
  fitness_data_path = rlang::syms(c("lrs_data_path",
                                    "lrs_data_path",
                                    "lrs_data_path",
                                    "lrs_data_path",
                                    "nestling_data_path")),
  fitdat_func = rlang::syms(c("make_data_adult",
                              "make_data_adult",
                              "make_data_parent",
                              "make_data_parent",
                              "make_data_nest")),
  trait = c("annual reproductive success",
            "annual survival",
            "annual reproductive success",
            "annual survival",
            "nestling survival"),
  xlab_start = c("B", "B", "Parental b", "Parental b", "Parental b"),
  make_sim = rlang::syms(c("make_sim_ars_adult",
                           "make_sim_surv_adult",
                           "make_sim_ars_parent",
                           "make_sim_surv_parent",
                           "make_sim_nest")),
  sim_par_vec = list(list("alpha" = 0.26,
                          "beta_bv" = -0.16,
                          "beta_bv2" = -1.2,
                          "beta_age_q1" = 9.4,
                          "beta_age_q2" = -14,
                          "beta_f" = -2.3,
                          "alpha_zi" = -0.80,
                          "sigma_ye" = 0.28,
                          "sigma_ll" = 0.18,
                          "sigma_id" = 0.30,
                          "sigma_res" = 0),
                     list("alpha" = 0.63,
                          "beta_bv" = -0.067,
                          "beta_bv2" = -1.3,
                          "beta_age_q1" = 24,
                          "beta_age_q2" = -2.6,
                          "beta_f" = -0.90,
                          "sigma_ye" = 0.33,
                          "sigma_ll" = 0.19,
                          "sigma_id" = 0.16,
                          "sigma_res" = 0),
                     list("alpha" = 0.081,
                          "beta_bv" = -0.037,
                          "beta_bv2" = -0.51,
                          "beta_age_q1" = 5.4,
                          "beta_age_q2" = -11,
                          "beta_f" = 0.66,
                          "alpha_zi" = -0.70,
                          "sigma_ye" = 0.25,
                          "sigma_ll" = 0.086,
                          "sigma_id" = 0.42,
                          "sigma_par" = 0.13,
                          "sigma_res" = 0),
                     list("alpha" = 0.63,
                          "beta_bv" = -0.067,
                          "beta_bv2" = -1.3,
                          "beta_age_q1" = -24,
                          "beta_age_q2" = 2.6,
                          "beta_f" = -0.90,
                          "sigma_ye" = 0.33,
                          "sigma_ll" = 0.19,
                          "sigma_id" = 0.16,
                          "sigma_par" = 0.084,
                          "sigma_res" = 0),
                     list("alpha" = -1.3,
                          "beta_bv" = -0.39,
                          "beta_bv2" = -0.29,
                          "beta_f" = 1.1,
                          "sigma_hy" = 0.44,
                          "sigma_hi" = 0.36,
                          "sigma_par" = 0.40,
                          "sigma_res" = 0)),
  stan_data_func = rlang::syms(c("make_stan_data_adult",
                                 "make_stan_data_adult",
                                 "make_stan_data_parent",
                                 "make_stan_data_parent",
                                 "make_stan_data_nest")),
  stan_pars = list(c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_age_q1",
                     "beta_age_q2",
                     "beta_f",
                     "alpha_std",
                     "beta_bv_std",
                     "beta_bv2_std",
                     "beta_age_q1_std",
                     "beta_age_q2_std",
                     "beta_f_std",
                     "ye",
                     "ll",
                     "id",
                     "bv_lat",
                     "sigma_ll",
                     "sigma_ye",
                     "sigma_id",
                     "alpha_zi",
                     "alpha_zi_std",
                     "theta",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_age_q1",
                     "beta_age_q2",
                     "beta_f",
                     "alpha_std",
                     "beta_bv_std",
                     "beta_bv2_std",
                     "beta_age_q1_std",
                     "beta_age_q2_std",
                     "beta_f_std",
                     "ye",
                     "ll",
                     "id",
                     "bv_lat",
                     "sigma_ll",
                     "sigma_ye",
                     "sigma_id",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_age_q1",
                     "beta_age_q2",
                     "beta_f",
                     "alpha_std",
                     "beta_bv_std",
                     "beta_bv2_std",
                     "beta_age_q1_std",
                     "beta_age_q2_std",
                     "beta_f_std",
                     "ye",
                     "ll",
                     "id",
                     "par",
                     "bv_lat",
                     "sigma_ll",
                     "sigma_ye",
                     "sigma_id",
                     "sigma_par",
                     "alpha_zi",
                     "alpha_zi_std",
                     "theta",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_age_q1",
                     "beta_age_q2",
                     "beta_f",
                     "alpha_std",
                     "beta_bv_std",
                     "beta_bv2_std",
                     "beta_age_q1_std",
                     "beta_age_q2_std",
                     "beta_f_std",
                     "ye",
                     "ll",
                     "id",
                     "par",
                     "bv_lat",
                     "sigma_ll",
                     "sigma_ye",
                     "sigma_id",
                     "sigma_par",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_f",
                     "alpha_std",
                     "beta_bv_std",
                     "beta_bv2_std",
                     "beta_f_std",
                     "hy",
                     "hi",
                     "par",
                     "bv_lat",
                     "sigma_hi",
                     "sigma_hy",
                     "sigma_par",
                     "y_rep")),
  y_col = c("sum_recruit", "survival", "sum_recruit", "survival", "recruit"),
  stan_file_name = c("r/zinf_ars_covmat.stan",
                     "r/adult_surv_covmat.stan",
                     "r/parent_zinf_ars_covmat.stan",
                     "r/parent_surv_covmat.stan",
                     "r/nestling_surv_covmat.stan"),
  stan_file_name_co_n = c("r/zinf_ars_co_n.stan",
                          "r/adult_surv_co_n.stan",
                          "r/parent_zinf_ars_co_n.stan",
                          "r/parent_surv_co_n.stan",
                          "r/nestling_surv_co_n.stan"),
  pred_marg_func = rlang::syms(c("make_zip_preds_and_marg",
                                 "make_logit_preds_and_marg",
                                 "make_zip_preds_and_marg",
                                 "make_logit_preds_and_marg",
                                 "make_logit_preds_and_marg")),
  ppc_fun = rlang::syms(c("ppc_ars",
                          "ppc_surv",
                          "ppc_ars",
                          "ppc_surv",
                          "ppc_nest"))
)

fitmod_map <- tar_map(
  values = values_fitmod,
  names = "mod",
  tar_target(
    co_data_gp,
    gp_data_func(pheno_data = co_data,
                 lrs_path = lrs_data_path,
                 lrs_path2 = lrs_data_path2,
                 nestling_path = nestling_data_path,
                 fam_path = geno_data_paths[3],
                 ped_path = pedigree_path,
                 sex_num = sex_num_lrs,
                 sex_keep = sex_keep)
  ),
  tar_target(
    co_grm_files,
    make_grm(analysis_inds = unique(getElement(co_data_gp, "ringnr")),
             bfile = gsub(x = geno_data_paths[1], ".bed", ""),
             ncores = 4,
             mem = 4 * 6000,
             # Path to plink program:
             plink_path = plink_path,
             # Where to save result:
             dir = paste0("data/co_", mod, "_", sex_lc, "_grm")),
  ),
  tar_target(
    co_grm_obj,
    load_grm(dir = gsub(x = `[`(co_grm_files, 1), "/grm.rel.bin",  ""),
             pheno_data = co_data_gp)
  ),
  tar_target(
    co_gp,
    run_gp(pheno_data = co_data_gp,
           inverse_relatedness_matrix = getElement(co_grm_obj, "inv_grm"),
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc),
           comp_conf = TRUE)
  ),
  tar_target(
    fitness_data,
    fitdat_func(gp_data = co_data_gp,
                gp_model = co_gp,
                fitness_data_path,
                sex_num = sex_num_lrs,
                inbreeding = inbreeding,
                ped_path = pedigree_path)
  ),
  tar_target(
    # We generally find more extreme estimated bvs with the more measurements?
    bv_vs_n_meas_plot,
    dplyr::filter(fitness_data, !duplicated(ringnr)) %>%
      ggplot(aes(y = abs(bv_mean), x = ifelse(is.na(co_n), 0, co_n))) +
      geom_point() +
      labs(y = "Absolute estimated breeding value",
           x = "Number of crossover count measurements") +
      geom_smooth(formula = y ~ x, method = "loess"),
    deployment = "main"
  ),
  tar_target(
    # We generally find more extreme estimated bvs with the more measurements?
    n_vs_fitness,
    ggplot(fitness_data,
           aes(y = get(y_col), x = ifelse(is.na(co_n), 0, co_n))) +
      geom_point() +
      labs(y = "Fitness",
           x = "Number of crossover count measurements") +
      geom_smooth(formula = y ~ x, method = "loess"),
    deployment = "main"
  ),
  tar_target(
    # We generally find more extreme estimated bvs with the more measurements?
    n_vs_fitness2,
    fitness_data %>%
      dplyr::mutate(grp = cut(co_n,
                              unique(quantile(co_n, seq(0, 1, length = 5))),
                              include.lowest = TRUE)) %>%
      ggplot(aes(y = get(y_col), x = grp, group = grp)) +
      geom_boxplot() +
      labs(y = "Fitness",
           x = "Number of crossover count measurements"),
    deployment = "main"
  ),
  tar_rep(
    sim_data,
    make_sim(data = fitness_data,
             stan_data = stan_data,
             pars = sim_par_vec),
    batches = 20
  ),
  tar_rep(
    sim_data_null,
    sim_par_vec %>%
      (function(x) {x$beta_bv <- x$beta_bv2 <- 0; x}) %>%
      make_sim(data = fitness_data,
               stan_data = stan_data,
               pars = .),
    batches = 20
  ),
  tar_target(
    bv_covmat,
    inla_bv_covmat(model = co_gp, n_samp = 1e4, ncores = 16)
  ),
  tar_target(
    stan_data,
    stan_data_func(data = fitness_data, gp_data = co_data_gp, bv_covmat),
    deployment = "main"
  ),
  tar_target(
    stan_file,
    stan_file_name,
    deployment = "main",
    format = "file"
  ),
  tar_target(
    stan_model,
    stan(file = stan_file,
         data = c(stan_data, list(Y = getElement(stan_data, y_col))),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         thin = 1.6e2, # to keep final object reasonably small
         pars = stan_pars,
         model_name = paste0("stan_", mod, "_", sex_lc),
         control = list(adapt_delta = 0.96)) # up to 0.99 did not help parent
  ),
  tar_target(
    stan_samps,
    get_samps(model = stan_model,
              pars = stan_pars)
  ),
  tar_target(
    stan_post_stats,
    summary(stan_model)$summary
  ),
  tar_target(
    stan_bv_out_vs_in_plot,
    plot_bv_out_vs_in(stats = stan_post_stats, dat = stan_data)
  ),
  tar_target(
    stan_bv_out_vs_in_plot_pdf,
    ggsave_path(paste0("figs/stan_bv_in_vs_out_", mod, "_", sex_lc, ".pdf"),
                plot = stan_bv_out_vs_in_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_bv_pred_marg,
    pred_marg_func(samp = stan_samps,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f"),
                                    action = list(avg_fun_alpha,
                                                  pred_fun_bv,
                                                  pred_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  avg_fun_f),
                                    x_axis_fun = x_axis_fun_bv,
                                    marg_eff_fun = marg_eff_fun_bv))
  ),
  tar_target(
    stan_bv_pred_plot,
    plot_lines_posterior(df = getElement(stan_bv_pred_marg, "df_pred"),
                         xlab = paste0(xlab_start,
                                       "reeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex, " ", trait),
                         title = "")
  ),
  tar_target(
    stan_bv_marg_plot,
    plot_lines_posterior(df = getElement(stan_bv_pred_marg, "df_marg"),
                         xlab = paste0(xlab_start,
                                       "reeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Marginal effect on ",
                                       sex,
                                       " ",
                                       trait),
                         title = "")
  ),
  tar_target(
    stan_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/stan_", mod, "_bv_pred_", sex_lc, ".pdf"),
                plot = stan_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_bv_marg_plot_pdf,
    ggsave_path(paste0("figs/stan_", mod, "_bv_marg_", sex_lc, ".pdf"),
                plot = stan_bv_marg_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_age_pred_marg,
    pred_marg_func(samp = stan_samps,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f"),
                                    action = list(avg_fun_alpha,
                                                  avg_fun_bv,
                                                  avg_fun_bv2,
                                                  pred_fun_age_q1,
                                                  pred_fun_age_q2,
                                                  avg_fun_f),
                                    x_axis_fun = x_axis_fun_age,
                                    marg_eff_fun = marg_eff_fun_age))
  ),
  tar_target(
    stan_age_pred_plot,
    plot_lines_posterior(df = getElement(stan_age_pred_marg, "df_pred"),
                         xlab = paste0("Age"),
                         ylab = paste0("Predicted ", sex, " ", trait),
                         title = "")
  ),
  tar_target(
    stan_age_pred_plot_pdf,
    ggsave_path(paste0("figs/stan_", mod, "_age_pred_", sex_lc, ".pdf"),
                plot = stan_age_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_f_pred_marg,
    pred_marg_func(samp = stan_samps,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f"),
                                    action = list(avg_fun_alpha,
                                                  avg_fun_bv,
                                                  avg_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  pred_fun_f),
                                    x_axis_fun = x_axis_fun_f,
                                    marg_eff_fun = marg_eff_fun_f))
  ),
  tar_target(
    stan_f_pred_plot,
    plot_lines_posterior(df = getElement(stan_f_pred_marg, "df_pred"),
                         xlab = paste0("Inbreeding coefficient"),
                         ylab = paste0("Predicted ", sex, " ", trait),
                         title = "")
  ),
  tar_target(
    stan_f_pred_plot_pdf,
    ggsave_path(paste0("figs/stan_", mod, "_f_pred_", sex_lc, ".pdf"),
                plot = stan_f_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_ppc,
    ppc_fun(dat = stan_data,
            samp = stan_samps)
  ),
  ##### simulation models
  tar_target(
    stan_model_sim,
    sim_data %>%
      `[[`(., 1) %>%
      stan_data_func(data = ., gp_data = co_data_gp, bv_covmat) %>%
      stan(file = stan_file,
           data = c(., list(Y = getElement(., y_col))),
           iter = 4.8e4,
           warmup = 8e3,
           chains = 16,
           cores = 16,
           thin = 1.6e2, # to keep final object reasonably small
           pars = stan_pars,
           model_name = paste0("stan_", mod, "_", sex_lc, "_sim"),
           control = list(adapt_delta = 0.96)), # up to 0.99 did not help parent
    pattern = map(sim_data)
  ),
  tar_target(
    stan_model_sim_null,
    sim_data_null %>%
      `[[`(., 1) %>%
      stan_data_func(data = ., gp_data = co_data_gp, bv_covmat) %>%
      stan(file = stan_file,
           data = c(., list(Y = getElement(., y_col))),
           iter = 4.8e4,
           warmup = 8e3,
           chains = 16,
           cores = 16,
           thin = 1.6e2, # to keep final object reasonably small
           pars = stan_pars,
           model_name = paste0("stan_", mod, "_", sex_lc, "_sim_null"),
           control = list(adapt_delta = 0.96)), # up to 0.99 did not help parent
    pattern = map(sim_data_null)
  ),
  tar_target(
    stan_sim_summ,
    summary(stan_model_sim)$summary,
    pattern = map(stan_model_sim)
  ),
  tar_target(
    stan_sim_samps,
    get_samps(model = stan_model_sim,
              pars = stan_pars),
    pattern = map(stan_model_sim)
  ),
  tar_target(
    stan_sim_error,
    sapply(names(sim_par_vec),
           function(par) {
             x <- getElement(stan_sim_samps, par) - getElement(sim_par_vec, par)
             if (length(x) == 0) {
               x <- rep(0, 10)
             } else
               x <- x / abs(getElement(sim_par_vec, par))
             c("mean" = mean(x),
               "mode" = suppressWarnings(posterior.mode(x)),
               "median" = median(x),
               "sd" = sd(x),
               "var" = var(x),
               hdi(x, credMass = 0.95)["lower"],
               hdi(x, credMass = 0.95)["upper"])
           }),
    pattern = map(stan_sim_samps)
  ),
  tar_target(
    stan_sim_miss,
    stan_sim_error %>%
      sapply(function(mat) mat["lower", ] * mat["upper", ] > 0) %>%
      rowMeans()
  ),
  tar_target(
    stan_sim_null_summ,
    summary(stan_model_sim_null)$summary,
    pattern = map(stan_model_sim_null)
  ),
  tar_target(
    stan_sim_null_samps,
    get_samps(model = stan_model_sim_null,
              pars = stan_pars),
    pattern = map(stan_model_sim_null)
  ),
  tar_target(
    stan_sim_null_error,
    sapply(names(sim_par_vec),
           function(par) {
             if (!par %in% c("beta_bv", "beta_bv2")) {
               x <- getElement(stan_sim_null_samps, par) -
                 getElement(sim_par_vec, par)
             } else {
               x <- getElement(stan_sim_null_samps, par)
             }
             if (length(x) == 0)
               x <- rep(0, 10)
             else
               x <- x / abs(getElement(sim_par_vec, par))
             c("mean" = mean(x),
               "mode" = suppressWarnings(posterior.mode(x)),
               "median" = median(x),
               "sd" = sd(x),
               "var" = var(x),
               hdi(x, credMass = 0.95)["lower"],
               hdi(x, credMass = 0.95)["upper"])
           }),
    pattern = map(stan_sim_null_samps)
  ),
  tar_target(
    stan_sim_null_falsepos_or_miss,
    stan_sim_null_error %>%
      sapply(function(mat) mat["lower", ] * mat["upper", ] > 0) %>%
      rowMeans()
  ),
  tar_target(
    stan_sim_bv_reliability,
    (stan_sim_summ %>%
       `[`(grepl(x = rownames(.), pattern = "bv_lat"), "mean") %>%
       var()
    ) / (sim_data %>%
           (function(dat) {
             dat <- dat[[1]]
             dat %>%
               `$`("bv_true") %>%
               `[`(order(dat$ringnr_num)) %>%
               `[`(match(unique(dat[order(dat$ringnr_num)]$ringnr),
                         dat[order(dat$ringnr_num)]$ringnr)) %>%
               var()
           })),
    pattern = map(sim_data, stan_sim_summ)
  ),
  tar_target(
    stan_sim_bv_plot,
    make_sim_bv_plot(summ = stan_sim_summ, sim_data = sim_data),
    pattern = map(sim_data, stan_sim_summ)
  ),
  tar_target(
    stan_sim_bv_plot_pdf,
    ggsave_path(paste0("figs/", mod, "_sim_bv_plot_", sex_lc, ".pdf"),
                plot = stan_sim_bv_plot[[1]],
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  ################ co_n modeller
  tar_target(
    stan_file_co_n,
    stan_file_name_co_n,
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_model_co_n,
    stan(file = stan_file_co_n,
         data = c(stan_data, list(Y = getElement(stan_data, y_col))),
         iter = 1.92e6,
         warmup = 3.2e5,
         chains = 16,
         cores = 16,
         # Remove random effects in zero-inflation component
         pars = c(stan_pars, "beta_co_n", "beta_co_n_std"),
         control = list(adapt_delta = 0.96),
         model_name = paste0("stan_", mod, "_", sex_lc, "_co_n"),
         thin = 6.4e3) # to keep final object reasonably small
  ),
  tar_target(
    stan_samps_co_n,
    get_samps(model = stan_model_co_n,
              pars = c(stan_pars, "beta_co_n", "beta_co_n_std"))
  ),
  tar_target(
    stan_post_stats_co_n,
    summary(stan_model_co_n)$summary
  ),
  tar_target(
    stan_bv_pred_marg_co_n,
    pred_marg_func(samp = stan_samps_co_n,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f",
                                                  "beta_co_n"),
                                    action = list(avg_fun_alpha,
                                                  pred_fun_bv,
                                                  pred_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  avg_fun_f,
                                                  avg_fun_co_n),
                                    x_axis_fun = x_axis_fun_bv,
                                    marg_eff_fun = marg_eff_fun_bv))
  ),
  tar_target(
    stan_bv_pred_plot_co_n,
    plot_lines_posterior(df = getElement(stan_bv_pred_marg_co_n, "df_pred"),
                         xlab = paste0(xlab_start,
                                       "reeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex, " ", trait),
                         title = "")
  ),
  tar_target(
    stan_bv_pred_plot_co_n_pdf,
    ggsave_path(paste0("figs/stan_", mod, "_bv_pred_", sex_lc, "_co_n.pdf"),
                plot = stan_bv_pred_plot_co_n,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_co_n_pred_marg_co_n,
    pred_marg_func(samp = stan_samps_co_n,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f",
                                                  "beta_co_n"),
                                    action = list(avg_fun_alpha,
                                                  avg_fun_bv,
                                                  avg_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  avg_fun_f,
                                                  pred_fun_co_n),
                                    x_axis_fun = x_axis_fun_co_n,
                                    marg_eff_fun = marg_eff_fun_co_n))
  ),
  tar_target(
    stan_co_n_pred_plot_co_n,
    plot_lines_posterior(df = getElement(stan_co_n_pred_marg_co_n, "df_pred"),
                         xlab = "Number of crossover count measurements",
                         ylab = paste0("Predicted ", sex, " ", trait),
                         title = "")
  ),
  tar_target(
    stan_co_n_pred_plot_co_n_pdf,
    ggsave_path(paste0("figs/stan_", mod, "_co_n_pred_", sex_lc, "_co_n.pdf"),
                plot = stan_co_n_pred_plot_co_n,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  )
)


values_sex <- tibble(
  sex = c("female", "male"),
  sex_lc = c("f", "m"),
  sex_keep = c("F", "M"),
  sex_num_lrs = c(2, 1)
)

sex_map <- tar_map(
  values = values_sex,
  names = "sex_lc",
  fitmod_map,
  tar_target(
    co_data, # All measurements of co count
    prep_co_data(recomb_data_path2,
                 lrs_data_path,
                 lrs_data_path2,
                 sex_keep = sex_keep)
  ),
  tar_target(
    # Test sets of n-fold CV in within-population genomic prediction
    cv_test_sets,
    make_cv_test_sets(analysis_inds = unique(co_data$ringnr),
                      num_folds = 10),
    deployment = "main"
  ),
  tar_target(
    co_data_cv,
    make_co_data_cv(co_data,
                    test_set = cv_test_sets,
                    sex_lc = sex_lc),
    pattern = map(cv_test_sets)
  ),
  tar_target(
    co_data_mean, # Measurements of co count averaged over inds
    co_data %>%
      group_by(id, sex) %>%
      summarise(co_count = mean(co_count),
                total_coverage = mean(total_coverage),
                hatch_year = mean(hatch_year),
                first_locality = mean(first_locality),
                n = n())
  ),
  tar_target(
    co_data_rand, # 1 random measurement of co count
    co_data %>%
      group_by(id, sex) %>%
      summarise(co_count = co_count[sample(n(), 1)],
                total_coverage = total_coverage[sample(n(), 1)],
                n = n())
  ),
  tar_target(
    num_meas_vs_co_count,
    ggplot(data = co_data_rand, aes(x = co_count, y = n)) +
      geom_point() +
      geom_smooth(method = "lm")
  )
  # tar_target(
  #   co_gam, # genomic animal model for sex-specific crossover rate
  #   run_gp(pheno_data = co_data,
  #          inverse_relatedness_matrix = co_grm_obj$inv_grm,
  #          effects_vec = inla_effects_gp_vector_grm_all,
  #          y = paste0("co_count_", sex_lc))
  # ),
  # tar_target(
  #   co_10fcv, # genomic animal model for sex-specific crossover rate
  #   run_gp(pheno_data = co_data_cv,
  #          inverse_relatedness_matrix = co_grm_obj$inv_grm,
  #          effects_vec = inla_effects_gp_vector_grm_all,
  #          y = paste0("co_count_", sex_lc, "_test")),
  #   pattern = map(co_data_cv)
  # ),
  # tar_target(
  #   co_10fcv_bvacc, # genomic animal model for sex-specific crossover rate
  #   get_10fcv_acc(data = co_data_cv,
  #                 model = co_10fcv,
  #                 y_str = paste0("co_count_", sex_lc)),
  #   pattern = map(co_data_cv, co_10fcv),
  #   iteration = "vector"
  # ),
  # tar_target(
  #   co_gam_cv, # genomic animal model for sex-specific crossover rate
  #   inla_cv(model = co_gam,
  #           pheno_data = co_data)
  # ),
  # tar_target(
  #   surv_samp_pairs_plot,
  #   as.data.frame(surv_samps_covmat)[, c("energy", "ll.1", "ye.1", "id.1",
  #                                        surv_pars[c(1:6, 17:19)])] %>%
  #     dplyr::mutate(sigma_ll = log(sigma_ll),
  #                   sigma_ye = log(sigma_ye),
  #                   sigma_id = log(sigma_id)) %>%
  #     ggpairs() # aes(color = factor(surv_samps_f$divergent)))
  # ),
)

list(
  sex_map,
  tar_target(
    recomb_data_path,
    "data/20240910_Sparrow_Recomb_Data.txt",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    recomb_data_path2,
    co_data_rename_cols("data/20250706_Sparrow_YAPP/2_recsumm_Crossover_Count_per_individual_post_QC.txt"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    geno_data_paths,
    paste0("data/70K_200K_maf_geno_mind_v5.", c("bed", "bim", "fam")),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    nestling_data_path,
    "data/Nestlingdata_20240312_fix.csv",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    lrs_data_path,
    "data/LRS_data_20230803_Helgeland_fix.csv",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    lrs_data_path2, # only to be used for hatch/year info, not fitness data
    "data/Missing_LRS_Sparrows_revised_WithInfo_fix.csv",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    morph_data_path,
    "data/AdultMorphology_20240201_fix.csv",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    pedigree_path,
    "data/20230317_Sparrow_Pedigree.txt",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    inla_effects_gp_vector_grm_all,
    c(################## Fixed effects:
      # Intercept
      "1",
      # Total coverage effect
      "total_coverage_scaled",
      ################# Random effects:
      # Year effect (iid random effect)
      "f(hatch_year, model = \"iid\", hyper = prior$hyperpar_var)",
      # Island effect (iid random effect)
      "f(first_locality, model = \"iid\", hyper = prior$hyperpar_var)",
      # Breeding values (using a GRM)
      "f(id1,
      values = as.numeric(colnames(inverse_relatedness_matrix)),
      model = \"generic0\",
      hyper = prior$hyperpar_var,
      constr = FALSE,
      Cmatrix = inverse_relatedness_matrix)",
      # ID effect (use this when using single measurements as the response)
      "f(id2,
      model = \"iid\",
      values = unique(pheno_data$id2),
      hyper = prior$hyperpar_var)"
    ),
    deployment = "main"
  ),
  tar_target(
    figureuruurue,
    list(post_stats_sim_rel_error_ars_adult_f,
         post_stats_sim_rel_error_ars_adult_m,
         post_stats_sim_rel_error_surv_adult_f,
         post_stats_sim_rel_error_surv_adult_m,
         post_stats_sim_rel_error_ars_parent_f,
         post_stats_sim_rel_error_ars_parent_m,
         post_stats_sim_rel_error_surv_parent_f,
         post_stats_sim_rel_error_surv_parent_m,
         post_stats_sim_rel_error_nest_f,
         post_stats_sim_rel_error_nest_m) %>%
      bind_rows(.id = "model") %>%
      ggplot(data = .) +
      geom_boxplot(aes(y = value, group = name)) +
      facet_grid(samp_type ~ model),
    deployment = "main"
  ),
  # tar_target(
  #   stan_file_adult_surv_covmat_co_n,
  #   "r/adult_surv_covmat_co_n.stan",
  #   format = "file",
  #   deployment = "main"
  # ),
  # tar_target(
  #   stan_file_adult_ars_zinf_covmat_co_n,
  #   "r/zinf_ars_covmat_co_n.stan",
  #   format = "file",
  #   deployment = "main"
  # ),
  tar_target(
    inbreeding,
    find_inbreeding(gsub(x = geno_data_paths[1], ".bed", ""),
                    ncores = 4,
                    mem = 4 * 6000,
                    # Path to plink program:
                    plink_path = plink_path),
    format = "file"
  ),
  tar_target(
    stan_model_sim_arstest,
    sim_data_ars_adult_m %>%
      `[[`(., 1) %>%
      (function(df) {df$bv_mean <- df$bv_true; df}) %>%
      make_stan_data_adult(data = ., gp_data = co_data_gp_ars_adult_m, bv_covmat_ars_adult_m) %>%
      stan(file = "r/zinf_ars_truetest.stan",
           data = c(., list(Y = getElement(., "sum_recruit"))),
           iter = 4.8e4,
           warmup = 8e3,
           chains = 16,
           cores = 16,
           thin = 1.6e2, # to keep final object reasonably small
           pars = c("alpha",
                    "beta_bv",
                    "beta_bv2",
                    "beta_age_q1",
                    "beta_age_q2",
                    "beta_f",
                    "alpha_std",
                    "beta_bv_std",
                    "beta_bv2_std",
                    "beta_age_q1_std",
                    "beta_age_q2_std",
                    "beta_f_std",
                    "ye",
                    "ll",
                    "id",
                    "sigma_ll",
                    "sigma_ye",
                    "sigma_id",
                    "alpha_zi",
                    "alpha_zi_std",
                    "theta",
                    "y_rep"),
           model_name = paste0("stan_ars_adult_m_sim_arstest"),
           control = list(adapt_delta = 0.96)), # up to 0.99 did not help parent
    pattern = map(sim_data_ars_adult_m)
  ),
  tar_target(
    stan_sim_summ_arstest,
    summary(stan_model_sim_arstest)$summary,
    pattern = map(stan_model_sim_arstest)
  ),
  tar_target(
    stan_sim_samps_arstest,
    get_samps(model = stan_model_sim_arstest,
              pars = c("alpha",
                       "beta_bv",
                       "beta_bv2",
                       "beta_age_q1",
                       "beta_age_q2",
                       "beta_f",
                       "alpha_std",
                       "beta_bv_std",
                       "beta_bv2_std",
                       "beta_age_q1_std",
                       "beta_age_q2_std",
                       "beta_f_std",
                       "ye",
                       "ll",
                       "id",
                       "sigma_ll",
                       "sigma_ye",
                       "sigma_id",
                       "alpha_zi",
                       "alpha_zi_std",
                       "theta",
                       "y_rep")),
    pattern = map(stan_model_sim_arstest)
  ),
  tar_target(
    stan_sim_error_arstest,
    sapply(names(list("alpha" = 0.26,
                      "beta_bv" = -0.16,
                      "beta_bv2" = -1.2,
                      "beta_age_q1" = 9.4,
                      "beta_age_q2" = -14,
                      "beta_f" = -2.3,
                      "alpha_zi" = -0.80,
                      "sigma_ye" = 0.28,
                      "sigma_ll" = 0.18,
                      "sigma_id" = 0.30,
                      "sigma_res" = 0)),
           function(par) {
             x <- getElement(stan_sim_samps_arstest, par) - getElement(list("alpha" = 0.26,
                                                                            "beta_bv" = -0.16,
                                                                            "beta_bv2" = -1.2,
                                                                            "beta_age_q1" = 9.4,
                                                                            "beta_age_q2" = -14,
                                                                            "beta_f" = -2.3,
                                                                            "alpha_zi" = -0.80,
                                                                            "sigma_ye" = 0.28,
                                                                            "sigma_ll" = 0.18,
                                                                            "sigma_id" = 0.30,
                                                                            "sigma_res" = 0), par)
             if (length(x) == 0)
               x <- rep(0, 10)
             else
               x <- x / abs(getElement(list("alpha" = 0.26,
                                            "beta_bv" = -0.16,
                                            "beta_bv2" = -1.2,
                                            "beta_age_q1" = 9.4,
                                            "beta_age_q2" = -14,
                                            "beta_f" = -2.3,
                                            "alpha_zi" = -0.80,
                                            "sigma_ye" = 0.28,
                                            "sigma_ll" = 0.18,
                                            "sigma_id" = 0.30,
                                            "sigma_res" = 0), par))
             c("mean" = mean(x),
               "mode" = suppressWarnings(posterior.mode(x)),
               "median" = median(x),
               "sd" = sd(x),
               "var" = var(x),
               hdi(x, credMass = 0.95)["lower"],
               hdi(x, credMass = 0.95)["upper"])
           }),
    pattern = map(stan_sim_samps_arstest)
  ),
  tar_target(
    stan_sim_miss_arstest,
    stan_sim_error_arstest %>%
      sapply(function(mat) mat["lower", ] * mat["upper", ] > 0) %>%
      rowMeans()
  )
)
