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
  workers = 40,
  seconds_idle = 120,
  tasks_max = 30,
  options_cluster = crew_options_slurm(
    script_lines = paste("#SBATCH --account=share-nv-bio \n module load",
                         "R/4.4.2-gfbf-2024a R-bundle-CRAN/2024.11-foss-2024a",
                         "CMake/3.29.3-GCCcore-13.3.0"),
    log_output = "Jobs/%A_%a.log",
    memory_gigabytes_per_cpu = 6,
    cpus_per_task = 16,
    time_minutes = 60 * 1 * 1, # minutes * hours * days
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
               "grid",
               "MCMCglmm",
               "HDInterval",
               "bayesplot",
               "patchwork"),
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
          "nest", # parental
          "n2"), # self
  gp_data_func = rlang::syms(c("make_gp_data_adult",
                               "make_gp_data_adult",
                               "make_gp_data_parent",
                               "make_gp_data_parent",
                               "make_gp_data_nest",
                               "make_gp_data_n2")),
  fitness_data_path = rlang::syms(c("lrs_data_path",
                                    "lrs_data_path",
                                    "lrs_data_path",
                                    "lrs_data_path",
                                    "nestling_data_path",
                                    "nestling_data_path")),
  fitdat_func = rlang::syms(c("make_data_adult",
                              "make_data_adult",
                              "make_data_parent",
                              "make_data_parent",
                              "make_data_nest",
                              "make_data_n2")),
  trait = c("annual reproductive success",
            "annual survival",
            "annual reproductive success",
            "annual survival",
            "nestling survival",
            "nestling survival"),
  model_label = c("GV-F", "GV-F", "PGV-F", "PGV-F", "PGV-F", "GV-F"),
  trait_short = c("ARS", "AS", "ARS", "AS", "NS", "NS"),
  xlab_start = c("G", "G", "Parental g", "Parental g", "Parental g", "G"),
  stan_data_func = rlang::syms(c("make_stan_data_adult",
                                 "make_stan_data_adult",
                                 "make_stan_data_parent",
                                 "make_stan_data_parent",
                                 "make_stan_data_nest",
                                 "make_stan_data_n2")),
  stan_pars = list(c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_age_q1",
                     "beta_age_q2",
                     "beta_f",
                     "ye",
                     "ll",
                     "id",
                     "bv_lat",
                     "sigma_ll",
                     "sigma_ye",
                     "sigma_id",
                     "alpha_zi",
                     "theta",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_age_q1",
                     "beta_age_q2",
                     "beta_f",
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
                     "ye",
                     "ll",
                     "id",
                     "dam",
                     "sire",
                     "bv_lat",
                     "sigma_ll",
                     "sigma_ye",
                     "sigma_id",
                     "sigma_dam",
                     "sigma_sire",
                     "alpha_zi",
                     "theta",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_age_q1",
                     "beta_age_q2",
                     "beta_f",
                     "ye",
                     "ll",
                     "id",
                     "dam",
                     "sire",
                     "bv_lat",
                     "sigma_ll",
                     "sigma_ye",
                     "sigma_id",
                     "sigma_dam",
                     "sigma_sire",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_f",
                     "beta_hatch_doy",
                     "beta_hatch_doy2",
                     "beta_first_dna_age",
                     "hy",
                     "hi",
                     "dam",
                     "sire",
                     "clutch",
                     "bv_lat",
                     "sigma_hi",
                     "sigma_hy",
                     "sigma_dam",
                     "sigma_sire",
                     "sigma_clutch",
                     "y_rep"),
                   c("alpha",
                     "beta_bv",
                     "beta_bv2",
                     "beta_f",
                     "beta_hatch_doy",
                     "beta_hatch_doy2",
                     "beta_first_dna_age",
                     "hy",
                     "hi",
                     "dam",
                     "sire",
                     "clutch",
                     "bv_lat",
                     "sigma_hi",
                     "sigma_hy",
                     "sigma_clutch",
                     "sigma_dam",
                     "sigma_sire",
                     "y_rep")),
  y_col = c("sum_recruit", "survival", "sum_recruit", "survival", rep("recruit", 2)),
  stan_file_name = c("r/zinf_ars_covmat.stan",
                     "r/adult_surv_covmat.stan",
                     "r/parent_zinf_ars_covmat.stan",
                     "r/parent_surv_covmat.stan",
                     "r/nestling_surv_covmat.stan",
                     "r/n2.stan"),
  stan_file_name_co_n = c("r/zinf_ars_co_n.stan",
                          "r/adult_surv_co_n.stan",
                          "r/parent_zinf_ars_co_n.stan",
                          "r/parent_surv_co_n.stan",
                          "r/nestling_surv_co_n.stan",
                          "r/n2_co_n.stan"),
  pred_marg_func = rlang::syms(c("make_zip_preds_and_marg",
                                 "make_logit_preds_and_marg",
                                 "make_zip_preds_and_marg",
                                 "make_logit_preds_and_marg",
                                 "make_logit_preds_and_marg",
                                 "make_logit_preds_and_marg")),
  ppc_fun = rlang::syms(c("ppc_ars",
                          "ppc_surv",
                          "ppc_ars",
                          "ppc_surv",
                          "ppc_nest",
                          "ppc_nest"))
)

fitmod_map <- tar_map(
  values = values_fitmod,
  names = "mod",
  ########## GP ##########
  tar_target(
    co_data_gp,
    gp_data_func(pheno_data = co_data,
                 lrs_path = lrs_data_path,
                 lrs_path2 = lrs_data_path2,
                 nestling_path = nestling_data_path,
                 fam_path = geno_data_paths[3],
                 ped_path = pedigree_path,
                 sex_num = sex_num_lrs,
                 sex_keep = sex_keep,
                 froh_file = froh_file)
  ),
  tar_target(
    cv_test_sets,
    co_data_gp %>%
      dplyr::filter(!is.na(n)) %>%
      `$`(id_red) %>%
      unique() %>%
      make_cv_test_sets(num_folds = 10)
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
    co_gp_cv,
    run_gp_cv(pheno_data = co_data_gp,
              inverse_relatedness_matrix = getElement(co_grm_obj, "inv_grm"),
              effects_vec = inla_effects_gp_vector_grm_all,
              y = paste0("co_count_", sex_lc),
              comp_conf = FALSE,
              test_set = cv_test_sets),
    pattern = map(cv_test_sets)
  ),
  tar_target(
    co_gp_cv_acc,
    cv_acc_fun(model = co_gp_cv,
               pheno_data = co_data_gp %>% dplyr::filter(!is.na(n)),
               test_set = cv_test_sets,
               y = paste0("co_count_", sex_lc)),
    pattern = map(co_gp_cv, cv_test_sets)
  ),
  tar_target(
    co_gp_cv_mean_acc,
    co_gp_cv_acc %>% as.data.frame() %>% t() %>% apply(2, mean),
    deployment = "main"
  ),
  ########## (P)GV-F ##########
  tar_target(
    fitness_data,
    fitdat_func(gp_data = co_data_gp,
                gp_model = co_gp,
                fitness_data_path,
                sex_num = sex_num_lrs,
                froh_file = froh_file,
                ped_path = pedigree_path)
  ),
  tar_target(
    num_co_meas_plot,
    dplyr::filter(fitness_data, !duplicated(ringnr)) %>%
      ggplot(aes(x = co_n)) +
      geom_histogram(bins = 42) +
      labs(y = "Count",
           x = "Number of ACC measurements"),
    deployment = "main"
  ),
  tar_target(
    # We generally find more extreme estimated bvs with the more measurements?
    bv_absmean_vs_co_meas_plot,
    dplyr::filter(fitness_data, !duplicated(ringnr)) %>%
      ggplot(aes(y = abs(bv_mean), group = (co_n != 0), x = (co_n != 0))) +
      geom_boxplot() +
      labs(y = "Absolute estimated genetic value",
           x = "Phenotyped for ACC or not"),
    deployment = "main"
  ),
  tar_target(
    # We generally find more extreme estimated bvs with the more measurements?
    bv_sd_vs_co_meas_plot,
    dplyr::filter(fitness_data, !duplicated(ringnr)) %>%
      ggplot(aes(y = bv_sd, group = (co_n != 0), x = (co_n != 0))) +
      geom_boxplot() +
      labs(y = "S.d. of genetic value",
           x = "Phenotyped for ACC or not"),
    deployment = "main"
  ),
  tar_target(
    # We generally find more extreme estimated bvs with the more measurements?
    bv_vs_n_meas_plot,
    dplyr::filter(fitness_data, !duplicated(ringnr)) %>%
      ggplot(aes(y = abs(bv_mean), x = ifelse(is.na(co_n), 0, co_n))) +
      geom_point() +
      labs(y = "Absolute estimated genetic value",
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
  tar_target(
    mvn_test,
    test_mvnorm(model = co_gp)
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
         control = list(adapt_delta = 0.96))
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
    stan_bv_pred_marg,
    pred_marg_func(samp = stan_samps,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f",
                                                  "beta_hatch_doy",
                                                  "beta_hatch_doy2",
                                                  "beta_first_dna_age"),
                                    action = list(avg_fun_alpha,
                                                  pred_fun_bv,
                                                  pred_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  avg_fun_f,
                                                  avg_fun_hatch_doy,
                                                  avg_fun_hatch_doy2,
                                                  avg_fun_first_dna_age),
                                    x_axis_fun = x_axis_fun_bv,
                                    marg_eff_fun = marg_eff_fun_bv))
  ),
  ########## (P)GV-F result plots ##########
  tar_target(
    stan_bv_pred_plot,
    plot_lines_posterior(df = getElement(stan_bv_pred_marg, "df_pred"),
                         xlab = paste0(xlab_start, "enetic value for ACC"),
                         ylab = paste0("Predicted ", trait_short),
                         title = "")
  ),
  tar_target(
    stan_bv_marg_plot,
    plot_lines_posterior(df = getElement(stan_bv_pred_marg, "df_marg"),
                         xlab = paste0(xlab_start, "enetic value for ACC"),
                         ylab = paste0("Marginal effect on ", trait_short),
                         title = "")
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
                                                  "beta_f",
                                                  "beta_hatch_doy",
                                                  "beta_hatch_doy2",
                                                  "beta_first_dna_age"),
                                    action = list(avg_fun_alpha,
                                                  avg_fun_bv,
                                                  avg_fun_bv2,
                                                  pred_fun_age_q1,
                                                  pred_fun_age_q2,
                                                  avg_fun_f,
                                                  avg_fun_hatch_doy,
                                                  avg_fun_hatch_doy2,
                                                  avg_fun_first_dna_age),
                                    x_axis_fun = x_axis_fun_age,
                                    marg_eff_fun = marg_eff_fun_age))
  ),
  tar_target(
    stan_age_pred_plot,
    plot_lines_posterior(df = getElement(stan_age_pred_marg, "df_pred"),
                         xlab = paste0("Age"),
                         ylab = paste0("Predicted ", trait),
                         title = "")
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
                                                  "beta_f",
                                                  "beta_hatch_doy",
                                                  "beta_hatch_doy2",
                                                  "beta_first_dna_age"),
                                    action = list(avg_fun_alpha,
                                                  avg_fun_bv,
                                                  avg_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  pred_fun_f,
                                                  avg_fun_hatch_doy,
                                                  avg_fun_hatch_doy2,
                                                  avg_fun_first_dna_age),
                                    x_axis_fun = x_axis_fun_f,
                                    marg_eff_fun = marg_eff_fun_f))
  ),
  tar_target(
    stan_f_pred_plot,
    plot_lines_posterior(df = getElement(stan_f_pred_marg, "df_pred"),
                         xlab = paste0("Inbreeding coefficient"),
                         ylab = paste0("Predicted ", trait),
                         title = "")
  ),
  tar_target(
    stan_hatch_doy_pred_marg,
    pred_marg_func(samp = stan_samps,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f",
                                                  "beta_hatch_doy",
                                                  "beta_hatch_doy2",
                                                  "beta_first_dna_age"),
                                    action = list(avg_fun_alpha,
                                                  avg_fun_bv,
                                                  avg_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  avg_fun_f,
                                                  pred_fun_hatch_doy,
                                                  pred_fun_hatch_doy2,
                                                  avg_fun_first_dna_age),
                                    x_axis_fun = x_axis_fun_hatch_doy,
                                    marg_eff_fun = marg_eff_fun_hatch_doy))
  ),
  tar_target(
    stan_hatch_doy_pred_plot,
    plot_lines_posterior(df = getElement(stan_hatch_doy_pred_marg, "df_pred"),
                         xlab = paste0("Hatch date (day of year)"),
                         ylab = paste0("Predicted ", trait),
                         title = "")
  ),
  tar_target(
    stan_first_dna_age_pred_marg,
    pred_marg_func(samp = stan_samps,
                   data = fitness_data,
                   pred_info = list(coef_name = c("alpha",
                                                  "beta_bv",
                                                  "beta_bv2",
                                                  "beta_age_q1",
                                                  "beta_age_q2",
                                                  "beta_f",
                                                  "beta_hatch_doy",
                                                  "beta_hatch_doy2",
                                                  "beta_first_dna_age"),
                                    action = list(avg_fun_alpha,
                                                  avg_fun_bv,
                                                  avg_fun_bv2,
                                                  avg_fun_age_q1,
                                                  avg_fun_age_q2,
                                                  avg_fun_f,
                                                  avg_fun_hatch_doy,
                                                  avg_fun_hatch_doy2,
                                                  pred_fun_first_dna_age),
                                    x_axis_fun = x_axis_fun_first_dna_age,
                                    marg_eff_fun = marg_eff_fun_first_dna_age))
  ),
  tar_target(
    stan_first_dna_age_pred_plot,
    plot_lines_posterior(df = getElement(stan_first_dna_age_pred_marg, "df_pred"),
                         xlab = paste0("Age at first DNA sampling (days)"),
                         ylab = paste0("Predicted ", trait),
                         title = "")
  ),
  tar_target(
    stan_ppc,
    ppc_fun(dat = fitness_data, samp = stan_samps)
  ),
  tar_target(
    stan_fixed_eff_coefs_plot,
    plot_fixed_eff_coefs(
      samp_df = (stan_samps %>%
                   as.data.frame() %>%
                   select(matches("^(alph|bet)"))),
      subtit = paste0(model_label, ", ", sex, " ", trait_short))
  ),
  tar_target(
    stan_fixed_eff_coefs_plot_png,
    ggsave_path(paste0("figs/coef_figs/",
                       trait_short, "_", model_label, "_", sex, ".png"),
                plot = stan_fixed_eff_coefs_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_rand_eff_sd_plot,
    plot_rand_eff_sd(
      samp_df = (stan_samps %>%
                   as.data.frame() %>%
                   select(matches("^sigma_"))),
      subtit = paste0(model_label, ", ", sex, " ", trait_short))
  ),
  tar_target(
    stan_rand_eff_sd_plot_png,
    ggsave_path(paste0("figs/sigma_figs/",
                       trait_short, "_", model_label, "_", sex, ".png"),
                plot = stan_rand_eff_sd_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_year_levels_plot,
    plot_levels(
      samp_df = (stan_samps %>%
                   as.data.frame() %>%
                   select(matches("^(ye|hy)"))),
      tit = paste0(c("Y", "Hatch y"),
                   "ear effect levels (linear predictor scale)"),
      subtit = paste0(model_label, ", ", sex, " ", trait_short),
      dat = fitness_data,
      samp_colname = c("ye", "hy"),
      num_col = c("y_num", "hy_num"),
      dat_col = c("year", "hatch_year"),
      isl = FALSE)
  ),
  tar_target(
    stan_year_levels_plot_png,
    ggsave_path(paste0("figs/year_level_figs/",
                       trait_short, "_", model_label, "_", sex, ".png"),
                plot = stan_year_levels_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_island_levels_plot,
    plot_levels(
      samp_df = (stan_samps %>%
                   as.data.frame() %>%
                   select(matches("^(ll|hi)"))),
      tit = paste0(c("I", "Hatch i"),
                   "sland effect levels (linear predictor scale)"),
      subtit = paste0(model_label, ", ", sex, " ", trait_short),
      dat = fitness_data,
      samp_colname = c("ll", "hi"),
      num_col = c("ll_num", "hi_num"),
      dat_col = c("last_locality", "hatch_island"),
      isl = TRUE)
  ),
  tar_target(
    stan_island_levels_plot_png,
    ggsave_path(paste0("figs/isl_level_figs/",
                       trait_short, "_", model_label, "_", sex, ".png"),
                plot = stan_island_levels_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  ########## co_n models ##########
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
                         xlab = paste0(xlab_start, "enetic value for ACC"),
                         ylab = paste0("Predicted ", trait),
                         title = "")
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
                         xlab = "Number of ACC measurements",
                         ylab = paste0("Predicted ", trait),
                         title = "")
  )
)

################# Models for direct impact of co on offspring fitness ########

values_dirfit <- tibble(
  dirfit = c("ars", "as", "ns"),
  fitness_data_path_dir = rlang::syms(c("lrs_data_path",
                                        "lrs_data_path",
                                        "nestling_data_path")),
  fitness_data_dir_func = rlang::syms(c("make_data_parent_dir",
                                        "make_data_parent_dir",
                                        "make_data_nest_dir")),
  dirfit_func = rlang::syms(c("dirfit_func_ars",
                              "dirfit_func_as",
                              "dirfit_func_ns")),
  dirfit_samps_func = rlang::syms(c("get_dirfit_samps_ars_as",
                                    "get_dirfit_samps_ars_as",
                                    "get_dirfit_samps_ns")),
  eval_func = rlang::syms(c("eval_func_ars",
                            "eval_func_as",
                            "eval_func_ns")),
  dirfit_func_parsum = rlang::syms(c("dirfit_func_parsum_ars",
                                     "dirfit_func_parsum_as",
                                     "dirfit_func_parsum_ns")),
  dirfit_samps_parsum_func = rlang::syms(c("get_dirfit_samps_parsum_ars_as",
                                           "get_dirfit_samps_parsum_ars_as",
                                           "get_dirfit_samps_parsum_ns")),
  eval_func_parsum = rlang::syms(c("eval_func_parsum_ars",
                                   "eval_func_parsum_as",
                                   "eval_func_parsum_ns")),
  pred_marg_func_dirfit = rlang::syms(c("make_zip_preds_and_marg",
                                        "make_logit_preds_and_marg",
                                        "make_logit_preds_and_marg")),
  trait = c("annual reproductive success",
            "annual survival",
            "nestling survival"),
  trait_short = c("ARS", "AS", "NS"))

dirfit_map <- tar_map(
  values = values_dirfit,
  names = "dirfit",
  tar_target(
    fitness_data_dir,
    fitness_data_dir_func(fitness_data_path_dir,
                          froh_file = froh_file,
                          sex_num = sex_num_lrs,
                          sex_lc = sex_lc,
                          co_dat_path = recomb_data_path2)
  ),
  tar_target(
    direct_fitness,
    dirfit_func(data = fitness_data_dir)
  ),
  tar_target(
    direct_fitness_parsum,
    dirfit_func_parsum(data = fitness_data_dir)
  ),
  tar_target(
    dirfit_samps,
    dirfit_samps_func(model = direct_fitness,
                      data = fitness_data_dir,
                      eval_func = eval_func)
  ),
  tar_target(
    dirfit_samps_parsum,
    dirfit_samps_parsum_func(model = direct_fitness_parsum,
                             data = fitness_data_dir,
                             eval_func = eval_func_parsum)
  ),
  tar_target(
    dirfit_co_count_sire_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_sire",
                                     "beta_co_count_sire2",
                                     "beta_co_count_dam",
                                     "beta_co_count_dam2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     pred_fun_co_count_sire,
                                     pred_fun_co_count_sire2,
                                     avg_fun_co_count_dam,
                                     avg_fun_co_count_dam2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     avg_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_co_count_sire,
                       marg_eff_fun = marg_eff_fun_co_count_sire))
  ),
  tar_target(
    dirfit_co_count_sire_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_co_count_sire_pred_marg,
                                         "df_pred"),
                         xlab = paste0("ACC in gamete from sire"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_co_count_dam_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_sire",
                                     "beta_co_count_sire2",
                                     "beta_co_count_dam",
                                     "beta_co_count_dam2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_sire,
                                     avg_fun_co_count_sire2,
                                     pred_fun_co_count_dam,
                                     pred_fun_co_count_dam2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     avg_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_co_count_dam,
                       marg_eff_fun = marg_eff_fun_co_count_dam))
  ),
  tar_target(
    dirfit_co_count_dam_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_co_count_dam_pred_marg,
                                         "df_pred"),
                         xlab = paste0("ACC in gamete from dam"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_age_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_sire",
                                     "beta_co_count_sire2",
                                     "beta_co_count_dam",
                                     "beta_co_count_dam2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_sire,
                                     avg_fun_co_count_sire2,
                                     avg_fun_co_count_dam,
                                     avg_fun_co_count_dam2,
                                     pred_fun_age_q1,
                                     pred_fun_age_q2,
                                     avg_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_age,
                       marg_eff_fun = marg_eff_fun_age))
  ),
  tar_target(
    dirfit_age_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_age_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Age"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_f_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_sire",
                                     "beta_co_count_sire2",
                                     "beta_co_count_dam",
                                     "beta_co_count_dam2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_sire,
                                     avg_fun_co_count_sire2,
                                     avg_fun_co_count_dam,
                                     avg_fun_co_count_dam2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     pred_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_f,
                       marg_eff_fun = marg_eff_fun_f))
  ),
  tar_target(
    dirfit_f_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_f_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Inbreeding coefficient"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_hatch_doy_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_sire",
                                     "beta_co_count_sire2",
                                     "beta_co_count_dam",
                                     "beta_co_count_dam2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_sire,
                                     avg_fun_co_count_sire2,
                                     avg_fun_co_count_dam,
                                     avg_fun_co_count_dam2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     avg_fun_f,
                                     pred_fun_hatch_doy,
                                     pred_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_hatch_doy,
                       marg_eff_fun = marg_eff_fun_hatch_doy))
  ),
  tar_target(
    dirfit_hatch_doy_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_hatch_doy_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Hatch date (day of year)"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_first_dna_age_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_sire",
                                     "beta_co_count_sire2",
                                     "beta_co_count_dam",
                                     "beta_co_count_dam2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_sire,
                                     avg_fun_co_count_sire2,
                                     avg_fun_co_count_dam,
                                     avg_fun_co_count_dam2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     avg_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     pred_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_first_dna_age,
                       marg_eff_fun = marg_eff_fun_first_dna_age))
  ),
  tar_target(
    dirfit_first_dna_age_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_first_dna_age_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Age at first DNA sampling (days)"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_co_count_parsum_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps_parsum,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_parsum",
                                     "beta_co_count_parsum2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     pred_fun_co_count_parsum,
                                     pred_fun_co_count_parsum2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     avg_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_co_count_parsum,
                       marg_eff_fun = marg_eff_fun_co_count_parsum))
  ),
  tar_target(
    dirfit_co_count_parsum_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_co_count_parsum_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Sum of ACCs in gametes from each parent"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_ps_age_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps_parsum,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_parsum",
                                     "beta_co_count_parsum2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_parsum,
                                     avg_fun_co_count_parsum2,
                                     pred_fun_age_q1,
                                     pred_fun_age_q2,
                                     avg_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_age,
                       marg_eff_fun = marg_eff_fun_age))
  ),
  tar_target(
    dirfit_ps_age_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_ps_age_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Age"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_ps_f_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps_parsum,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_parsum",
                                     "beta_co_count_parsum2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_parsum,
                                     avg_fun_co_count_parsum2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     pred_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_f,
                       marg_eff_fun = marg_eff_fun_f))
  ),
  tar_target(
    dirfit_ps_f_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_ps_f_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Inbreeding coefficient"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_ps_hatch_doy_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps_parsum,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_parsum",
                                     "beta_co_count_parsum2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_parsum,
                                     avg_fun_co_count_parsum2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     avg_fun_f,
                                     pred_fun_hatch_doy,
                                     pred_fun_hatch_doy2,
                                     avg_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_hatch_doy,
                       marg_eff_fun = marg_eff_fun_hatch_doy))
  ),
  tar_target(
    dirfit_ps_hatch_doy_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_ps_hatch_doy_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Hatch date (day of year)"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_ps_first_dna_age_pred_marg,
    pred_marg_func_dirfit(
      samp = dirfit_samps_parsum,
      data = fitness_data_dir,
      pred_info = list(coef_name = c("alpha",
                                     "beta_co_count_parsum",
                                     "beta_co_count_parsum2",
                                     "beta_age_q1",
                                     "beta_age_q2",
                                     "beta_f",
                                     "beta_hatch_doy",
                                     "beta_hatch_doy2",
                                     "beta_first_dna_age"),
                       action = list(avg_fun_alpha,
                                     avg_fun_co_count_parsum,
                                     avg_fun_co_count_parsum2,
                                     avg_fun_age_q1,
                                     avg_fun_age_q2,
                                     avg_fun_f,
                                     avg_fun_hatch_doy,
                                     avg_fun_hatch_doy2,
                                     pred_fun_first_dna_age),
                       x_axis_fun = x_axis_fun_first_dna_age,
                       marg_eff_fun = marg_eff_fun_first_dna_age))
  ),
  tar_target(
    dirfit_ps_first_dna_age_pred_plot,
    plot_lines_posterior(df = getElement(dirfit_ps_first_dna_age_pred_marg,
                                         "df_pred"),
                         xlab = paste0("Age at first DNA sampling (days)"),
                         ylab = paste0("Predicted ", trait_short),
                         title = toTitleCase(sex))
  ),
  tar_target(
    dirfit_fixed_eff_coefs_plot,
    plot_fixed_eff_coefs(
      samp_df = (dirfit_samps %>%
                   as.data.frame() %>%
                   select(matches("^(alph|bet)"))),
      subtit = paste0("S&D, ", sex, " ", trait_short))
  ),
  tar_target(
    dirfit_fixed_eff_coefs_plot_png,
    ggsave_path(paste0("figs/coef_figs/",
                       trait_short, "_S&D_", sex, ".png"),
                plot = dirfit_fixed_eff_coefs_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_rand_eff_sd_plot,
    plot_rand_eff_sd(
      samp_df = (dirfit_samps %>%
                   as.data.frame() %>%
                   select(matches("^sigma_"))),
      subtit = paste0("S&D, ", sex, " ", trait_short))
  ),
  tar_target(
    dirfit_rand_eff_sd_plot_png,
    ggsave_path(paste0("figs/sigma_figs/",
                       trait_short, "_S&D_", sex, ".png"),
                plot = dirfit_rand_eff_sd_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_year_levels_plot,
    plot_levels(
      samp_df = (dirfit_samps %>%
                   as.data.frame() %>%
                   select(matches("^(ye|hy)"))),
      tit = paste0(c("Y", "Hatch y"),
                   "ear effect levels (linear predictor scale)"),
      subtit = paste0("S&D, ", sex, " ", trait_short),
      dat = fitness_data_dir,
      samp_colname = c("ye", "hy"),
      num_col = c("y_num", "hy_num"),
      dat_col = c("year", "hatch_year"),
      isl = FALSE)
  ),
  tar_target(
    dirfit_year_levels_plot_png,
    ggsave_path(paste0("figs/year_level_figs/",
                       trait_short, "_S&D_", sex, ".png"),
                plot = dirfit_year_levels_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_island_levels_plot,
    plot_levels(
      samp_df = (dirfit_samps %>%
                   as.data.frame() %>%
                   select(matches("^(ll|hi)"))),
      tit = paste0(c("I", "Hatch i"),
                   "sland effect levels (linear predictor scale)"),
      subtit = paste0("S&D, ", sex, " ", trait_short),
      dat = fitness_data_dir,
      samp_colname = c("ll", "hi"),
      num_col = c("ll_num", "hi_num"),
      dat_col = c("last_locality", "hatch_island"),
      isl = TRUE)
  ),
  tar_target(
    dirfit_island_levels_plot_png,
    ggsave_path(paste0("figs/isl_level_figs/",
                       trait_short, "_S&D_", sex, ".png"),
                plot = dirfit_island_levels_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_ps_fixed_eff_coefs_plot,
    plot_fixed_eff_coefs(
      samp_df = (dirfit_samps_parsum %>%
                   as.data.frame() %>%
                   select(matches("^(alph|bet)"))),
      subtit = paste0("PS, ", sex, " ", trait_short))
  ),
  tar_target(
    dirfit_ps_fixed_eff_coefs_plot_png,
    ggsave_path(paste0("figs/coef_figs/",
                       trait_short, "_PS_", sex, ".png"),
                plot = dirfit_ps_fixed_eff_coefs_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_ps_rand_eff_sd_plot,
    plot_rand_eff_sd(
      samp_df = (dirfit_samps_parsum %>%
                   as.data.frame() %>%
                   select(matches("^sigma_"))),
      subtit = paste0("PS, ", sex, " ", trait_short))
  ),
  tar_target(
    dirfit_ps_rand_eff_sd_plot_png,
    ggsave_path(paste0("figs/sigma_figs/",
                       trait_short, "_PS_", sex, ".png"),
                plot = dirfit_ps_rand_eff_sd_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_ps_year_levels_plot,
    plot_levels(
      samp_df = (dirfit_samps_parsum %>%
                   as.data.frame() %>%
                   select(matches("^(ye|hy)"))),
      tit = paste0(c("Y", "Hatch y"),
                   "ear effect levels (linear predictor scale)"),
      subtit = paste0("PS, ", sex, " ", trait_short),
      dat = fitness_data_dir,
      samp_colname = c("ye", "hy"),
      num_col = c("y_num", "hy_num"),
      dat_col = c("year", "hatch_year"),
      isl = FALSE)
  ),
  tar_target(
    dirfit_ps_year_levels_plot_png,
    ggsave_path(paste0("figs/year_level_figs/",
                       trait_short, "_PS_", sex, ".png"),
                plot = dirfit_ps_year_levels_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_ps_island_levels_plot,
    plot_levels(
      samp_df = (dirfit_samps_parsum %>%
                   as.data.frame() %>%
                   select(matches("^(ll|hi)"))),
      tit = paste0(c("I", "Hatch i"),
                   "sland effect levels (linear predictor scale)"),
      subtit = paste0("PS, ", sex, " ", trait_short),
      dat = fitness_data_dir,
      samp_colname = c("ll", "hi"),
      num_col = c("ll_num", "hi_num"),
      dat_col = c("last_locality", "hatch_island"),
      isl = TRUE)
  ),
  tar_target(
    dirfit_ps_island_levels_plot_png,
    ggsave_path(paste0("figs/isl_level_figs/",
                       trait_short, "_PS_", sex, ".png"),
                plot = dirfit_ps_island_levels_plot,
                width = 9,
                height = 9,
                device = "png"),
    format = "file"
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
  dirfit_map,
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
  ),
  tar_combine(
    co_gp_cv_mean_acc_fitmod,
    fitmod_map[["co_gp_cv_mean_acc"]],
    command = dplyr::bind_rows(!!!.x, .id = "fitmod") %>%
      dplyr::mutate(fitmod = gsub(pattern = "co_gp_cv_mean_acc_",
                                  replacement = "",
                                  x = fitmod)),
    deployment = "main"
  )
)

list(
  sex_map,
  ########## Data files ##########
  tar_target(
    recomb_data_path2,
    co_data_rename_cols(
      paste0("data/20260427_Sparrow_YAPP/",
             "2_recsumm_Crossover_Count_per_individual_post_QC.txt")),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    geno_data_paths,
    paste0("data/20260427_Sparrow_57K_SUPER/20260427_Sparrow_57K_SUPER.",
           c("bed", "bim", "fam")),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    nestling_data_path,
    c("data/nestling_data_for_kenneth_20260707.csv"),
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
    paste0("data/combined_200k_70k_sparrow_genotype_data/pedigrees/",
           "2024-03-04_pedigree/helge_ped_err0.0027_adj_04-03-2024.txt"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    inla_effects_gp_vector_grm_all,
    c(##### Fixed effects:
      # Intercept
      "1",
      # Total coverage effect
      "total_coverage_scaled",
      # # Inbreeding coefficient
      "froh",
      ################# Random effects:
      # Year effect (iid random effect)
      "f(hatch_year, model = \"iid\", hyper = prior$hyperpar_var)",
      # Island effect (iid random effect)
      "f(first_locality, model = \"iid\", hyper = prior$hyperpar_var)",
      # Genetic values (using a GRM)
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
    froh_file,
    "data/20260121-FROH2.5_helgeland.txt",
    format = "file"
  ),
  ########## Figures ##########
  tar_target(
    stan_bv_pred_plot_ars2x2,
    layout_2x2_fig(p1 = stan_bv_pred_plot_ars_adult_f,
                   p2 = stan_bv_pred_plot_ars_adult_m,
                   p3 = stan_bv_pred_plot_ars_parent_f,
                   p4 = stan_bv_pred_plot_ars_parent_m,
                   tit_str = "Effect of ACC breeding value on ARS",
                   labs = paste(LETTERS[1:4],
                                c("GV-F", "GV-F", "PGV-F", "PGV-F"),
                                sep = " - "))
  ),
  tar_target(
    stan_bv_pred_plot_ars2x2_png,
    ggsave_path("figs/stan_bv_pred_plot_ars2x2.png",
                plot = stan_bv_pred_plot_ars2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_age_pred_plot_ars2x2,
    layout_2x2_fig(p1 = stan_age_pred_plot_ars_adult_f,
                   p2 = stan_age_pred_plot_ars_adult_m,
                   p3 = stan_age_pred_plot_ars_parent_f,
                   p4 = stan_age_pred_plot_ars_parent_m,
                   tit_str = "Effect of age on ARS")
  ),
  tar_target(
    stan_f_pred_plot_ars2x2,
    layout_2x2_fig(p1 = stan_f_pred_plot_ars_adult_f,
                   p2 = stan_f_pred_plot_ars_adult_m,
                   p3 = stan_f_pred_plot_ars_parent_f,
                   p4 = stan_f_pred_plot_ars_parent_m,
                   tit_str = "Effect of inbreeding on ARS")
  ),
  tar_target(
    stan_bv_pred_plot_surv2x2,
    layout_2x2_fig(p1 = stan_bv_pred_plot_surv_adult_f,
                   p2 = stan_bv_pred_plot_surv_adult_m,
                   p3 = stan_bv_pred_plot_surv_parent_f,
                   p4 = stan_bv_pred_plot_surv_parent_m,
                   tit_str = "Effect of ACC breeding value on AS",
                   labs = paste(LETTERS[1:4],
                                c("GV-F", "GV-F", "PGV-F", "PGV-F"),
                                sep = " - "))
  ),
  tar_target(
    stan_bv_pred_plot_surv2x2_png,
    ggsave_path("figs/stan_bv_pred_plot_surv2x2.png",
                plot = stan_bv_pred_plot_surv2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_age_pred_plot_surv2x2,
    layout_2x2_fig(p1 = stan_age_pred_plot_surv_adult_f,
                   p2 = stan_age_pred_plot_surv_adult_m,
                   p3 = stan_age_pred_plot_surv_parent_f,
                   p4 = stan_age_pred_plot_surv_parent_m,
                   tit_str = "Effect of age on AS")
  ),
  tar_target(
    stan_f_pred_plot_surv2x2,
    layout_2x2_fig(p1 = stan_f_pred_plot_surv_adult_f,
                   p2 = stan_f_pred_plot_surv_adult_m,
                   p3 = stan_f_pred_plot_surv_parent_f,
                   p4 = stan_f_pred_plot_surv_parent_m,
                   tit_str = "Effect of inbreeding on AS")
  ),
  tar_target(
    stan_bv_pred_plot_nest2x2,
    layout_2x2_fig(p1 = stan_bv_pred_plot_n2_f,
                   p2 = stan_bv_pred_plot_n2_m,
                   p3 = stan_bv_pred_plot_nest_f,
                   p4 = stan_bv_pred_plot_nest_m,
                   tit_str = "Effect of ACC genetic value on NS",
                   labs = paste(LETTERS[1:4],
                                c("GV-F", "GV-F", "PGV-F", "PGV-F"),
                                sep = " - "))
  ),
  tar_target(
    stan_bv_pred_plot_nest2x2_png,
    ggsave_path("figs/stan_bv_pred_plot_nest2x2.png",
                plot = stan_bv_pred_plot_nest2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_f_pred_plot_nest2x2,
    layout_2x2_fig(p1 = stan_f_pred_plot_n2_f,
                   p2 = stan_f_pred_plot_n2_m,
                   p3 = stan_f_pred_plot_nest_f,
                   p4 = stan_f_pred_plot_nest_m,
                   tit_str = "Inbreeding effect on NS")
  ),
  tar_target(
    dirfit_pred_plot_ars_3x2,
    layout_3x2_fig(p1 = dirfit_co_count_sire_pred_plot_ars_f,
                   p2 = dirfit_co_count_sire_pred_plot_ars_m,
                   p3 = dirfit_co_count_dam_pred_plot_ars_f,
                   p4 = dirfit_co_count_dam_pred_plot_ars_m,
                   p5 = dirfit_co_count_parsum_pred_plot_ars_f,
                   p6 = dirfit_co_count_parsum_pred_plot_ars_m,
                   tit_str = "S&D/PS models for ARS",
                   labs = paste(LETTERS[1:6],
                                c(rep("S&D", 4), "PS", "PS"),
                                sep = " - "))
  ),
  tar_target(
    dirfit_pred_plot_ars_3x2_png,
    ggsave_path("figs/dirfit_pred_plot_ars_3x2.png",
                plot = dirfit_pred_plot_ars_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_pred_plot_as_3x2,
    layout_3x2_fig(p1 = dirfit_co_count_sire_pred_plot_as_f,
                   p2 = dirfit_co_count_sire_pred_plot_as_m,
                   p3 = dirfit_co_count_dam_pred_plot_as_f,
                   p4 = dirfit_co_count_dam_pred_plot_as_m,
                   p5 = dirfit_co_count_parsum_pred_plot_as_f,
                   p6 = dirfit_co_count_parsum_pred_plot_as_m,
                   tit_str = "S&D/PS models for AS",
                   labs = paste(LETTERS[1:6],
                                c(rep("S&D", 4), "PS", "PS"),
                                sep = " - "))
  ),
  tar_target(
    dirfit_pred_plot_as_3x2_png,
    ggsave_path("figs/dirfit_pred_plot_as_3x2.png",
                plot = dirfit_pred_plot_as_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_pred_plot_ns_3x2,
    layout_3x2_fig(p1 = dirfit_co_count_sire_pred_plot_ns_f,
                   p2 = dirfit_co_count_sire_pred_plot_ns_m,
                   p3 = dirfit_co_count_dam_pred_plot_ns_f,
                   p4 = dirfit_co_count_dam_pred_plot_ns_m,
                   p5 = dirfit_co_count_parsum_pred_plot_ns_f,
                   p6 = dirfit_co_count_parsum_pred_plot_ns_m,
                   tit_str = "S&D/PS models for NS",
                   labs = paste(LETTERS[1:6],
                                c(rep("S&D", 4), "PS", "PS"),
                                sep = " - "))
  ),
  tar_target(
    dirfit_pred_plot_ns_3x2_png,
    ggsave_path("figs/dirfit_pred_plot_ns_3x2.png",
                plot = dirfit_pred_plot_ns_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_age_pred_plot_2x2,
    layout_2x2_fig(p1 = dirfit_age_pred_plot_ars_f,
                   p2 = dirfit_age_pred_plot_ars_m,
                   p3 = dirfit_age_pred_plot_as_f,
                   p4 = dirfit_age_pred_plot_as_m,
                   tit_str = "Effect of age in S&D models")
  ),
  tar_target(
    dirfit_age_pred_plot_2x2_png,
    ggsave_path("figs/dirfit_age_pred_plot_2x2.png",
                plot = dirfit_age_pred_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_f_pred_plot_3x2,
    layout_3x2_fig(p1 = dirfit_f_pred_plot_ars_f,
                   p2 = dirfit_f_pred_plot_ars_m,
                   p3 = dirfit_f_pred_plot_as_f,
                   p4 = dirfit_f_pred_plot_as_m,
                   p5 = dirfit_f_pred_plot_ns_f,
                   p6 = dirfit_f_pred_plot_ns_m,
                   tit_str = "Effect of inbreeding in S&D models")
  ),
  tar_target(
    dirfit_f_pred_plot_3x2_png,
    ggsave_path("figs/dirfit_f_pred_plot_3x2.png",
                plot = dirfit_f_pred_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_nfix_pred_plot_2x2,
    layout_2x2_fig(p1 = dirfit_hatch_doy_pred_plot_ns_f,
                   p2 = dirfit_hatch_doy_pred_plot_ns_m,
                   p3 = dirfit_first_dna_age_pred_plot_ns_f,
                   p4 = dirfit_first_dna_age_pred_plot_ns_m,
                   tit_str = "Nestling fixed effects in S&D models")
  ),
  tar_target(
    dirfit_nfix_pred_plot_2x2_png,
    ggsave_path("figs/dirfit_nfix_pred_plot_2x2.png",
                plot = dirfit_nfix_pred_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_ps_age_pred_plot_2x2,
    layout_2x2_fig(p1 = dirfit_ps_age_pred_plot_ars_f,
                   p2 = dirfit_ps_age_pred_plot_ars_m,
                   p3 = dirfit_ps_age_pred_plot_as_f,
                   p4 = dirfit_ps_age_pred_plot_as_m,
                   tit_str = "Effect of age in PS models")
  ),
  tar_target(
    dirfit_ps_age_pred_plot_2x2_png,
    ggsave_path("figs/dirfit_ps_age_pred_plot_2x2.png",
                plot = dirfit_ps_age_pred_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_ps_f_pred_plot_3x2,
    layout_3x2_fig(p1 = dirfit_ps_f_pred_plot_ars_f,
                   p2 = dirfit_ps_f_pred_plot_ars_m,
                   p3 = dirfit_ps_f_pred_plot_as_f,
                   p4 = dirfit_ps_f_pred_plot_as_m,
                   p5 = dirfit_ps_f_pred_plot_ns_f,
                   p6 = dirfit_ps_f_pred_plot_ns_m,
                   tit_str = "Effect of inbreeding in PS models")
  ),
  tar_target(
    dirfit_ps_f_pred_plot_3x2_png,
    ggsave_path("figs/dirfit_ps_f_pred_plot_3x2.png",
                plot = dirfit_ps_f_pred_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    dirfit_ps_nfix_pred_plot_2x2,
    layout_2x2_fig(p1 = dirfit_ps_hatch_doy_pred_plot_ns_f,
                   p2 = dirfit_ps_hatch_doy_pred_plot_ns_m,
                   p3 = dirfit_ps_first_dna_age_pred_plot_ns_f,
                   p4 = dirfit_ps_first_dna_age_pred_plot_ns_m,
                   tit_str = "Nestling fixed effects in PS models")
  ),
  tar_target(
    dirfit_ps_nfix_pred_plot_2x2_png,
    ggsave_path("figs/dirfit_ps_nfix_pred_plot_2x2.png",
                plot = dirfit_ps_nfix_pred_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_gvf_bv_plot_3x2,
    layout_3x2_fig(p1 = stan_bv_pred_plot_ars_adult_f,
                   p2 = stan_bv_pred_plot_ars_adult_m,
                   p3 = stan_bv_pred_plot_surv_adult_f,
                   p4 = stan_bv_pred_plot_surv_adult_m,
                   p5 = stan_bv_pred_plot_n2_f,
                   p6 = stan_bv_pred_plot_n2_m,
                   tit_str = "GV-F models")
  ),
  tar_target(
    stan_gvf_bv_plot_3x2_png,
    ggsave_path("figs/stan_gvf_bv_plot_3x2.png",
                plot = stan_gvf_bv_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_gvf_bv_marg_plot_3x2,
    layout_3x2_fig(p1 = stan_bv_marg_plot_ars_adult_f,
                   p2 = stan_bv_marg_plot_ars_adult_m,
                   p3 = stan_bv_marg_plot_surv_adult_f,
                   p4 = stan_bv_marg_plot_surv_adult_m,
                   p5 = stan_bv_marg_plot_n2_f,
                   p6 = stan_bv_marg_plot_n2_m,
                   tit_str = "GV-F models")
  ),
  tar_target(
    stan_gvf_bv_marg_plot_3x2_png,
    ggsave_path("figs/stan_gvf_bv_marg_plot_3x2.png",
                plot = stan_gvf_bv_marg_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_gvf_bv_out_vs_in_plot_3x2,
    layout_3x2_fig(p1 = stan_bv_out_vs_in_plot_ars_adult_f,
                   p2 = stan_bv_out_vs_in_plot_ars_adult_m,
                   p3 = stan_bv_out_vs_in_plot_surv_adult_f,
                   p4 = stan_bv_out_vs_in_plot_surv_adult_m,
                   p5 = stan_bv_out_vs_in_plot_n2_f,
                   p6 = stan_bv_out_vs_in_plot_n2_m,
                   tit_str = "GV-F models")
  ),
  tar_target(
    stan_gvf_bv_out_vs_in_plot_3x2_png,
    ggsave_path("figs/stan_gvf_bv_out_vs_in_plot_3x2.png",
                plot = stan_gvf_bv_out_vs_in_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_gvf_f_plot_3x2,
    layout_3x2_fig(p1 = stan_f_pred_plot_ars_adult_f,
                   p2 = stan_f_pred_plot_ars_adult_m,
                   p3 = stan_f_pred_plot_surv_adult_f,
                   p4 = stan_f_pred_plot_surv_adult_m,
                   p5 = stan_f_pred_plot_n2_f,
                   p6 = stan_f_pred_plot_n2_m,
                   tit_str = "GV-F models")
  ),
  tar_target(
    stan_gvf_f_plot_3x2_png,
    ggsave_path("figs/stan_gvf_f_plot_3x2.png",
                plot = stan_gvf_f_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_gvf_age_plot_2x2,
    layout_2x2_fig(p1 = stan_age_pred_plot_ars_adult_f,
                   p2 = stan_age_pred_plot_ars_adult_m,
                   p3 = stan_age_pred_plot_surv_adult_f,
                   p4 = stan_age_pred_plot_surv_adult_m,
                   tit_str = "GV-F models")
  ),
  tar_target(
    stan_gvf_age_plot_2x2_png,
    ggsave_path("figs/stan_gvf_age_plot_2x2.png",
                plot = stan_gvf_age_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_gvf_nfix_plot_2x2,
    layout_2x2_fig(p1 = stan_hatch_doy_pred_plot_n2_f,
                   p2 = stan_hatch_doy_pred_plot_n2_m,
                   p3 = stan_first_dna_age_pred_plot_n2_f,
                   p4 = stan_first_dna_age_pred_plot_n2_m,
                   tit_str = "GV-F models")
  ),
  tar_target(
    stan_gvf_nfix_plot_2x2_png,
    ggsave_path("figs/stan_gvf_nfix_plot_2x2.png",
                plot = stan_gvf_nfix_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_pgvf_bv_plot_3x2,
    layout_3x2_fig(p1 = stan_bv_pred_plot_ars_parent_f,
                   p2 = stan_bv_pred_plot_ars_parent_m,
                   p3 = stan_bv_pred_plot_surv_parent_f,
                   p4 = stan_bv_pred_plot_surv_parent_m,
                   p5 = stan_bv_pred_plot_nest_f,
                   p6 = stan_bv_pred_plot_nest_m,
                   tit_str = "PGV-F models")
  ),
  tar_target(
    stan_pgvf_bv_plot_3x2_png,
    ggsave_path("figs/stan_pgvf_bv_plot_3x2.png",
                plot = stan_pgvf_bv_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_pgvf_bv_marg_plot_3x2,
    layout_3x2_fig(p1 = stan_bv_marg_plot_ars_parent_f,
                   p2 = stan_bv_marg_plot_ars_parent_m,
                   p3 = stan_bv_marg_plot_surv_parent_f,
                   p4 = stan_bv_marg_plot_surv_parent_m,
                   p5 = stan_bv_marg_plot_nest_f,
                   p6 = stan_bv_marg_plot_nest_m,
                   tit_str = "PGV-F models")
  ),
  tar_target(
    stan_pgvf_bv_marg_plot_3x2_png,
    ggsave_path("figs/stan_pgvf_bv_marg_plot_3x2.png",
                plot = stan_pgvf_bv_marg_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_pgvf_bv_out_vs_in_plot_3x2,
    layout_3x2_fig(p1 = stan_bv_out_vs_in_plot_ars_parent_f,
                   p2 = stan_bv_out_vs_in_plot_ars_parent_m,
                   p3 = stan_bv_out_vs_in_plot_surv_parent_f,
                   p4 = stan_bv_out_vs_in_plot_surv_parent_m,
                   p5 = stan_bv_out_vs_in_plot_nest_f,
                   p6 = stan_bv_out_vs_in_plot_nest_m,
                   tit_str = "PGV-F models")
  ),
  tar_target(
    stan_pgvf_bv_out_vs_in_plot_3x2_png,
    ggsave_path("figs/stan_pgvf_bv_out_vs_in_plot_3x2.png",
                plot = stan_pgvf_bv_out_vs_in_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_pgvf_f_plot_3x2,
    layout_3x2_fig(p1 = stan_f_pred_plot_ars_parent_f,
                   p2 = stan_f_pred_plot_ars_parent_m,
                   p3 = stan_f_pred_plot_surv_parent_f,
                   p4 = stan_f_pred_plot_surv_parent_m,
                   p5 = stan_f_pred_plot_nest_f,
                   p6 = stan_f_pred_plot_nest_m,
                   tit_str = "PGV-F models")
  ),
  tar_target(
    stan_pgvf_f_plot_3x2_png,
    ggsave_path("figs/stan_pgvf_f_plot_3x2.png",
                plot = stan_pgvf_f_plot_3x2,
                width = 7.5,
                height = 9,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_pgvf_age_plot_2x2,
    layout_2x2_fig(p1 = stan_age_pred_plot_ars_parent_f,
                   p2 = stan_age_pred_plot_ars_parent_m,
                   p3 = stan_age_pred_plot_surv_parent_f,
                   p4 = stan_age_pred_plot_surv_parent_m,
                   tit_str = "PGV-F models")
  ),
  tar_target(
    stan_pgvf_age_plot_2x2_png,
    ggsave_path("figs/stan_pgvf_age_plot_2x2.png",
                plot = stan_pgvf_age_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_target(
    stan_pgvf_nfix_plot_2x2,
    layout_2x2_fig(p1 = stan_hatch_doy_pred_plot_nest_f,
                   p2 = stan_hatch_doy_pred_plot_nest_m,
                   p3 = stan_first_dna_age_pred_plot_nest_f,
                   p4 = stan_first_dna_age_pred_plot_nest_m,
                   tit_str = "PGV-F models")
  ),
  tar_target(
    stan_pgvf_nfix_plot_2x2_png,
    ggsave_path("figs/stan_pgvf_nfix_plot_2x2.png",
                plot = stan_pgvf_nfix_plot_2x2,
                width = 7.5,
                height = 6.5,
                device = "png"),
    format = "file"
  ),
  tar_combine(
    co_gp_cv_mean_acc_fitmod_sex,
    sex_map[["co_gp_cv_mean_acc_fitmod"]],
    command = dplyr::bind_rows(!!!.x, .id = "sex") %>%
      dplyr::mutate(sex = gsub(pattern = "co_gp_cv_mean_acc_fitmod_",
                               replacement = "",
                               x = sex)),
    deployment = "main"
  ),
  ############################## Permutation models ############################
  tar_target(
    perm_idx, # Indices of permutation branching
    seq_len(250),
    deployment = "main"
  ),
  tar_target(
    stan_data_perm,
    permute_data(dat = stan_data_surv_adult_m),
    pattern = map(perm_idx),
    deployment = "main"
  ),
  tar_target(
    stan_model_perm,
    stan(file = stan_file_surv_adult_m,
         data = c(stan_data_perm,
                  list(Y = getElement(stan_data_perm, "survival"))),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         thin = 1.6e2, # to keep final object reasonably small
         pars = c("beta_bv_std", "beta_bv2_std"),
         model_name = "stan_surv_adult_m_perm",
         control = list(adapt_delta = 0.96)),
    pattern = map(stan_data_perm)
  ),
  tar_target(
    # For both linear and quadratic coefficient, check if 0 is contained in the
    # 95% middle percentile (if not, consider it a false positive)
    # Acceptable range for 250 permutations: (0.024 to 0.08)
    stan_perm_falsepos_rate,
    stan_model_perm %>%
      lapply(FUN = function(x) {
        summary(x) %>%
          getElement("summary") %>%
          {((.)[1:2, "2.5%"] * (.)[1:2, "97.5%"]) > 0}
      }) %>%
      do.call(what = rbind) %>%
      colMeans(),
    deployment = "main"
  )
)
