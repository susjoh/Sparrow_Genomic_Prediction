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
  workers = 16,
  seconds_idle = 120,
  tasks_max = 1,
  launch_max = 5,
  script_lines = paste("#SBATCH --account=share-nv-bio \n module load",
                       "R/4.3.2-gfbf-2023a R-bundle-CRAN/2023.12-foss-2023a",
                       "CMake/3.26.3-GCCcore-12.3.0"),
  slurm_log_output = "Jobs/%A.log",
  slurm_memory_gigabytes_per_cpu = 6,
  slurm_cpus_per_task = 16,
  slurm_time_minutes = 60 * 24 * 1, # minutes * hours * days
  slurm_partition = "CPUQ",
  verbose = TRUE
)

# Set target options:
tar_option_set(
  packages = c("tibble",
               "data.table",
               "stats",
               "parallel",
               "ggplot2",
               "GGally",
               "dplyr",
               "qs",
               "magrittr",
               "tools",
               "INLA",
               "rstan",
               "bayesplot"), # Packages that your targets need for their tasks.
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

values_sex <- tibble(
  sex = c("female", "male"),
  sex_lc = c("f", "m"),
  sex_keep = c("F", "M")
)

sex_map <- tar_map(
  values = values_sex,
  names = "sex_lc",
  tar_target(
    co_data, # All measurements of co count
    prep_co_data(recomb_data_path,
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
                intra_shuff = mean(intra_shuff),
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
                intra_shuff = intra_shuff[sample(n(), 1)],
                total_coverage = total_coverage[sample(n(), 1)],
                n = n())
  ),
  tar_target(
    num_meas_vs_co_count,
    ggplot(data = co_data_rand, aes(x = co_count, y = n)) +
      geom_point() +
      geom_smooth(method = "lm")
  ),
  tar_target(
    co_grm_files, # GRM file for inds with sex-specific crossover measurement
    make_grm(analysis_inds = unique(co_data$ringnr), # Vector of inds.
             bfile = gsub(x = geno_data_paths[1], ".bed", ""),
             ncores = 4,
             mem = 4 * 6000,
             # Path to plink program:
             plink_path = plink_path,
             # Where to save result:
             dir = paste0("data/co_", sex_lc, "_grm")),
    format = "file"
  ),
  tar_target(
    co_grm_obj, # GRM file for inds with sex-specific crossover measurement
    load_grm(dir = gsub(x = co_grm_files[1],  "/grm.rel.bin",  ""),
             pheno_data = co_data)
  ),
  tar_target(
    co_gam, # genomic animal model for sex-specific crossover rate
    run_gp(pheno_data = co_data,
           inverse_relatedness_matrix = co_grm_obj$inv_grm,
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc))
  ),
  tar_target(
    co_10fcv, # genomic animal model for sex-specific crossover rate
    run_gp(pheno_data = co_data_cv,
           inverse_relatedness_matrix = co_grm_obj$inv_grm,
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc, "_test")),
    pattern = map(co_data_cv)
  ),
  tar_target(
    co_10fcv_bvacc, # genomic animal model for sex-specific crossover rate
    get_10fcv_acc(data = co_data_cv,
                  model = co_10fcv,
                  y_str = paste0("co_count_", sex_lc)),
    pattern = map(co_data_cv, co_10fcv),
    iteration = "vector"
  ),
  tar_target(
    co_gam_cv, # genomic animal model for sex-specific crossover rate
    inla_cv(model = co_gam,
            pheno_data = co_data)
  ),
  tar_target(
    co_data_adult_ars,
    make_adult_ars_gp_data(pheno_data = co_data,
                           lrs_path = lrs_data_path,
                           fam_path = geno_data_paths[3])
  ),
  tar_target(
    co_adult_ars_grm_files,
    make_grm(analysis_inds = unique(co_data_adult_ars$ringnr),
             bfile = gsub(x = geno_data_paths[1], ".bed", ""),
             ncores = 4,
             mem = 4 * 6000,
             # Path to plink program:
             plink_path = plink_path,
             # Where to save result:
             dir = paste0("data/co_adult_ars_", sex_lc, "_grm")),
    format = "file"
  ),
  tar_target(
    co_adult_ars_grm_obj, # GRM file for females with crossover measurement
    load_grm(dir = gsub(x = co_adult_ars_grm_files[1],  "/grm.rel.bin",  ""),
             pheno_data = co_data_adult_ars)
  ),
  tar_target(
    co_gp_adult_ars, # genomic animal model for crossover rate in females
    run_gp(pheno_data = co_data_adult_ars,
           inverse_relatedness_matrix = co_adult_ars_grm_obj$inv_grm,
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc))
  ),
  tar_target(
    data_adult_ars,
    make_data_adult_ars(gp_data = co_data_adult_ars,
                        gp_model = co_gp_adult_ars,
                        lrs_data_path,
                        sex_lc = sex_lc,
                        inbreeding = inbreeding)
  ),
  tar_target(
    data_adult_ars_bvpost,
    make_data_adult_ars_bvpost(gp_data = co_data_adult_ars,
                               gp_model = co_gp_adult_ars,
                               lrs_data_path,
                               sex_lc = sex_lc,
                               inbreeding = inbreeding)
  ),
  tar_target(
    bv_std_vec,
    sapply(X = 1:1e3,
           FUN = function(i, N, means, sds, id_idx) rnorm(N, means, sds)[id_idx],
           N = nrow(data_adult_ars_bvpost$lrs),
           means = data_adult_ars_bvpost$lrs$bv_mean,
           sds = data_adult_ars_bvpost$lrs$bv_sd,
           id_idx = data_adult_ars_bvpost$lrs$ringnr_num)
  ),
  tar_target(
    stan_data_adult_ars_bvpost,
    data_adult_ars_bvpost$lrs$sum_recruit %>%
      list(N = nrow(data_adult_ars_bvpost$lrs),
           N_ll = max(data_adult_ars_bvpost$lrs$ll_num),
           N_ye = max(data_adult_ars_bvpost$lrs$y_num),
           N_id = max(data_adult_ars_bvpost$lrs$ringnr_num),
           ye_idx = data_adult_ars_bvpost$lrs$y_num,
           ll_idx = data_adult_ars_bvpost$lrs$ll_num,
           id_idx = data_adult_ars_bvpost$lrs$ringnr_num,
           bv_mean = unique(data_adult_ars_bvpost$lrs$bv_mean[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           bv_sd = unique(data_adult_ars_bvpost$lrs$bv_sd[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           Y = .,
           age = data_adult_ars_bvpost$lrs$age,
           age_q1 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 1],
           age_q2 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 2],
           f = data_adult_ars_bvpost$lrs$fhat3,
           phi_inv_rate = pc_rate(mean(.)^2 / (var(.) - mean(.))),
           exp_rate = 40,
           bv_mean_std = mean(bv_std_vec),
           bv_sd_std = sd(bv_std_vec),
           alpha_prior_mean = log(mean(.)),
           beta_prior_sd = 0.05),
    deployment = "main"
  ),
  tar_target(
    stan_adult_ars_bvpost_3e4,
    stan(file = stan_file_adult_ars_bvpost,
         data = stan_data_adult_ars_bvpost,
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         pars = c("alpha",
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
                  "phi",
                  "y_rep"),
         model_name = paste0("stan_adult_ars_bvpost_", sex_lc),
         thin = 1.2e2) # to keep final object reasonably small
  ),
  tar_rep(
    stan_data_adult_ars_bvpost_sim,
    ars_bvpost_sim(stan_data_adult_ars_bvpost) %>%
      list(N = nrow(data_adult_ars_bvpost$lrs),
           N_ll = max(data_adult_ars_bvpost$lrs$ll_num),
           N_ye = max(data_adult_ars_bvpost$lrs$y_num),
           N_id = max(data_adult_ars_bvpost$lrs$ringnr_num),
           ye_idx = data_adult_ars_bvpost$lrs$y_num,
           ll_idx = data_adult_ars_bvpost$lrs$ll_num,
           id_idx = data_adult_ars_bvpost$lrs$ringnr_num,
           bv_mean = unique(data_adult_ars_bvpost$lrs$bv_mean[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           bv_sd = unique(data_adult_ars_bvpost$lrs$bv_sd[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           Y = .,
           age = data_adult_ars_bvpost$lrs$age,
           age_q1 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 1],
           age_q2 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 2],
           f = data_adult_ars_bvpost$lrs$fhat3,
           phi_inv_rate = pc_rate(mean(.)^2 / (var(.) - mean(.))),
           exp_rate = 40,
           bv_mean_std = mean(bv_std_vec),
           bv_sd_std = sd(bv_std_vec),
           alpha_prior_mean = log(mean(.)),
           beta_prior_sd = 0.05),
    batches = 3,
    deployment = "main"
  ),
  tar_target(
    stan_adult_ars_bvpost_sim_3e4,
    stan(file = stan_file_adult_ars_bvpost,
         data = `[[`(stan_data_adult_ars_bvpost_sim, 1),
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         pars = c("alpha",
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
                  "phi",
                  "y_rep"),
         model_name = paste0("stan_adult_ars_bvpost_", sex_lc, "_sim"),
         thin = 1.2e2), # to keep final object reasonably small
    pattern = map(stan_data_adult_ars_bvpost_sim)
  ),
  tar_rep(
    stan_data_adult_surv_bvpost_sim,
    surv_bvpost_sim(stan_data_adult_surv_bvpost) %>%
      list(N = nrow(data_adult_ars_bvpost$lrs),
           N_ll = max(data_adult_ars_bvpost$lrs$ll_num),
           N_ye = max(data_adult_ars_bvpost$lrs$y_num),
           N_id = max(data_adult_ars_bvpost$lrs$ringnr_num),
           ye_idx = data_adult_ars_bvpost$lrs$y_num,
           ll_idx = data_adult_ars_bvpost$lrs$ll_num,
           id_idx = data_adult_ars_bvpost$lrs$ringnr_num,
           bv_mean = unique(data_adult_ars_bvpost$lrs$bv_mean[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           bv_sd = unique(data_adult_ars_bvpost$lrs$bv_sd[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           Y = .,
           age = data_adult_ars_bvpost$lrs$age,
           age_q1 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 1],
           age_q2 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 2],
           f = data_adult_ars_bvpost$lrs$fhat3,
           exp_rate = 40,
           bv_mean_std = mean(bv_std_vec),
           bv_sd_std = sd(bv_std_vec),
           alpha_prior_mean = log(mean(.) / (1 - mean(.))),
           beta_prior_sd = 0.05),
    batches = 3,
    deployment = "main"
  ),
  tar_target(
    stan_adult_surv_bvpost_sim_3e4,
    stan(file = stan_file_adult_surv_bvpost,
         data = `[[`(stan_data_adult_surv_bvpost_sim, 1),
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         pars = c("alpha",
                  "beta_bv",
                  "beta_bv2",
                  "beta_age_q1",
                  "beta_age_q2",
                  "beta_f",
                  "ye",
                  "ll",
                  "id",
                  "res",
                  "bv_lat",
                  "sigma_ll",
                  "sigma_ye",
                  "sigma_id",
                  "sigma_res",
                  "y_rep"),
         model_name = paste0("stan_adult_surv_bvpost_", sex_lc, "_sim"),
         thin = 1.2e2), # to keep final object reasonably small
    pattern = map(stan_data_adult_surv_bvpost_sim)
  ),
  tar_target(
    stan_data_adult_surv_bvpost,
    data_adult_ars$lrs$survival %>%
      list(N = nrow(data_adult_ars_bvpost$lrs),
           N_ll = max(data_adult_ars_bvpost$lrs$ll_num),
           N_ye = max(data_adult_ars_bvpost$lrs$y_num),
           N_id = max(data_adult_ars_bvpost$lrs$ringnr_num),
           ye_idx = data_adult_ars_bvpost$lrs$y_num,
           ll_idx = data_adult_ars_bvpost$lrs$ll_num,
           id_idx = data_adult_ars_bvpost$lrs$ringnr_num,
           bv_mean = unique(data_adult_ars_bvpost$lrs$bv_mean[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           bv_sd = unique(data_adult_ars_bvpost$lrs$bv_sd[
             order(data_adult_ars_bvpost$lrs$ringnr_num)]),
           Y = .,
           age = data_adult_ars_bvpost$lrs$age,
           age_q1 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 1],
           age_q2 = poly(data_adult_ars_bvpost$lrs$age, degree = 2)[, 2],
           f = data_adult_ars_bvpost$lrs$fhat3,
           exp_rate = 40,
           bv_mean_std = mean(bv_std_vec),
           bv_sd_std = sd(bv_std_vec),
           alpha_prior_mean = log(mean(.) / (1 - mean(.))),
           beta_prior_sd = 0.05),
    deployment = "main"
  ),
  tar_target(
    stan_adult_surv_bvpost_3e4,
    stan(file = stan_file_adult_surv_bvpost,
         data = stan_data_adult_surv_bvpost,
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         pars = c("alpha",
                  "beta_bv",
                  "beta_bv2",
                  "beta_age_q1",
                  "beta_age_q2",
                  "beta_f",
                  "ye",
                  "ll",
                  "id",
                  "res",
                  "bv_lat",
                  "sigma_ll",
                  "sigma_ye",
                  "sigma_id",
                  "sigma_res",
                  "y_rep"),
         model_name = paste0("stan_adult_surv_bvpost_", sex_lc),
         thin = 1.2e2) # to keep final object reasonably small
  ),
  tar_target(
    ars_samps,
    get_ars_samps(model = stan_adult_ars_bvpost_3e4)
  ),
  tar_target(
    ars_samp_pairs_plot,
    ggpairs(as.data.frame(ars_samps)[, c(
      "energy",
      "sigma_ll",
      "sigma_ye",
      "sigma_id",
      "phi",
      "alpha",
      "beta_bv",
      "beta_bv2",
      "beta_age_q1",
      "beta_age_q2",
      "beta_f",
      "ye.2",
      "ll.2",
      "id.3",
      "bv_lat.245"
    )])
  ),
  tar_target(
    surv_samps,
    get_surv_samps(model = stan_adult_surv_bvpost_3e4)
  ),
  tar_target(
    surv_samp_pairs_plot,
    ggpairs(as.data.frame(surv_samps)[, c(
      "energy",
      "sigma_ll",
      "sigma_ye",
      "sigma_id",
      "sigma_res",
      "alpha",
      "beta_bv",
      "beta_bv2",
      "beta_age_q1",
      "beta_age_q2",
      "beta_f",
      "ye.2",
      "ll.2",
      "id.3",
      "bv_lat.245",
      "res.2314"
    )])
  ),
  tar_target(
    ars_bv_preds_and_marg,
    bv_preds_and_marg(samps = ars_samps,
                      data = stan_data_adult_ars_bvpost,
                      inv_link = exp)
  ),
  tar_target(
    surv_bv_preds_and_marg,
    bv_preds_and_marg(samps = surv_samps,
                      data = stan_data_adult_surv_bvpost,
                      inv_link = function(x) 1 / (1 + exp(-x)))
  ),
  tar_target(
    ars_bv_pred_plot,
    plot_lines_posterior(df = ars_bv_preds_and_marg$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = "Predicted annual recruits (both sexes)",
                         title = "")
  ),
  tar_target(
    ars_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/ars_bv_pred_", sex_lc, ".pdf"),
                plot = ars_bv_pred_plot,
                width = 7,
                height = 7,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    ars_bv_marg_plot,
    plot_lines_posterior(df = ars_bv_preds_and_marg$df_marg,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = "Marginal effect on annual recruits (both sexes)",
                         title = "")
  ),
  tar_target(
    surv_bv_pred_plot,
    plot_lines_posterior(df = surv_bv_preds_and_marg$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = "Predicted annual survival (both sexes)",
                         title = "")
  ),
  tar_target(
    surv_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/surv_bv_pred_", sex_lc, ".pdf"),
                plot = surv_bv_pred_plot,
                width = 7,
                height = 7,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    surv_bv_marg_plot,
    plot_lines_posterior(df = surv_bv_preds_and_marg$df_marg,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = "Marginal effect on annual survival (both sexes)",
                         title = "")
  ),
  tar_target(
    ars_ppc_mean,
    ppc_stat(y = stan_data_adult_ars_bvpost$Y, yrep = ars_samps$y_rep),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_sd,
    ppc_stat(y = stan_data_adult_ars_bvpost$Y, yrep = ars_samps$y_rep, stat = "sd"),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_zeros,
    ppc_stat(y = stan_data_adult_ars_bvpost$Y, yrep = ars_samps$y_rep, stat = function(y) mean(y == 0)),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_ones,
    ppc_stat(y = stan_data_adult_ars_bvpost$Y, yrep = ars_samps$y_rep, stat = function(y) mean(y == 1)),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_twos,
    ppc_stat(y = stan_data_adult_ars_bvpost$Y, yrep = ars_samps$y_rep, stat = function(y) mean(y == 2)),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_bar,
    ppc_bars(y = stan_data_adult_ars_bvpost$Y, yrep = ars_samps$y_rep),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_bar,
    ppc_bars(y = stan_data_adult_surv_bvpost$Y, yrep = surv_samps$y_rep),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_sd,
    ppc_stat(y = stan_data_adult_surv_bvpost$Y, yrep = surv_samps$y_rep, stat = "sd"),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_bar_ll,
    ppc_bars_grouped(y = stan_data_adult_surv_bvpost$Y, yrep = surv_samps$y_rep, group = stan_data_adult_surv_bvpost$ll_idx),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_sd_ll,
    ppc_stat_grouped(y = stan_data_adult_surv_bvpost$Y, yrep = surv_samps$y_rep, stat = "sd", group = stan_data_adult_surv_bvpost$ll_idx),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_bar_ye,
    ppc_bars_grouped(y = stan_data_adult_surv_bvpost$Y, yrep = surv_samps$y_rep, group = stan_data_adult_surv_bvpost$ye_idx),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_sd_ye,
    ppc_stat_grouped(y = stan_data_adult_surv_bvpost$Y, yrep = surv_samps$y_rep, stat = "sd", group = stan_data_adult_surv_bvpost$ye_idx),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_p_mean,
    mean(colMeans(ars_samps$y_rep) > mean(stan_data_adult_ars_bvpost$Y)),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_p_sd,
    mean(apply(ars_samps$y_rep, 2, sd) > sd(stan_data_adult_ars_bvpost$Y)),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_p_zeros,
    mean(apply(ars_samps$y_rep, 2, function(y) mean(y == 0)) > mean(stan_data_adult_ars_bvpost$Y == 0)),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_p_ones,
    mean(apply(ars_samps$y_rep, 2, function(y) mean(y == 1)) > mean(stan_data_adult_ars_bvpost$Y == 1)),
    deployment = "main"
  ),
  tar_target(
    ars_ppc_p_twos,
    mean(apply(ars_samps$y_rep, 2, function(y) mean(y == 2)) > mean(stan_data_adult_ars_bvpost$Y == 2)),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_p_mean,
    mean(colMeans(surv_samps$y_rep) > mean(stan_data_adult_surv_bvpost$Y)),
    deployment = "main"
  ),
  tar_target(
    surv_ppc_p_sd,
    mean(apply(surv_samps$y_rep, 2, sd) > sd(stan_data_adult_surv_bvpost$Y)),
    deployment = "main"
  )
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
    lrs_data_path2,
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
    stan_file_adult_ars_bvpost,
    "r/adult_ars_bvpost.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_surv_bvpost,
    "r/adult_surv_bvpost.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    inbreeding,
    find_inbreeding(gsub(x = geno_data_paths[1], ".bed", ""),
                    ncores = 4,
                    mem = 4 * 6000,
                    # Path to plink program:
                    plink_path = plink_path),
    format = "file"
  )
)
