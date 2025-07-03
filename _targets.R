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
               "bayesplot",
               "tidyr",
               "mvtnorm"), # Packages that your targets need for their tasks.
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
  sex_keep = c("F", "M"),
  sex_num_lrs = c(2, 1)
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
  tar_target(
    co_data_adult_ars,
    make_adult_ars_gp_data(pheno_data = co_data,
                           lrs_path = lrs_data_path,
                           fam_path = geno_data_paths[3],
                           sex_num = sex_num_lrs)
  ),
  tar_target(
    co_data_parent_gp,
    make_parent_gp_data(pheno_data = co_data,
                        lrs_path = lrs_data_path,
                        lrs_path2 = lrs_data_path2,
                        fam_path = geno_data_paths[3],
                        ped_path = pedigree_path,
                        sex_keep = sex_keep)
  ),
  tar_target(
    co_data_nestling_gp,
    make_nestling_gp_data(pheno_data = co_data,
                          lrs_path = lrs_data_path,
                          lrs_path2 = lrs_data_path2,
                          nestling_path = nestling_data_path,
                          fam_path = geno_data_paths[3],
                          ped_path = pedigree_path,
                          sex_keep = sex_keep)
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
    co_parent_grm_files,
    make_grm(analysis_inds = unique(co_data_parent_gp$ringnr),
             bfile = gsub(x = geno_data_paths[1], ".bed", ""),
             ncores = 4,
             mem = 4 * 6000,
             # Path to plink program:
             plink_path = plink_path,
             # Where to save result:
             dir = paste0("data/co_parent_", sex_lc, "_grm")),
    format = "file"
  ),
  tar_target(
    co_nestling_grm_files,
    make_grm(analysis_inds = unique(co_data_nestling_gp$ringnr),
             bfile = gsub(x = geno_data_paths[1], ".bed", ""),
             ncores = 4,
             mem = 4 * 6000,
             # Path to plink program:
             plink_path = plink_path,
             # Where to save result:
             dir = paste0("data/co_nestling_", sex_lc, "_grm")),
    format = "file"
  ),
  tar_target(
    co_adult_ars_grm_obj,
    load_grm(dir = gsub(x = co_adult_ars_grm_files[1], "/grm.rel.bin",  ""),
             pheno_data = co_data_adult_ars)
  ),
  tar_target(
    co_parent_grm_obj,
    load_grm(dir = gsub(x = co_parent_grm_files[1], "/grm.rel.bin",  ""),
             pheno_data = co_data_parent_gp)
  ),
  tar_target(
    co_nestling_grm_obj,
    load_grm(dir = gsub(x = co_nestling_grm_files[1], "/grm.rel.bin",  ""),
             pheno_data = co_data_nestling_gp)
  ),
  tar_target(
    co_gp_adult_ars,
    run_gp(pheno_data = co_data_adult_ars,
           inverse_relatedness_matrix = co_adult_ars_grm_obj$inv_grm,
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc))
  ),
  tar_target(
    co_gp_adult_ars2,
    run_gp(pheno_data = co_data_adult_ars,
           inverse_relatedness_matrix = co_adult_ars_grm_obj$inv_grm,
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc),
           comp_conf = TRUE)
  ),
  tar_target(
    co_gp_parent,
    run_gp(pheno_data = co_data_parent_gp,
           inverse_relatedness_matrix = co_parent_grm_obj$inv_grm,
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc),
           comp_conf = TRUE)
  ),
  tar_target(
    co_gp_nestling,
    run_gp(pheno_data = co_data_nestling_gp,
           inverse_relatedness_matrix = co_nestling_grm_obj$inv_grm,
           effects_vec = inla_effects_gp_vector_grm_all,
           y = paste0("co_count_", sex_lc),
           comp_conf = TRUE)
  ),
  tar_target(
    bv_covmat,
    inla_bv_covmat(model = co_gp_adult_ars2,
                   ncores = 16)
  ),
  tar_target(
    bv_covmat_parent,
    inla_bv_covmat(model = co_gp_parent,
                   ncores = 16)
  ),
  tar_target(
    bv_covmat_nestling,
    inla_bv_covmat(model = co_gp_nestling,
                   ncores = 16)
  ),
  tar_target(
    data_adult_ss, # Single sex
    make_data_adult(gp_data = co_data_adult_ars,
                    gp_model = co_gp_adult_ars,
                    lrs_data_path,
                    sex_num = sex_num_lrs,
                    inbreeding = inbreeding)
  ),
  tar_target(
    data_parent_adult_ss,
    make_data_parent_adult(gp_data = co_data_parent_gp,
                           gp_model = co_gp_parent,
                           lrs_data_path,
                           sex_num = sex_num_lrs,
                           inbreeding = inbreeding,
                           ped_path = pedigree_path)
  ),
  tar_target(
    data_nestling,
    make_data_nestling(gp_data = co_data_nestling_gp,
                       gp_model = co_gp_nestling,
                       nestling_data_path = nestling_data_path,
                       sex_num = sex_num_lrs,
                       inbreeding = inbreeding,
                       ped_path = pedigree_path)
  ),
  tar_target(
    # We generally find more extreme estimated bvs with the more measurements
    bv_vs_n_meas_plot,
    dplyr::filter(data_adult_ss, !duplicated(ringnr)) %>%
      ggplot(aes(y = abs(bv_mean), x = ifelse(is.na(co_n), 0, co_n))) +
      geom_point() +
      labs(y = "Absolute estimated breeding value",
           x = "Number of crossover count measurements") +
      geom_smooth(method = "loess")
  ),
  tar_target(
    stan_data_adult_ss,
    make_stan_data_adult(data_adult_ss),
    deployment = "main"
  ),
  tar_target(
    stan_data_adult_ss_covmat,
    make_stan_data_adult_covmat(data_adult_ss, co_data_adult_ars, bv_covmat),
    deployment = "main"
  ),
  tar_target(
    stan_data_parent_adult_ss_covmat,
    make_stan_data_parent_covmat(data = data_parent_adult_ss,
                                 gp_data = co_data_parent_gp,
                                 covmat = bv_covmat_parent),
    deployment = "main"
  ),
  tar_target(
    stan_data_nestling_covmat,
    make_stan_data_nestling_covmat(data = data_nestling,
                                   gp_data = co_data_nestling_gp,
                                   covmat = bv_covmat_nestling),
    deployment = "main"
  ),
  tar_target(
    stan_adult_ars_ss,
    stan(file = stan_file_adult_ars_zinf,
         data = c(stan_data_adult_ss,
                  list(Y = stan_data_adult_ss$sum_recruit,
                       alpha_prior_mean = log(mean(c(stan_data_adult_ss$sum_recruit[stan_data_adult_ss$sum_recruit != 0],
                                                     stan_data_adult_ss$sum_recruit))),
                       beta_prior_sd = 0.2,
                       exp_rate = 1 / 0.2,
                       alpha_zi_prior_mean = log(mean(stan_data_adult_ss$sum_recruit == 0) /
                                                   (1 - mean(stan_data_adult_ss$sum_recruit == 0))),
                       beta_zi_prior_sd = 0.5,
                       exp_rate_zi = 1 / 0.5)),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         pars = ars_pars,
         control = list(adapt_delta = 0.9),
         model_name = paste0("stan_adult_ars_ss_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  tar_target(
    stan_adult_ars_ss_covmat,
    stan(file = stan_file_adult_ars_zinf_covmat,
         data = c(stan_data_adult_ss_covmat,
                  list(Y = stan_data_adult_ss_covmat$sum_recruit,
                       alpha_prior_mean = log(mean(c(stan_data_adult_ss_covmat$sum_recruit[stan_data_adult_ss_covmat$sum_recruit != 0],
                                                     stan_data_adult_ss_covmat$sum_recruit))),
                       beta_prior_sd = 0.2,
                       exp_rate = 1 / 0.2,
                       alpha_zi_prior_mean = log(mean(stan_data_adult_ss_covmat$sum_recruit == 0) /
                                                   (1 - mean(stan_data_adult_ss_covmat$sum_recruit == 0))),
                       beta_zi_prior_sd = 0.5,
                       exp_rate_zi = 1 / 0.5)),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         pars = ars_pars,
         control = list(adapt_delta = 0.9),
         model_name = paste0("stan_adult_ars_ss_covmat_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  tar_target(
    stan_parent_adult_ars_ss_covmat,
    stan(file = stan_file_parent_ars_zinf_covmat,
         data = c(stan_data_parent_adult_ss_covmat,
                  list(Y = stan_data_parent_adult_ss_covmat$sum_recruit,
                       alpha_prior_mean = log(mean(c(stan_data_parent_adult_ss_covmat$sum_recruit[stan_data_parent_adult_ss_covmat$sum_recruit != 0],
                                                     stan_data_parent_adult_ss_covmat$sum_recruit))),
                       beta_prior_sd = 0.2,
                       exp_rate = 1 / 0.2,
                       alpha_zi_prior_mean = log(mean(stan_data_parent_adult_ss_covmat$sum_recruit == 0) /
                                                   (1 - mean(stan_data_parent_adult_ss_covmat$sum_recruit == 0))),
                       beta_zi_prior_sd = 0.5,
                       exp_rate_zi = 1 / 0.5)),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         pars = ars_pars,
         control = list(adapt_delta = 0.9),
         model_name = paste0("stan_parent_adult_ars_ss_covmat_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  tar_target(
    ars_sim,
    lapply(1:3, function(i) make_ars_sim(stan_data_adult_ss)),
    deployment = "main"
  ),
  tar_target(
    stan_adult_surv_ss,
    stan(file = stan_file_adult_surv,
         data = c(stan_data_adult_ss,
                  list(Y = stan_data_adult_ss$survival,
                       alpha_prior_mean = stan_data_adult_ss$survival_logit_mean,
                       beta_prior_sd = 0.5,
                       exp_rate = sqrt(1 / 0.5^2))),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.95),
         pars = surv_pars,
         model_name = paste0("stan_adult_surv_ss_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  tar_target(
    stan_adult_surv_ss_covmat,
    stan(file = stan_file_adult_surv_covmat,
         data = c(stan_data_adult_ss_covmat,
                  list(Y = stan_data_adult_ss_covmat$survival,
                       alpha_prior_mean = stan_data_adult_ss_covmat$survival_logit_mean,
                       beta_prior_sd = 0.5,
                       exp_rate = sqrt(1 / 0.5^2))),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.95),
         pars = surv_pars,
         model_name = paste0("stan_adult_surv_ss_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  # Also with co_n covariate
  tar_target(
    stan_adult_surv_ss_co_n,
    stan(file = stan_file_adult_surv_co_n,
         data = c(stan_data_adult_ss,
                  list(Y = stan_data_adult_ss$survival,
                       alpha_prior_mean = stan_data_adult_ss$survival_logit_mean,
                       beta_prior_sd = 0.5,
                       exp_rate = sqrt(1 / 0.5^2))),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.95),
         pars = c(surv_pars, "beta_co_n", "beta_co_n_std"),
         model_name = paste0("stan_adult_surv_ss_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  tar_target(
    stan_adult_surv_ss_covmat_co_n,
    stan(file = stan_file_adult_surv_covmat_co_n,
         data = c(stan_data_adult_ss_covmat,
                  list(Y = stan_data_adult_ss_covmat$survival,
                       alpha_prior_mean = stan_data_adult_ss_covmat$survival_logit_mean,
                       beta_prior_sd = 0.5,
                       exp_rate = sqrt(1 / 0.5^2))),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.95),
         pars = c(surv_pars, "beta_co_n", "beta_co_n_std"),
         model_name = paste0("stan_adult_surv_ss_covmat_co_n_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  tar_target(
    stan_parent_adult_surv_ss_covmat,
    stan(file = stan_file_parent_surv_covmat,
         data = c(stan_data_parent_adult_ss_covmat,
                  list(Y = stan_data_parent_adult_ss_covmat$survival,
                       alpha_prior_mean = stan_data_parent_adult_ss_covmat$survival_logit_mean,
                       beta_prior_sd = 0.5,
                       exp_rate = sqrt(1 / 0.5^2))),
         iter = 4.8e5,
         warmup = 8e4,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         pars = surv_pars,
         model_name = paste0("stan_parent_adult_surv_ss_covmat_", sex_lc),
         thin = 1.6e3) # to keep final object reasonably small
  ),
  tar_target(
    stan_nestling_surv_covmat,
    stan(file = stan_file_nestling_surv_covmat,
         data = c(stan_data_nestling_covmat,
                  list(Y = stan_data_nestling_covmat$recruit,
                       alpha_prior_mean = stan_data_nestling_covmat$recruit_logit_mean,
                       beta_prior_sd = 0.5,
                       exp_rate = sqrt(1 / 0.5^2))),
         iter = 4.8e4,
         warmup = 8e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.95),
         pars = nest_pars,
         model_name = paste0("stan_nestling_surv_", sex_lc),
         thin = 1.6e2) # to keep final object reasonably small
  ),
  tar_target(
    surv_sim,
    lapply(1:3, function(i) make_surv_sim(stan_data_adult_ss)),
    deployment = "main"
  ),
  # tar_target(
  #   stan_adult_surv_sim,
  #   stan(file = stan_file_adult_surv,
  #        data = c(stan_data_adult,
  #                 list(Y = surv_sim,
  #                      alpha_prior_mean = log(1 / (1 / mean(surv_sim) - 1)),
  #                      beta_prior_sd = 0.5,
  #                      exp_rate = sqrt(1 / 0.5^2))),
  #        iter = 3.6e4,
  #        warmup = 6e3,
  #        chains = 16,
  #        cores = 16,
  #        control = list(adapt_delta = 0.9),
  #        pars = surv_pars,
  #        model_name = paste0("stan_adult_surv_", sex_lc, "_sim"),
  #        thin = 1.2e2), # to keep final object reasonably small
  #   pattern = map(surv_sim)
  # ),
  tar_target(
    ars_samps,
    get_samps(model = stan_adult_ars_ss,
              pars = ars_pars)
  ),
  tar_target(
    ars_samps_covmat,
    get_samps(model = stan_adult_ars_ss_covmat,
              pars = ars_pars)
  ),
  tar_target(
    ars_samps_parent,
    get_samps(model = stan_parent_adult_ars_ss_covmat,
              pars = ars_pars)
  ),
  # tar_target(
  #   ars_samp_pairs_plot,
  #   ggpairs(as.data.frame(ars_samps)[, c("energy", ars_pars[c(1:6, 17:20)])])
  # ),
  tar_target(
    surv_samps,
    get_samps(model = stan_adult_surv_ss,
              pars = surv_pars)
  ),
  tar_target(
    surv_samps_covmat,
    get_samps(model = stan_adult_surv_ss_covmat,
              pars = surv_pars)
  ),
  tar_target(
    surv_samps_covmat_co_n,
    get_samps(model = stan_adult_surv_ss_covmat_co_n,
              pars = c(surv_pars, "beta_co_n", "beta_co_n_std"))
  ),
  tar_target(
    surv_samps_parent,
    get_samps(model = stan_parent_adult_surv_ss_covmat,
              pars = surv_pars)
  ),
  tar_target(
    nest_samps,
    get_samps(model = stan_nestling_surv_covmat,
              pars = nest_pars)
  ),
  tar_target(
    surv_samp_pairs_plot,
    as.data.frame(surv_samps)[, c("energy", surv_pars[c(1:6, 17:19)])] %>%
      dplyr::mutate(sigma_ll = log(sigma_ll),
                    sigma_ye = log(sigma_ye),
                    sigma_id = log(sigma_id)) %>%
      ggpairs() # aes(color = factor(surv_samps_f$divergent)))
  ),
  tar_target(
    ars_bv_preds_and_marg,
    make_ars_bv_preds_and_marg(samps = ars_samps,
                               data = stan_data_adult_ss)
  ),
  tar_target(
    surv_bv_preds_and_marg,
    make_surv_bv_preds_and_marg(samps = surv_samps,
                                data = stan_data_adult_ss)
  ),
  tar_target(
    ars_parent_bv_preds_and_marg,
    make_ars_bv_preds_and_marg(samps = ars_samps_parent,
                               data = stan_data_parent_adult_ss_covmat)
  ),
  tar_target(
    surv_parent_bv_preds_and_marg,
    make_surv_bv_preds_and_marg(samps = surv_samps_parent,
                                data = stan_data_parent_adult_ss_covmat)
  ),
  tar_target(
    ars_covmat_bv_preds_and_margt,
    make_ars_bv_preds_and_marg(samps = ars_samps_covmat,
                               data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    surv_covmat_bv_preds_and_marg,
    make_surv_bv_preds_and_marg(samps = surv_samps_covmat,
                                data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    surv_covmat_co_n_bv_preds_and_marg,
    make_surv_bv_preds_and_marg_co_n(samps = surv_samps_covmat_co_n,
                                     data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    nest_bv_preds_and_marg,
    make_nest_bv_preds_and_marg(samps = nest_samps,
                                data = stan_data_nestling_covmat)
  ),
  tar_target(
    ars_bv_pred_plot,
    plot_lines_posterior(df = ars_bv_preds_and_marg$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ",
                                       sex,
                                       " annual reproductive success"),
                         title = "",
                         data = stan_data_adult_ss)
  ),
  tar_target(
    ars_bv_marg_plot,
    plot_lines_posterior(df = ars_bv_preds_and_marg$df_marg_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Marginal effect on ",
                                       sex,
                                       " annual reproductive success"),
                         title = "",
                         data = stan_data_adult_ss)
  ),
  tar_target(
    ars_parent_bv_pred_plot,
    plot_lines_posterior(df = ars_parent_bv_preds_and_marg$df_pred,
                         xlab = paste0("Parental breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ",
                                       sex,
                                       " annual reproductive success"),
                         title = "")
  ),
  tar_target(
    ars_covmat_bv_pred_plot,
    plot_lines_posterior(df = ars_covmat_bv_preds_and_margt$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ",
                                       sex,
                                       " annual reproductive success"),
                         title = "",
                         data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    ars_covmat_bv_marg_plot,
    plot_lines_posterior(df = ars_covmat_bv_preds_and_margt$df_marg_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ",
                                       sex,
                                       " annual reproductive success"),
                         title = "",
                         data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    surv_bv_pred_plot,
    plot_lines_posterior(df = surv_bv_preds_and_marg$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex," annual survival"),
                         title = "",
                         data = stan_data_adult_ss)
  ),
  tar_target(
    surv_bv_marg_plot,
    plot_lines_posterior(df = surv_bv_preds_and_marg$df_marg,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Marginal effect on ",
                                       sex,
                                       " annual survival"),
                         title = "",
                         data = stan_data_adult_ss)
  ),
  tar_target(
    surv_parent_bv_pred_plot,
    plot_lines_posterior(df = surv_parent_bv_preds_and_marg$df_pred,
                         xlab = paste0("Parental breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex," annual survival"),
                         title = "")
  ),
  tar_target(
    surv_covmat_bv_pred_plot,
    plot_lines_posterior(df = surv_covmat_bv_preds_and_marg$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex," annual survival"),
                         title = "",
                         data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    surv_covmat_bv_marg_plot,
    plot_lines_posterior(df = surv_covmat_bv_preds_and_marg$df_marg,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex," annual survival"),
                         title = "",
                         data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    surv_covmat_co_n_bv_pred_plot,
    plot_lines_posterior(df = surv_covmat_co_n_bv_preds_and_marg$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex," annual survival"),
                         title = "",
                         data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    surv_covmat_co_n_bv_marg_plot,
    plot_lines_posterior(df = surv_covmat_co_n_bv_preds_and_marg$df_marg,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex," annual survival"),
                         title = "",
                         data = stan_data_adult_ss_covmat)
  ),
  tar_target(
    nest_bv_pred_plot,
    plot_lines_posterior(df = nest_bv_preds_and_marg$df_pred,
                         xlab = paste0("Parental breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted nestling survival"),
                         title = "")
  ),
  tar_target(
    nest_bv_marg_plot,
    plot_lines_posterior(df = nest_bv_preds_and_marg$df_marg,
                         xlab = paste0("Parental breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted nestling survival"),
                         title = "")
  ),
  tar_target(
    ars_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/ars_bv_pred_", sex_lc, ".pdf"),
                plot = ars_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    ars_parent_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/ars_parent_bv_pred_", sex_lc, ".pdf"),
                plot = ars_parent_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    ars_covmat_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/ars_covmat_bv_pred_", sex_lc, ".pdf"),
                plot = ars_covmat_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    surv_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/surv_bv_pred_", sex_lc, ".pdf"),
                plot = surv_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    surv_parent_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/surv_parent_bv_pred_", sex_lc, ".pdf"),
                plot = surv_parent_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    surv_covmat_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/surv_covmat_bv_pred_", sex_lc, ".pdf"),
                plot = surv_covmat_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    surv_covmat_co_n_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/surv_covmat_co_n_bv_pred_", sex_lc, ".pdf"),
                plot = surv_covmat_co_n_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    nest_bv_pred_plot_pdf,
    ggsave_path(paste0("figs/nest_bv_pred_", sex_lc, ".pdf"),
                plot = nest_bv_pred_plot,
                width = 7,
                height = 5,
                device = "pdf"),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    ars_ppc,
    do_ars_ppc(y = stan_data_adult_ss$sum_recruit,
               yrep = ars_samps$y_rep,
               ll_i = stan_data_adult_ss$ll_idx,
               ye_i = stan_data_adult_ss$ye_idx,
               id_i = stan_data_adult_ss$id_idx)
  ),
  tar_target(
    surv_ppc,
    do_surv_ppc(y = stan_data_adult_ss$survival,
                yrep = surv_samps$y_rep,
                # co_n = replace_na(data_adult_ss$co_n, 0),
                # co_meas = data_adult_ss$co_meas,
                ll_i = stan_data_adult_ss$ll_idx,
                ye_i = stan_data_adult_ss$ye_idx,
                id_i = stan_data_adult_ss$id_idx)
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
    stan_file_adult_ars,
    "r/adult_ars.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_surv,
    "r/adult_surv.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_surv_covmat,
    "r/adult_surv_covmat.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_parent_surv_covmat,
    "r/parent_surv_covmat.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_nestling_surv_covmat,
    "r/nestling_surv_covmat.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_surv_covmat_co_n,
    "r/adult_surv_covmat_co_n.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_surv_co_n,
    "r/adult_surv_co_n.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_ars_zinf,
    "r/zinf_ars.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_ars_zinf_covmat,
    "r/zinf_ars_covmat.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_parent_ars_zinf_covmat,
    "r/parent_zinf_ars_covmat.stan",
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
  ),
  tar_target(
    ars_pars,
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
      "alpha_zi",
      "beta_zi_bv",
      "beta_zi_bv2",
      "beta_zi_age_q1",
      "beta_zi_age_q2",
      "beta_zi_f",
      "alpha_zi_std",
      "beta_zi_bv_std",
      "beta_zi_bv2_std",
      "beta_zi_age_q1_std",
      "beta_zi_age_q2_std",
      "beta_zi_f_std",
      "ye_zi",
      "ll_zi",
      "id_zi",
      # "res_zi",
      "bv_lat",
      "sigma_ll_zi",
      "sigma_ye_zi",
      "sigma_id_zi",
      # "sigma_res_zi",
      # "phi_inv",
      # "phi",
      "theta",
      "y_rep")
  ),
  tar_target(
    surv_pars,
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
      # "res",
      "bv_lat",
      "sigma_ll",
      "sigma_ye",
      "sigma_id",
      # "sigma_res",
      "y_rep")
  ),
  tar_target(
    nest_pars,
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
      "id",
      # "res",
      "bv_lat",
      "sigma_hi",
      "sigma_hy",
      "sigma_id",
      # "sigma_res",
      "y_rep")
  )
)
