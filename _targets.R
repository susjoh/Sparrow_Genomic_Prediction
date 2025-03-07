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
               "dplyr",
               "magrittr",
               "tools",
               "INLA",
               "rstan",
               "bayesplot",
               "tidyr"), # Packages that your targets need for their tasks.
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
    data_adult,
    make_data_adult(gp_data = co_data_adult_ars,
                    gp_model = co_gp_adult_ars,
                    lrs_data_path,
                    sex_num = "all",
                    inbreeding = inbreeding)
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
    stan_data_adult,
    make_stan_data_adult(data_adult),
    deployment = "main"
  ),
  tar_target(
    stan_data_adult_ss,
    make_stan_data_adult(data_adult_ss),
    deployment = "main"
  ),
  # tar_target(
  #   stan_adult_ars_3e4,
  #   stan(file = stan_file_adult_ars_zinf,
  #        data = c(stan_data_adult,
  #            list(Y = stan_data_adult$sum_recruit,
  #                 phi_inv_rate = mean(stan_data_adult$sum_recruit)^2 /
  #                   (var(stan_data_adult$sum_recruit) -
  #                      mean(stan_data_adult$sum_recruit)),
  #                 alpha_prior_mean = stan_data_adult$sum_recruit_log_mean,
  #                 beta_prior_sd = 0.05,
  #                 exp_rate = 40,
  #                 alpha_zi_prior_mean = 0,
  #                 beta_zi_prior_sd = 0.5,
  #                 exp_rate_zi = sqrt(1 / 0.5^2))),
  #        iter = 3.6e4,
  #        warmup = 6e3,
  #        chains = 16,
  #        cores = 16,
  #        control = list(adapt_delta = 0.9),
  #        model_name = paste0("stan_adult_ars_", sex_lc),
  #        thin = 1.2e2) # to keep final object reasonably small
  # ),
  tar_target(
    stan_adult_ars_3e4_ss,
    stan(file = stan_file_adult_ars_zinf,
         data = c(stan_data_adult_ss,
                  list(Y = stan_data_adult_ss$sum_recruit,
                       phi_inv_rate = mean(stan_data_adult_ss$sum_recruit)^2 /
                         (var(stan_data_adult_ss$sum_recruit) -
                            mean(stan_data_adult_ss$sum_recruit)),
                       alpha_prior_mean = stan_data_adult_ss$sum_recruit_log_mean,
                       beta_prior_sd = 0.05,
                       exp_rate = 40, # sqrt(1 / 0.05^2)
                       alpha_zi_prior_mean = 0,
                       beta_zi_prior_sd = 0.5,
                       exp_rate_zi = sqrt(1 / 0.5^2))),
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         pars = ars_pars,
         control = list(adapt_delta = 0.9),
         model_name = paste0("stan_adult_ars_ss_", sex_lc),
         thin = 1.2e2) # to keep final object reasonably small
  ),
  tar_target(
    ars_sim,
    lapply(1:3, function(i) make_ars_sim(stan_data_adult)),
    deployment = "main"
  ),
  # tar_target(
  #   stan_adult_ars_sim_3e4,
  #   stan(file = stan_file_adult_ars,
  #        data = c(stan_data_adult,
  #                 list(Y = ars_sim,
  #                      phi_inv_rate = mean(ars_sim)^2 /
  #                        (var(ars_sim) - mean(ars_sim)),
  #                      alpha_prior_mean = log(mean(ars_sim)),
  #                      beta_prior_sd = 0.05,
  #                      exp_rate = 40)),
  #        iter = 3.6e4,
  #        warmup = 6e3,
  #        chains = 16,
  #        cores = 16,
  #        control = list(adapt_delta = 0.9),
  #        pars = ars_pars,
  #        model_name = paste0("stan_adult_ars_", sex_lc, "_sim"),
  #        thin = 1.2e2), # to keep final object reasonably small
  #   pattern = map(ars_sim)
  # ),
  tar_target(
    stan_adult_surv_3e4_ss,
    stan(file = stan_file_adult_surv,
         data = c(stan_data_adult_ss,
                  list(Y = stan_data_adult_ss$survival,
                       alpha_prior_mean = stan_data_adult_ss$survival_logit_mean,
                       beta_prior_sd = 0.5,
                       exp_rate = sqrt(1 / 0.5^2))),
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         pars = surv_pars,
         model_name = paste0("stan_adult_surv_", sex_lc),
         thin = 1.2e2) # to keep final object reasonably small
  ),
  tar_target(
    surv_sim,
    lapply(1:3, function(i) make_surv_sim(stan_data_adult)),
    deployment = "main"
  ),
  # tar_target(
  #   stan_adult_surv_sim_3e4,
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
    get_samps(model = stan_adult_ars_3e4_ss,
              pars = ars_pars)
  ),
  # tar_target(
  #   ars_samp_pairs_plot,
  #   ggpairs(as.data.frame(ars_samps)[, c("energy", ars_pars[c(1:6, 17:20)])])
  # ),
  tar_target(
    surv_samps,
    get_samps(model = stan_adult_surv_3e4_ss,
              pars = surv_pars)
  ),
  # tar_target(
  #   surv_samp_pairs_plot,
  #   ggpairs(as.data.frame(surv_samps)[, c("energy", surv_pars[c(1:6, 18:21)])]),
  #   pattern = map(surv_models)
  # ),
  tar_target(
    ars_bv_preds_and_marg,
    bv_preds_and_marg(samps = ars_samps,
                      data = stan_data_adult_ss,
                      inv_link = exp)
  ),
  tar_target(
    surv_bv_preds_and_marg,
    bv_preds_and_marg(samps = surv_samps,
                      data = stan_data_adult_ss,
                      inv_link = function(x) 1 / (1 + exp(-x)))
  ),
  tar_target(
    ars_bv_pred_plot,
    plot_lines_posterior(df = ars_bv_preds_and_marg$df_pred,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ",
                                       sex ,
                                       " annual reproductive success"),
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
    ars_bv_marg_plot,
    plot_lines_posterior(df = ars_bv_preds_and_marg$df_marg,
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Marginal effect on ",
                                       sex,
                                       " annual reproductive success"),
                         title = "")
  ),
  tar_target(
    surv_bv_pred_plot,
    plot_lines_posterior(df = getElement(surv_bv_preds_and_marg, "df_pred"),
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Predicted ", sex," annual survival"),
                         title = "")
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
    surv_bv_marg_plot,
    plot_lines_posterior(df = getElement(surv_bv_preds_and_marg, "df_marg"),
                         xlab = paste0("Breeding value for ",
                                       sex,
                                       " crossover count"),
                         ylab = paste0("Marginal effect on ",
                                       sex,
                                       " annual survival"),
                         title = "")
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
    stan_file_adult_ars_zinf,
    "r/zinf_ars.stan",
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
      "res",
      "bv_lat",
      "sigma_ll",
      "sigma_ye",
      "sigma_id",
      "sigma_res",
      "y_rep")
  )
)
