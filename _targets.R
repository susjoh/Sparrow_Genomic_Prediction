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
  workers = 3*2*2,
  seconds_idle = 120,
  tasks_max = 1,
  launch_max = 5,
  script_lines = paste("#SBATCH --account=share-nv-bio \n module load",
                       "R/4.3.2-gfbf-2023a R-bundle-CRAN/2023.12-foss-2023a",
                       "CMake/3.26.3-GCCcore-12.3.0"),
  slurm_log_output = "Jobs/%A.log",
  slurm_memory_gigabytes_per_cpu = 6,
  slurm_cpus_per_task = 4,
  slurm_time_minutes = 60 * 24 * 0.5, # minutes * hours * days
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
               "INLA",
               "rstan"), # Packages that your targets need for their tasks.
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
tar_source(files = "r/gamma.parms.from.quantiles.R")

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
  # tar_target(
  #   berkson_error_sd_est1, # simply the differences of variances
  #   sqrt(var(data_adult_ars$lrs$pred_bv[data_adult_ars$lrs$co_meas]) -
  #          var(data_adult_ars$lrs$pred_bv[!data_adult_ars$lrs$co_meas]))
  # ),
  # tar_target(
  #   berkson_error_sd_est2, # using the accuracy, and fixing sigma_w^2
  #   sqrt(var(data_adult_ars$lrs$pred_bv[!data_adult_ars$lrs$co_meas]) *
  #          (1 / colMeans(co_10fcv_bvacc)^2 - 1))
  # ),
  # tar_target(
  #   berkson_error_sd_est, # using the accuracy, and fixing sigma_w^2+sigma_u^2
  #   sqrt(var(data_adult_ars$lrs$pred_bv[data_adult_ars$lrs$co_meas]) *
  #          (1 - colMeans(co_10fcv_bvacc)^2))
  # ),
  # tar_target(
  #   berkson_error_sd_priorpars,
  #   gamma.parms.from.quantiles(c(berkson_error_sd_est[3],
  #                                berkson_error_sd_est[2]))
  # ),
  tar_target(
    stan_data_adult_ars_bvpost_hyp,
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
           f = data_adult_ars_bvpost$lrs$fhat3,
           phi_inv_rate = pc_rate(mean(.)^2 / (var(.) - mean(.))),
           exp_rate = 10,
           hyperprior_shape = 2,
           hyperprior_rate = 20,
           quad_scale = sqrt(2),
           alpha_prior_mean = log(mean(.)))
  ),
  tar_target(
    stan_adult_ars_bvpost_hyp_3e4,
    stan(file = stan_file_adult_ars_bvpost_hyp,
         data = stan_data_adult_ars_bvpost_hyp,
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         model_name = paste0("stan_adult_ars_bvpost_hyp_", sex_lc),
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
           f = data_adult_ars_bvpost$lrs$fhat3,
           phi_inv_rate = pc_rate(mean(.)^2 / (var(.) - mean(.))),
           exp_rate = pc_rate(U = 0.12),
           hyperprior_shape = 2,
           hyperprior_rate = 40,
           quad_scale = 2),
    batches = 3
  ),
  tar_target(
    stan_adult_ars_bvpost_sim_1e4,
    stan(file = stan_file_adult_ars_bvpost_hyp,
         data = `[[`(stan_data_adult_ars_bvpost_sim, 1),
         iter = 1.2e4,
         warmup = 2e3,
         chains = 16,
         cores = 16,
         model_name = paste0("stan_adult_ars_", sex_lc),
         thin = 4e1), # to keep final object reasonably small
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
           f = data_adult_ars_bvpost$lrs$fhat3,
           exp_rate = 10,
           hyperprior_shape = 2,
           hyperprior_rate = 20,
           quad_scale = sqrt(2),
           alpha_prior_mean = log(mean(.) / (1 - mean(.)))),
    batches = 3
  ),
  tar_target(
    stan_adult_surv_bvpost_sim_1e4,
    stan(file = stan_file_adult_surv_bvpost_hyp,
         data = `[[`(stan_data_adult_surv_bvpost_sim, 1),
         iter = 1.2e4,
         warmup = 2e3,
         chains = 16,
         cores = 16,
         model_name = paste0("stan_adult_surv_bvpost_hyp_", sex_lc),
         thin = 4e1), # to keep final object reasonably small
    pattern = map(stan_data_adult_surv_bvpost_sim)
  ),
  tar_target(
    stan_data_adult_surv_bvpost_hyp,
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
           f = data_adult_ars_bvpost$lrs$fhat3,
           exp_rate = 10,
           hyperprior_shape = 2,
           hyperprior_rate = 20,
           quad_scale = sqrt(2),
           alpha_prior_mean = log(mean(.) / (1 - mean(.))))
  ),
  tar_target(
    stan_adult_surv_bvpost_hyp_3e4,
    stan(file = stan_file_adult_surv_bvpost_hyp,
         data = stan_data_adult_surv_bvpost_hyp,
         iter = 3.6e4,
         warmup = 6e3,
         chains = 16,
         cores = 16,
         control = list(adapt_delta = 0.9),
         model_name = paste0("stan_adult_surv_bvpost_hyp_", sex_lc),
         thin = 1.2e2) # to keep final object reasonably small
  ),
  tar_target(
    ars_samps,
    get_ars_samps(model = stan_adult_ars_bvpost_hyp_3e4)
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
      "beta_age",
      "beta_age2",
      "beta_f",
      "beta_prior_sd",
      "ye.2",
      "ll.2",
      "id.3",
      "bv_lat.245"
    )])
  ),
  tar_target(
    surv_samps,
    get_surv_samps(model = stan_adult_surv_bvpost_hyp_3e4)
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
      "beta_age",
      "beta_age2",
      "beta_f",
      "beta_prior_sd",
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
                      data = stan_data_adult_ars_bvpost_hyp,
                      inv_link = exp)
  ),
  tar_target(
    surv_bv_preds_and_marg,
    bv_preds_and_marg(samps = surv_samps,
                      data = stan_data_adult_surv_bvpost_hyp,
                      inv_link = function(x) 1 / (1 + exp(-x)))
  ),
  tar_target(
    ars_bv_pred_plot,
    plot_lines_posterior(df = ars_bv_preds_and_marg$df_pred,
                         xlab = "Breeding value",
                         ylab = "Predicted annual recruits",
                         title = "")
  ),
  tar_target(
    ars_bv_pred_plot_pdf,
    ggsave_path("figs/ars_bv_pred_plot.pdf",
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
                         xlab = "Breeding value",
                         ylab = "Marginal effect on annual recruits",
                         title = "")
  ),
  tar_target(
    surv_bv_pred_plot,
    plot_lines_posterior(df = surv_bv_preds_and_marg$df_pred,
                         xlab = "Breeding value",
                         ylab = "Predicted annual survival",
                         title = "")
  ),
  tar_target(
    surv_bv_pred_plot_pdf,
    ggsave_path("figs/surv_bv_pred_plot.pdf",
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
                         xlab = "Breeding value",
                         ylab = "Marginal effect on annual survival",
                         title = "")
  )
)

list(
  sex_map,
  tar_target(
    recomb_data_path,
    "data/20240910_Sparrow_Recomb_Data.txt",
    format = "file"
  ),
  tar_target(
    geno_data_paths,
    paste0("data/70K_200K_maf_geno_mind_v5.", c("bed", "bim", "fam")),
    format = "file"
  ),
  tar_target(
    nestling_data_path,
    "data/Nestlingdata_20240312_fix.csv",
    format = "file"
  ),
  tar_target(
    lrs_data_path,
    "data/LRS_data_20230803_Helgeland_fix.csv",
    format = "file"
  ),
  tar_target(
    lrs_data_path2,
    "data/Missing_LRS_Sparrows_revised_WithInfo_fix.csv",
    format = "file"
  ),
  tar_target(
    morph_data_path,
    "data/AdultMorphology_20240201_fix.csv",
    format = "file"
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
    )
  ),
  tar_target(
    stan_file_adult_ars_bvpost_hyp,
    "r/adult_ars_bvpost_hyp.stan",
    format = "file",
    deployment = "main"
  ),
  tar_target(
    stan_file_adult_surv_bvpost_hyp,
    "r/adult_surv_bvpost_hyp.stan",
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
