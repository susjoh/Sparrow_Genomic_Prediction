make_cv_test_sets <- function(analysis_inds,
                              num_folds = 10) {
  n <- length(analysis_inds)

  # Split individuals into n folds
  scrambled <- sample(n)
  folds <- cut(seq_len(n), num_folds)
  levels(folds) <- seq_len(num_folds)

  # Create list of test individuals in each fold
  lapply(1:num_folds, function(fo) {
    x <- analysis_inds[sort(scrambled[folds == fo])]
    attributes(x) <- list(fo = fo)
    x
  })
}

make_co_data_cv <- function(co_data,
                            test_set,
                            sex_lc) {
  co_data$test <- co_data$ringnr %in% test_set
  co_data[, paste0("co_count_", sex_lc, "_test")] <-
    co_data[, paste0("co_count_", sex_lc)]
  co_data[co_data$test, paste0("co_count_", sex_lc, "_test")] <- NA
  co_data
}

prep_co_data <- function(recomb_data_path,
                         lrs_data_path,
                         lrs_data_path2,
                         sex_keep) {
  ############## Data Wrangling ##################

  # Load data
  pheno_data_all <- fread(recomb_data_path, data.table = FALSE)

  pheno_data_all %<>% filter(., !is.na(co_count))

  # These two files contain all our inds with co_count measurements
  lrs <- fread(file = lrs_data_path)
  lrs2 <- fread(file = lrs_data_path2)
  pheno_data_all$id_red <- gsub(pattern = "_.+",
                                replacement = "",
                                x = pheno_data_all$id)
  missing_lrs <- with(pheno_data_all, id_red[which(!id_red %in% lrs$ringnr)])
  missing_lrs2 <- with(pheno_data_all, id_red[which(!id_red %in% lrs2$ringnr)])
  missing_all <- intersect(missing_lrs, missing_lrs2) # none missing in both

  # Extract hatch year, and first locality
  pheno_data_all$hatch_year <-
    lrs$hatch_year[match(pheno_data_all$id_red, lrs$ringnr)]
  pheno_data_all$hatch_year[is.na(pheno_data_all$hatch_year)] <-
    na.omit(lrs2$hatch_year[match(pheno_data_all$id_red, lrs2$ringnr)])

  pheno_data_all$first_locality <-
    lrs$first_locality[match(pheno_data_all$id_red, lrs$ringnr)]
  pheno_data_all$first_locality[is.na(pheno_data_all$first_locality)] <-
    na.omit(lrs2$first_locality[match(pheno_data_all$id_red, lrs2$ringnr)])


  pheno_data_mean <- pheno_data_all %>%
    group_by(id, sex) %>%
    summarise(co_count = mean(co_count),
              intra_shuff = mean(intra_shuff),
              total_coverage = mean(total_coverage),
              hatch_year = mean(hatch_year),
              first_locality = mean(first_locality),
              n = n())

  pheno_data_all$n <-
    pheno_data_mean$n[match(pheno_data_all$id, pheno_data_mean$id)]

  pheno_data <- pheno_data_all

  pheno_data$hatch_year <- as.factor(pheno_data$hatch_year)
  pheno_data$first_locality <- as.factor(pheno_data$first_locality)

  # Make columns for sex-specific crossover rates
  pheno_data$co_count_m <- pheno_data$co_count_f <- pheno_data$co_count
  pheno_data$co_count_m <- ifelse(pheno_data$sex == "M",
                                  pheno_data$co_count_m,
                                  NA)
  pheno_data$co_count_f <- ifelse(pheno_data$sex == "M",
                                  NA,
                                  pheno_data$co_count_f)

  pheno_data <- dplyr::filter(pheno_data, sex == sex_keep)

  # Subset to be able to test models quickly
  # n_obs_subset <- 700
  # set.seed(1)
  # pheno_data <- pheno_data[sort(sample(dim(pheno_data)[1], n_obs_subset)), ]

  # Center/scale these effects (INLA can't handle the huge numbers)
  pheno_data$total_coverage_scaled <- scale(pheno_data$total_coverage)
  # pheno_data$total_coverage2_scaled <- scale(pheno_data$total_coverage2)

  # Numerical ids
  ids_num <- match(pheno_data$id, unique(pheno_data$id))
  # Numerical ids, with only one entry for each individual
  ids_num_unique <- match(unique(pheno_data$id), pheno_data$id)

  pheno_data$ringnr <- pheno_data$id # just an alternate name

  # Numerical id columns needed for random effects in INLA
  pheno_data$id1 <- # For breeding values
    pheno_data$id2 <- # For ID effect
    pheno_data$id3 <- # For meas. error
    ids_num # numerical ids

  pheno_data
}

make_adult_ars_gp_data <- function(pheno_data,
                                   lrs_path,
                                   fam_path) {

  geno_inds <- fread(file = fam_path, select = 2)$V2
  geno_inds_red <- gsub(x = geno_inds, pattern = "_.+", "")
  dupes <- duplicated(geno_inds_red)

  lrs <- fread(file = lrs_path)
  # remove last loc. outside of main islands
  lrs <- lrs[lrs$last_locality %in%
               c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332),
             c("ringnr", "sex", "hatch_year", "first_locality")]
  # Only need a single row per ind now
  lrs <- lrs[!duplicated(lrs$ringnr), ]
  # Remove non-genotyped inds
  lrs <- lrs[lrs$ringnr %in% geno_inds_red, ]
  # For duplicated genotypes:
  lrs$id_red <- lrs$ringnr
  lrs$ringnr[which(!lrs$ringnr %in% geno_inds)] <-
    sapply(which(!lrs$ringnr %in% geno_inds),
           function(ind) {
             geno_inds[match(lrs$ringnr[ind], geno_inds_red)]
           })
  lrs$id <- lrs$ringnr
  # remove inds for which co_count is measured
  lrs <- lrs[!lrs$ringnr %in% pheno_data$id, ]
  # Recasting
  lrs$sex <- ifelse(lrs$sex > 1.5, "F", "M")
  lrs$hatch_year <- as.factor(lrs$hatch_year)
  lrs$id1 <- lrs$id2 <- lrs$id3 <- seq(from = max(pheno_data$id1) + 1,
                                       length = nrow(lrs),
                                       by = 1)
  # Merge
  plyr::rbind.fill(pheno_data, lrs)
}

make_data_adult_ars <- function(gp_data,
                                gp_model,
                                lrs_data_path,
                                sex_lc,
                                inbreeding) {

  lrs <- fread(file = lrs_data_path)
  lrs <- lrs[lrs$last_locality %in%
               c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332)]
  lrs <- lrs[lrs$ringnr %in% gp_data$id_red]

  # add predicted phenotypes and breeding values from GP model to data
  gp_data$pred_pheno <- gp_model$summary.fitted.values$mean
  gp_data$pred_bv <- gp_model$summary.random$id1$mean[
    order(gp_model$summary.random$id1$ID)][
      match(gp_data$id, unique(gp_data$id))]

  lrs$pred_bv <- gp_data$pred_bv[match(lrs$ringnr, gp_data$id_red)]
  lrs$age <- lrs$year - lrs$hatch_year
  # Remove inds with missing sex
  lrs <- lrs[!is.na(lrs$sex), ]
  lrs$ringnr_num <- match(lrs$ringnr, unique(lrs$ringnr))
  lrs$ll_num <- match(lrs$last_locality, unique(lrs$last_locality))
  lrs$y_num <- match(lrs$year, unique(lrs$year))
  lrs$hy_num <- match(lrs$hatch_year, unique(lrs$hatch_year))

  # Add inbreeding info
  ibc <- fread(file = inbreeding)
  gp_data$fhat3 <- ibc$Fhat3[match(gp_data$ringnr, ibc$IID)]
  lrs$fhat3 <- gp_data$fhat3[match(lrs$ringnr, gp_data$id_red)]

  # Number of crossover measurements
  lrs$co_n <- gp_data$n[match(lrs$ringnr, gp_data$id_red)]
  lrs$co_meas <- !is.na(lrs$co_n)
  lrs <- lrs[order(lrs$co_meas), ]

  # Add crossover count measurements, and phenotypic predictions, with repeats
  ph_dat <- lapply(1:nrow(lrs),
                   function(i) {
                     lrs_row <- lrs[i, ]
                     co_red <- gp_data[gp_data$id_red == lrs_row$ringnr, ]
                     lrs_df <- lrs_row[rep(1, nrow(co_red)), ]
                     cbind(lrs_df, co_red[, c(paste0("co_count_", sex_lc),
                                              "pred_pheno")])
                   }) %>%
    do.call(what = "rbind", .)

  lst(lrs, ph_dat)
}

make_data_adult_ars_bvpost <- function(gp_data,
                                       gp_model,
                                       lrs_data_path,
                                       sex_lc,
                                       inbreeding) {

  lrs <- fread(file = lrs_data_path)
  lrs <- lrs[lrs$last_locality %in%
               c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332)]
  lrs <- lrs[lrs$ringnr %in% gp_data$id_red]

  # add predicted phenotypes and breeding values from GP model to data
  gp_data$pred_pheno <- gp_model$summary.fitted.values$mean
  gp_data$bv_mean <- gp_model$summary.random$id1$mean[
    order(gp_model$summary.random$id1$ID)][gp_data$id1]
  gp_data$bv_sd <- gp_model$summary.random$id1$sd[
    order(gp_model$summary.random$id1$ID)][gp_data$id1]

  lrs$bv_mean <- gp_data$bv_mean[match(lrs$ringnr, gp_data$id_red)]
  lrs$bv_sd <- gp_data$bv_sd[match(lrs$ringnr, gp_data$id_red)]
  lrs$age <- lrs$year - lrs$hatch_year
  # Remove inds with missing sex
  lrs <- lrs[!is.na(lrs$sex), ]
  lrs$ringnr_num <- match(lrs$ringnr, unique(lrs$ringnr))
  lrs$ll_num <- match(lrs$last_locality, unique(lrs$last_locality))
  lrs$y_num <- match(lrs$year, unique(lrs$year))
  lrs$hy_num <- match(lrs$hatch_year, unique(lrs$hatch_year))

  # Add inbreeding info
  ibc <- fread(file = inbreeding)
  gp_data$fhat3 <- ibc$Fhat3[match(gp_data$ringnr, ibc$IID)]
  lrs$fhat3 <- gp_data$fhat3[match(lrs$ringnr, gp_data$id_red)]

  # Number of crossover measurements
  lrs$co_n <- gp_data$n[match(lrs$ringnr, gp_data$id_red)]
  lrs$co_meas <- !is.na(lrs$co_n)
  lrs <- lrs[order(lrs$co_meas), ]

  # Add crossover count measurements, and phenotypic predictions, with repeats
  ph_dat <- lapply(1:nrow(lrs),
                   function(i) {
                     lrs_row <- lrs[i, ]
                     co_red <- gp_data[gp_data$id_red == lrs_row$ringnr, ]
                     lrs_df <- lrs_row[rep(1, nrow(co_red)), ]
                     cbind(lrs_df, co_red[, c(paste0("co_count_", sex_lc),
                                              "pred_pheno")])
                   }) %>%
    do.call(what = "rbind", .)

  lst(lrs, ph_dat)
}

find_inbreeding <- function(bfile,
                            mem,
                            ncores,
                            plink_path) {

  # LD pruning
  exit_code <- system2(plink_path,
                       paste0("--bfile ", bfile, " ",
                              "--chr-set 32 ",
                              "--memory ", mem, " ",
                              "--indep 50 5 2 ",
                              "--threads ", ncores, " ",
                              "--out ", bfile))
  if (exit_code != 0) {
    stop("Error in plink")
  }
  # Calc inbreeding
  exit_code <- system2(plink_path,
                       paste0("--bfile ", bfile, " ",
                              "--chr-set 32 ",
                              "--memory ", mem, " ",
                              "--exclude ", bfile, ".prune.out ",
                              "--ibc ",
                              "--threads ", ncores, " ",
                              "--out ", bfile))
  if (exit_code != 0) {
    stop("Error in plink")
  }

  paste0(bfile, ".ibc")
}

#### Some functions for genomic prediction

make_grm <- function(analysis_inds,
                     bfile,
                     ncores,
                     mem,
                     plink_path,
                     dir) {

  dir.create(dir, showWarnings = FALSE)

  # Create a file listing which individuals to make GRM for
  write.table(data.frame(1, analysis_inds),
              file = paste0(dir, "/keep.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  # Create GRM files
  exit_code <- system2(plink_path,
                       paste0("--bfile ", bfile, " ",
                              "--keep ", dir, "/keep.txt ",
                              "--chr-set 32 ",
                              "--memory ", mem, " ",
                              "--maf 0.001 ",
                              "--freq ",
                              "--make-rel square bin ", # calculate raw GRM
                              "--make-just-bim ", # to know which SNPS were kept
                              "--threads ", ncores, " ",
                              "--out ", dir, "/grm"))

  if (exit_code != 0) {
    stop("Error in plink")
  }

  paste0(dir, "/grm.", c("rel.bin", "rel.id", "frq", "bim"))
}

load_grm <- function(dir,
                     pheno_data) {

  ids <- fread(file = paste0(dir, "/grm.rel.id"),
               select = 2,
               data.table = FALSE,
               header = FALSE)
  n <- dim(ids)[1]

  # Load into R
  grm <- matrix(readBin(paste0(dir, "/grm.rel.bin"),
                        "numeric",
                        size = 8,
                        n ^ 2),
                nrow = n)

  dimnames(grm)[[1]] <- dimnames(grm)[[2]] <-
    pheno_data$id1[match(ids$V2, pheno_data$ringnr)]

  # Add small value to diagonal to get positive definite matrix
  grm <- grm + diag(1e-2, nrow(grm))

  if (!isSymmetric(grm))
    stop("GRM not symmetric")

  # Check positive definiteness:
  e_vals <- eigen(grm, only.values = TRUE)$values
  if (!all(Im(e_vals) == 0 & Re(e_vals) > 0))
    stop("GRM not positive definite")

  # Return GRM and its inverse
  lst(grm = grm, inv_grm = solve(grm))
}

make_prior <- function(pc_matrix = NULL,
                       va_apriori = NULL,
                       pc_prec_upper_var,
                       var_init,
                       tau = 0.05) {

  rr_effect_var <- NULL
  if (!is.null(pc_matrix) && !is.null(va_apriori)) {
    # Assume each PC-"marker" has N(0, var), where var is the
    # PC-equivalent of va / 2sum(p(1-p))
    # See Aspheim et al. (2024)
    rr_effect_var <- list(
      prec = list(fixed = FALSE, # TRUE: we estimate the numerator
                  initial = log(1 / (va_apriori / sum(diag(var(pc_matrix))))))
    )
  }

  # PC priors for random effects
  hyperpar_var <- list(
    prec = list(initial = log(1 / var_init),
                prior = "pc.prec",
                param = c(sqrt(pc_prec_upper_var), tau),
                fixed = FALSE))

  lst(hyperpar_var, rr_effect_var)
}

run_gp <- function(pheno_data,
                   inverse_relatedness_matrix,
                   effects_vec,
                   y,
                   pc_matrix = NULL,
                   va_apriori = NULL) {

  pheno_data$y <- pheno_data[, y]

  # Use PC-priors for random effect variances,
  # and default (wide Gaussian) for fixed effects
  prior <- make_prior(pc_matrix,
                      va_apriori,
                      pc_prec_upper_var = var(pheno_data$y, na.rm = TRUE) / 2,
                      var_init = var(pheno_data$y, na.rm = TRUE) / 3,
                      tau = 0.05)

  inla_formula <- reformulate(effects_vec, response = "y")

  # Run model
  inla(inla_formula,
       family = "gaussian",
       data = pheno_data,
       verbose  = TRUE,
       # control.compute = list(config = TRUE),
       control.family = list(hyper = prior$hyperpar_var)) %>%
    inla.rerun()
}

inla_posterior_variances <- function(prec_marginal) {
  sigma_marg <- inla.tmarginal(function(x) 1 / x, prec_marginal)
  inla.zmarginal(sigma_marg, silent = TRUE)
}

run_gp_pois <- function(pheno_data,
                        inverse_relatedness_matrix,
                        effects_vec,
                        y,
                        E = 1,
                        pc_matrix = NULL,
                        va_apriori = NULL) {

  pheno_data$y <- pheno_data[, y]

  # Use PC-priors for random effect variances,
  # and default (wide Gaussian) for fixed effects
  prior <- make_prior(pc_matrix,
                      va_apriori,
                      pc_prec_upper_var =
                        var(log(pheno_data$y), na.rm = TRUE) / 2,
                      var_init = var(log(pheno_data$y), na.rm = TRUE) / 3,
                      tau = 0.05)

  inla_formula <- reformulate(effects_vec, response = "y")

  inla(inla_formula,
       data = pheno_data,
       verbose  = TRUE,
       # control.compute = list(config = TRUE),
       family = "poisson",
       E = E) %>%
    inla.rerun()
}

inla_cv <- function(model,
                    pheno_data) {
  ########### Leave-individual-out CV
  # INLA can do this without refitting the model
  groups <- lapply(pheno_data$id,
                   FUN = function(id) {
                     which(pheno_data$id == id)
                   })

  cvs <- inla.group.cv(result = model, groups = groups)

  data.frame(
    cv_mean = cvs$mean, # Prediction phenotype when left out
    cv_ls = log(cvs$cv) # Log-score (goodness of fit) for each obs
  )
}

# Does PCA
plink_pca <- function(analysis_inds,
                      bfile,
                      ncores,
                      mem,
                      plink_path,
                      dir) {

  max_num_pc <- length(analysis_inds)

  dir.create(dir, showWarnings = FALSE)
  # Create file containing ringnr of inds to include
  write.table(data.frame(1, analysis_inds),
              file = paste0(dir, "/keep.txt"),
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)

  exit_code <- system2(
    plink_path,
    paste0("--bfile ", bfile, " ",
           "--pca ", max_num_pc, " \'header\' ", # Do PCA
           "--keep ", dir, "/keep.txt ", # Include only the inds. in keep.txt
           "--chr-set 32 ", # Sparrow chromosomes
           "--memory ", mem, " ",
           "--threads ", ncores, " ",
           "--out ", dir, "/pca"))

  if (exit_code != 0) {
    stop("Error in plink")
  }
}

make_pc_matrix <- function(eigenvec_file,
                           analysis_inds) {

  # Load PCA results
  pca <- fread(file = eigenvec_file)

  # Create eigenvector matrix (V)
  zz <- as.matrix(pca[, -c("FID", "IID")])[match(analysis_inds, pca$IID), ]

  # Standardizing so column 1 has variance 1
  zz / sqrt(var(zz)[1])
}

find_num_pc <- function(eigenval_file,
                        min_var_explained) {

  # Find amount of variance explained by each PC cumulatively
  eigenvals <- fread(file = eigenval_file)
  variance_proportion <- cumsum(eigenvals$V1) / sum(eigenvals$V1)
  enough_variance_explained <- cumsum(variance_proportion > min_var_explained)

  # Number of PCs to explain a given amount of the genomic variation
  num_pc <- which(enough_variance_explained == 1)
  # To avoid crashing (?) INLA when using uniform PC scalings:
  if (num_pc == 1) num_pc <- 2

  var_explained_by_num_pc <- variance_proportion[num_pc]

  lst(min_var_explained, num_pc, var_explained_by_num_pc)
}

reduce_pc_matrix <- function(full_pc_matrix,
                             num_pc_obj) {
  full_pc_matrix[, 1:num_pc_obj$num_pc]
}

get_exp_acc <- function(p,
                        h2,
                        M,
                        N) {

  # sqrt(p / (1 + M / (N * h2 * p)))

  N * (h2 * p)^2 / (N * h2 * p + M)

}

inla_post_var <- function(prec_marg) {
  inla.tmarginal(function(x) 1 / x, prec_marg) %>%
    inla.zmarginal(., silent = TRUE) %>%
    as.data.frame()
}

get_10fcv_acc <- function(model, data, y_str) {
  bv <- model$summary.random$id1$mean[order(model$summary.random$id1$ID)]
  bv_rep <- bv[data$id1]

  y <- data[, y_str]

  bv_acc_unscaled <- cor.test(bv_rep[data$test], y[data$test])
  va_est <- 1 / model$summary.hyperpar["Precision for id1", "mean"]
  bv_acc <- bv_acc_unscaled$estimate / sqrt(va_est / var(y))
  bv_acc_lower <- bv_acc_unscaled$conf.int[1] / sqrt(va_est / var(y))
  bv_acc_upper <- bv_acc_unscaled$conf.int[2] / sqrt(va_est / var(y))

  data.frame(bv_acc, bv_acc_lower, bv_acc_upper)
}

pc_rate <- function(U, a = 0.05) {
  -log(a) / U
}

get_ars_samps <- function(model) {
  # Extract posterior samples
  samps <- rstan::extract(model,
                          pars = c("alpha",
                                   "beta_bv",
                                   "beta_bv2",
                                   "beta_age",
                                   "beta_age2",
                                   "beta_f",
                                   "ye",
                                   "ll",
                                   "id",
                                   "bv_lat",
                                   "sigma_ye",
                                   "sigma_ll",
                                   "sigma_id",
                                   "phi",
                                   "beta_prior_sd",
                                   "bv_lat_full_mean",
                                   "bv_lat_full_sd"))

  samp_params <- get_sampler_params(model, inc_warmup = FALSE)
  samps$energy <- unlist(lapply(samp_params, function(x) x[, "energy__"]))
  samps$divergent <- unlist(lapply(samp_params, function(x) x[, "divergent__"]))

  samps
}

get_surv_samps <- function(model) {
  # Extract posterior samples
  samps <- rstan::extract(model,
                          pars = c("alpha",
                                   "beta_bv",
                                   "beta_bv2",
                                   "beta_age",
                                   "beta_age2",
                                   "beta_f",
                                   "ye",
                                   "ll",
                                   "id",
                                   "res",
                                   "bv_lat",
                                   "sigma_ye",
                                   "sigma_ll",
                                   "sigma_id",
                                   "sigma_res",
                                   "beta_prior_sd",
                                   "bv_lat_full_mean",
                                   "bv_lat_full_sd"))

  samp_params <- get_sampler_params(model, inc_warmup = FALSE)
  samps$energy <- unlist(lapply(samp_params, function(x) x[, "energy__"]))
  samps$divergent <- unlist(lapply(samp_params, function(x) x[, "divergent__"]))

  samps
}

ars_bvpost_sim <- function(data) {

  alpha_std <- alpha <- -0.077466411
  beta_bv_std <- beta_bv <-  0.1
  beta_bv2_std <- beta_bv2 <- -0.5
  beta_age_std <- beta_age <- 0.182861276
  beta_age2_std <- beta_age2 <- -0.044908613
  beta_f_std <- beta_f <- -0.111374375

  sigma_ye <- 0.367186607
  sigma_ll <- 0.411729306
  sigma_id <- 0.599423570

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)

  phi <- 1.670098105

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = bv_lat[data$id_idx]

  eta1 <-
    # alpha_std +
    # beta_bv_std * scale(bv_lat_full) +
    # beta_bv2_std * scale(bv_lat_full)^2 +
    # beta_age_std * scale(data$age) +
    # beta_age2_std * scale(data$age)^2 +
    # beta_f_std * scale(data$f)
    alpha +
    beta_bv * bv_lat_full +
    beta_bv2 * (bv_lat_full)^2 +
    beta_age * data$age +
    beta_age2 * (data$age)^2 +
    beta_f * data$f
  eta2 <- ye_full +
    ll_full +
    id_full

  eta <- eta1 + eta2

  mu <- exp(eta)
  rnbinom(n = data$N, mu = mu, size = phi)
}

surv_bvpost_sim <- function(data) {

  alpha_std <- alpha <- -0.077466411
  beta_bv_std <- beta_bv <-  0.1
  beta_bv2_std <- beta_bv2 <- 0.5
  beta_age_std <- beta_age <- 0.182861276
  beta_age2_std <- beta_age2 <- -0.044908613
  beta_f_std <- beta_f <- -0.111374375

  sigma_ye <- 0.367186607
  sigma_ll <- 0.411729306
  sigma_id <- 0.599423570
  sigma_res <- 0.4

  ye <- rnorm(n = data$N_ye, 0, sigma_ye)
  ll <- rnorm(n = data$N_ll, 0, sigma_ll)
  id <- rnorm(n = data$N_id, 0, sigma_id)
  res <- rnorm(n = data$N, 0, sd = sigma_res)

  bv_lat <- rnorm(n = data$N_id, mean = data$bv_mean, sd = data$bv_sd)

  ye_full = ye[data$ye_idx]
  ll_full = ll[data$ll_idx]
  id_full = id[data$id_idx]
  bv_lat_full = bv_lat[data$id_idx]

  eta1 <-
    # alpha_std +
    # beta_bv_std * scale(bv_lat_full) +
    # beta_bv2_std * scale(bv_lat_full)^2 +
    # beta_age_std * scale(data$age) +
    # beta_age2_std * scale(data$age)^2 +
    # beta_f_std * scale(data$f)
    alpha +
    beta_bv * bv_lat_full +
    beta_bv2 * (bv_lat_full)^2 +
    beta_age * data$age +
    beta_age2 * (data$age)^2 +
    beta_f * data$f
  eta2 <- ye_full +
    ll_full +
    id_full +
    res

  eta <- eta1 + eta2

  p <- 1 / (1 + exp(-eta))
  # print(summary(p))
  rbinom(n = data$N, p = p, size = 1)
}

samp_plot_df <- function(x, y, n_samp) {
  y %>%
    as.data.frame() %>%
    mutate(sample = 1:n_samp) %>%
    tidyr::pivot_longer(-sample, names_to = "x_idx", values_to = "y_samp") %>%
    mutate(x = rep(x, n_samp),
           y_mean = rep(apply(y, 2, mean), n_samp),
           y_lower = rep(apply(y, 2, quantile, p = 0.025), n_samp),
           y_upper = rep(apply(y, 2, quantile, p = 0.975), n_samp))
}

bv_preds_and_marg <- function(data,
                              samps,
                              inv_link,
                              n_plot = 200) {
  n_samp <- length(samps$alpha)
  avg_age <- mean(data$age)
  avg_f <- mean(data$f)

  bv <- seq(min(c(data$bv_mean) - 2 * max(data$bv_sd)),
            max(c(data$bv_mean) + 2 * max(data$bv_sd)),
            length.out = n_plot)

  y_pred_samples <- array(NA, dim = c(n_samp, n_plot))

  for (i in seq_len(n_samp)) {
    for (j in seq_len(n_plot)) {
      y_pred_samples[i, j] <- samps$alpha[i] +
        samps$beta_bv[i] * bv[j] +
        samps$beta_bv2[i] * bv[j]^2 +
        samps$beta_age[i] * avg_age +
        samps$beta_age2[i] * avg_age^2 +
        samps$beta_f[i] * avg_f
    }
  }

  df_pred <- samp_plot_df(y = inv_link(y_pred_samples), x = bv, n_samp = n_samp)

  marg_effect <- sapply(bv, function(x) samps$beta_bv + 2 * samps$beta_bv2 * x)
  df_marg <- samp_plot_df(y = marg_effect, x = bv, n_samp = n_samp)

  lst(df_pred, df_marg)
}

plot_lines_posterior <- function(df, xlab, ylab, title) {
  ggplot(df, aes(x = x)) +
    # Transparent lines for posterior samples
    geom_line(aes(y = y_samp,
                  group = sample,
                  color = "Posterior samples",
                  linetype = "Posterior samples"),
              alpha = 0.01) +
    # Posterior mean
    geom_line(aes(y = y_mean,
                  color = "Posterior mean",
                  linetype = "Posterior mean"),
              linewidth = 1) +
    # Posterior 95% CI
    geom_line(aes(y = y_lower, color = "95% CI", linetype = "95% CI")) +
    geom_line(aes(y = y_upper, color = "95% CI", linetype = "95% CI")) +
    # Adding a title and axis labels
    labs(x = xlab,
         y = ylab,
         title = title,
         color = "",
         linetype = "") + # Combine legends
    # Minimal theme for clean appearance
    theme_minimal() +
    # Customize the color and linetype scales
    scale_color_manual(values = c("Posterior mean" = "blue",
                                  "95% CI" = "black",
                                  "Posterior samples" = "black")) +
    scale_linetype_manual(values = c("Posterior mean" = "solid",
                                     "95% CI" = "dashed",
                                     "Posterior samples" = "solid")) +
    # Ensure the correct legend appearance
    guides(color = guide_legend(
      override.aes = list(
        # Match transparency for samples, full for others
        alpha = c(1, 1, 0.01),
        # Make legend dashed looked like in plot
        linewidth = c(0.5, 1, 1))),
      linetype = guide_legend(override.aes = list(
        # Correct dashing in legend
        linetype = c("dashed", "solid", "solid"))))
}

ggsave_path <- function(filename,
                        ...) {
  ggsave(filename, ...)
  filename
}
