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
                   va_apriori = NULL,
                   scale = NULL) {

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
       scale = scale,
       # control.compute = list(config = TRUE),
       control.family = list(hyper = prior$hyperpar_var)) %>%
    inla.rerun() %>%
    inla.rerun()
}

inla_posterior_variances <- function(prec_marginal) {
  sigma_marg <- inla.tmarginal(function(x) 1 / x, prec_marginal)
  inla.zmarginal(sigma_marg, silent = TRUE)
}

run_gp_pois <- function(pheno_data,
                        inverse_relatedness_matrix,
                        effects_vec,
                        y) {

  pheno_data$y <- pheno_data[, ..y]

  # Use PC-priors for random effect variances,
  # and default (wide Gaussian) for fixed effects
  prior <- make_prior(pc_prec_upper_var = var(pheno_data$y, na.rm = TRUE),
                      var_init = var(pheno_data$y, na.rm = TRUE),
                      tau = 0.05)

  inla_formula <- reformulate(effects_vec, response = "y")

  inla(inla_formula,
       data = pheno_data,
       verbose  = TRUE,
       control.compute = list(config = TRUE),
       control.predictor = list(link = 1),
       family = "poisson") %>%
    inla.rerun() %>%
    inla.rerun()
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
