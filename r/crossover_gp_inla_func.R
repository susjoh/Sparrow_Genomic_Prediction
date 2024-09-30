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

make_prior <- function(pc_prec_upper_var,
                       var_init,
                       tau = 0.05) {

  # PC priors for random effects
  hyperpar_var <- list(
    prec = list(initial = log(1 / var_init),
                prior = "pc.prec",
                param = c(sqrt(pc_prec_upper_var), tau),
                fixed = FALSE))

  lst(hyperpar_var)
}

run_gp <- function(pheno_data,
                   inverse_relatedness_matrix,
                   effects_vec,
                   y,
                   scale = NULL) {

  pheno_data$y <- pheno_data[, y]

  # Use PC-priors for random effect variances,
  # and default (wide Gaussian) for fixed effects
  prior <- make_prior(pc_prec_upper_var = var(pheno_data$y, na.rm = TRUE) / 2,
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
