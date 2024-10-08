###############################################################################
#
# crossover_gp.R
#
# Author: Kenneth Aase
#
# Purpose:
#     Genomic Prediction of sex-specific house sparrow crossover rates in the
#     Helgeland house sparrow system, using a GRM and R-INLA
#
# Accompanying scripts:
#     crossover_gp_func.R - contains R functions used by this script
###############################################################################

# Load libraries
library(data.table)
library(tibble)
library(stats)
library(ggplot2)
library(dplyr)
library(INLA) # I used version INLA_23.04.24
# INLA is not on CRAN, you can install it with:
# install.packages(
#   "INLA",
#   repos = c(getOption("repos"),
#             INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)

# Load custom functions
source("r/crossover_gp_inla_func.R")

# Filepaths that might need to be changed
recomb_data_path <- "data/20240910_Sparrow_Recomb_Data.txt"
geno_data_path <- "data/70K_200K_maf_geno_mind_v5"
plink_path <- "PLINK/plink_linux" # path to plink program

############## Data Wrangling ##################

# Load data
pheno_data_all <- fread(recomb_data_path, data.table = FALSE)

pheno_data_mean <- pheno_data_all %>%
  group_by(id, sex) %>%
  summarise(co_count = mean(co_count),
            intra_shuff = mean(intra_shuff),
            total_coverage = mean(total_coverage),
            var_co_count = var(co_count),
            var_intra_shuff = var(intra_shuff),
            var_total_coverage = var(total_coverage),
            n = n())

pheno_data_rand <- pheno_data_all %>%
  group_by(id, sex) %>%
  summarise(co_count = co_count[1],
            intra_shuff = intra_shuff[1],
            total_coverage = total_coverage[1],
            n = n())

pheno_data <- pheno_data_mean

# Make columns for sex-specific crossover rates
pheno_data$co_count_m <- pheno_data$co_count_f <- pheno_data$co_count
pheno_data$co_count_m <- ifelse(pheno_data$sex == "M",
                                pheno_data$co_count_m,
                                NA)
pheno_data$co_count_f <- ifelse(pheno_data$sex == "M",
                                NA,
                                pheno_data$co_count_f)

# For now: only female
pheno_data <- dplyr::filter(pheno_data, sex == "F")

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

# Test set
test_inds_num <- sort(sample(max(ids_num), round(max(ids_num) * 0.3)))
pheno_data$test <- ifelse(ids_num %in% test_inds_num, TRUE, FALSE)
pheno_data$co_count_f_test <- pheno_data$co_count_f
pheno_data$co_count_f_test[pheno_data$test] <- NA
pheno_data$n_obs <- c(diff(ids_num_unique),
                      1 + length(ids_num) - tail(ids_num_unique, 1))[ids_num]

############## Modelling ##################

pca_dir <- "data/pca_crossover_helgeland"

plink_pca(analysis_inds = unique(pheno_data$ringnr), # Vector of inds.,
          bfile = geno_data_path,
          ncores = 4,
          mem = 20 * 6000,
          plink_path = plink_path,
          dir = pca_dir)

npc <- lapply(seq(from = 0.01, by = 0.01, to = 0.99),
              function(i) {
                find_num_pc(eigenval_file = paste0(pca_dir, "/pca.eigenval"),
                            min_var_explained = i) %>% as.data.frame()
              }) %>%
  do.call(what = rbind, .)

npc$exp_acc <- apply(
  npc,
  1,
  function(row) {
    get_exp_acc(p = row[3],
                h2 = 0.2963445,
                N = length(unique(pheno_data$id[!pheno_data$test])),
                M = row[2])
  })
# Use best number of PCs according to expected accuracy formula
num_pc <- npc[which.max(npc$exp_acc), "num_pc"]

pc_matrix <- make_pc_matrix(
  eigenvec_file = "data/pca_crossover_helgeland/pca.eigenvec",
  analysis_inds = unique(pheno_data$ringnr))[, 1:num_pc]

# Create GRM using PLINK
make_grm(analysis_inds = unique(pheno_data$ringnr), # Vector of inds.
         bfile = geno_data_path,
         ncores = 4,
         mem = 20 * 6000,
         # Path to plink program:
         plink_path = plink_path,
         # Where to save result:
         dir = "data/grm_crossover_helgeland")

grm_obj <- load_grm(dir = "data/grm_crossover_helgeland",
                    pheno_data = pheno_data)

effects <- c(
  # Fixed effects:
  "1", # intercept
  "total_coverage_scaled",
  # "total_coverage2_scaled",
  # Random effects:
  # Hatch year
  # "f(hatch_year, model = \"iid\", hyper = prior$hyperpar_var)",
  # Breeding values
  "f(id1,
     values = as.numeric(colnames(inverse_relatedness_matrix)),
     model = \"generic0\",
     hyper = prior$hyperpar_var,
     constr = FALSE,
     Cmatrix = inverse_relatedness_matrix)",
  # Breeding values (BPCRR) version
  # "f(id1, model = \"z\", Z = pc_matrix, hyper = prior$rr_effect_var)",
  # ID effect (use when using repeated measurements)
  # "f(id2,
  #   model = \"iid\",
  #   values = unique(pheno_data$id2),
  #   hyper = prior$hyperpar_var)",
  # Measurement error
  "f(id3,
    model = \"iid\",
    values = unique(pheno_data$id3),
    scale = pheno_data$n[match(unique(pheno_data$id), pheno_data$id)],
    hyper = prior$hyperpar_var)"
)

# Model female-specific crossover rates
model <- run_gp(pheno_data = pheno_data,
                inverse_relatedness_matrix = grm_obj$inv_grm,
                effects_vec = effects,
                # scale = sqrt(pheno_data$n),
                pc_matrix = pc_matrix,
                va_apriori = 6,
                y = "co_count_f_test")

######################### Predictions ####################

# Predicted breeding values and phenotypes (mean posterior predictions)
# Female-specific crossover
pheno_data$pred_bv <-
  model$summary.random$id1$mean[order(model$summary.random$id1$ID)][ids_num]
pheno_data$pred_pheno <- model$summary.fitted.values$mean

ggplot(data = pheno_data[pheno_data$test, ],
       aes(x = pred_bv, y = co_count_f, col = n > 2)) +
  geom_point()

ggplot(data = pheno_data[pheno_data$test & pheno_data$n <= 2, ],
       aes(x = pred_bv, y = co_count_f)) +
  geom_point()



