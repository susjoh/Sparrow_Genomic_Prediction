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
library(qs)
library(magrittr)
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
nestling_data_path <- "data/Nestlingdata_20240312_fix.csv"
lrs_data_path <- "data/LRS_data_20230803_Helgeland_fix.csv"
lrs_data_path2 <- "data/Missing_LRS_Sparrows_revised_WithInfo_fix.csv"
morph_data_path <- "data/AdultMorphology_20240201_fix.csv"

############## Data Wrangling ##################

# Load data
pheno_data_all <- fread(recomb_data_path, data.table = FALSE)

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

pheno_data_rand <- pheno_data_all %>%
  group_by(id, sex) %>%
  summarise(co_count = co_count[sample(n(), 1)],
            intra_shuff = intra_shuff[sample(n(), 1)],
            total_coverage = total_coverage[sample(n(), 1)],
            n = n())

# # No U (or upside-down U) shape
# ggplot(data = pheno_data_rand, aes (x = co_count, y = n)) +
#   geom_point() +
#   geom_smooth(method = "lm")

pheno_data_all$n <-
  pheno_data_mean$n[match(pheno_data_all$id, pheno_data_mean$id)]

pheno_data <- pheno_data_all

pheno_data$hatch_year <- as.factor(pheno_data$hatch_year)
pheno_data$first_locality <- as.factor(pheno_data$first_locality)

#########################

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

lrs_red <- rbind(
  lrs[!lrs$ringnr %in% pheno_data$id,
      c("ringnr", "sex", "hatch_year", "first_locality")],
  lrs2[!lrs2$ringnr %in% pheno_data$id,
       c("ringnr", "sex", "hatch_year", "first_locality")])
lrs_red$sex <- ifelse(lrs_red$sex > 1.5, "F", "M")
lrs_red$hatch_year <- as.factor(lrs_red$hatch_year)
lrs_red$first_locality <- as.factor(lrs_red$first_locality)

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
                # N = length(unique(pheno_data$id[!pheno_data$test])),
                N = length(unique(pheno_data$id)) - 1,
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

# Here I specify the effects to include in the model
effects <- c(
  ################## Fixed effects:
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

  # Breeding values (alternative to GRM: ridge regression on PCs of SNP matrix)
  # "f(id1, model = \"z\", Z = pc_matrix, hyper = prior$rr_effect_var)",

  # ID effect (use this when using single measurements as the response)
  "f(id2,
    model = \"iid\",
    values = unique(pheno_data$id2),
    hyper = prior$hyperpar_var)"

  # If using means as response: include this to account for different variances
  # "f(id3,
  #   model = \"iid\",
  #   values = unique(pheno_data$id3),
  #   scale = pheno_data$n[match(unique(pheno_data$id), pheno_data$id)],
  #   hyper = prior$hyperpar_var)"
)

# Full genomic animal model for female-specific crossover rates
# model <- run_gp(pheno_data = pheno_data,
#                 inverse_relatedness_matrix = grm_obj$inv_grm,
#                 effects_vec = effects,
#                 pc_matrix = pc_matrix,
#                 va_apriori = 6,
#                 # y = "co_count_f_test"
#                 y = "co_count_f"
# )

# qsave(model, file = "results/fullmodel")

model <- qread(file = "results/fullmodel")

# Alternatively: Same model with Poisson likelihood
# - only a slight improvement in accuracy though, so don't bother.
model_pois <- run_gp_pois(pheno_data = pheno_data,
                     inverse_relatedness_matrix = grm_obj$inv_grm,
                     effects_vec = effects,
                     pc_matrix = pc_matrix,
                     va_apriori = 6,
                     y = "co_count_f",
                     E = 1
)


qsave(model_pois, file = "results/fullmodel_pois")

model_pred <- run_gp(pheno_data = pheno_data,
                     inverse_relatedness_matrix = grm_obj$inv_grm,
                     effects_vec = effects,
                     va_apriori = 6,
                     y = "co_count_f_test"
)

######################### Variance components ###########

# Transform precisions to variances
lapply(model$marginals.hyperpar, inla_post_var) %>%
  do.call(rbind, .) ->
  var_comps
rownames(var_comps) %<>% gsub("Precision", "Variance", .) # Fix rownames
# Posterior statistics for the variance components
var_comps

lambda <- sqrt(var_comps["Variance for id1", "mean"] /
                 var(pheno_data$co_count_f))

######################### Predictions ####################

# Predicted breeding values and phenotypes
# (we use the posterior mean as a point estimate)
pheno_data$pred_bv <-
  model_pred$summary.random$id1$mean[order(model$summary.random$id1$ID)][ids_num]
pheno_data$pred_pheno <- model_pred$summary.fitted.values$mean

pheno_data$est_bv <-
  model$summary.random$id1$mean[order(model$summary.random$id1$ID)][ids_num]
pheno_data$est_pheno <- model$summary.fitted.values$mean

########### Leave-individual-out CV
# INLA can do this without refitting the model
groups <- lapply(pheno_data$id,
                 FUN = function(id) {
                   which(pheno_data$id == id)
                 })

cvs <- inla.group.cv(result = model, groups = groups)

pheno_data$cv_mean <- cvs$mean # Prediction phenotype when left out
pheno_data$cv_ls <- log(cvs$cv) # Log-score (goodness of fit) for each obs

################ CV plots
ggplot(data = pheno_data,
       aes(x = cv_mean, y = co_count_f)) +
  geom_point() +
  stat_smooth(method = "lm") +
  geom_abline()

# Residuals
ggplot(data = pheno_data, aes(cv_mean, co_count_f - cv_mean)) +
  geom_point() + stat_smooth(method = "lm")

# Better fit for individuals with intermediate co_count values
ggplot(data = pheno_data, aes(co_count_f, cv_ls, color = n > 2)) +
  geom_point() +
  stat_smooth()

# No difference in goodness of fit inds with different #obs
ggplot(data = pheno_data,
       aes(x = n, y = cv_ls)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = "loess")

ggplot(
  mapping = aes(x = pheno_data$cv_mean[unique(pheno_data$id1)],
                y = pheno_data_mean$co_count[pheno_data_mean$sex == "F"])) +
  geom_point() +
  stat_smooth(method = "lm") +
  geom_abline()

# Leave-individual-out prediction accuracy (of phenotype, not breeding value)
# Unscaled:
liocv_pheno_acc_unscaled <- with(pheno_data, cor(x = cv_mean, y = co_count_f))
liocv_pheno_acc_unscaled
# Properly scaled:
l <- sqrt(1 - var_comps$mean[1] / var(pheno_data$co_count)) # scaling factor
liocv_pheno_acc <- liocv_pheno_acc_unscaled / l


##### Test-set model plots
ggplot(data = pheno_data[pheno_data$test, ],
       aes(x = pred_bv, y = co_count_f, col = n > 2)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(data = pheno_data[pheno_data$test, ],
       aes(x = pred_bv, y = n)) +
  geom_point() +
  stat_smooth(method = "lm")

# BERKSON!?
ggplot(data = pheno_data[pheno_data$test, ],
       aes(x = pred_pheno, y = co_count_f)) +
  geom_point() +
  stat_smooth(method = "lm") +
  coord_equal() +
  geom_abline()

# u and w independent?
lm(pred_pheno ~ I(co_count_f - pred_pheno),
   data = pheno_data[pheno_data$test, ]) %>%
  summary()

# but normality mainly holds
ggplot(data = pheno_data[pheno_data$test, ],
       aes(sample = co_count_f - pred_pheno)) +
  geom_qq() +
  geom_qq_line()
ggplot(data = pheno_data[pheno_data$test, ], aes(x = co_count_f - pred_pheno)) +
  geom_density()

# Berkson bv?
ggplot(data = pheno_data[pheno_data$test, ],
       aes(x = pred_bv, y = est_bv)) + # 0.62 corr between est_bv and co_count_f
  geom_point() +
  stat_smooth(method = "lm") +
  coord_equal() +
  geom_abline()

# Indep?
lm(pred_bv ~ I(est_bv - pred_bv),
   data = pheno_data[pheno_data$test, ]) %>%
  summary()

# normality?
ggplot(data = pheno_data[pheno_data$test, ],
       aes(sample = est_bv - pred_bv)) +
  geom_qq() +
  geom_qq_line()
ggplot(data = pheno_data[pheno_data$test, ], aes(x = est_bv - pred_bv)) +
  geom_density()

ggplot(data = pheno_data[pheno_data$test, ],
       aes(x = pred_pheno, y = co_count_f - pred_pheno)) +
  geom_point() +
  stat_smooth(method = "loess") +
  geom_hline(yintercept = 0)


with(pheno_data, qqplot(pred_bv[test & n > 2], pred_bv[test & n <= 2]))
abline(0, 1)

with(pheno_data, qqplot(cv_m[test & n > 2], pred_pheno[test & n <= 2]))
abline(0, 1)
