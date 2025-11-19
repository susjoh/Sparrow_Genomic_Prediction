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

# Rename columns so that existing code remains compatible with both data sets
co_data_rename_cols <- function(recomb_data_path) {

  dat <- fread(recomb_data_path, data.table = FALSE)

  dat$co_count <- dat$total_CO_count
  dat$id <- dat$parent
  dat$total_coverage <- dat$Total_Coverage

  path <- gsub(x = recomb_data_path, pattern = ".txt", "")
  new_path <- paste0(path, "_newcols.txt")

  write.table(dat, file = new_path, quote = FALSE, row.names = FALSE)
  new_path
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

make_adult_gp_data <- function(pheno_data,
                               lrs_path,
                               lrs_path2,
                               nestling_path,
                               ped_path,
                               fam_path,
                               sex_num,
                               sex_keep) {

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
  # Remove inds with missing sex
  lrs %<>% dplyr::filter(!is.na(sex))
  # Sex filtering
  if (sex_num != "all") {
    lrs %<>%
      dplyr::filter(sex != 1.5) %>%
      mutate(sex = round(sex)) %>%
      filter(sex %in% sex_num)
  }
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

# Sex-GP on inds with offspring with fitness measures
make_parent_gp_data <- function(pheno_data,
                                lrs_path,
                                lrs_path2,
                                nestling_path,
                                fam_path,
                                ped_path,
                                sex_num,
                                sex_keep) {

  geno_inds <- fread(file = fam_path, select = 2)$V2
  geno_inds_red <- gsub(x = geno_inds, pattern = "_.+", "")

  pedigree <- read.table(file = ped_path, header = TRUE)
  pedigree$id_red <- gsub(x = pedigree$id, pattern = "_.+", "")
  pedigree$dam_red <- gsub(x = pedigree$dam, pattern = "_.+", "")
  pedigree$sire_red <- gsub(x = pedigree$sire, pattern = "_.+", "")

  lrs <- fread(file = lrs_path)
  # remove last loc. outside of main islands
  lrs <- lrs[lrs$last_locality %in%
               c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332),
             c("ringnr", "sex", "hatch_year", "first_locality")]
  # Only need a single row per ind now
  lrs <- lrs[!duplicated(lrs$ringnr), ]

  lrs2 <- fread(file = lrs_path2)
  lrs2 <- lrs2[!duplicated(lrs2$ringnr), ]
  lrs2 <- lrs2[lrs2$last_locality %in%
                 c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332),
               c("ringnr", "sex", "hatch_year", "first_locality")]

  # Keep only adults with fitness data
  pedigree %<>% dplyr::filter(id_red %in% lrs$ringnr)
  # Inds we want to do GP on (genotyped, and have child with fitness data)
  if (sex_keep == "F") {
    gp_inds <- unique(pedigree$dam[pedigree$dam %in% geno_inds])
  } else if (sex_keep == "M") {
    gp_inds <- unique(pedigree$sire[pedigree$sire %in% geno_inds])
  } else {
    stop("sex error")
  }
  # Don't need to do GP for inds already phenotyped
  gp_inds <- gp_inds[!gp_inds %in% pheno_data$id]
  gp_inds_red <- gsub(x = gp_inds, pattern = "_.+", "")

  # Data frame to add to pheno_data
  gp_inds_df <- lrs[lrs$ringnr %in% gp_inds_red, ]
  # Add inds missing from lrs, but present in lrs2
  gp_inds_df <- rbind(gp_inds_df,
                      lrs2[(!lrs2$ringnr %in% gp_inds_df$ringnr) &
                             (lrs2$ringnr %in% gp_inds_red), ])
  # For duplicated genotypes:
  gp_inds_df$id_red <- gp_inds_df$ringnr
  gp_inds_df$ringnr[which(!gp_inds_df$ringnr %in% geno_inds)] <-
    sapply(which(!gp_inds_df$ringnr %in% geno_inds),
           function(ind) {
             geno_inds[match(gp_inds_df$ringnr[ind], geno_inds_red)]
           })
  gp_inds_df$id <- gp_inds_df$ringnr
  # Recasting
  gp_inds_df %<>% dplyr::filter(!is.na(sex))
  gp_inds_df$sex <- ifelse(gp_inds_df$sex > 1.5, "F", "M")
  gp_inds_df$hatch_year <- as.factor(gp_inds_df$hatch_year)
  gp_inds_df$id1 <- gp_inds_df$id2 <- gp_inds_df$id3 <-
    seq(from = max(pheno_data$id1) + 1,  length = nrow(gp_inds_df), by = 1)
  # Merge
  plyr::rbind.fill(pheno_data, gp_inds_df)
}

# Sex-GP on inds with offspring with fitness measures
make_nestling_gp_data <- function(pheno_data,
                                  lrs_path,
                                  lrs_path2,
                                  nestling_path,
                                  fam_path,
                                  ped_path,
                                  sex_num,
                                  sex_keep) {

  geno_inds <- fread(file = fam_path, select = 2)$V2
  geno_inds_red <- gsub(x = geno_inds, pattern = "_.+", "")

  pedigree <- read.table(file = ped_path, header = TRUE)
  pedigree$id_red <- gsub(x = pedigree$id, pattern = "_.+", "")
  pedigree$dam_red <- gsub(x = pedigree$dam, pattern = "_.+", "")
  pedigree$sire_red <- gsub(x = pedigree$sire, pattern = "_.+", "")

  lrs <- fread(file = lrs_path)
  # remove last loc. outside of main islands
  lrs <- lrs[lrs$last_locality %in%
               c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332),
             c("ringnr", "sex", "hatch_year", "first_locality")]
  # Only need a single row per ind now
  lrs <- lrs[!duplicated(lrs$ringnr), ]

  lrs2 <- fread(file = lrs_path2)
  lrs2 <- lrs2[!duplicated(lrs2$ringnr), ]
  lrs2 <- lrs2[lrs2$last_locality %in%
                 c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332),
               c("ringnr", "sex", "hatch_year", "first_locality")]

  nest <- fread(file = nestling_path)

  # Keep only nestlings with fitness data
  pedigree %<>% dplyr::filter(id_red %in% nest$ringnr)
  # Inds we want to do GP on (genotyped, and have child with fitness data)
  if (sex_keep == "F") {
    gp_inds <- unique(pedigree$dam[pedigree$dam %in% geno_inds])
  } else if (sex_keep == "M") {
    gp_inds <- unique(pedigree$sire[pedigree$sire %in% geno_inds])
  } else {
    stop("sex error")
  }
  # Don't need to do GP for inds already phenotyped
  gp_inds <- gp_inds[!gp_inds %in% pheno_data$id]
  gp_inds_red <- gsub(x = gp_inds, pattern = "_.+", "")

  # Data frame to add to pheno_data
  gp_inds_df <- lrs[lrs$ringnr %in% gp_inds_red, ]
  # Add inds missing from lrs, but present in lrs2
  gp_inds_df <- rbind(gp_inds_df,
                      lrs2[(!lrs2$ringnr %in% gp_inds_df$ringnr) &
                             (lrs2$ringnr %in% gp_inds_red), ])
  # For duplicated genotypes:
  gp_inds_df$id_red <- gp_inds_df$ringnr
  gp_inds_df$ringnr[which(!gp_inds_df$ringnr %in% geno_inds)] <-
    sapply(which(!gp_inds_df$ringnr %in% geno_inds),
           function(ind) {
             geno_inds[match(gp_inds_df$ringnr[ind], geno_inds_red)]
           })
  gp_inds_df$id <- gp_inds_df$ringnr
  # Recasting
  gp_inds_df %<>% dplyr::filter(!is.na(sex))
  gp_inds_df$sex <- ifelse(gp_inds_df$sex > 1.5, "F", "M")
  gp_inds_df$hatch_year <- as.factor(gp_inds_df$hatch_year)
  gp_inds_df$id1 <- gp_inds_df$id2 <- gp_inds_df$id3 <-
    seq(from = max(pheno_data$id1) + 1,  length = nrow(gp_inds_df), by = 1)
  # Merge
  plyr::rbind.fill(pheno_data, gp_inds_df)
}

make_data_adult <- function(gp_data,
                            gp_model,
                            lrs_data_path,
                            inbreeding,
                            sex_num,
                            ped_path = NULL) {

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

  # Remove inds with missing sex
  lrs %<>% dplyr::filter(!is.na(sex))
  # Sex filtering
  if (sex_num != "all") {
    lrs %<>%
      dplyr::filter(sex != 1.5) %>%
      mutate(sex = round(sex)) %>%
      filter(sex %in% sex_num)
  }

  lrs$bv_mean <- gp_data$bv_mean[match(lrs$ringnr, gp_data$id_red)]
  lrs$bv_sd <- gp_data$bv_sd[match(lrs$ringnr, gp_data$id_red)]
  lrs$age <- lrs$year - lrs$hatch_year

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
  lrs %<>% mutate(co_n = ifelse(is.na(co_n), 0, co_n))
  lrs$co_meas <- (lrs$co_n > 0)

  age_poly <- poly(lrs$age, degree = 2)
  lrs <- cbind(lrs, age_q1 = age_poly[, 1], age_q2 = age_poly[, 2])

  lrs$idx <- seq_len(nrow(lrs))
  lrs
}

make_data_parent <- function(gp_data,
                             gp_model,
                             lrs_data_path,
                             inbreeding,
                             sex_num,
                             ped_path) {

  lrs <- fread(file = lrs_data_path)
  lrs <- lrs[lrs$last_locality %in%
               c(20, 22, 23, 24, 26, 27, 28, 34, 35, 38, 331, 332)]

  pedigree <- read.table(file = ped_path, header = TRUE)
  pedigree$id_red <- gsub(x = pedigree$id, pattern = "_.+", "")
  pedigree$dam_red <- gsub(x = pedigree$dam, pattern = "_.+", "")
  pedigree$sire_red <- gsub(x = pedigree$sire, pattern = "_.+", "")

  # add predicted phenotypes and breeding values from GP model to data
  gp_data$pred_pheno <- gp_model$summary.fitted.values$mean
  gp_data$bv_mean <- gp_model$summary.random$id1$mean[
    order(gp_model$summary.random$id1$ID)][gp_data$id1]
  gp_data$bv_sd <- gp_model$summary.random$id1$sd[
    order(gp_model$summary.random$id1$ID)][gp_data$id1]

  # Remove inds with missing sex
  lrs %<>% dplyr::filter(!is.na(sex))
  # Sex filtering
  if (sex_num != "all") {
    lrs %<>%
      dplyr::filter(sex != 1.5) %>%
      mutate(sex = round(sex)) %>%
      filter(sex %in% sex_num)
  }

  # Keep lrs entries where gp_data ringnrs are in pedigree
  if (sex_num == 2) {
    pedigree %<>% dplyr::filter(dam %in% gp_data$id_red)
  } else if (sex_num == 1) {
    pedigree %<>% dplyr::filter(sire %in% gp_data$id_red)
  } else {
    stop("sex mistake")
  }
  lrs %<>% dplyr::filter(ringnr %in% pedigree$id_red)

  # add parent's BV stats
  if (sex_num == 2) {
    lrs$parent <- pedigree$dam_red[match(lrs$ringnr, pedigree$id_red)]
  } else if (sex_num == 1) {
    lrs$parent <- pedigree$sire_red[match(lrs$ringnr, pedigree$id_red)]
  }
  lrs$bv_mean <- gp_data$bv_mean[match(lrs$parent, gp_data$id_red)]
  lrs$bv_sd <- gp_data$bv_sd[match(lrs$parent, gp_data$id_red)]

  # Add age
  lrs$age <- lrs$year - lrs$hatch_year

  lrs$ringnr_num <- match(lrs$ringnr, unique(lrs$ringnr))
  lrs$parent_num <- match(lrs$parent, unique(lrs$parent))
  lrs$ll_num <- match(lrs$last_locality, unique(lrs$last_locality))
  lrs$y_num <- match(lrs$year, unique(lrs$year))
  lrs$hy_num <- match(lrs$hatch_year, unique(lrs$hatch_year))

  # Add inbreeding info
  ibc <- fread(file = inbreeding)
  gp_data$fhat3 <- ibc$Fhat3[match(gp_data$ringnr, ibc$IID)]
  lrs$fhat3 <- gp_data$fhat3[match(lrs$parent, gp_data$id_red)]

  # # Number of crossover measurements
  lrs$co_n <- gp_data$n[match(lrs$parent, gp_data$id_red)]
  lrs %<>% mutate(co_n = ifelse(is.na(co_n), 0, co_n))
  lrs$co_meas <- (lrs$co_n > 0)

  age_poly <- poly(lrs$age, degree = 2)
  lrs <- cbind(lrs, age_q1 = age_poly[, 1], age_q2 = age_poly[, 2])

  lrs$idx <- seq_len(nrow(lrs))
  lrs
}

make_data_nest <- function(gp_data,
                           gp_model,
                           nestling_data_path,
                           inbreeding,
                           sex_num,
                           ped_path) {

  nest <- fread(file = nestling_data_path)

  pedigree <- read.table(file = ped_path, header = TRUE)
  pedigree$id_red <- gsub(x = pedigree$id, pattern = "_.+", "")
  pedigree$dam_red <- gsub(x = pedigree$dam, pattern = "_.+", "")
  pedigree$sire_red <- gsub(x = pedigree$sire, pattern = "_.+", "")

  # add predicted phenotypes and breeding values from GP model to data
  gp_data$pred_pheno <- gp_model$summary.fitted.values$mean
  gp_data$bv_mean <- gp_model$summary.random$id1$mean[
    order(gp_model$summary.random$id1$ID)][gp_data$id1]
  gp_data$bv_sd <- gp_model$summary.random$id1$sd[
    order(gp_model$summary.random$id1$ID)][gp_data$id1]

  # Keep nest entries where gp_data ringnrs are in pedigree
  if (sex_num == 2) {
    pedigree %<>% dplyr::filter(dam %in% gp_data$id_red)
  } else if (sex_num == 1) {
    pedigree %<>% dplyr::filter(sire %in% gp_data$id_red)
  } else {
    stop("sex mistake")
  }
  nest %<>% dplyr::filter(ringnr %in% pedigree$id_red)

  # add parent's BV stats
  if (sex_num == 2) {
    nest$parent <- pedigree$dam_red[match(nest$ringnr, pedigree$id_red)]
  } else if (sex_num == 1) {
    nest$parent <- pedigree$sire_red[match(nest$ringnr, pedigree$id_red)]
  }
  nest$bv_mean <- gp_data$bv_mean[match(nest$parent, gp_data$id_red)]
  nest$bv_sd <- gp_data$bv_sd[match(nest$parent, gp_data$id_red)]

  nest$ringnr_num <- match(nest$ringnr, unique(nest$ringnr))
  nest$parent_num <- match(nest$parent, unique(nest$parent))
  nest$hi_num <- match(nest$hatch_island, unique(nest$hatch_island))
  nest$hy_num <- match(nest$hatch_year, unique(nest$hatch_year))

  # Add inbreeding info
  ibc <- fread(file = inbreeding)
  gp_data$fhat3 <- ibc$Fhat3[match(gp_data$ringnr, ibc$IID)]
  nest$fhat3 <- gp_data$fhat3[match(nest$parent, gp_data$id_red)]

  # # Number of crossover measurements
  nest$co_n <- gp_data$n[match(nest$parent, gp_data$id_red)]
  nest %<>% mutate(co_n = ifelse(is.na(co_n), 0, co_n))
  nest$co_meas <- (nest$co_n > 0)

  nest$idx <- seq_len(nrow(nest))
  nest
}

ars_adult_mod_func <- function(data,
                               bv_colname,
                               bv_samp = NULL,
                               gp_data = NULL) {
  hyperpar_var <- list(
    prec = list(initial = log(1),
                prior = "pc.prec",
                # sd, prob larger than sd
                param = c(sqrt(1), 0.05),
                fixed = FALSE))

  if (bv_colname == "bv_samp") {
    if ("parent" %in% colnames(data))
      data$bv_samp <- bv_samp[match(data$parent, gp_data$id_red)]
    else
      data$bv_samp <- bv_samp[match(data$ringnr, gp_data$id_red)]
  }

  data$bv <- getElement(data, bv_colname)
  data$bv2 <- (data$bv)^2

  effects_vec <- c("bv",
                   "bv2",
                   "age_q1",
                   "age_q2",
                   "fhat3",
                   "f(y_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ll_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ringnr_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(idx, model = \"iid\", hyper = hyperpar_var)")

  inla_formula <- reformulate(effects_vec, response = "sum_recruit")

  inla(inla_formula,
       family = "zeroinflatedpoisson1",
       control.compute = list(config = TRUE),
       control.family = list(hyper = list(theta = list(param = c(0, 1 / 1.75^2)))),
       data = data, verbose = TRUE) %>%
    INLA::inla.rerun()
}

surv_adult_mod_func <- function(data,
                                bv_colname,
                                bv_samp = NULL,
                                gp_data = NULL) {
  hyperpar_var <- list(
    prec = list(initial = log(1),
                prior = "pc.prec",
                # sd, prob larger than sd
                param = c(sqrt(1), 0.05),
                fixed = FALSE))

  if (bv_colname == "bv_samp") {
    if ("parent" %in% colnames(data))
      data$bv_samp <- bv_samp[match(data$parent, gp_data$id_red)]
    else
      data$bv_samp <- bv_samp[match(data$ringnr, gp_data$id_red)]
  }

  data$bv <- getElement(data, bv_colname)
  data$bv2 <- (data$bv)^2

  effects_vec <- c("bv",
                   "bv2",
                   "age_q1",
                   "age_q2",
                   "fhat3",
                   "f(y_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ll_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ringnr_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(idx, model = \"iid\", hyper = hyperpar_var)")

  inla_formula <- reformulate(effects_vec, response = "survival")

  inla(inla_formula,
       family = "binomial",
       control.compute = list(config = TRUE),
       control.family = list(link = "logit"),
       data = data,
       verbose = TRUE) %>%
    INLA::inla.rerun()
}

ars_parent_mod_func <- function(data,
                                bv_colname,
                                bv_samp = NULL,
                                gp_data = NULL) {
  hyperpar_var <- list(
    prec = list(initial = log(1),
                prior = "pc.prec",
                # sd, prob larger than sd
                param = c(sqrt(1), 0.05),
                fixed = FALSE))

  if (bv_colname == "bv_samp") {
    if ("parent" %in% colnames(data))
      data$bv_samp <- bv_samp[match(data$parent, gp_data$id_red)]
    else
      data$bv_samp <- bv_samp[match(data$ringnr, gp_data$id_red)]
  }

  data$bv <- getElement(data, bv_colname)
  data$bv2 <- (data$bv)^2

  effects_vec <- c("bv",
                   "bv2",
                   "age_q1",
                   "age_q2",
                   "fhat3",
                   "f(y_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ll_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ringnr_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(parent_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(idx, model = \"iid\", hyper = hyperpar_var)")

  inla_formula <- reformulate(effects_vec, response = "sum_recruit")

  inla(inla_formula,
       family = "zeroinflatedpoisson1",
       control.compute = list(config = TRUE),
       control.family = list(hyper = list(theta = list(param = c(0, 1 / 1.75^2)))),
       data = data, verbose = TRUE) %>%
    INLA::inla.rerun()
}

surv_parent_mod_func <- function(data,
                                 bv_colname,
                                 bv_samp = NULL,
                                 gp_data = NULL) {
  hyperpar_var <- list(
    prec = list(initial = log(1),
                prior = "pc.prec",
                # sd, prob larger than sd
                param = c(sqrt(1), 0.05),
                fixed = FALSE))

  if (bv_colname == "bv_samp") {
    if ("parent" %in% colnames(data))
      data$bv_samp <- bv_samp[match(data$parent, gp_data$id_red)]
    else
      data$bv_samp <- bv_samp[match(data$ringnr, gp_data$id_red)]
  }

  data$bv <- getElement(data, bv_colname)
  data$bv2 <- (data$bv)^2

  effects_vec <- c("bv",
                   "bv2",
                   "age_q1",
                   "age_q2",
                   "fhat3",
                   "f(y_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ll_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(ringnr_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(parent_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(idx, model = \"iid\", hyper = hyperpar_var)")

  inla_formula <- reformulate(effects_vec, response = "survival")

  inla(inla_formula,
       family = "binomial",
       control.compute = list(config = TRUE),
       control.family = list(link = "logit"),
       data = data,
       verbose = TRUE) %>%
    INLA::inla.rerun()
}

nest_mod_func <- function(data,
                          bv_colname,
                          bv_samp = NULL,
                          gp_data = NULL) {
  hyperpar_var <- list(
    prec = list(initial = log(1),
                prior = "pc.prec",
                # sd, prob larger than sd
                param = c(sqrt(1), 0.05),
                fixed = FALSE))

  if (bv_colname == "bv_samp") {
    if ("parent" %in% colnames(data))
      data$bv_samp <- bv_samp[match(data$parent, gp_data$id_red)]
    else
      data$bv_samp <- bv_samp[match(data$ringnr, gp_data$id_red)]
  }

  data$bv <- getElement(data, bv_colname)
  data$bv2 <- (data$bv)^2

  effects_vec <- c("bv",
                   "bv2",
                   "fhat3",
                   "f(hy_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(hi_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(parent_num, model = \"iid\", hyper = hyperpar_var)",
                   "f(idx, model = \"iid\", hyper = hyperpar_var)")

  inla_formula <- reformulate(effects_vec, response = "recruit")

  inla(inla_formula,
       family = "binomial",
       control.compute = list(config = TRUE),
       control.family = list(link = "logit"),
       data = data,
       verbose = TRUE) %>%
    INLA::inla.rerun()
}

ars_adult_parsamp <- function(s) {
  c("alpha"= (Intercept),
    "beta_bv" = bv,
    "beta_bv2" = bv2,
    "beta_age_q1" = age_q1,
    "beta_age_q2" = age_q2,
    "beta_f" = fhat3,
    "alpha_zi" = log(theta[1] / (1 - theta[1])),
    "sigma_ye" = sqrt(1 / theta[2]),
    "sigma_ll" = sqrt(1 / theta[3]),
    "sigma_id" = sqrt(1 / theta[4]),
    "sigma_res" = sqrt(1 / theta[5]))
}

surv_adult_parsamp <- function(s) {
  c("alpha"= (Intercept),
    "beta_bv" = bv,
    "beta_bv2" = bv2,
    "beta_age_q1" = age_q1,
    "beta_age_q2" = age_q2,
    "beta_f" = fhat3,
    "sigma_ye" = sqrt(1 / theta[1]),
    "sigma_ll" = sqrt(1 / theta[2]),
    "sigma_id" = sqrt(1 / theta[3]),
    "sigma_res" = sqrt(1 / theta[4]))
}

ars_parent_parsamp <- function(s) {
  c("alpha"= (Intercept),
    "beta_bv" = bv,
    "beta_bv2" = bv2,
    "beta_age_q1" = age_q1,
    "beta_age_q2" = age_q2,
    "beta_f" = fhat3,
    "alpha_zi" = log(theta[1] / (1 - theta[1])),
    "sigma_ye" = sqrt(1 / theta[2]),
    "sigma_ll" = sqrt(1 / theta[3]),
    "sigma_id" = sqrt(1 / theta[4]),
    "sigma_par" = sqrt(1 / theta[5]),
    "sigma_res" = sqrt(1 / theta[6]))
}

surv_parent_parsamp <- function(s) {
  c("alpha"= (Intercept),
    "beta_bv" = bv,
    "beta_bv2" = bv2,
    "beta_age_q1" = age_q1,
    "beta_age_q2" = age_q2,
    "beta_f" = fhat3,
    "sigma_ye" = sqrt(1 / theta[1]),
    "sigma_ll" = sqrt(1 / theta[2]),
    "sigma_id" = sqrt(1 / theta[3]),
    "sigma_par" = sqrt(1 / theta[4]),
    "sigma_res" = sqrt(1 / theta[5]))
}

nest_parsamp <- function(s) {
  c("alpha"= (Intercept),
    "beta_bv" = bv,
    "beta_bv2" = bv2,
    "beta_f" = fhat3,
    "sigma_hy" = sqrt(1 / theta[1]),
    "sigma_hi" = sqrt(1 / theta[2]),
    "sigma_par" = sqrt(1 / theta[3]),
    "sigma_res" = sqrt(1 / theta[4]))
}

samp_sampmod <- function(model, func, n_samp = 3) {

  inla.posterior.sample(n = n_samp, result = model, add.names = FALSE) %>%
    inla.posterior.sample.eval(fun = func, samples = .)
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
                   va_apriori = NULL,
                   comp_conf = FALSE) {

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
       control.compute = list(config = comp_conf),
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

get_samps <- function(model, pars) {
  # Extract posterior samples
  samps <- rstan::extract(model, pars = pars)

  samp_params <- get_sampler_params(model, inc_warmup = FALSE)
  samps$energy <- unlist(lapply(samp_params, function(x) x[, "energy__"]))
  samps$divergent <- unlist(lapply(samp_params, function(x) x[, "divergent__"]))

  samps
}

get_bfmi <- function(model) {
  samp_params <- get_sampler_params(model, inc_warmup = FALSE)

  sapply(samp_params, function(x) {
    E <- x[, "energy__"]
    sum((diff(E))^2) / sum((E - mean(E))^2)
  })
}

make_sim_ars_adult <- function(data,
                               stan_data,
                               pars) {

  bv <- rmvnorm(1, stan_data$bv_mean, stan_data$bv_covmat)[1, ]
  data$bv_true <- bv[stan_data$id_idx]

  ye <- rnorm(n = max(data$y_num), 0, sd = pars$sigma_ye)
  ll <- rnorm(n = max(data$ll_num), 0, sd = pars$sigma_ll)
  id <- rnorm(n = max(data$ringnr_num), 0, sd = pars$sigma_id)
  idx <- rnorm(n = max(data$idx), 0, sd = pars$sigma_res)

  ye_full <- ye[data$y_num]
  ll_full <- ll[data$ll_num]
  id_full <- id[data$ringnr_num]
  idx_full <- idx[data$idx]

  eta <- pars$alpha +
    pars$beta_bv * data$bv_true +
    pars$beta_bv2 * (data$bv_true)^2 +
    pars$beta_age_q1 * data$age_q1 +
    pars$beta_age_q2 * data$age_q2 +
    pars$beta_f * data$fhat3 +
    ye_full +
    ll_full +
    id_full +
    idx_full

  pois <- rpois(n = nrow(data), lambda = exp(eta))
  zeros <- rbinom(n = nrow(data), size = 1, prob = inv_logit(pars$alpha_zi))
  data$sum_recruit_orig <- data$sum_recruit
  data$sum_recruit <- pois * (1 - zeros)
  data
}

make_sim_surv_adult <- function(data,
                                stan_data,
                                pars) {

  bv <- rmvnorm(1, stan_data$bv_mean, stan_data$bv_covmat)[1, ]
  data$bv_true <- bv[stan_data$id_idx]

  ye <- rnorm(n = max(data$y_num), 0, pars$sigma_ye)
  ll <- rnorm(n = max(data$ll_num), 0, pars$sigma_ll)
  id <- rnorm(n = max(data$ringnr_num), 0, pars$sigma_id)
  idx <- rnorm(n = max(data$idx), 0, pars$sigma_res)

  ye_full <- ye[data$y_num]
  ll_full <- ll[data$ll_num]
  id_full <- id[data$ringnr_num]
  idx_full <- idx[data$idx]

  eta <- pars$alpha +
    pars$beta_bv * data$bv_true +
    pars$beta_bv2 * (data$bv_true)^2 +
    pars$beta_age_q1 * data$age_q1 +
    pars$beta_age_q2 * data$age_q2 +
    pars$beta_f * data$fhat3 +
    ye_full +
    ll_full +
    id_full +
    idx_full

  data$survival_orig <- data$survival
  data$survival <- rbinom(n = nrow(data), size = 1, prob = inv_logit(eta))
  data
}

make_sim_ars_parent <- function(data,
                                stan_data,
                                pars) {

  bv <- rmvnorm(1, stan_data$bv_mean, stan_data$bv_covmat)[1, ]
  data$bv_true <- bv[stan_data$par_idx]

  ye <- rnorm(n = max(data$y_num), 0, pars$sigma_ye)
  ll <- rnorm(n = max(data$ll_num), 0, pars$sigma_ll)
  id <- rnorm(n = max(data$ringnr_num), 0, pars$sigma_id)
  parent <- rnorm(n = max(data$parent_num), 0, pars$sigma_par)
  idx <- rnorm(n = max(data$idx), 0, pars$sigma_res)

  ye_full <- ye[data$y_num]
  ll_full <- ll[data$ll_num]
  id_full <- id[data$ringnr_num]
  parent_full <- parent[data$parent_num]
  idx_full <- idx[data$idx]

  eta <- pars$alpha +
    pars$beta_bv * data$bv_true +
    pars$beta_bv2 * (data$bv_true)^2 +
    pars$beta_age_q1 * data$age_q1 +
    pars$beta_age_q2 * data$age_q2 +
    pars$beta_f * data$fhat3 +
    ye_full +
    ll_full +
    id_full +
    parent_full +
    idx_full

  pois <- rpois(n = nrow(data), lambda = exp(eta))
  zeros <- rbinom(n = nrow(data), size = 1, prob = inv_logit(pars$alpha_zi))
  data$sum_recruit_orig <- data$sum_recruit
  data$sum_recruit <- pois * (1 - zeros)
  data
}

make_sim_surv_parent <- function(data,
                                 stan_data,
                                 pars) {

  bv <- rmvnorm(1, stan_data$bv_mean, stan_data$bv_covmat)[1, ]
  data$bv_true <- bv[stan_data$par_idx]

  ye <- rnorm(n = max(data$y_num), 0, pars$sigma_ye)
  ll <- rnorm(n = max(data$ll_num), 0, pars$sigma_ll)
  id <- rnorm(n = max(data$ringnr_num), 0, pars$sigma_id)
  parent <- rnorm(n = max(data$parent_num), 0, pars$sigma_par)
  idx <- rnorm(n = max(data$idx), 0, pars$sigma_res)

  ye_full <- ye[data$y_num]
  ll_full <- ll[data$ll_num]
  id_full <- id[data$ringnr_num]
  parent_full <- parent[data$parent_num]
  idx_full <- idx[data$idx]

  eta <- pars$alpha +
    pars$beta_bv * data$bv_true +
    pars$beta_bv2 * (data$bv_true)^2 +
    pars$beta_age_q1 * data$age_q1 +
    pars$beta_age_q2 * data$age_q2 +
    pars$beta_f * data$fhat3 +
    ye_full +
    ll_full +
    id_full +
    parent_full +
    idx_full

  data$survival_orig <- data$survival
  data$survival <- rbinom(n = nrow(data), size = 1, prob = inv_logit(eta))
  data
}

make_sim_nest <- function(data,
                          stan_data,
                          pars) {

  bv <- rmvnorm(1, stan_data$bv_mean, stan_data$bv_covmat)[1, ]
  data$bv_true <- bv[stan_data$par_idx]

  hy <- rnorm(n = max(data$hy_num), 0, pars$sigma_hy)
  hi <- rnorm(n = max(data$hi_num), 0, pars$sigma_hi)
  parent <- rnorm(n = max(data$parent_num), 0, pars$sigma_par)
  idx <- rnorm(n = max(data$idx), 0, pars$sigma_res)

  hy_full <- hy[data$hy_num]
  hi_full <- hi[data$hi_num]
  parent_full <- parent[data$parent_num]
  idx_full <- idx[data$idx]

  eta <- pars$alpha +
    pars$beta_bv * data$bv_true +
    pars$beta_bv2 * (data$bv_true)^2 +
    pars$beta_f * data$fhat3 +
    hi_full +
    hy_full +
    parent_full +
    idx_full

  data$recruit_orig <- data$recruit
  data$recruit <- rbinom(n = nrow(data), size = 1, prob = inv_logit(eta))
  data
}

samp_plot_df <- function(x, y, n_samp) {
  y %>%
    as.data.frame() %>%
    mutate(sample = seq_len(n_samp)) %>%
    tidyr::pivot_longer(-sample, names_to = "x_idx", values_to = "y_samp") %>%
    mutate(x = rep(x, n_samp),
           y_mean = rep(apply(y, 2, mean), n_samp),
           y_lower = rep(apply(y, 2, quantile, p = 0.025), n_samp),
           y_upper = rep(apply(y, 2, quantile, p = 0.975), n_samp))
}

inv_logit <- function(x) 1 / (1 + exp(-x))

make_ars_bv_preds_and_marg <- function(data,
                                       samp,
                                       n_plot = 200) {

  n_samp <- length(samp[[1]])
  avg_age <- mean(data$age)
  avg_f <- mean(data$fhat3)

  age_poly <- poly(data$age, degree = 2)
  avg_age_q <- predict(age_poly, newdata = mean(data$age))

  bv <- seq(min(c(data$bv_mean)), max(c(data$bv_mean)), length.out = n_plot)

  count_pred <- zinf_pred <- array(NA, dim = c(n_samp, n_plot))
  for (i in seq_len(n_samp)) {
    for (j in seq_len(n_plot)) {
      (samp$alpha[i] +
         samp$beta_bv[i] * bv[j] +
         samp$beta_bv2[i] * bv[j]^2 +
         samp$beta_age_q1[i] * avg_age_q[, 1] +
         samp$beta_age_q2[i] * avg_age_q[, 2] +
         samp$beta_f[i] * avg_f) %>%
        exp() ->
        count_pred[i, j]

      zinf_pred[i, j] <- inv_logit(samp$alpha_zi[i])

    }
  }
  y_pred <- count_pred * (1 - zinf_pred)

  df_count <- samp_plot_df(y = count_pred, x = bv, n_samp = n_samp)
  df_zinf <- samp_plot_df(y = zinf_pred, x = bv, n_samp = n_samp)
  df_pred <- samp_plot_df(y = y_pred, x = bv, n_samp = n_samp)

  marg_count <-
    sapply(bv, function(x) samp$beta_bv + 2 * samp$beta_bv2 * x) *
    count_pred
  marg_zinf <- matrix(0, nrow = n_samp, ncol = n_plot)
  # marg_zinf <-  sapply(bv, function(x) samp$beta_zi_bv + 2 * samp$beta_zi_bv2 * x) *
  # zinf_pred *
  # (1 - zinf_pred)
  marg_pred <- sapply(bv, function(x) {
    samp$beta_bv + # - samp$beta_zi_bv +
      2 * x * samp$beta_bv2 # - samp$beta_zi_bv2)
  }) * count_pred * (1 - zinf_pred)

  df_marg_count <- samp_plot_df(y = marg_count, x = bv, n_samp = n_samp)
  df_marg_zinf <- samp_plot_df(y = marg_zinf, x = bv, n_samp = n_samp)
  df_marg <- samp_plot_df(y = marg_pred, x = bv, n_samp = n_samp)

  lst(df_count, df_zinf, df_pred, df_marg_count, df_marg_zinf, df_marg)
}

# make_ars_bv_preds_and_marg_co_n <- function(data,
#                                             samps,
#                                             n_plot = 200) {
#   n_samp <- length(samps$alpha)
#   avg_age <- mean(data$age)
#   avg_f <- mean(data$f)
#   avg_co_n <- mean(data$co_n)
#
#   age_poly <- poly(data$age, degree = 2)
#   avg_age_q <- predict(age_poly, newdata = mean(data$age))
#
#   bv <- seq(min(c(data$bv_mean)), max(c(data$bv_mean)), length.out = n_plot)
#
#   count_pred <- zinf_pred <- array(NA, dim = c(n_samp, n_plot))
#   for (i in seq_len(n_samp)) {
#     for (j in seq_len(n_plot)) {
#       (samps$alpha[i] +
#          samps$beta_bv[i] * bv[j] +
#          samps$beta_bv2[i] * bv[j]^2 +
#          samps$beta_age_q1[i] * avg_age_q[, 1] +
#          samps$beta_age_q2[i] * avg_age_q[, 2] +
#          samps$beta_f[i] * avg_f +
#          samps$beta_co_n[i] * avg_co_n) %>%
#         exp() ->
#         count_pred[i, j]
#
#       (samps$alpha_zi[i]# +
#         # samps$beta_zi_bv[i] * bv[j] +
#         # samps$beta_zi_bv2[i] * bv[j]^2 +
#         # samps$beta_zi_age_q1[i] * avg_age_q[, 1] +
#         # samps$beta_zi_age_q2[i] * avg_age_q[, 2] +
#         # samps$beta_zi_f[i] * avg_f
#       ) %>%
#         inv_logit() ->
#         zinf_pred[i, j]
#     }
#   }
#   y_pred <- count_pred * (1 - zinf_pred)
#
#   df_count <- samp_plot_df(y = count_pred, x = bv, n_samp = n_samp)
#   df_zinf <- samp_plot_df(y = zinf_pred, x = bv, n_samp = n_samp)
#   df_pred <- samp_plot_df(y = y_pred, x = bv, n_samp = n_samp)
#
#   marg_count <- sapply(bv, function(x) samps$beta_bv + 2 * samps$beta_bv2 * x) *
#     count_pred
#   marg_zinf <- matrix(0, nrow = n_samp, ncol = n_plot)
#   # marg_zinf <- sapply(bv, function(x) samps$beta_zi_bv + 2 * samps$beta_zi_bv2 * x) *
#   #   zinf_pred *
#   #   (1 - zinf_pred)
#   marg_pred <- sapply(bv, function(x) {
#     samps$beta_bv + # - samps$beta_zi_bv +
#       2 * x * (samps$beta_bv2)# - samps$beta_zi_bv2)
#   }) * count_pred * (1 - zinf_pred)
#
#   df_marg_count <- samp_plot_df(y = marg_count, x = bv, n_samp = n_samp)
#   df_marg_zinf <- samp_plot_df(y = marg_zinf, x = bv, n_samp = n_samp)
#   df_marg_pred <- samp_plot_df(y = marg_pred, x = bv, n_samp = n_samp)
#
#   lst(df_count, df_zinf, df_pred, df_marg_count, df_marg_zinf, df_marg_pred)
# }

make_surv_bv_preds_and_marg <- function(data,
                                        samp,
                                        n_plot = 200) {
  n_samp <- length(samp[[1]])
  avg_age <- mean(data$age)
  avg_f <- mean(data$fhat3)

  age_poly <- poly(data$age, degree = 2)
  avg_age_q <- predict(age_poly, newdata = mean(data$age))

  bv <- seq(min(c(data$bv_mean)), max(c(data$bv_mean)), length.out = n_plot)
  predictor_samples <- array(NA, dim = c(n_samp, n_plot))

  for (i in seq_len(n_samp)) {
    for (j in seq_len(n_plot)) {
      predictor_samples[i, j] <- samp$alpha[i] +
        samp$beta_bv[i] * bv[j] +
        samp$beta_bv2[i] * bv[j]^2 +
        samp$beta_age_q1[i] * avg_age_q[, 1] +
        samp$beta_age_q2[i] * avg_age_q[, 2] +
        samp$beta_f[i] * avg_f
    }
  }

  pred_prob_samples <- inv_logit(predictor_samples)
  df_pred <- samp_plot_df(y = pred_prob_samples, x = bv, n_samp = n_samp)

  marg_effect <-
    sapply(bv, function(x) samp$beta_bv + 2 * samp$beta_bv2 * x) *
    pred_prob_samples *
    (1 - pred_prob_samples)

  df_marg <- samp_plot_df(y = marg_effect, x = bv, n_samp = n_samp)

  lst(df_pred, df_marg)
}

# make_surv_bv_preds_and_marg_co_n <- function(data,
#                                              samps,
#                                              n_plot = 200) {
#   n_samp <- length(samps$alpha)
#   avg_age <- mean(data$age)
#   avg_f <- mean(data$f)
#   avg_co_n <- mean(data$co_n)
#
#   age_poly <- poly(data$age, degree = 2)
#   avg_age_q <- predict(age_poly, newdata = mean(data$age))
#
#   bv <- seq(min(c(data$bv_mean)), max(c(data$bv_mean)), length.out = n_plot)
#
#   predictor_samples <- array(NA, dim = c(n_samp, n_plot))
#   for (i in seq_len(n_samp)) {
#     for (j in seq_len(n_plot)) {
#       predictor_samples[i, j] <- samps$alpha[i] +
#         samps$beta_bv[i] * bv[j] +
#         samps$beta_bv2[i] * bv[j]^2 +
#         samps$beta_age_q1[i] * avg_age_q[, 1] +
#         samps$beta_age_q2[i] * avg_age_q[, 2] +
#         samps$beta_f[i] * avg_f +
#         samps$beta_co_n[i] * avg_co_n
#     }
#   }
#   pred_prob_samples <- inv_logit(predictor_samples)
#   df_pred <- samp_plot_df(y = pred_prob_samples, x = bv, n_samp = n_samp)
#
#   marg_effect <-
#     sapply(bv, function(x) samps$beta_bv + 2 * samps$beta_bv2 * x) *
#     pred_prob_samples *
#     (1 - pred_prob_samples)
#   df_marg <- samp_plot_df(y = marg_effect, x = bv, n_samp = n_samp)
#
#   lst(df_pred, df_marg)
# }

make_nest_bv_preds_and_marg <- function(data,
                                        samp,
                                        n_plot = 200) {
  n_samp <- length(samp[[1]])
  avg_f <- mean(data$fhat3)

  bv <- seq(min(c(data$bv_mean)), max(c(data$bv_mean)), length.out = n_plot)

  predictor_samples <- array(NA, dim = c(n_samp, n_plot))
  for (i in seq_len(n_samp)) {
    for (j in seq_len(n_plot)) {
      predictor_samples[i, j] <- samp$alpha[i] +
        samp$beta_bv[i] * bv[j] +
        samp$beta_bv2[i] * bv[j]^2 +
        samp$beta_f[i] * avg_f
    }
  }
  pred_prob_samples <- inv_logit(predictor_samples)
  df_pred <- samp_plot_df(y = pred_prob_samples, x = bv, n_samp = n_samp)

  marg_effect <-
    sapply(bv, function(x) samp$beta_bv + 2 * samp$beta_bv2 * x) *
    pred_prob_samples *
    (1 - pred_prob_samples)
  df_marg <- samp_plot_df(y = marg_effect, x = bv, n_samp = n_samp)

  lst(df_pred, df_marg)
}

plot_lines_posterior <- function(df,
                                 xlab,
                                 ylab,
                                 title,
                                 data = NULL,
                                 n_bins = 25,
                                 bs = 11,
                                 lw = 1,
                                 leg.x = 0.9) {

  if (!is.null(data)) {
    # Stats for sample size histogram hist
    breaks <- seq(from = range(data$bv_mean)[1],
                  range(data$bv_mean)[2],
                  length = n_bins)
    bv_int <- findInterval(sort(data$bv_mean), vec = breaks)
    bv_count <- sapply(1:n_bins, function(i) sum(bv_int == i))
    x_int <- findInterval(df$x, vec = breaks)
    scaling <- max(bv_count) / max(df$y_samp)
    df$hist_height <- bv_count[x_int] / scaling
  }

  plot <- ggplot(df, aes(x = x)) +
    # Transparent lines for posterior samples
    geom_line(aes(y = y_samp,
                  group = sample,
                  color = "Posterior samples",
                  linetype = "Posterior samples"),
              alpha = 0.01,
              linewidth = 0.5 * lw) +
    # Posterior mean
    geom_line(aes(y = y_mean,
                  color = "Posterior mean",
                  linetype = "Posterior mean"),
              linewidth = 1 * lw) +
    # Posterior 95% CI
    geom_line(aes(y = y_lower, color = "95% CI", linetype = "95% CI"),
              linewidth = 0.5 * lw) +
    geom_line(aes(y = y_upper, color = "95% CI", linetype = "95% CI"),
              linewidth = 0.5 * lw) +
    # Adding a title and axis labels
    labs(x = xlab,
         y = ylab,
         title = title,
         color = "",
         linetype = "") + # Combine legends
    # Minimal theme for clean appearance
    theme_minimal(base_size = bs) +
    # Customize the color and linetype scales
    scale_color_manual(values = c("Posterior mean" = "blue",
                                  "95% CI" = "darkred",
                                  "Posterior samples" = "black")) +
    scale_linetype_manual(values = c("Posterior mean" = "solid",
                                     "95% CI" = "dashed",
                                     "Posterior samples" = "solid")) +
    # Ensure the correct legend appearance
    guides(color = guide_legend(
      override.aes = list(
        # Match transparency for samples, full for others
        alpha = c(1, 1, 0.03),
        # Make legend dashed looked like in plot
        linewidth = lw * c(0.5, 1, 0.5))),
      linetype = guide_legend(override.aes = list(
        # Correct dashing in legend
        linetype = c("dashed", "solid", "solid")))) +
    theme(legend.position = c(leg.x, 0.9),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA)) +
    scale_x_continuous(expand = c(0, 0))

  if (!is.null(data)) {
    plot <- plot +
      geom_ribbon(aes(ymin = 0, ymax = hist_height),
                  alpha = 0.25,
                  fill = "red") +
      scale_y_continuous(
        name = ylab,
        # Add a second y-axis
        sec.axis = sec_axis(trans = ~ . * scaling,
                            name = "Number of individuals",
                            labels = function(x) {
                              ifelse(x >= 0, as.character(x), "")
                            })) # +
    # guides(alpha = guide_legend(title = waiver(),
    #                             label = "Number of individuals"))
  }
  plot
}

ggsave_path <- function(filename,
                        ...) {
  ggsave(filename, ...)
  filename
}

make_stan_data_adult <- function(data, gp_data, covmat) {

  # Make bv stats same order as ringnr, and give one entry per ind.
  bv_mean <- data$bv_mean[order(data$ringnr_num)][
    match(unique(data[order(data$ringnr_num)]$ringnr),
          data[order(data$ringnr_num)]$ringnr)]

  # All the same reordering that was done for bv_mean
  bv_covmat <- covmat[
    gp_data$id1, gp_data$id1][
      match(data$ringnr, gp_data$id_red),
      match(data$ringnr, gp_data$id_red)][
        order(data$ringnr_num), order(data$ringnr_num)][
          match(unique(data[order(data$ringnr_num)]$ringnr),
                data[order(data$ringnr_num)]$ringnr),
          match(unique(data[order(data$ringnr_num)]$ringnr),
                data[order(data$ringnr_num)]$ringnr)]

  # sample, and add repeats
  bv_std_vec <- rmvnorm(1e4, bv_mean, bv_covmat)[, data$ringnr_num]

  list(N = nrow(data),
       sex = data$sex,
       N_ll = max(data$ll_num),
       N_ye = max(data$y_num),
       N_id = max(data$ringnr_num),
       ye_idx = data$y_num,
       ll_idx = data$ll_num,
       id_idx = data$ringnr_num,
       bv_mean = bv_mean,
       bv_covmat = bv_covmat,
       bv_covmat_chol = t(chol(bv_covmat)), # pre-multiply this with z-vec
       sum_recruit = data$sum_recruit,
       sum_recruit_log_mean = log(mean(data$sum_recruit)),
       survival = data$survival,
       age = data$age,
       age_q1 = poly(data$age, degree = 2)[, 1],
       age_q2 = poly(data$age, degree = 2)[, 2],
       f = data$fhat3,
       bv_mean_std = mean(apply(bv_std_vec, 2, mean)),
       bv_sd_std = mean(apply(bv_std_vec, 2, sd)),
       co_n = data$co_n,
       co_meas = data$co_meas,
       beta_prior_sd_ars = 0.2,
       exp_rate_ars = 1 / 0.2,
       beta_prior_sd_surv = 0.5,
       exp_rate_surv = 1 / 0.5,
       beta_zi_prior_sd = 0.5,
       exp_rate_zi = 1 / 0.5,
       alpha_zi_prior_mean = log(mean(data$sum_recruit == 0) /
                                   (1 - mean(data$sum_recruit == 0))),
       alpha_prior_mean_surv = log(1 / (1 / mean(data$survival) - 1)),
       alpha_prior_mean_ars = log(mean(data$sum_recruit[data$sum_recruit != 0]))
  )
}

make_stan_data_parent <- function(data, gp_data, covmat) {

  # Make bv stats same order as ringnr, and give one entry per ind.
  bv_mean <- data$bv_mean[order(data$parent_num)][
    match(unique(data[order(data$parent_num)]$parent),
          data[order(data$parent_num)]$parent)]

  # All the same reordering that was done for bv_mean
  bv_covmat <- covmat[
    gp_data$id1, gp_data$id1][
      match(data$parent, gp_data$id_red),
      match(data$parent, gp_data$id_red)][
        order(data$parent_num), order(data$parent_num)][
          match(unique(data[order(data$parent_num)]$parent),
                data[order(data$parent_num)]$parent),
          match(unique(data[order(data$parent_num)]$parent),
                data[order(data$parent_num)]$parent)]

  # sample, and add repeats
  bv_std_vec <- rmvnorm(1e4, bv_mean, bv_covmat)[, data$parent_num]

  idx <-

    list(N = nrow(data),
         sex = data$sex,
         N_ll = max(data$ll_num),
         N_ye = max(data$y_num),
         N_id = max(data$ringnr_num),
         N_par = max(data$parent_num),
         ye_idx = data$y_num,
         ll_idx = data$ll_num,
         id_idx = data$ringnr_num,
         par_idx = data$parent_num,
         id_to_par_idx = data$parent_num[match(seq_len(max(data$ringnr_num)),
                                               data$ringnr_num)],
         bv_mean = bv_mean,
         bv_covmat = bv_covmat,
         bv_covmat_chol = t(chol(bv_covmat)), # pre-multiply this with z-vec
         sum_recruit = data$sum_recruit,
         sum_recruit_log_mean = log(mean(data$sum_recruit)),
         survival = data$survival,
         survival_logit_mean = log(1 / (1 / mean(data$survival) - 1)),
         age = data$age,
         age_q1 = poly(data$age, degree = 2)[, 1],
         age_q2 = poly(data$age, degree = 2)[, 2],
         f = data$fhat3,
         bv_mean_std = mean(apply(bv_std_vec, 2, mean)),
         bv_sd_std = mean(apply(bv_std_vec, 2, sd)),
         co_n = data$co_n,
         co_meas = data$co_meas,
         beta_prior_sd_ars = 0.2,
         exp_rate_ars = 1 / 0.2,
         beta_prior_sd_surv = 0.5,
         exp_rate_surv = 1 / 0.5,
         beta_zi_prior_sd = 0.5,
         exp_rate_zi = 1 / 0.5,
         alpha_zi_prior_mean = log(mean(data$sum_recruit == 0) /
                                     (1 - mean(data$sum_recruit == 0))),
         alpha_prior_mean_surv = log(1 / (1 / mean(data$survival) - 1)),
         alpha_prior_mean_ars = log(mean(data$sum_recruit[data$sum_recruit != 0])))
}

make_stan_data_nest <- function(data, gp_data, covmat) {

  # Make bv stats same order as ringnr, and give one entry per parent
  bv_mean <- data$bv_mean[order(data$parent_num)][
    match(unique(data[order(data$parent_num)]$parent),
          data[order(data$parent_num)]$parent)]

  # All the same reordering that was done for bv_mean
  bv_covmat <- covmat[
    gp_data$id1, gp_data$id1][
      match(data$parent, gp_data$id_red),
      match(data$parent, gp_data$id_red)][
        order(data$parent_num), order(data$parent_num)][
          match(unique(data[order(data$parent_num)]$parent),
                data[order(data$parent_num)]$parent),
          match(unique(data[order(data$parent_num)]$parent),
                data[order(data$parent_num)]$parent)]

  # sample, and add repeats
  bv_std_vec <- rmvnorm(1e4, bv_mean, bv_covmat)[, data$parent_num]

  list(N = nrow(data),
       N_hi = max(data$hi_num),
       N_hy = max(data$hy_num),
       N_id = max(data$ringnr_num),
       N_par = max(data$parent_num),
       hi_idx = data$hi_num,
       hy_idx = data$hy_num,
       id_idx = data$ringnr_num,
       par_idx = data$parent_num,
       bv_mean = bv_mean,
       bv_covmat = bv_covmat,
       bv_covmat_chol = t(chol(bv_covmat)), # pre-multiply this with z-vec
       recruit = data$recruit,
       alpha_prior_mean_nestling = log(1 / (1 / mean(data$recruit) - 1)),
       f = data$fhat3,
       bv_mean_std = mean(apply(bv_std_vec, 2, mean)),
       bv_sd_std = mean(apply(bv_std_vec, 2, sd)),
       co_n = data$co_n,
       co_meas = data$co_meas,
       beta_prior_sd_ars = 0.2,
       exp_rate_ars = 1 / 0.2,
       beta_prior_sd_surv = 0.5,
       exp_rate_surv = 1 / 0.5,
       beta_zi_prior_sd = 0.5,
       exp_rate_zi = 1 / 0.5)
}

do_ars_ppc <- function(y, yrep, ll_i, ye_i, id_i) {

  lst(mean = ppc_stat(y, yrep),
      sd = ppc_stat(y, yrep, stat = "sd"),
      zeros = ppc_stat(y, yrep, stat = function(y) mean(y == 0)),
      ones = ppc_stat(y, yrep, stat = function(y) mean(y == 1)),
      twos = ppc_stat(y, yrep, stat = function(y) mean(y == 2)),
      threes = ppc_stat(y, yrep, stat = function(y) mean(y == 3)),
      fours = ppc_stat(y, yrep, stat = function(y) mean(y == 4)),
      fives = ppc_stat(y, yrep, stat = function(y) mean(y == 5)),
      sixes = ppc_stat(y, yrep, stat = function(y) mean(y == 6)),
      sevens = ppc_stat(y, yrep, stat = function(y) mean(y == 7)),
      eights = ppc_stat(y, yrep, stat = function(y) mean(y == 8)),
      bar = ppc_bars(y, yrep),
      zeros_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 0), group = ll_i),
      ones_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 1), group = ll_i),
      twos_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 2), group = ll_i),
      threes_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 3), group = ll_i),
      fours_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 4), group = ll_i),
      fives_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 5), group = ll_i),
      sixes_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 6), group = ll_i),
      sevens_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 7), group = ll_i),
      eights_ll = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 8), group = ll_i),
      sd_ll = ppc_stat_grouped(y, yrep, stat = "sd", group = ll_i),
      zeros_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 0), group = ye_i),
      ones_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 1), group = ye_i),
      twos_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 2), group = ye_i),
      threes_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 3), group = ye_i),
      fours_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 4), group = ye_i),
      fives_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 5), group = ye_i),
      sixes_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 6), group = ye_i),
      sevens_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 7), group = ye_i),
      eights_ye = ppc_stat_grouped(y, yrep, stat = function(y) mean(y == 8), group = ye_i),
      sd_ye = ppc_stat_grouped(y, yrep, stat = "sd", group = ye_i),
      p_mean = mean(colMeans(yrep) > mean(y)),
      p_sd = mean(apply(yrep, 2, sd) > sd(y)),
      p_zeros = mean((apply(yrep, 2, function(y) mean(y == 0))) > mean(y == 0)),
      p_ones = mean((apply(yrep, 2, function(y) mean(y == 1))) > mean(y == 1)),
      p_twos = mean((apply(yrep, 2, function(y) mean(y == 2))) > mean(y == 2)))
}

do_surv_ppc <- function(y, yrep, co_n, co_meas, ll_i, ye_i, id_i) {

  lst(bar = ppc_bars(y, yrep),
      bar_ll = ppc_bars_grouped(y, yrep, group = ll_i),
      bar_ye = ppc_bars_grouped(y, yrep, group = ye_i),
      # bar_co_meas = ppc_bars_grouped(y, yrep, group = co_meas),
      # bar_co_n = ppc_bars_grouped(y, yrep, group = co_n),
      bar_longlife = ppc_bars_grouped(
        y, yrep, group = factor(tapply(y, factor(id_i), length)[id_i] > 2)),
      sd = ppc_stat(y, yrep, stat = "sd"),
      sd_ll = ppc_stat_grouped(y, yrep, group = ll_i, stat = "sd"),
      sd_ye = ppc_stat_grouped(y, yrep, group = ye_i, stat = "sd"),
      # sd_co_meas = ppc_stat_grouped(y, yrep, group = co_meas, stat = "sd"),
      # sd_co_n = ppc_stat_grouped(y, yrep, group = co_n, stat = "sd"),
      p_mean = mean(colMeans(yrep) > mean(y)),
      p_sd = mean(apply(yrep, 2, sd) > sd(y)))
}

inla_bv_covmat <- function(model, n_samp = 1e4, ncores) {

  covmat <- INLA::inla.posterior.sample(n = n_samp,
                                        result = model,
                                        add.names = FALSE,
                                        num.threads = ncores) %>%
    INLA::inla.posterior.sample.eval(fun = function() id1) %>%
    t() %>%
    cov()
  # Same order as pheno_data
  covmat[order(model$summary.random$id1$ID),
         order(model$summary.random$id1$ID)]
}

make_sim_bv_plot <- function(summ,
                             sim_data) {

  y <- summ %>%
    `[`(grepl(x = rownames(.), pattern = "bv_lat"), "mean")

  if ("parent" %in% colnames(sim_data[[1]])) {
    x <- sim_data %>%
      (function(dat) {
        dat <- dat[[1]]
        dat %>%
          `$`("bv_true") %>%
          `[`(order(dat$parent_num)) %>%
          `[`(match(unique(dat[order(dat$parent_num)]$parent),
                    dat[order(dat$parent_num)]$parent))
      })
  } else {
    x <- sim_data %>%
      (function(dat) {
        dat <- dat[[1]]
        dat %>%
          `$`("bv_true") %>%
          `[`(order(dat$ringnr_num)) %>%
          `[`(match(unique(dat[order(dat$ringnr_num)]$ringnr),
                    dat[order(dat$ringnr_num)]$ringnr))
      })
  }

  ggplot(mapping = aes(x = x, y = y)) +
    geom_point(pch = 16) +
    xlab("true bv") +
    ylab ("posterior mean bv") +
    stat_smooth(formula = y ~ x, method = "lm") +
    geom_abline(slope = 1, intercept = 0) +
    theme_minimal() +
    coord_fixed(ratio = 1) +
    scale_x_continuous(limits = range(x, y)) +
    scale_y_continuous(limits = range(x, y)) +
    theme(panel.border = element_rect(fill = NA),
          panel.grid.minor = element_blank())
}
