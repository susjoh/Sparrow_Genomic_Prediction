#
# Genomic Prediction of Sparrow Recombination Rate
# Susan Johnston
# 2024-09-10
#
# NB. This analysis focusses on females only!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up the working environment & load data #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(hibayes)
library(dplyr)
library(ggplot2)

# Read in pedigree

pedigree <- read.table("data/20230317_Sparrow_Pedigree.txt", header = T)

# Read in genomic data

bin <- read_plink(bfile = "data/70K_200K_maf_geno_mind_v5", mode = "A", threads = 4)

fam <- bin[["fam"]]
geno <- bin[["geno"]]
map <- bin[["map"]]
geno.id <- fam[, 2]

# Read in phenotype data and subset to just females (FULL data)

recfull <- read.table("data/20240910_Sparrow_Recomb_Data.txt", header = T)
recfull <- na.omit(recfull)
recfull <- subset(recfull, sex == "F")
recfull$total_coverage <- recfull$total_coverage/1e6
recfull$meiosis <- NULL

# Make individual mean phenotype data (MEAN data)

recmeans <- recfull %>%
  group_by(id, sex) %>%
  summarise(mean_co_count = mean(co_count),
            mean_intra_shuff = mean(intra_shuff),
            mean_total_coverage = mean(total_coverage),
            var_co_count = var(co_count),
            var_intra_shuff = var(intra_shuff),
            var_total_coverage = var(total_coverage),
            n = n())

# Population struction PCA

# system("plink --bfile data/70K_200K_maf_geno_mind_v5 --cow --autosome --pca --out data/70K_200K_maf_geno_mind_v5")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Basic demo model in BayesCpi               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 2000 iterations is ~98 sec

full_gebvs <- list()

# FULL dataset ###################################################

fitCpi_full <- ibrm(co_count ~ total_coverage,
                    data = recfull,
                    M = geno,
                    M.id = geno.id,
                    method = "BayesCpi",
                    Pi = c(0.98, 0.02), niter = 20000, nburn = 16000, thin = 5,
                    printfreq = 100, seed = 666666, verbose = TRUE)

summary(fitCpi_full)
plot(fitCpi_full$MCMCsamples$h2[1,], type = "l")

#~~ match gebvs to phenos

gebv <- fitCpi_full[["g"]]
names(gebv)[2] <- "fitCpi_full_gebv"
gebv$id <- as.character(gebv$id)

full_gebvs[[length(full_gebvs)+1]] <- gebv

recfull <- left_join(recfull, gebv)

ggplot(recfull, aes(fitCpi_full_gebv, co_count)) +
  geom_point(alpha = 0.2) +
  stat_smooth()

cor.test(recfull$fitCpi_full_gebv, recfull$co_count)

# save(fitCpi_full, file = "results/2_F_fitCpi_full_gebv.RData")

# MEAN dataset ###################################################

fitCpi_mean <- ibrm(mean_co_count ~ mean_total_coverage + n,
               data = recmeans,
               M = geno,
               M.id = geno.id,
               method = "BayesCpi",
               Pi = c(0.98, 0.02), niter = 20000, nburn = 16000, thin = 5,
               printfreq = 100, seed = 666666, verbose = TRUE)

summary(fitCpi_mean)
plot(fitCpi_mean$MCMCsamples$h2[1,], type = "l")

#~~ match gebvs to phenos

gebv <- fitCpi_mean[["g"]]
names(gebv)[2] <- "fitCpi_mean_gebv"
gebv$id <- as.character(gebv$id)

full_gebvs[[length(full_gebvs)+1]] <- gebv

recmeans <- left_join(recmeans, gebv)

ggplot(recmeans, aes(fitCpi_mean_gebv, mean_co_count, col = n)) +
  geom_point(alpha = 0.2) +
  stat_smooth() +
  scale_color_gradient(low = "red", high = "white")

cor.test(recmeans$fitCpi_mean_gebv, recmeans$mean_co_count)

# save(fitCpi_mean, file = "results/2_F_fitCpi_mean_gebv.RData")

# MEAN dataset with > 2 measures ###################################################

fitCpi_mean_g2 <- ibrm(mean_co_count ~ mean_total_coverage + n,
                    data = subset(recmeans, n > 2),
                    M = geno,
                    M.id = geno.id,
                    method = "BayesCpi",
                    Pi = c(0.98, 0.02), niter = 20000, nburn = 16000, thin = 5,
                    printfreq = 100, seed = 666666, verbose = TRUE)

summary(fitCpi_mean_g2)
plot(fitCpi_mean_g2$MCMCsamples$h2[1,], type = "l")


#~~ match gebvs to phenos

gebv <- fitCpi_mean_g2[["g"]]
names(gebv)[2] <- "fitCpi_mean_g2_gebv"
gebv$id <- as.character(gebv$id)

full_gebvs[[length(full_gebvs)+1]] <- gebv


recmeans <- left_join(recmeans, gebv)

ggplot(recmeans, aes(fitCpi_mean_g2_gebv, mean_co_count, col = n > 3)) +
  geom_point(alpha = 0.2) +
  stat_smooth()

cor.test(recmeans$fitCpi_mean_g2_gebv, recmeans$mean_co_count)

# save(fitCpi_mean, file = "results/2_F_fitCpi_mean_g2_gebv.RData")


# SAMPLE ONE GAMETE dataset ###################################################

recsamp1 <- recfull %>%
  group_by(id) %>%
  summarise(co_count = sample(co_count, 1),
            total_coverage = mean(total_coverage))   # slight hack


fitCpi_samp1 <- ibrm(co_count ~ total_coverage,
                    data = recsamp1,
                    M = geno,
                    M.id = geno.id,
                    method = "BayesCpi",
                    Pi = c(0.98, 0.02), niter = 20000, nburn = 16000, thin = 5,
                    printfreq = 100, seed = 666666, verbose = TRUE)

summary(fitCpi_samp1)
plot(fitCpi_samp1$MCMCsamples$h2[1,], type = "l")

#~~ match gebvs to phenos

gebv <- fitCpi_samp1[["g"]]
names(gebv)[2] <- "fitCpi_samp1_gebv"
gebv$id <- as.character(gebv$id)

full_gebvs[[length(full_gebvs)+1]] <- gebv

recmeans <- left_join(recmeans, gebv)

ggplot(recmeans, aes(fitCpi_samp1_gebv, mean_co_count)) +
  geom_point(alpha = 0.2) +
  stat_smooth(method = "lm")

cor.test(recmeans$fitCpi_samp1_gebv, recmeans$mean_co_count)

# save(fitCpi_samp1, file = "results/2_F_fitCpi_samp1_gebv.RData")

rm(bin, fam, geno, map, geno.id, gebv, pedigree, recsumm)

save.image("results/2_F_fitCpi_results.RData")

