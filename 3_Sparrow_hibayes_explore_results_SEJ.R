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

load("results/2_F_fitCpi_results.RData")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Plot GEBVs against known phenotypes        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# FULL model

ggplot(recfull, aes(fitCpi_full_gebv, co_count)) +
  geom_point(alpha = 0.2) +
  stat_smooth()

cor.test(recfull$fitCpi_full_gebv, recfull$co_count)

# MEANS model

ggplot(recmeans, aes(fitCpi_mean_gebv, mean_co_count, col = n)) +
  geom_point(alpha = 0.2) +
  stat_smooth()

cor.test(recmeans$fitCpi_mean_gebv, recmeans$mean_co_count)

# MEANS >2 samples

ggplot(recmeans, aes(fitCpi_mean_g2_gebv, mean_co_count)) +
  geom_point(alpha = 0.2) +
  stat_smooth() +
  facet_wrap(~n > 3)

cor.test(recmeans$fitCpi_mean_g2_gebv, recmeans$mean_co_count)

# FULL vs MEANS

full_gebvs_tab <- left_join(full_gebvs[[1]], full_gebvs[[2]])
full_gebvs_tab <- left_join(full_gebvs_tab, full_gebvs[[3]])
full_gebvs_tab <- left_join(full_gebvs_tab, recmeans[,1:9])
full_gebvs_tab$n[which(is.na(full_gebvs_tab$n))] <- 0

ggplot(full_gebvs_tab, aes(fitCpi_full_gebv, fitCpi_mean_gebv)) +
  geom_point(alpha = 0.2) +
  stat_smooth() +
  facet_wrap(~id %in% recmeans$id)

cor.test(full_gebvs_tab$fitCpi_full_gebv, full_gebvs_tab$fitCpi_mean_gebv)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Sorting out disruptive selection on GEBVs  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Look at relationship between gebv and n:

ggplot(full_gebvs_tab, aes(fitCpi_full_gebv, `n`)) +
  geom_point(alpha = 0.2) +
  stat_smooth()

# Are zeros less related to the population? unconnected in the pedigree? CLUES?
# What is the population structure like?

# system("plink --bfile data/70K_200K_maf_geno_mind_v5 --cow --autosome --pca --out data/70K_200K_maf_geno_mind_v5")

sparrowpca <- read.table("data/70K_200K_maf_geno_mind_v5.eigenvec")
names(sparrowpca) <- c("fam", "id", paste0("PC", 1:20))


full_gebvs_tab <- left_join(full_gebvs_tab, sparrowpca[,2:8])

ggplot(full_gebvs_tab, aes(PC1, PC2, col = n>0)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~n>0)

# island effect?

islands <- read.csv("data/SNPtypedind_PhenotypicData.csv", header = T)
islands$id <- islands$ringnr
full_gebvs_tab <- left_join(full_gebvs_tab, subset(islands, select = c(id, first_locality, hatchyear)))


ggplot(subset(full_gebvs_tab, !is.na(first_locality)), aes(PC1, PC2, col = factor(first_locality))) +
  geom_point(alpha = 0.2)



