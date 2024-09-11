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
library(pedtricks)

load("results/2_F_fitCpi_results.RData")

pedsumm <- pedtricks::ped_stats(pedigree)
summary(pedsumm)

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


