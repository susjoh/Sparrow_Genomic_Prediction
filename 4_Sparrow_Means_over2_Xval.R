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

subdata <- subset(recmeans, n > 2)
subdata <- subdata[sample(nrow(subdata), replace = F)]
subdata$order <- 1:nrow(subdata)
ntoremove <- 79
full_gebvs <- list()


for(i in seq(1, nrow(subdata)-ntoremove, ntoremove)){
  message(paste0("Run ", i))

  fitCpi_mean_g2 <- ibrm(mean_co_count ~ mean_total_coverage + n,
                         data = subset(subdata, !order %in% i:(i + ntoremove-1)),
                         M = geno,
                         M.id = geno.id,
                         method = "BayesCpi",
                         Pi = c(0.98, 0.02), niter = 2000, nburn = 1600, thin = 5,
                         printfreq = 100, seed = 666666, verbose = TRUE)

  gebv <- fitCpi_mean_g2[["g"]]
  gebv$id <- as.character(gebv$id)
  gebv$run <- i
  gebv$indata <- "no"
  gebv$indata <- ifelse(gebv$id %in% subset(subdata, !order %in% i:(i + ntoremove-1))$id, "yes", gebv$indata)
  gebv$indata <- ifelse(gebv$id %in% subset(subdata, order %in% i:(i + ntoremove-1))$id, "xval", gebv$indata)

  full_gebvs[[length(full_gebvs)+1]] <- gebv

}

full_gebvs_hold <- full_gebvs

full_gebvs <- bind_rows(full_gebvs)
full_gebvs <- left_join(full_gebvs, subdata[,c("id", "mean_co_count")])
full_gebvs <- subset(full_gebvs, indata != "no")

ggplot(full_gebvs, aes(gebv, mean_co_count, col = indata)) +
  geom_point(alpha = 0.4) +
  facet_wrap(~factor(run)) +
  scale_colour_manual(values = c("red", "black"))

ggplot(subset(full_gebvs, indata == "xval"),
       aes(gebv, mean_co_count)) +
  geom_point(alpha = 0.4) +
  stat_smooth()

cor.test(subset(full_gebvs, indata == "xval")$gebv, subset(full_gebvs, indata == "xval")$mean_co_count)


