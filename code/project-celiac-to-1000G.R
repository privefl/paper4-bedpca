library(bigsnpr)
library(ggplot2)

obj.bed <- bed("~/Bureau/Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC.bed")
bed.ref <- bed(download_1000G("tmp-data"))

system.time(
  test <- bed_projectPCA(bed.ref, obj.bed, k = 20, ncores = nb_cores(),
                         build.new = "hg18", liftOver = "tmp-data/liftOver")
)

PC.ref <- predict(test$obj.svd.ref)
proj1 <- test$simple_proj
proj2 <- test$OADP_proj

# shrinkage coefficients
shrinkage <- unname(sapply(1:20, function(k) {
  MASS::rlm(proj2[, k] ~ proj1[, k] + 0)$coef
}))
round(shrinkage, 2)
#  [1] 1.01 1.01 1.04 1.05 1.30 1.33 1.51 1.72 1.93 2.08 2.30 2.31 2.36
# [14] 2.39 2.54 2.87 2.97 3.25 3.65 3.91

plot_grid(plotlist = lapply(1:8, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj1[, k1], proj1[, k2]), color = "red") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj2[, k1], proj2[, k2]), color = "blue") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    scale = 0.95
  )
}), nrow = 4)

# ggsave("figures/proj1000G-Celiac.png", width = 11, height = 7)

#### Ancestry estimation ####

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")

plot_grid(plotlist = lapply(5:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = fam2$Population) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}))

unAsIs <- function(X) {
  class(X) <- setdiff(class(X), "AsIs")
  X
}

library(dplyr)
seq_PC <- c(1:19)
all_covRob <- cbind(fam2, PC = I(PC.ref)) %>%
  group_by(`Super Population`, `Population Description`, Population) %>%
  summarize(maha = list(bigutilsr::covRob(unAsIs(PC[, seq_PC]), estim = "pairwiseGK")))
POP <- paste(all_covRob$`Super Population`, all_covRob$Population, sep = "_")

all_dist <- sapply(all_covRob$maha, function(maha) {
  mahalanobis(proj2[, seq_PC], center = maha$center, cov = maha$cov)
})
colnames(all_dist) <- POP

all_prob <- pchisq(all_dist, df = length(seq_PC), lower.tail = FALSE)
choose_pop <- apply(all_prob, 1, which.max)
pval_max <- all_prob[cbind(rows_along(all_prob), choose_pop)]
choose_pop2 <- ifelse(pval_max < 0.05, NA, POP[choose_pop])

ggplot() +
  geom_histogram(aes(pval_max), breaks = seq(0, 1, by = 0.05),
                 color = "#FFFFFF", fill = "#000000", alpha = 0.5) +
  geom_vline(xintercept = 0.05, color = "red") +
  theme_bigstatsr() +
  labs(x = "Maximum p-value testing population belonging")

# ggsave("figures/hist-pval-max.pdf", width = 8.5, height = 5.5)

# Get population from external files
pop.files <- list.files(path = "~/Bureau/Dubois2010_data/",
                        pattern = "cluster_*", full.names = TRUE)
pop <- snp_getSampleInfos(obj.bed, pop.files)[[1]]
pop_celiac <- c("Netherlands", "Italy", "UK", "UK", "Finland")[pop]

mean(pval_max < 0.05) # 25.8%
tapply(pval_max, pop_celiac, function(x) mean(x < 0.05))
#   Finland       Italy Netherlands          UK
# 0.3469305   0.7028846   0.2359005   0.1934177


# https://www.internationalgenome.org/category/population/
table(pop_celiac, substr(choose_pop2, 1, 3))
# pop_celiac     EUR
# Finland       1617
# Italy          309
# Netherlands   1260
# UK            8161

table(choose_pop2, pop_celiac) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Self-reported ancestry (top) of UKBB individuals and their matching to 1000G populations (left) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred",
                 align = "|l|c|c|c|c|") %>%
  print(caption.placement = "top")

#         Finland Italy Netherlands   UK
# EUR_CEU       5     0         930 4753
# EUR_FIN    1612     0           0    2
# EUR_GBR       0     0         329 3390
# EUR_IBS       0    37           1   10
# EUR_TSI       0   272           0    6
