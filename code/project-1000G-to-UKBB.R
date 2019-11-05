library(bigsnpr)
library(bigreadr)
rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
ind.row <- which(!fam$V2 %in% rel$ID2)
length(ind.row)  # 406545
bed.ref <- bed("data/ukbb_bed/ukbb_488282.bed")
obj.bed <- bed(download_1000G("data"))

system.time(
  test <- bed_projectPCA(bed.ref, obj.bed, ind.row.ref = ind.row,
                         k = 20, ncores = nb_cores())
) # 4.4 H
length(attr(test$obj.svd.ref, "subset")) # 171,987

PC.ref <- predict(test$obj.svd.ref)
proj1 <- test$simple_proj
proj2 <- test$OADP_proj

# shrinkage coefficients
shrinkage <- unname(sapply(1:20, function(k) {
  MASS::rlm(proj2[, k] ~ proj1[, k] + 0)$coef
}))
round(shrinkage, 2)
#  [1] 1.00 1.00 1.00 1.01 1.01 1.02 1.03 1.03 1.04 1.04 1.04 1.05 1.05
# [14] 1.06 1.07 1.07 1.08 1.08 1.08 1.08


library(ggplot2)
ind <- sample(nrow(PC.ref), 20e3)

plot_grid(plotlist = lapply(1:8, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(
    qplot(PC.ref[ind, k1], PC.ref[ind, k2], size = I(2)) +
      geom_point(aes(proj1[, k1], proj1[, k2]), color = "red") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    qplot(PC.ref[ind, k1], PC.ref[ind, k2], size = I(2)) +
      geom_point(aes(proj2[, k1], proj2[, k2]), color = "blue") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    scale = 0.95
  )
}), nrow = 4)

# ggsave("figures/proj1000G-UKBB2.png", width = 10, height = 7)
