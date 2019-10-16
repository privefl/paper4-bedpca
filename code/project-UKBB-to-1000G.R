library(bigsnpr)
library(ggplot2)

obj.bed <- bed("data/ukbb_488282.bed")
bed.ref <- bed(download_1000G("data"))

system.time(
  test <- bed_projectPCA(bed.ref, obj.bed, k = 20, ncores = nb_cores())
) # 20 min (including 12 min for projection)

PC.ref <- predict(test$obj.svd.ref)
proj1 <- test$simple_proj
proj2 <- test$OADP_proj

# shrinkage coefficients
shrinkage <- unname(sapply(1:20, function(k) {
  MASS::rlm(proj2[, k] ~ proj1[, k] + 0)$coef
}))
round(shrinkage, 2)
#  [1] 1.01 1.02 1.06 1.08 1.35 1.82 2.32 2.35 2.78 2.84 2.98 3.47 4.37
# [14] 4.64 4.94 5.30 5.74 6.49 6.70 6.74

ind <- sample(nrow(proj1), 20e3)

plot_grid(plotlist = lapply(1:8, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj1[ind, k1], proj1[ind, k2]), color = "red") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj2[ind, k1], proj2[ind, k2]), color = "blue") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    scale = 0.95
  )
}), nrow = 4)

# ggsave("figures/proj1000G-UKBB.png", width = 11, height = 7)


