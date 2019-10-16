library(bigsnpr)
library(ggplot2)

bedfile.ref <- download_1000G("tmp-data")
obj.bed <- bed(bedfile.ref)

set.seed(1); train <- sample(nrow(obj.bed), 0.6 * nrow(obj.bed))

# obj.svd <- bed_autoSVD2(obj.bed, ind.row = train, k = 20, ncores = nb_cores())
# saveRDS(obj.svd, "tmp-data/SVD_1000G.rds")
obj.svd <- readRDS("tmp-data/SVD_1000G.rds")

test <- bed_projectSelfPCA(obj.svd, obj.bed,
                           ind.row = rows_along(obj.bed)[-train],
                           ncores = nb_cores())

PC.ref <- predict(test$obj.svd.ref)
proj1 <- test$simple_proj
proj2 <- test$OADP_proj

# shrinkage coefficients
shrinkage <- unname(sapply(1:20, function(k) {
  MASS::rlm(proj2[, k] ~ proj1[, k] + 0)$coef
}))
round(shrinkage, 2)
#  [1] 1.01 1.02 1.07 1.09 1.51 1.68 1.94 1.40 2.88 3.17 2.90 2.92 3.23
# [14] 5.13 5.25 5.04 4.58 5.69 6.26 6.32

plot_grid(plotlist = lapply(1:4, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj1[, k1], proj1[, k2]), color = "red") +
      theme_bigstatsr() +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj2[, k1], proj2[, k2]), color = "blue") +
      theme_bigstatsr() +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    scale = 0.95
  )
}))

# ggsave("figures/proj1000G-PC1-8.pdf", width = 15, height = 7.5)

plot_grid(plotlist = lapply(5:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj1[, k1], proj1[, k2]), color = "red") +
      theme_bigstatsr() +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj2[, k1], proj2[, k2]), color = "blue") +
      theme_bigstatsr() +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    scale = 0.95
  )
}), nrow = 3)

# ggsave("figures/proj1000G-PC9-20.pdf", width = 13.5, height = 8)


ind.keep <- attr(obj.svd, "subset")
K <- bed_tcrossprodSelf(obj.bed, ind.row = train, ind.col = ind.keep)
eigval <- eigen(K[], symmetric = TRUE, only.values = TRUE)$values
stopifnot(all.equal(sqrt(eigval[1:20]), obj.svd$d))
proj3 <- hdpca::pc_adjust(eigval, p = length(ind.keep), n = length(train),
                          test.scores = proj1, n.spikes.max = 50)
round(shrinkage2 <- (proj3 / proj1)[1, ], 2)
#  [1] 1.01 1.02 1.07 1.10 1.43 1.54 1.73 1.82 2.51 2.78 2.88 2.91 3.11
# [14] 1.00 1.00 1.00 1.00 1.00 1.00 1.00

plot_grid(plotlist = lapply(3:6, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj3[, k1], proj3[, k2]), color = "limegreen") +
      theme_bigstatsr() +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj2[, k1], proj2[, k2]), color = "blue") +
      theme_bigstatsr() +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    scale = 0.95
  )
}))

# ggsave("figures/proj1000G-PC5-12.pdf", width = 15, height = 7.5)

test2 <- bed_projectSelfPCA(obj.svd, obj.bed,
                            ind.row = train,
                            ncores = nb_cores())
proj4 <- test2$OADP_proj

plot_grid(plotlist = lapply(3:6, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(proj4[, k1], proj4[, k2], color = I("blue")) +
    geom_point(aes(PC.ref[, k1], PC.ref[, k2])) +
    theme_bigstatsr() +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}), scale = 0.95)

# ggsave("figures/proj1000G-related.pdf", width = 8.5, height = 6.5)
