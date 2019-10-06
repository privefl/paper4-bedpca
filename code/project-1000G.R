library(bigsnpr)
bedfile.ref <- download_1000G("tmp-data")
bed.ref <- bed(bedfile.ref)

set.seed(1)
train <- sample(nrow(bed.ref), 0.6 * nrow(bed.ref))
test <- bed_projectPCA(
  bed.new = bed.ref,
  bed.ref = bed.ref,
  ind.row.new = rows_along(bed.ref)[-train],
  ind.row.ref = train,
  k = 20,
  kNN = 6,
  ncores = nb_cores()
)

PC.ref <- predict(test$obj.svd.ref)
proj1 <- test$simple_proj
proj2 <- test$OADP_proj

library(ggplot2)
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
