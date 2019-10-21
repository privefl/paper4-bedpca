library(bigsnpr)
library(bigreadr)
rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
ind.row <- which(!fam$V2 %in% rel$ID2)
obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")

system.time(
  obj.svd <- bed_autoSVD2(obj.bed, ind.row = ind.row, k = 20, ncores = nb_cores())
) # 5h
# Phase of clumping (on MAC) at r^2 > 0.2.. keep 261307 variants.
# Discarding 0 variant with MAC < 10.
#
# Iteration 1:
# Computing SVD..
# 3105 outlier variants detected..
# 8 long-range LD regions detected..
#
# Iteration 2:
# Computing SVD..
# 4550 outlier variants detected..
# 11 long-range LD regions detected..
#
# Iteration 3:
# Computing SVD..
# 3882 outlier variants detected..
# 11 long-range LD regions detected..
#
# Iteration 4:
# Computing SVD..
# 3832 outlier variants detected..
# 15 long-range LD regions detected..
#
# Iteration 5:
# Computing SVD..
# 5500 outlier variants detected..
# 17 long-range LD regions detected..
#
# Iteration 6:
# Computing SVD..
# Maximum number of iterations reached.


system.time(
  ind.keep2 <- bed_clumping(obj.bed, ind.row = ind.row, ncores = nb_cores())
) # 22 min
length(ind.keep2)  # 261307

system.time(
  obj.svd2 <- bed_randomSVD(obj.bed, ind.row = ind.row, ind.col = ind.keep2,
                            k = 20, ncores = nb_cores())
) # 34 min

library(ggplot2)
plot(obj.svd)
# ggsave("figures/UKBB-screeplot.pdf", width = 12, height = 7)

plot(obj.svd, type = "loadings", loadings = 1:20, coeff = 0.4)
# ggsave("figures/UKBB-loadings.pdf", width = 13, height = 7)

plot(obj.svd, type = "scores", scores = 1:20, coeff = 0.4, ncol = 5)
# ggsave("figures/UKBB-PC1-20.png", width = 9, height = 7)


PCs <- predict(obj.svd)
system.time(
  knn <- bigutilsr::knn_parallel(PCs, k = 31, ncores = nb_cores())
) # 6 min
