library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
ind.row <- which(!fam$V2 %in% rel$ID2)
obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")

system.time(
  obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.row, k = 20, ncores = nb_cores())
) # 5h
# Phase of clumping (on MAC) at r^2 > 0.2.. keep 261306 variants.
# Discarding 0 variant with MAC < 10.
#
# Iteration 1:
# Computing SVD..
# 3105 outlier variants detected..
# 8 long-range LD regions detected..
#
# Iteration 2:
# Computing SVD..
# 4549 outlier variants detected..
# 11 long-range LD regions detected..
#
# Iteration 3:
# Computing SVD..
# 3881 outlier variants detected..
# 11 long-range LD regions detected..
#
# Iteration 4:
# Computing SVD..
# 3839 outlier variants detected..
# 15 long-range LD regions detected..
#
# Iteration 5:
# Computing SVD..
# 5488 outlier variants detected..
# 17 long-range LD regions detected..
#
# Iteration 6:
# Computing SVD..
# Maximum number of iterations reached.

# saveRDS(obj.svd, "tmp-results/SVD_UKBB.rds")

system.time(
  ind.keep2 <- bed_clumping(obj.bed, ind.row = ind.row, ncores = nb_cores())
) # 22 min
length(ind.keep2)  # 261,306

system.time(
  obj.svd2 <- bed_randomSVD(obj.bed, ind.row = ind.row, ind.col = ind.keep2,
                            k = 20, ncores = nb_cores())
) # 34 min

library(ggplot2)
plot(obj.svd)
# ggsave("figures/UKBB-screeplot.pdf", width = 12, height = 7)

plot(obj.svd, type = "loadings", loadings = 1:20, coeff = 0.4)
# ggsave("figures/UKBB-loadings.pdf", width = 13, height = 7)

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("data/coding1001.tsv")
csv <- "data/ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "21000-0.0", "22006-0.0"),
              col.names = c("eid", "pop", "is_caucasian")) %>%
  mutate(
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning),
    is_caucasian = as.logical(is_caucasian)
  )
pop_UKBB <- df0$pop[match(obj.bed$fam$sample.ID, df0$eid)]
table(pop_UKBB)

plot_grid(plotlist = lapply(1:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.5) +
    aes(color = pop_UKBB[ind.row]) +
    theme(legend.position = "none")
}), ncol = 5)
# ggsave("figures/UKBB-PC1-20.png", width = 11, height = 8)


# Benchmark KNN
PCs <- predict(obj.svd)
system.time(
  knn <- bigutilsr::knn_parallel(PCs, k = 31, ncores = nb_cores())
) # 6 min


#### Reported by the UK Biobank ####

snp_qc <- bigreadr::fread2("https://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.txt")
str(loadings <- snp_qc[grep("_loading$", names(snp_qc))])

plot_grid(plotlist = lapply(1:40, function(i) {
  l_i <- loadings[[i]]
  ind <- which(!is.na(l_i))
  ggplot() +
    geom_hex(aes(ind, l_i[ind])) +
    theme_bigstatsr(0.4) +
    labs(title = paste0("Loadings of PC", i), x = "Column index", y = NULL) +
    scale_fill_viridis_c()
}), align = "hv", ncol = 5, scale = 0.95)

# ggsave("figures/UKBB-loadings1-40.pdf", width = 15, height = 15)

nPC <- 20
PC <- fread2(
  "data/ukb22544.csv",
  select = paste0("22009-0.", 1:nPC),
  col.names = paste0("PC", 1:nPC)
)

plot_grid(plotlist = lapply(1:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC[ind, k1], PC[ind, k2]) +
    theme_bigstatsr(0.4) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    coord_equal()
}), scale = 0.95)

# ggsave("figures/UKBB-scores1-20.pdf", width = 15, height = 15)
