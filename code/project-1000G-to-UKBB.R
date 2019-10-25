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
length(attr(test$obj.svd.ref, "subset")) # 171977

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

library(bigreadr)
## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("data/coding1001.tsv")
csv <- "data/ukb22544.csv"
library(dplyr)
df0 <- fread2(csv, select = c("eid", "21000-0.0", "22006-0.0"),
              col.names = c("eid", "pop", "is_caucasian")) %>%
  mutate(
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning),
    is_caucasian = as.logical(is_caucasian)
  )

eid <- bed.ref$fam$sample.ID[ind.row[PC.ref[, 6] > -50 & PC.ref[, 5] < -10]]
pop_UKBB <- df0$pop[match(eid, df0$eid)]
table(pop_UKBB)


fam2 <- bigreadr::fread2("data/1000G_phase3_common_norel.fam2")
fam2$Population[proj2[, 6] > -50 & proj2[, 5] < -10]
