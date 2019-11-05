library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
is_rel <- (fam$V2 %in% rel$ID1 | fam$V2 %in% rel$ID2)
obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")

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

ind_used <- c(sample(which(pop_UKBB == "British"), 10e3),
              sample(which(pop_UKBB == "Irish"), 5e3),
              which(!pop_UKBB %in% c("British", "Irish")))
ind.train <- intersect(ind_used, which(!is_rel))
table(pop_UKBB[ind.train])

system.time(
  obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.train,
                         k = 50, ncores = nb_cores())
) # 2.9H

plot(obj.svd)
# ggsave("figures/UKBB-screeplot-restricted.pdf", width = 7, height = 5)

plot_grid(plotlist = lapply(14:25, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.4) +
    aes(color = pop_UKBB[ind.train]) +
    theme(legend.position = "none")
}))
# ggsave("figures/UKBB-scores-restricted.png", width = 11, height = 6)


system.time(
  proj <- bed_projectSelfPCA(obj.svd, obj.bed,
                             ind.row = rows_along(obj.bed)[-ind.train],
                             ncores = nb_cores())
) # 21 min
proj1 <- proj$simple_proj
proj2 <- proj$OADP_proj

# shrinkage coefficients
shrinkage <- unname(sapply(1:50, function(k) {
  MASS::rlm(proj2[, k] ~ proj1[, k] + 0)$coef
}))
round(shrinkage, 2)
#  [1] 1.00 1.00 1.01 1.01 1.04 1.05 1.07 1.07 1.11 1.12 1.13 1.20 1.23 1.25 1.30 1.31 1.35
# [18] 1.37 1.40 1.43 1.44 1.50 1.52 1.54 1.56 1.56 1.59 1.65 1.66 1.70 1.73 1.74 1.74 1.77
# [35] 1.77 1.78 1.78 1.79 1.79 1.79 1.80 1.80 1.80 1.80 1.80 1.80 1.80 1.80 1.80 1.80


PC.ref <- predict(obj.svd)
ind <- rows_along(PC.ref) # sample(nrow(PC.ref), 10e3)
ind2 <- sample(nrow(proj1), 20e3)
coeff <- 0.5
plot_grid(plotlist = lapply(14:25, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(

    qplot(PC.ref[ind, k1], PC.ref[ind, k2]) +
      geom_point(aes(proj1[ind2, k1], proj1[ind2, k2]), color = "red") +
      theme_bigstatsr(coeff) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
      coord_equal(),

    qplot(PC.ref[ind, k1], PC.ref[ind, k2]) +
      geom_point(aes(proj2[ind2, k1], proj2[ind2, k2]), color = "blue") +
      theme_bigstatsr(coeff) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
      coord_equal(),

    scale = 0.95, ncol = 2
  )
}), scale = 0.95, ncol = 3)
# ggsave("figures/projUKBB-related.png", width = 13, height = 9)
