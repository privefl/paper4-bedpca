library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
ind.row <- which(!fam$V2 %in% rel$ID2)
obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")

obj.svd <- readRDS("tmp-results/SVD_UKBB.rds")
PC <- predict(obj.svd)

dist <- bigutilsr::covRob(PC, estim = "pairwiseGK")$dist

ggplot() +
  geom_histogram(aes(log(dist)), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  theme_bigstatsr() +
  labs(x = "log Mahalanobis distance", y = "Frequency")

# ggsave("figures/hist-Maha-dist.pdf", width = 9, height = 6)

plot_grid(plotlist = lapply(1:10, function(k) {
  k1 <- 2 * k - 1; k2 <- 2 * k
  qplot(PC[, k1], PC[, k2]) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme_bigstatsr(0.5) +
    aes(color = dist) +
    theme(legend.position = "none") +
    scale_color_viridis_c(trans = "log") +
    coord_equal()
}), ncol = 5)

# ggsave("figures/UKBB-Maha-dist.png", width = 11, height = 8)

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

pval <- pchisq(dist, df = ncol(PC), lower.tail = FALSE)

ggplot() +
  geom_histogram(aes(pval), breaks = 0:20 / 20,
                 color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_y_sqrt(breaks = 10^(1:6)) +
  theme_bigstatsr() +
  labs(x = "P-value", y = "Frequency (sqrt-scale)")

# ggsave("figures/hist-Maha-pval.pdf", width = 9, height = 6)


sapply(setNames(nm = c(0, 1e-20, 1e-8, 1e-5, 1e-4, 1e-3, 0.01, 0.02, 0.05, 0.1, 0.2)), function(thr) {
  tab <- table(pop_UKBB[ind.row][pval >= thr])
  c(tab, All = sum(tab))
})
#                                 0  1e-20  1e-08  1e-05  1e-04  0.001   0.01   0.02   0.05    0.1    0.2
# Prefer not to answer         1347    867    844    830    823    808    777    762    719    673    581
# Do not know                   178     65     61     58     57     53     51     50     47     44     39
# White                         456    372    356    343    338    332    311    302    289    275    241
# Mixed                          38      4      4      4      4      4      2      2      2      2      2
# Asian or Asian British         39      0      0      0      0      0      0      0      0      0      0
# Black or Black British         20      1      1      1      1      1      1      1      1      1      1
# Chinese                      1436      1      1      1      1      1      1      1      1      1      1
# Other ethnic group           4102    220    203    199    194    172    133    122    109     95     76
# British                    354557 340933 334013 329188 326376 321659 310963 304622 290116 270424 236302
# Irish                       10628  10488  10435  10243   9809   8604   6230   5353   4088   3128   2230
# Any other white background  14954   6390   5079   4535   4274   3980   3493   3324   2973   2646   2177
# White and Black Caribbean     543      4      4      4      4      4      4      4      4      2      2
# White and Black African       364      3      3      3      3      3      3      3      2      2      1
# White and Asian               716     12     10      6      6      6      5      5      5      5      5
# Any other mixed background    891     73     66     58     57     56     49     48     44     41     33
# Indian                       5242      1      1      1      1      1      1      1      1      0      0
# Pakistani                    1606      0      0      0      0      0      0      0      0      0      0
# Bangladeshi                   213      0      0      0      0      0      0      0      0      0      0
# Any other Asian background   1658      0      0      0      0      0      0      0      0      0      0
# Caribbean                    3793      0      0      0      0      0      0      0      0      0      0
# African                      3103      0      0      0      0      0      0      0      0      0      0
# Any other Black background    106      0      0      0      0      0      0      0      0      0      0
# All                        405990 359434 351081 345474 341948 335684 322024 314600 298401 277339 241691

plot_grid(plotlist = lapply(1:10, function(k) {
  k1 <- 2 * k - 1; k2 <- 2 * k
  qplot(PC[, k1], PC[, k2]) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme_bigstatsr(0.7) +
    aes(color = log(dist) < 5) +
    theme(legend.position = "none") +
    scale_color_viridis_d(direction = -1) +
    coord_equal()
}), ncol = 5)

# ggsave("figures/UKBB-Maha-outlier.png", width = 11, height = 8)


# Reported by the UK Biobank

is_WB <- df0$is_caucasian[match(obj.bed$fam$sample.ID, df0$eid)]
length(ind1 <- na.omit(ind.row[is_WB[ind.row]]))
table(pop_UKBB[ind1])  ## All British

plot_grid(plotlist = lapply(1:10, function(k) {
  k1 <- 2 * k - 1; k2 <- 2 * k
  qplot(PC[, k1], PC[, k2]) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2), color = "Is White British?") +
    theme_bigstatsr(0.5) +
    aes(color = ind.row %in% ind1) +
    theme(legend.position = "none") +
    scale_color_viridis_d(direction = -1) +
    coord_equal()
}), ncol = 5)

# ggsave("figures/UKBB-White-British.png", width = 11, height = 8)
