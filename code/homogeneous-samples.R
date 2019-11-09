library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("data/coding1001.tsv")
df0 <- fread2(
  "data/ukb22544.csv",
  select = c("eid", "21000-0.0", "22006-0.0", paste0("22009-0.", 1:20)),
  col.names = c("eid", "pop", "is_caucasian", paste0("PC", 1:20))
) %>%
  mutate(
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning),
    is_caucasian = as.logical(is_caucasian)
  )

PC <- na.omit(as.matrix(df0[-(1:3)]))
ind_na <- attr(PC, "na.action")
pop_UKBB <- df0$pop[-ind_na]


dist <- bigutilsr::covRob(PC, estim = "pairwiseGK")$dist

ggplot() +
  geom_histogram(aes(log(dist)), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  theme_bigstatsr() +
  labs(x = "log Mahalanobis distance", y = "Frequency")

# ggsave("figures/hist-Maha-dist.pdf", width = 10, height = 6)


plot_grid(plotlist = lapply(1:9, function(k) {
  k1 <- 2 * k - 1; k2 <- 2 * k
  qplot(PC[, k1], PC[, k2]) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme_bigstatsr(0.5) +
    aes(color = dist) +
    theme(legend.position = "none") +
    scale_color_viridis_c(trans = "log") +
    coord_equal()
}), ncol = 3)

# ggsave("figures/UKBB-Maha-dist.png", width = 10, height = 9)


sapply(setNames(3:12, paste("<", 3:12)), function(thr) {
  tab <- table(pop_UKBB[dist < exp(thr)])
  c(tab, All = sum(tab))
}) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Number of UKBB individuals with (log) Mahalanobis distance lower than some threshold (top), grouped by self-reported ancestry (left).",
                 label = "tab:homogeneous",
                 align = "|l|r|r|r|r|r|r|r|r|r|r|") %>%
  print(caption.placement = "top")
#                               < 3    < 4    < 5    < 6    < 7    < 8    < 9   < 10   < 11   < 12
# Prefer not to answer          484   1013   1062   1099   1139   1177   1279   1405   1471   1583
# Do not know                    36     68     76     84     92    118    155    188    196    204
# White                         186    422    457    483    513    533    543    543    545    546
# Mixed                           2      6      6      7      8     15     26     42     46     46
# Asian or Asian British          0      0      0      0      0      3     20     40     42     42
# Black or Black British          1      2      2      2      2      2      2      4      6     26
# Chinese                         0      1      1      1      1      2      5     21   1423   1504
# Other ethnic group             57    230    261    314    469    885   1939   2761   3681   4356
# British                    191713 400516 416492 424490 427769 429172 431026 431082 431089 431090
# Irish                        1416  12039  12620  12700  12734  12743  12759  12759  12759  12759
# Any other white background   1468   4747   6953   9341  12979  14613  15741  15810  15820  15820
# White and Black Caribbean       1      4      4      4      9     35    142    537    589    597
# White and Black African         1      3      3      4      6     29     99    333    400    402
# White and Asian                 4      7     13     23     79    350    651    790    802    802
# Any other mixed background     24     66     87    155    274    391    595    884    990    996
# Indian                          0      2      2      5      6     29   1682   5700   5716   5716
# Pakistani                       0      0      0      0      1     13    532   1747   1748   1748
# Bangladeshi                     0      0      0      0      0      0      2    220    221    221
# Any other Asian background      0      0      0      1      6     66    427   1364   1730   1747
# Caribbean                       0      0      0      0      0      0      3    113   1323   4299
# African                         0      1      1      1      1      1      3     58    350   3205
# Any other Black background      0      0      0      0      0      1      3     22     49    118
# All                        195393 419127 438040 448714 456088 460178 467634 476423 480996 487827

plot_grid(plotlist = lapply(1:9, function(k) {
  k1 <- 2 * k - 1; k2 <- 2 * k
  qplot(PC[, k1], PC[, k2]) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme_bigstatsr(0.7) +
    aes(color = log(dist) < 5) +
    theme(legend.position = "none") +
    scale_color_viridis_d(direction = -1) +
    coord_equal()
}), ncol = 3)

# ggsave("figures/UKBB-Maha-outlier.png", width = 10, height = 9)


# Reported by the UK Biobank
is_WB <- df0$is_caucasian[-ind_na]; is_WB[is.na(is_WB)] <- FALSE
table(pop_UKBB[is_WB])  ## All British

plot_grid(plotlist = lapply(1:9, function(k) {
  k1 <- 2 * k - 1; k2 <- 2 * k
  qplot(PC[, k1], PC[, k2]) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme_bigstatsr(0.5) +
    aes(color = is_WB) +
    theme(legend.position = "none") +
    scale_color_viridis_d(direction = -1) +
    coord_equal()
}), ncol = 3)

# ggsave("figures/UKBB-White-British.png", width = 10, height = 9)
