library(bigsnpr)
library(ggplot2)

obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")
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
#  [1] 1.01 1.02 1.06 1.08 1.36 1.82 2.33 2.36 2.78 2.84 2.99 3.51 4.38
# [14] 4.67 4.99 5.31 5.74 6.55 6.71 6.75

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

#### Ancestry estimation ####

fam2 <- bigreadr::fread2("data/1000G_phase3_common_norel.fam2")

plot(test$obj.svd.ref)
plot_grid(plotlist = lapply(5:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = fam2$Population) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    theme(legend.position = "none")
}))

# ggsave("figures/PC-1000G-UKBB.pdf", width = 10, height = 6)

seq_PC <- c(1:14, 16)
pop <- paste(fam2$`Super Population`, fam2$Population, sep = "_")
pop_PCs <- vctrs::vec_split(PC.ref, pop)

all_pval <- sapply(pop_PCs$val, function(PC) {
  maha <- bigutilsr::covRob(PC[, seq_PC], estim = "pairwiseGK")
  dist <- mahalanobis(proj2[, seq_PC], center = maha$center, cov = maha$cov)
  pval <- pchisq(dist, df = length(seq_PC), lower.tail = FALSE)
})

choose_pop <- apply(all_pval, 1, which.max)
pval_max <- all_pval[cbind(rows_along(all_pval), choose_pop)]
choose_pop2 <- ifelse(pval_max < 0.05, NA, pop_PCs$key[choose_pop])

ggplot() +
  geom_histogram(aes(pval_max), breaks = seq(0, 1, by = 0.05),
                 color = "#FFFFFF", fill = "#000000", alpha = 0.5) +
  geom_vline(xintercept = 0.05, color = "red") +
  theme_bigstatsr() +
  labs(x = "Maximum p-value testing population belonging")

# ggsave("figures/hist-pval-max.pdf", width = 8.5, height = 5.5)

library(bigreadr)
## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("data/coding1001.tsv")
unknown <- c("Prefer not to answer", "Do not know", "Mixed", "Other ethnic group",
             "Asian or Asian British", "Black or Black British",
             "White and Black Caribbean", "White and Black African",
             "White and Asian", "Any other mixed background")
csv <- "data/ukb22544.csv"
library(dplyr)
df0 <- fread2(csv, select = c("eid", "21000-0.0", "22006-0.0"),
              col.names = c("eid", "pop", "is_caucasian")) %>%
  mutate(
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning) %>%
      forcats::fct_recode("Other White" = "Any other white background",
                          "Other Asian" = "Any other Asian background",
                          "Other Black" = "Any other Black background") %>%
      forcats::fct_other(drop = unknown, other_level = NA) %>%
      droplevels(pop, exclude = NA) %>%
      forcats::fct_relevel(c(
        "British", "Irish", "White", "Other White",
        "Indian", "Pakistani", "Bangladeshi", "Chinese", "Other Asian",
        "Caribbean", "African", "Other Black")),
    is_caucasian = as.logical(is_caucasian)
  )

pop_UKBB <- df0$pop[match(obj.bed$fam$sample.ID, df0$eid)]

mean(pval_max < 0.05) # 15.9%
tapply(pval_max, pop_UKBB, function(x) mean(x < 0.05))
#     British       Irish       White Other White      Indian   Pakistani
#   0.1107768   0.2033712   0.2018349   0.6136191   0.5720784   0.5606407
# Bangladeshi     Chinese Other Asian   Caribbean     African Other Black
#   0.7420814   0.4448138   0.6485404   0.4980219   0.5948814   0.6186441


# https://www.internationalgenome.org/category/population/
table(pop_UKBB, substr(choose_pop2, 1, 3))
#                AFR    AMR    EAS    EUR    SAS
# British          0      1      1 383277      2
# Irish            0      0      0  10161      0
# White            0      1      0    434      0
# Other White      0     55      0   6056      0
# Indian           0      0      0      0   2446
# Pakistani        0      0      0      0    768
# Bangladeshi      0      0      0      0     57
# Chinese          0      0    834      1      0
# Other Asian      0      0     60      0    554
# Caribbean     2152      0      0      0      5
# African       1297      0      0      1      0
# Other Black     45      0      0      0      0

table(choose_pop2, pop_UKBB, exclude = NULL) %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Self-reported ancestry (top) of UKBB individuals and their matching to 1000G populations (left) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred",
                 align = "|l|c|c|c|c|c|c|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 7, 10, 15, 20, 25, 26))
