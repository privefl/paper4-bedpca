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

#### Ancestry estimation ####

fam2 <- bigreadr::fread2("data/1000G_phase3_common_norel.fam2")

plot_grid(plotlist = lapply(5:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = fam2$Population) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}))

unAsIs <- function(X) {
  class(X) <- setdiff(class(X), "AsIs")
  X
}

library(dplyr)
seq_PC <- c(1:14, 16)
all_covRob <- cbind(fam2, PC = I(PC.ref)) %>%
  group_by(`Super Population`, `Population Description`, Population) %>%
  summarize(maha = list(bigutilsr::covRob(unAsIs(PC[, seq_PC]), estim = "pairwiseGK")))
POP <- paste(all_covRob$`Super Population`, all_covRob$Population, sep = "_")

all_dist <- sapply(all_covRob$maha, function(maha) {
  mahalanobis(proj2[, seq_PC], center = maha$center, cov = maha$cov)
})
colnames(all_dist) <- POP

round(all_prob <- pchisq(all_dist, df = length(seq_PC), lower.tail = FALSE), 3)
round(corr <- 100 * cor(all_prob), 1)
ind <- which(corr > 10 & lower.tri(corr), arr.ind = TRUE)
data.frame(pop1 = POP[ind[, 1]], pop2 = POP[ind[, 2]], cor = corr[ind])
# 1  ACB  ASW 26.57052
# 2  YRI  ESN 55.01614
# 3  KHV  CDX 39.74580
# 4  CHS  CHB 65.06130
# 5  CEU  GBR 65.59252
# 6  STU  ITU 70.99639

choose_pop <- apply(all_prob, 1, which.max)
pval_max <- all_prob[cbind(rows_along(all_prob), choose_pop)]
choose_pop2 <- ifelse(pval_max < 0.05, NA, POP[choose_pop])

ggplot() +
  geom_histogram(aes(pval_max), breaks = seq(0, 1, by = 0.05),
                 color = "#FFFFFF", fill = "#000000", alpha = 0.5) +
  geom_vline(xintercept = 0.05, color = "red") +
  theme_bigstatsr() +
  labs(x = "Maximum p-value testing population belonging")

# ggsave("figures/hist-pval-max.pdf", width = 8.5, height = 5.5)

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("data/coding1001.tsv")
unknown <- c("Prefer not to answer", "Do not know", "Mixed", "Other ethnic group",
             "Asian or Asian British", "Black or Black British",
             "White and Black Caribbean", "White and Black African",
             "White and Asian", "Any other mixed background")
csv <- "data/ukb22544.csv"
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

mean(pval_max < 0.05) # 15.5%
tapply(pval_max, pop_UKBB, function(x) mean(x < 0.05))
#                      White                    Chinese                    British                      Irish
#                  0.2036697                  0.4853723                  0.1046357                  0.1980400
# Any other white background                     Indian                  Pakistani                Bangladeshi
#                  0.6078655                  0.6471309                  0.6327231                  0.8144796
# Any other Asian background                  Caribbean                    African Any other Black background
#                  0.6800229                  0.5205958                  0.6357678                  0.6355932


# https://www.internationalgenome.org/category/population/
table(pop_UKBB, substr(choose_pop2, 1, 3))
#                               AFR    AMR    EAS    EUR    SAS
# White                           0      1      0    433      0
# Chinese                         0      0    773      1      0
# British                         0      1      1 385925      1
# Irish                           0      0      0  10229      0
# Any other white background      0     51      0   6151      0
# Indian                          0      0      0      0   2017
# Pakistani                       0      0      0      0    642
# Bangladeshi                     0      0      0      0     41
# Any other Asian background      0      0     55      0    504
# Caribbean                    2056      0      0      0      4
# African                      1166      0      0      1      0
# Any other Black background     43      0      0      0      0

table(choose_pop2, pop_UKBB) %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Self-reported ancestry (top) of UKBB individuals and their matching to 1000G populations (left) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred",
                 align = "|l|c|c|c|c|c|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top", hline.after = c(-1, 0, 7, 11, 16, 21, 26))
