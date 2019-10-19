library(bigsnpr)
library(ggplot2)

snp_plinkQC(
  plink.path = download_plink("tmp-data"),
  prefix.in = "~/Bureau/Dubois2010_data/FinnuncorrNLITUK3hap550",
  geno = 0.01, mind = 0.01,
  autosome.only = TRUE
)

obj.bed <- bed("~/Bureau/Dubois2010_data/FinnuncorrNLITUK3hap550_QC.bed")
bed.ref <- bed(download_1000G("tmp-data"))

system.time(
  test <- bed_projectPCA(bed.ref, obj.bed, k = 20, ncores = nb_cores(),
                         build.new = "hg18", liftOver = "tmp-data/liftOver")
)

PC.ref <- predict(test$obj.svd.ref)
proj1 <- test$simple_proj
proj2 <- test$OADP_proj

plot_grid(plotlist = lapply(1:8, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  plot_grid(
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj1[, k1], proj1[, k2]), color = "red") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
      geom_point(aes(proj2[, k1], proj2[, k2]), color = "blue") +
      theme_bigstatsr(0.5) +
      labs(x = paste0("PC", k1), y = paste0("PC", k2)),
    scale = 0.95
  )
}), nrow = 4)

# ggsave("figures/proj1000G-Celiac.png", width = 11, height = 7)


#### Ancestry estimation ####

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")

plot_grid(plotlist = lapply(5:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = fam2$Population) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}))


seq_PC <- 1:19
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


# Get population from external files
pop.files <- list.files(path = "~/Bureau/Dubois2010_data/",
                        pattern = "cluster_*", full.names = TRUE)
pop <- snp_getSampleInfos(obj.bed, pop.files)[[1]]
pop_celiac <- c("Netherlands", "Italy", "UK", "UK", "Finland")[pop]

mean(pval_max < 0.05) # 29.6%
tapply(pval_max, pop_celiac, function(x) mean(x < 0.05))
#   Finland       Italy Netherlands          UK
# 0.3731283   0.7651588   0.2809466   0.1992893

# https://www.internationalgenome.org/category/population/
table(pop_celiac, substr(choose_pop2, 1, 3))
# pop_celiac     AMR  EUR
# Finland          0 1549
# Italy            0  244
# Netherlands      0 1185
# UK               1 5407

library(magrittr)
table(pop_celiac, choose_pop2, exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Ancestry (left) of Celiac individuals and their matching to 1000G populations (top) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred-celiac",
                 align = "|l|c|c|c|c|c|c|c|") %>%
  print(caption.placement = "top")

#             AMR_MXL EUR_CEU EUR_FIN EUR_GBR EUR_IBS EUR_TSI <NA>
# Finland           0       4    1543       1       0       1  922
# Italy             0       0       0       0      15     229  795
# Netherlands       0     802       0     382       1       0  463
# UK                1    2865       1    2532       4       5 1346
