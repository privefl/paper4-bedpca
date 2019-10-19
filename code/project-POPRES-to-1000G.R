library(bigsnpr)
library(ggplot2)

snp_plinkQC(
  plink.path = download_plink("tmp-data"),
  prefix.in = "~/Bureau/POPRES_data/POPRES_allchr",
  geno = 0.01, mind = 0.01,
  autosome.only = TRUE
)
obj.bed <- bed("~/Bureau/POPRES_data/POPRES_allchr_QC.bed")
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


#### Ancestry estimation ####

fam2 <- bigreadr::fread2("tmp-data/1000G_phase3_common_norel.fam2")

plot_grid(plotlist = lapply(5:10, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = fam2$Population) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}))


seq_PC <- c(1:18, 20)
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


pop <- obj.bed$fam$family.ID
pop2 <- dplyr::case_when(
  pop %in% c("Portugal", "Spain") ~ "SW Europe",
  pop %in% c("?orway", "Sweden", "Finland", "Denmark") ~ "Scandinavia",
  pop %in% c("Hungary", "Slovakia", "Czech", "Austria", "Slovenia",
             "Croatia") ~ "Central Europe",
  pop %in% c("Russian", "Latvia", "Ukraine", "Poland") ~ "Eastern Europe",
  pop %in% c("Greece", "Turkey", "Serbia", "Cyprus", "Albania", "Kosovo",
             "Bosnia", "Macedonia,",  "Romania", "Bulgaria") ~ "SE Europe",
  pop %in% c("Swiss-German", "Swiss-French", "Swiss-Italian") ~ "Switzerland",
  pop == "?etherlands" ~ "Netherlands",
  pop %in% c("United", "Scotland", "Ireland") ~ "Anglo-Irish Isles",
  TRUE ~ pop
)

mean(pval_max < 0.05) # 56.5%
tapply(pval_max, pop2, function(x) mean(x < 0.05))
# Anglo-Irish Isles           Belgium    Central Europe    Eastern Europe
#         0.1879699         0.2790698         0.9272727         0.9333333
#            France           Germany             Italy       Netherlands
#         0.6966292         0.4788732         0.8584475         0.2941176
#       Scandinavia         SE Europe         SW Europe       Switzerland
#         0.6666667         0.9468085         0.2689394         0.8198198


# https://www.internationalgenome.org/category/population/
table(substr(choose_pop2, 1, 3), pop2)
# All EUR

table(pop2, choose_pop2, exclude = NULL) %>%
  print() %>%
  { ifelse(. == 0, "", .) } %>%
  xtable::xtable(caption = "Ancestry (left) of POPRES individuals and their matching to 1000G populations (top) by our method. See the description of 1000G populations at \\url{https://www.internationalgenome.org/category/population/}.",
                 label = "tab:ancestry-pred-popres",
                 align = "|l|c|c|c|c|c|") %>%
  print(caption.placement = "top")

#                   EUR_CEU EUR_GBR EUR_IBS EUR_TSI <NA>
# Anglo-Irish Isles     159      57       0       0   50
# Belgium                28       3       0       0   12
# Central Europe          3       0       0       1   51
# Eastern Europe          2       0       0       0   28
# France                 14       0      13       0   62
# Germany                37       0       0       0   34
# Italy                   0       0      12      19  188
# Netherlands            10       2       0       0    5
# Scandinavia             5       0       0       0   10
# SE Europe               0       0       0       5   89
# SW Europe               0       0     193       0   71
# Switzerland            33       1       6       0  182
