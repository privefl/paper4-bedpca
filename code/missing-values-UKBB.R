library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)
library(runonce)

rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
ind.row <- which(!fam$V2 %in% rel$ID2)
obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")


#### IMPUTATION BY THE MEAN -> see pca-UKBB.R ####

obj.svd <- readRDS("tmp-results/SVD_UKBB.rds")
ind.keep <- attr(obj.svd, "subset")

counts <- bed_counts(obj.bed, ind.row, ind.keep, byrow = TRUE, ncores = nb_cores())
het <- counts[2, ] / colSums(counts[1:3, ])
nbNA <- counts[4, ]
PC <- predict(obj.svd)
summary(mylm <- lm(as.formula(paste0("het ~ (", paste0("PC[ ,", 1:6, "]", collapse = " + "), ")^2"))))
het_adj <- mylm$coefficients[1] + mylm$residuals
qplot(nbNA, het_adj, alpha = I(0.5), col = PC[, 16]) +
  scale_x_log10() +
  scale_color_viridis_c(direction = -1) +
  theme_bigstatsr() +
  labs(x = "Number of missing values (log-scale)",
       y = "PC-corrected heterozygosity",  color = "PC16")

ggsave("figures/het-col-PC16.png", width = 9, height = 6)


#### USING DOSAGES ####

info_snp_UKBB <- rbind_df(lapply(1:22, function(chr) {
  file <- paste0("data/ukb_imp_mfi/ukb_mfi_chr", chr, "_v3.txt")
  df <- fread2(file, select = c(3:5, 8), col.names = c("pos", "a0", "a1", "info"))
  cbind.data.frame(chr = chr, df)
}))

# subset samples
sample <- fread2("~/UKBiobank/ukb25589_imp_chr1_v3_s487327.sample")
str(sample)
(N <- readBin("data/ukb_imp_chr1_v3.bgen", what = 1L, size = 4, n = 4)[4])
sample <- sample[-1, ]
stopifnot(nrow(sample) == N)
ind.row2 <- na.omit(match(obj.bed$fam$sample.ID[ind.row], sample$ID_2))

info_snp <- bigsnpr::snp_match(
  setNames(obj.bed$map[ind.keep, ], c("chr", "rsid", "beta", "pos", "a1", "a0")),
  info_snp_UKBB, strand_flip = FALSE)
# 240,444 variants to be matched.
# 240,332 variants have been matched; 0 were flipped and 31,017 were reversed.

list_snp_id <- with(info_snp, split(paste(chr, pos, a0, a1, sep = "_"),
                                    factor(chr, levels = 1:22)))

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("data/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_PCA_BGEN",
    ind_row = ind.row2,
    bgi_dir = "data/ukb_imp_bgi",
    ncores = 12
  )
) # 42 min

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_PCA_BGEN.rds")
G <- ukbb$genotypes
file.size(G$backingfile) / 1024^3  # 91 GB
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

obj.svd2 <- save_run(
  big_randomSVD(G, fun.scaling = snp_scaleBinom(), k = 20, ncores = nb_cores()),
  file = "tmp-results/SVD_UKBB2.rds"
) # 76 min

round(100 * cor(obj.svd2$u, obj.svd$u[-attr(ind.row2, "na.action"), ]), 1)


#### USING RANDOM IMPUTATION ACCORDING TO FREQUENCIES ####

system.time(
  snp_readBed2(obj.bed$bedfile, ind.row = ind.row, ind.col = ind.keep,
               backingfile = "data/UKBB_PCA")
) # 6 min
ukbb <- snp_attach("data/UKBB_PCA.rds")
G <- ukbb$genotypes

system.time(
  G2 <- snp_fastImputeSimple(G, method = "random", ncores = nb_cores())
) # 2 min

obj.svd3 <- save_run(
  big_randomSVD(G2, fun.scaling = snp_scaleBinom(), k = 20, ncores = nb_cores()),
  file = "tmp-results/SVD_UKBB3.rds"
) # 91 min

round(100 * cor(obj.svd3$u, obj.svd$u), 1)

#### USING IMPUTATION WITH MACHINE LEARNING ####

system.time(
  snp_readBed2(obj.bed$bedfile, ind.row = ind.row, ind.col = ind.keep,
               backingfile = "data/UKBB_PCA2")
) # 6 min
ukbb <- snp_attach("data/UKBB_PCA2.rds")
G <- ukbb$genotypes

system.time(
  infos <- snp_fastImpute(G, infos.chr = ukbb$map$chromosome,
                          n.cor = 20e3, ncores = nb_cores())
) # 97H (4 days)

pvals <- c(0.01, 0.005, 0.002, 0.001); colvals <- 2:5
df <- data.frame(pNA = infos[1, ], pError = infos[2, ])
plot(subset(df, pNA > 0.001), pch = 20)
idc <- lapply(seq_along(pvals), function(i) {
  curve(pvals[i] / x, from = 0, lwd = 2,
        col = colvals[i], add = TRUE)
})
legend("topright", legend = pvals, title = "p(NA & Error)",
       col = colvals, lty = 1, lwd = 2)

G2 <- G$copy(CODE_IMPUTE_PRED)
obj.svd4 <- save_run(
  big_randomSVD(G2, fun.scaling = snp_scaleBinom(), k = 20, ncores = nb_cores()),
  file = "tmp-results/SVD_UKBB4.rds"
) # 53 min

round(100 * cor(obj.svd4$u, obj.svd$u), 1)

obj.svd5 <- save_run(
  big_randomSVD(G2, fun.scaling = snp_scaleBinom(),
                ind.row = rows_along(G2)[-attr(ind.row2, "na.action")],
                k = 20, ncores = nb_cores()),
  file = "tmp-results/SVD_UKBB5.rds"
) # 53 min

round(100 * cor(obj.svd5$u, obj.svd$u[-attr(ind.row2, "na.action"), ]), 1)


obj.svd6 <- save_run(
  bed_randomSVD(obj.bed, ind.row = ind.row[-attr(ind.row2, "na.action")],
                ind.col = ind.keep, k = 20, ncores = nb_cores()),
  file = "tmp-results/SVD_UKBB6.rds"
) # 50 min

round(100 * cor(obj.svd6$u, obj.svd$u[-attr(ind.row2, "na.action"), ]), 1) %>%
  xtable::xtable(digits = 1)

round(100 * cor(obj.svd6$u, obj.svd2$u), 1) %>%
  xtable::xtable(digits = 1)
