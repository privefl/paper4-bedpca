library(bigsnpr)
library(ggplot2)

plink <- download_plink("tmp-data")
mygenome <- "~/Bureau/thesis/backingfiles/genome_Florian_Prive_v5_Full_20190212023418.txt"

system(glue::glue(
  "{plink} --23file {mygenome}",
  " --autosome --geno 0 --snps-only",
  " --a2-allele {sub_bed(download_1000G('tmp-data'), '.bim')} 5 2",
  " --make-bed --out tmp-data/mygenome"
))

obj.bed <- bed("tmp-data/mygenome.bed")
bed.ref <- bed("~/Bureau/POPRES_data/POPRES_allchr_QC.bed")

system.time(
  test <- bed_projectPCA(bed.ref, obj.bed, k = 10, ncores = nb_cores(),
                         build.ref = "hg18", liftOver = "tmp-data/liftOver",
                         match.min.prop = 0, roll.size = 20)
) # using 18,764 variants only

PC.ref <- predict(test$obj.svd.ref)
proj2 <- test$OADP_proj


#### Ancestry estimation ####

pop <- bed.ref$fam$family.ID
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

plot_grid(plotlist = lapply(1:2, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = pop2) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}))

unAsIs <- function(X) {
  class(X) <- setdiff(class(X), "AsIs")
  X
}

library(dplyr)
seq_PC <- c(1:3)
all_covRob <- cbind.data.frame(pop2, PC = I(PC.ref)) %>%
  group_by(pop2) %>%
  summarize(maha = list(bigutilsr::covRob(unAsIs(PC[, seq_PC]), estim = "pairwiseGK")))
POP <- all_covRob$pop2

all_dist <- sapply(all_covRob$maha, function(maha) {
  mahalanobis(proj2[, seq_PC], center = maha$center, cov = maha$cov)
})
names(all_dist) <- POP

all_prob <- pchisq(all_dist, df = length(seq_PC), lower.tail = FALSE)
round(100 * all_prob, 2)
# Anglo-Irish Isles           Belgium    Central Europe    Eastern Europe
#              0.00              0.53              0.42              0.00
#            France           Germany             Italy       Netherlands
#             43.37              1.70              1.11              0.00
#       Scandinavia         SE Europe         SW Europe       Switzerland
#              0.00              0.00              0.00             25.19

qplot(PC.ref[, 2], PC.ref[, 1], size = I(2), color = pop2) +
  stat_ellipse(geom = "polygon", alpha = 0.5, aes(fill = pop2)) +
  theme_bigstatsr() +
  labs(x = paste0("PC", 2), y = paste0("PC", 1)) +
  geom_point(aes(proj2[2], proj2[1]), size = 5, color = "black") +
  labs(color = "Population", fill = "Population") +
  scale_x_reverse() +
  scale_y_reverse() +
  coord_equal()

ggsave("figures/me.pdf", width = 7, height = 5)
