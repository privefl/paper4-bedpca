library(bigsnpr)
library(dplyr)
library(ggplot2)

bedfiles <-
  sapply(c("~/Bureau/thesis/backingfiles/genome_Florian_Prive_v5_Full_20190212023418.txt",
           "tmp-data/genome_B_V_v5_Full_20191017230843.txt"),
         function(genome) {
           name <- sub("\\.txt$", "", basename(genome))
           genome <- bigreadr::fread2(genome)
           bim <- bigreadr::fread2(sub_bed(download_1000G('tmp-data'), '.bim'))
           genome2 <- genome %>%
             tidyr::separate(genotype, c("g1", "g2"), sep = 1) %>%
             mutate(chromosome = as.integer(chromosome)) %>%
             left_join(bim, by = c(chromosome = "V1", position = "V4")) %>%
             mutate_at(c("g1", "g2"), ~ {
               ifelse(. %in% c("A", "T", "C", "G"), . == V5, NA)
             })
           big_snp <- snp_fake(1, nrow(genome2))
           big_snp$map[c(1:2, 4:6)] <- genome2[c(2, 1, 3, 8:9)]
           big_snp$fam$sample.ID <- name
           g <- genome2$g1 + genome2$g2
           big_snp$genotypes[] <- ifelse(is.na(g), 3, g)
           snp_writeBed(big_snp, tempfile(fileext = ".bed"))
         })

write(sapply(bedfiles, sub_bed), tmp <- tempfile(), ncolumns = 1)
plink <- download_plink("tmp-data")

system(glue::glue(
  "{plink} --merge-list {tmp}",
  " --autosome --geno 0",
  " --make-bed --out tmp-data/our_genome"
))

obj.bed <- bed("tmp-data/our_genome.bed")
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

plot_grid(plotlist = lapply(1:3, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[, k1], PC.ref[, k2], size = I(2), color = pop2) +
    theme_bigstatsr(0.7) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2))
}))

seq_PC <- 1:4
pop_PCs <- vctrs::vec_split(PC.ref, pop2)

all_pval <- sapply(pop_PCs$val, function(PC) {
  maha <- bigutilsr::covRob(PC[, seq_PC], estim = "pairwiseGK")
  dist <- mahalanobis(proj2[, seq_PC], center = maha$center, cov = maha$cov)
  pval <- pchisq(dist, df = length(seq_PC), lower.tail = FALSE)
})

colnames(all_pval) <- pop_PCs$key
round(100 * all_pval, 2)
#      Central Europe SE Europe Eastern Europe Anglo-Irish Isles Scandinavia
# [1,]           0.05         0              0             45.63       85.82
# [2,]           1.24         0              0              0.22        0.18
#      Italy Germany SW Europe Switzerland Belgium France Netherlands
# [1,]  0.00    1.64         0        0.00    0.01   0.06        0.00
# [2,]  0.01    9.63         0       28.37   30.90  24.60        0.79

qplot(PC.ref[, 2], PC.ref[, 1], size = I(2), color = pop2) +
  stat_ellipse(geom = "polygon", alpha = 0.5, aes(fill = pop2)) +
  theme_bigstatsr() +
  labs(x = paste0("PC", 2), y = paste0("PC", 1)) +
  geom_point(aes(proj2[, 2], proj2[, 1]), size = 3, color = "black") +
  labs(color = "Population", fill = "Population") +
  scale_x_reverse() +
  scale_y_reverse() +
  coord_equal()

# ggsave("figures/us.pdf", width = 8, height = 8)
