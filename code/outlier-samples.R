#### Detection of individual outliers ####

library(bigsnpr)
library(ggplot2)

obj.bed <- bed(download_1000G("tmp-data"))
obj.svd <- bed_autoSVD(obj.bed, k = 20, ncores = nb_cores())

prob <- bigutilsr::prob_dist(obj.svd$u)
hist(S <- prob$dist.self / sqrt(prob$dist.nn))

plot_grid(plotlist = lapply(7:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S) +
    scale_colour_viridis_c()
}), scale = 0.95)

# ggsave("figures/outliers-1000G.pdf", width = 7, height = 5)


#### Restricting to homogeneous individuals ####

obj.bed <- bed("~/Bureau/Dubois2010_data/FinnuncorrNLITUK3hap550_QC.bed")
obj.svd <- bed_autoSVD(obj.bed, k = 10, ncores = nb_cores())

plot(obj.svd, type = "scores", scores = 1:10, coeff = 0.7)

dist <- bigutilsr::covRob(obj.svd$u, estim = "pairwiseGK")$dist
pval <- pchisq(dist, df = 10, lower.tail = FALSE)

# Get population from external files
pop.files <- list.files(path = "~/Bureau/Dubois2010_data/",
                        pattern = "cluster_*", full.names = TRUE)
pop <- snp_getSampleInfos(obj.bed, pop.files)[[1]]
pop_celiac <- c("Netherlands", "Italy", "UK", "UK", "Finland")[pop]

library(ggplot2)
plot_grid(plotlist = lapply(1:2, function(k) {
  plot_grid(
    plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.7) +
      aes(color = pop_celiac) +
      labs(color = "Population"),
    plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.7) +
      aes(color = pval < 0.001) +
      scale_colour_viridis_d(),
    ncol = 1
  )
}), scale = 0.95)

# ggsave("figures/homogeneous.pdf", width = 8, height = 6)
