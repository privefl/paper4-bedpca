#### Detection of individual outliers ####

library(bigsnpr)
library(ggplot2)

obj.bed <- bed(download_1000G("tmp-data"))
obj.svd <- bed_autoSVD(obj.bed, k = 20, ncores = nb_cores())

prob <- bigutilsr::prob_dist(obj.svd$u, ncores = nb_cores())
S <- prob$dist.self / sqrt(prob$dist.nn)

ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(breaks = 0:5 / 5, limits = c(0, NA)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)")

# ggsave("figures/hist-outliers-1000G.pdf", width = 9, height = 6)

plot_grid(plotlist = lapply(7:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S) +
    scale_colour_viridis_c()
}), scale = 0.95)

# ggsave("figures/outliers-1000G.pdf", width = 7, height = 5)

plot_grid(plotlist = lapply(7:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S > 0.4) +
    scale_colour_viridis_d()
}), scale = 0.95)

# ggsave("figures/be-outlier-1000G.pdf", width = 7, height = 5)

library(magrittr)
S2 <- apply(obj.svd$u, 2, function(x) abs(x - mean(x)) / sd(x)) %>%
  apply(1, max)

ggplot() +
  geom_histogram(aes(S2), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)") +
  geom_vline(xintercept = 6, color = "red")

# ggsave("figures/hist-outliers2-1000G.pdf", width = 9, height = 6)

plot_grid(plotlist = lapply(7:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S2) +
    scale_colour_viridis_c(trans = "log", breaks = c(1, 3, 6, 20))
}), scale = 0.95)

# ggsave("figures/outliers2-1000G.pdf", width = 7, height = 5)

plot_grid(plotlist = lapply(7:10, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S2 > 6) +
    scale_colour_viridis_d()
}), scale = 0.95)

# ggsave("figures/be-outlier2-1000G.pdf", width = 7, height = 5)
