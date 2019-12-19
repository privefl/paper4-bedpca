plot_grid2 <- function(plotlist, ..., title_ratio = 1/8, legend_ratio = 1/6) {

  main_grid <- plot_grid(plotlist = lapply(plotlist, function(p) {
    p + theme(legend.position = "none") + ggtitle(NULL)
  }), ...)

  plot_grid(
    cowplot::get_title(plotlist[[1]]),  ggplot() + theme_minimal(),
    main_grid,                          cowplot::get_legend(plotlist[[1]]),
    rel_heights = c(title_ratio, 1), ncol = 2,
    rel_widths = c(1, legend_ratio), nrow = 2
  )
}
