library(magick)

"figures/us.pdf" %>% {
  svg <- sub("\\.pdf$", ".svg", .)
  img <- image_read_pdf(.)
  image_write(img, svg, format = "SVG")
  print(svg)
}

"figures/outliers-1000G.pdf" %>% {
  svg <- sub("\\.pdf$", ".png", .)
  img <- image_read_pdf(.)
  image_write(img, svg, format = "PNG")
  print(svg)
}
