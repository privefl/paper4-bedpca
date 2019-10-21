library(magick)

"figures/proj1000G-related.pdf" %>% {
  svg <- sub("\\.pdf$", ".svg", .)
  img <- image_read_pdf(.)
  image_write(img, svg)
  print(svg)
}
