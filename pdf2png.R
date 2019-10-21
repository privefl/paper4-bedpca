library(magick)

"figures/us.pdf" %>% {
  svg <- sub("\\.pdf$", ".svg", .)
  img <- image_read_pdf(.)
  image_write(img, svg)
  print(svg)
}
