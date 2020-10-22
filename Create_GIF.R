library(magick)
library(purrr)
library(tidyverse)

list.files(path = "~/Desktop/Bioavailability_MS/remin_gif/", pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps = 0.5) %>% # animates, can opt for number of loops
  image_write("~/Desktop/Bioavailability_MS/remin_gif/remin.gif") 


