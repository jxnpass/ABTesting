library(tidyverse)

my_theme <-
  theme_classic() +
  theme(
    legend.title = element_text(family = "Arial", face = "bold"),
    text = element_text(family = "Arial", size = 10), 
    plot.title = element_text(size = 12, face = "bold", hjust = .5), 
    plot.subtitle = element_text(size = 11, hjust = .5), 
    plot.caption = element_text(size = 10), 
    legend.position = 'bottom',
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) 
