
library(tidyverse)
library(FME)
library(MEIGOR)
library(GGally)
library(lamW)
library(cowplot)
library(ggthemes)
library(wesanderson)

source("biodesigner.R")

p_Rtkw <- list( 
    Xmin = -1.425, 
    Xmax = 44.36, 
    b = .027, 
    c = .4)

p_CPM <- list(  # Silva et al. (https://doi.org/10.1016/j.foodres.2020.109476)
    Xmin = -1.425, 
    Xopt = 38.17,
    Xmax = 44.36, 
    n = 2,
    mu_opt = 0.976)


## sup. figure 1

my_cols <- wes_palette("Cavalcanti1", 4)

p <- tibble(x = seq(-1.5, 45.5, length = 1000),
       Ratkowsky = Ratkowsky(p_Rtkw, x)$sq_y,
       CPM = CPM_model(p_CPM, x)$sq_y) %>%
    pivot_longer(-x) %>%
    ggplot() +
    geom_line(aes(x = x, y = value, colour = name),
              linewidth = 1) +
    xlab("Temperature (ºC)") + ylab("Square root of mu (log CFU/h)^0.5") +
    theme_hc(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "top") +
    scale_colour_manual(values = my_cols)


ggsave(p, filename = "suppFigure1.png",
       width = 9, height = 6)

