library(tidyverse)

# Hg + 3/(1+exp(-4*Ni)) - (Sb^2) + 0.5*Sb + 1.5/(1+exp(-4*Sn)) + 

p1 <- ggplot(NULL) +
  geom_function(fun = function(x) x, 
                color = "darkorchid1") +
  xlim(-2, 2) +
  labs(x = "Hg")

p2 <- ggplot(NULL) +
  geom_function(fun = function(x) 3/(1+exp(-4*x)), 
                color = "deepskyblue3") +
  xlim(-2, 2) +
  labs(x = "Ni")

p3 <- ggplot(NULL) +
  geom_function(fun = function(x) 1.5/(1+exp(-4*x)), 
                color = "darkorange") +
  xlim(-2, 2) +
  labs(x = "Sn")

p4 <- ggplot(NULL) +
  geom_function(fun = function(x) (x^2) + 0.5*x, 
                color = "coral1") +
  xlim(-2, 2) +
  labs(x = "Sb")

cowplot::plot_grid(p1, p2, p3, p4)
ggsave("index/figures/univariatelines.png", width = 6, height = 4)
