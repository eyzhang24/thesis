library(plotly)
if(!require(reticulate)) {
  install.packages('reticulate')
  reticulate::install_miniconda()
  reticulate::conda_install('r-reticulate', 'python-kaleido')
  reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
  Sys.setenv(RETICULATE_PYTHON = 
               '/Users/elizabethzhang/Library/r-miniconda-arm64/envs/r-reticulate/bin/python')
  reticulate::use_miniconda('r-reticulate')
}
library(tidyverse)

# load in observed data
comb_small <- read_csv("madres_data/base_data.csv")

# log-transform and scale target data
comb_scale <- comb_small |> 
  mutate(across(10:19, ~scale(log(.)))) 

# check ranges of scaled predictors
range(comb_scale$Hg)
range(comb_scale$Ni)
range(comb_scale$Cd)
range(comb_scale$As)
range(comb_scale$Co)

# generate data covering 2d predictor surface
data <- expand.grid(x1 = seq(-3, 3, by = 0.1), 
                    x2 = seq(-3, 3, by = 0.1))
x1 <- data$x1
x2 <- data$x2

# function to create plot
create_plot <- function(xax, yax, Y, xname = NA, yname = NA) {
  plot_ly(x = ~xax, y = ~yax, z = ~Y, intensity = ~Y) |> 
    add_trace(type = "mesh3d") |> 
    layout(scene = list(
      xaxis = list(rangemode = "normal",
                   showgrid = FALSE,
                   showline = TRUE, 
                   mirror = TRUE, 
                   ticks = "outside", 
                   title = xname), 
      yaxis = list(rangemode = "normal", 
                   showgrid = FALSE,
                   showline = TRUE, 
                   mirror = TRUE, 
                   ticks = "outside", 
                   title = yname), 
      zaxis = list(rangemode = "normal", 
                   showgrid = FALSE,
                   showline = TRUE, 
                   mirror = TRUE, 
                   ticks = "outside", 
                   title = "Y"), 
      aspectmode = "cube"
    ))
}

##########
# marginally significant (Hg and Ni)
##########

# no interaction
y00 <- with(data, x1 + 3/(1+exp(-4*x2))) 
fig00 <- create_plot(x1, x2, y00, "Hg", "Ni")  |> 
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = 0.1))))
fig00
save_image(fig00, "index/figures/surfaces/p00.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, smaller effect size
yam1 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.35*x1*x2) 
figam1 <- create_plot(x1, x2, yam1, "Hg", "Ni") |> 
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = .1))))
figam1
save_image(figam1, "index/figures/surfaces/am1.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, larger effect size
yam2 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.75*x1*x2) 
figam2 <- create_plot(x1, x2, yam2, "Hg", "Ni") |> 
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = .1))))
figam2
save_image(figam2, "index/figures/surfaces/am2.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, smaller effect size
yap1 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.13*x1*((x2-1)^2))
figap1 <- create_plot(x1, x2, yap1, "Hg", "Ni")|> 
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = .1))))
figap1
save_image(figap1, "index/figures/surfaces/ap1.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, larger effect size
yap2 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.26*x1*((x2-1)^2))
figap2 <- create_plot(x1, x2, yap2, "Hg", "Ni") |> 
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = .1))))
figap2
save_image(figap2, "index/figures/surfaces/ap2.png", 
           width = 720, height = 480, scale = 3)

##########
# marginally insignificant (Cd and As)
##########

# multiplicative interaction, smaller effect size
ybm1 <- with(data, 0.35*x1*x2) 
figbm1 <- create_plot(x1, x2, ybm1, "Cd", "As") |> 
  layout(scene = list(camera = list(eye = list(x = 1.4, y = 1.4, z = 1.2))))
figbm1
save_image(figbm1, "index/figures/surfaces/bm1.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, larger effect size
ybm2 <- with(data, 0.7*x1*x2) 
figbm2 <- create_plot(x1, x2, ybm2, "Cd", "As") |> 
  layout(scene = list(camera = list(eye = list(x = 1.4, y = 1.4, z = 1.2))))
figbm2
save_image(figbm2, "index/figures/surfaces/bm2.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, smaller effect size
ybp1 <- with(data, 0.125*x1*((x2-1)^2)) 
figbp1 <- create_plot(x1, x2, ybp1, "Cd", "As") |> 
  layout(scene = list(camera = list(eye = list(x = 1.4, y = 1.4, z = 1.2))))
figbp1
save_image(figbp1, "index/figures/surfaces/bp1.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, larger effect size
ybp2 <- with(data, 0.25*x1*((x2-1)^2))
figbp2 <- create_plot(x1, x2, ybp2, "Cd", "As") |> 
  layout(scene = list(camera = list(eye = list(x = 1.4, y = 1.4, z = 1.2))))
figbp2
save_image(figbp2, "index/figures/surfaces/bp2.png", 
           width = 720, height = 480, scale = 3)

##########
# highly correlated (Ni and Co)
##########

# multiplicative interaction, smaller effect size
ycm1 <- with(data, 3/(1+exp(-4*x2)) + 0.3*x1*x2) 
figcm1 <- create_plot(x1, x2, ycm1, "Ni", "Co") |> 
  layout(scene = list(camera = list(eye = list(x = 1.2, y = 1.2, z = 1.5))))
figcm1
save_image(figcm1, "index/figures/surfaces/cm1.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, larger effect size
ycm2 <- with(data, 3/(1+exp(-4*x2)) + 0.6*x1*x2) 
figcm2 <- create_plot(x1, x2, ycm2, "Ni", "Co") |> 
  layout(scene = list(camera = list(eye = list(x = 1.2, y = 1.2, z = 1.5))))
figcm2
save_image(figcm2, "index/figures/surfaces/cm2.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, smaller effect size
ycp1 <- with(data, 3/(1+exp(-4*x2)) + 0.1*x1*((x2-1)^2))
figcp1 <- create_plot(x1, x2, ycp1, "Ni", "Co") |> 
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = .1))))
figcp1
save_image(figcp1, "index/figures/surfaces/cp1.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, larger effect size
ycp2 <- with(data, 3/(1+exp(-4*x2)) + 0.2*x1*((x2-1)^2))
figcp2 <- create_plot(x1, x2, ycp2, "Ni", "Co") |> 
  layout(scene = list(camera = list(eye = list(x = 1.5, y = 1.5, z = .1))))
figcp2
save_image(figcp2, "index/figures/surfaces/cp2.png", 
           width = 720, height = 480, scale = 3)

##########
# for presentation
##########

# no interaction
y98 <- with(data, x1 + x2) 
fig98 <- create_plot(x1, x2, y98, "x1", "x2")  |> 
  layout(scene = list(camera = list(eye = list(x = -1, y = 2, z = .1))))
fig98
save_image(fig98, "misc_code/p98.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, smaller effect size
y99 <- with(data, x1 + x2 + 0.5*x1*x2) 
fig99 <- create_plot(x1, x2, y99, "x1", "x2") |> 
  layout(scene = list(camera = list(eye = list(x = -1, y = 2, z = 0.1))))
fig99
save_image(fig99, "misc_code/p99.png", 
           width = 720, height = 480, scale = 3)


