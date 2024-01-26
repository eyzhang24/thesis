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
create_plot <- function(xax, yax, Y) {
  plot_ly(x = ~xax, y = ~yax, z = ~Y, intensity = ~Y) |> 
    add_trace(type = "mesh3d") |> 
    layout(scene = list(
      xaxis = list(rangemode = "normal",
                   showgrid = FALSE,
                   showline = TRUE, 
                   mirror = TRUE, 
                   ticks = "outside", 
                   title = "x1"), 
      yaxis = list(rangemode = "normal", 
                   showgrid = FALSE,
                   showline = TRUE, 
                   mirror = TRUE, 
                   ticks = "outside", 
                   title = "x2"), 
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
fig00 <- create_plot(x1, x2, y00)
fig00
save_image(fig00, "index/figures/surfaces/p00.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, smaller effect size
y01 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.3*x1*x2) 
fig01 <- create_plot(x1, x2, y01)
fig01
save_image(fig01, "index/figures/surfaces/p01.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, larger effect size
y02 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.6*x1*x2) 
fig02 <- create_plot(x1, x2, y02)
fig02
save_image(fig02, "index/figures/surfaces/p02.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, smaller effect size
y03 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.1*x1*((x2-1)^2))
fig03 <- create_plot(x1, x2, y03)
fig03
save_image(fig03, "index/figures/surfaces/p03.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, larger effect size
y04 <- with(data, x1 + 3/(1+exp(-4*x2)) + 0.2*x1*((x2-1)^2))
fig04 <- create_plot(x1, x2, y04)
fig04
save_image(fig04, "index/figures/surfaces/p04.png", 
           width = 720, height = 480, scale = 3)

# sinusoidal interaction, smaller effect size
y05 <- with(data, x1 + 3/(1+exp(-4*x2)) + sin(0.25*pi*x1*x2))
fig05 <- create_plot(x1, x2, y05)
fig05
save_image(fig05, "index/figures/surfaces/p05.png", 
           width = 720, height = 480, scale = 3)

# sinusoidal interaction, larger effect size
y06 <- with(data, x1 + 3/(1+exp(-4*x2)) + 2*sin(0.25*pi*x1*x2))
fig06 <- create_plot(x1, x2, y06)
fig06
save_image(fig06, "index/figures/surfaces/p06.png", 
           width = 720, height = 480, scale = 3)

##########
# marginally insignificant (Cd and As)
##########

# multiplicative interaction, smaller effect size
y11 <- with(data, 0.5*x1*x2) 
fig11 <- create_plot(x1, x2, y11)
fig11
save_image(fig11, "index/figures/surfaces/p11.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, larger effect size
y12 <- with(data, x1*x2) 
fig12 <- create_plot(x1, x2, y12)
fig12
save_image(fig12, "index/figures/surfaces/p12.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, smaller effect size
y13 <- with(data, 0.1*x1*((x2-1)^2))
fig13 <- create_plot(x1, x2, y13)
fig13
save_image(fig13, "index/figures/surfaces/p13.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, larger effect size
y14 <- with(data, 0.2*x1*((x2-1)^2))
fig14 <- create_plot(x1, x2, y14)
fig14
save_image(fig14, "index/figures/surfaces/p14.png", 
           width = 720, height = 480, scale = 3)

# sinusoidal interaction, smaller effect size
y15 <- with(data, sin(0.25*pi*x1*x2))
fig15 <- create_plot(x1, x2, y15)
fig15
save_image(fig15, "index/figures/surfaces/p15.png", 
           width = 720, height = 480, scale = 3)

# sinusoidal interaction, larger effect size
y16 <- with(data, 2*sin(0.25*pi*x1*x2))
fig16 <- create_plot(x1, x2, y16)
fig16
save_image(fig16, "index/figures/surfaces/p16.png", 
           width = 720, height = 480, scale = 3)

##########
# highly correlated (Ni and Co)
##########

# multiplicative interaction, smaller effect size
y21 <- with(data, 3/(1+exp(-4*x2)) + 0.3*x1*x2) 
fig21 <- create_plot(x1, x2, y21)
fig21
save_image(fig21, "index/figures/surfaces/p21.png", 
           width = 720, height = 480, scale = 3)

# multiplicative interaction, larger effect size
y22 <- with(data, 3/(1+exp(-4*x2)) + 0.6*x1*x2) 
fig22 <- create_plot(x1, x2, y22)
fig22
save_image(fig22, "index/figures/surfaces/p22.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, smaller effect size
y23 <- with(data, 3/(1+exp(-4*x2)) + 0.1*x1*((x2-1)^2))
fig23 <- create_plot(x1, x2, y23)
fig23
save_image(fig23, "index/figures/surfaces/p23.png", 
           width = 720, height = 480, scale = 3)

# polynomial interaction, larger effect size
y24 <- with(data, 3/(1+exp(-4*x2)) + 0.2*x1*((x2-1)^2))
fig24 <- create_plot(x1, x2, y24)
fig24
save_image(fig24, "index/figures/surfaces/p24.png", 
           width = 720, height = 480, scale = 3)

# sinusoidal interaction, smaller effect size
y25 <- with(data, 3/(1+exp(-4*x2)) + sin(0.25*pi*x1*x2))
fig25 <- create_plot(x1, x2, y25)
fig25
save_image(fig25, "index/figures/surfaces/p25.png", 
           width = 720, height = 480, scale = 3)

# sinusoidal interaction, larger effect size
y26 <- with(data, 3/(1+exp(-4*x2)) + 2*sin(0.25*pi*x1*x2))
fig26 <- create_plot(x1, x2, y26)
fig26
save_image(fig26, "index/figures/surfaces/p26.png", 
           width = 720, height = 480, scale = 3)

