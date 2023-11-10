library(tidyverse)
library(bkmr)
library(stats)
library(splines)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

########
# generate simulated points
########

# generate data from distribution
set.seed(0) # reproducibility
x <- seq(0, 25, length.out = 50)
Y <- exp(x/10) + 2*sin(x/2) + rnorm(50, mean = 0, sd = 0.5)
df <- data.frame(x, Y)

# plot data and linear regression line
ggplot(df, aes(x, Y)) +
  geom_point() +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), color = "darkorange") + 
  geom_smooth(method = "lm", formula = "y~x", 
              color = "deepskyblue3", fill = "gray70", 
              linetype = "dashed", linewidth = 0.5)

########
# kernel regression
########

# get normal distribution of weights around query points
df$Weight <- dnorm(df$x, mean = 12.5, sd = 1)

# plot points colored by their weights
p1 <- ggplot(df, aes(x, Y)) +
  geom_point(aes(color = Weight)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), color = "darkorange") + 
  geom_vline(xintercept = 12.5, linetype = "dotted") +
  theme(legend.position = c(0.95, 0.3))#c(0.1, 0.7))

# plot a curve of weights
normcurv <- data.frame(x = seq(0, 25, length.out = 250)) 
normcurv$Y <- dnorm(normcurv$x, mean = 12.5, sd = 1)
p2 <- ggplot(normcurv, aes(x, Y, color = Y)) +
  geom_line() +
  theme(legend.position = "none")

# stitch plots together
cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.7, 0.3))

# fit kernel regression
kmr_toy <- ksmooth(df$x, df$Y, kernel = "normal", bandwidth = 8/3, x.points = df$x)
df <- df |> 
  left_join(as.data.frame(kmr_toy), by = "x") |> 
  rename(Yhat = y)

# plot kernel regression estimation
ggplot(df) +
  geom_point(aes(x, Y)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), color = "darkorange") +
  geom_line(aes(x, Yhat), color = "deepskyblue3", linetype = "dashed") 

########
# spline regression
########

# fit spline regression
kn <- c(25/3, 50/3) # 3 knots of equal width
spline_toy <- lm(Y ~ bs(x, knots = kn), data = df)
p <- predict(spline_toy, se = T)
df$Yhats <- p$fit
df$Ylows <- df$Yhats - 1.96 * p$se.fit
df$Yups <- df$Yhats + 1.96 * p$se.fit

# plot spline regression estimation
ggplot(df) +
  geom_point(aes(x, Y)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), color = "darkorange") +
  geom_line(aes(x, Yhats), color = "deepskyblue3", linetype = "dashed") +
  geom_ribbon(aes(x, ymin = Ylows, ymax = Yups), fill = "gray70", alpha = 0.3) +
  geom_vline(xintercept = kn, linetype = "dotted")