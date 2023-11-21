library(tidyverse)
library(stats)
library(splines)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

########
# generate simulated points
########

# generate data from distribution
set.seed(0) # reproducibility
x <- seq(0, 25, length.out = 51)
Y <- exp(x/10) + 2*sin(x/2) + rnorm(51, mean = 0, sd = 0.5)
df <- data.frame(x, Y)

# plot data and linear regression line
q1 <- ggplot(df, aes(x, Y)) +
  geom_point() +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                linetype = "dashed", color = "darkorange") + 
  geom_smooth(method = "lm", formula = "y~x", 
              color = "deepskyblue3", fill = "gray70", linewidth = 0.5, se = F)

# save plot
ggsave("index/figures/ch3_toy1.png", plot = q1, device = "png", 
       width = 5, height = 3)

########
# kernel regression
########

# get normal distribution of weights around query points
df$Weight <- dnorm(df$x, mean = 12.5, sd = 1)

# plot points colored by their weights
p1 <- ggplot(df, aes(x, Y)) +
  geom_point(aes(color = Weight)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                linetype = "dashed", color = "darkorange") + 
  geom_vline(xintercept = 12.5, linetype = "dotted") +
  theme(legend.position = "none")
  # theme(legend.position = c(0.93, 0.29), 
  #       legend.key.size = unit(0.35, 'cm'))#c(0.1, 0.7))

# plot a curve of weights
normcurv <- data.frame(x = seq(0, 25, length.out = 250)) 
normcurv$Weight <- dnorm(normcurv$x, mean = 12.5, sd = 1)
p2 <- ggplot(normcurv, aes(x, Weight, color = Weight)) +
  geom_line() +
  scale_y_continuous(breaks = c(0, 0.2, 0.4)) +
  theme(legend.position = "none") 

# stitch plots together
q2 <- cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.7, 0.3))
q2

# save plot
ggsave("index/figures/ch3_toy2.png", plot = q2, device = "png", 
       width = 5, height = 4)

# fit kernel regression with sigma = 1, bandwidth = 8/3
kmr_toy <- ksmooth(df$x, df$Y, kernel = "normal", bandwidth = 8/3, x.points = df$x)
df <- df |> 
  left_join(as.data.frame(kmr_toy), by = "x") |> 
  rename(Yhat = y)

# plot kernel regression estimation
q3 <- ggplot(df) +
  geom_point(aes(x, Y)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                linetype = "dashed", color = "darkorange") +
  geom_line(aes(x, Yhat), color = "deepskyblue3") 
q3

# save plot
ggsave("index/figures/ch3_toy3.png", plot = q3, device = "png", 
       width = 5, height = 3)

# fit kernel regression with sigma = 5, bandwidth = 40/3
kmr_toy_5 <- ksmooth(df$x, df$Y, kernel = "normal", bandwidth = 40/3, x.points = df$x)

# fit kernel regression with sigma = 0.1, bandwith = 8/30
kmr_toy_1 <- ksmooth(df$x, df$Y, kernel = "normal", bandwidth = 8/30, x.points = df$x)

# re-join data
dfrho <- df |> 
  left_join(as.data.frame(kmr_toy_5), by = "x") |> 
  rename("rho = 50" = y) |> 
  left_join(as.data.frame(kmr_toy_1), by = "x") |> 
  rename("rho = 0.02" = y) |> 
  select(-Yhat) |> 
  pivot_longer(cols = c("rho = 50", "rho = 0.02"), values_to = "Yhat")

# plot kernel regression with two values of rho
qrho <- ggplot(dfrho) +
  geom_point(aes(x, Y)) +
  geom_line(aes(x, Yhat), color = "deepskyblue3") +
  facet_wrap(~name)
qrho

# save plot
ggsave("index/figures/ch3_toyrho.png", plot = qrho, device = "png", 
       width = 7, height = 3)

########
# spline regression
########

kn <- c(5, 10, 15, 20) # 3 knots of equal width

# fit linear spline regression
spline_toy_line <- lm(Y ~ bs(x, knots = kn, degree = 1), data = df)
p_line <- predict(spline_toy_line, se = T)
df$Yhats_line <- p_line$fit
# df$Ylows <- df$Yhats - 1.96 * p$se.fit
# df$Yups <- df$Yhats + 1.96 * p$se.fit

q4 <- ggplot(df) +
  geom_point(aes(x, Y)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                linetype = "dashed", color = "darkorange") +
  geom_line(aes(x, Yhats_line), color = "deepskyblue3") +
  # geom_ribbon(aes(x, ymin = Ylows_line, ymax = Yups_line), fill = "gray70", alpha = 0.3) +
  geom_vline(xintercept = kn, linetype = "dotted")
q4

# save plot
ggsave("index/figures/ch3_toy4.png", plot = q4, device = "png", 
       width = 5, height = 3)

# fit cubic spline regression
spline_toy_cub <- lm(Y ~ bs(x, knots = kn), data = df)
p_cub <- predict(spline_toy_cub, se = T)
df$Yhats_cub <- p_cub$fit

# plot spline regression estimation
q5 <- ggplot(df) +
  geom_point(aes(x, Y)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                linetype = "dashed", color = "darkorange") +
  geom_line(aes(x, Yhats_cub), color = "deepskyblue3") +
  geom_vline(xintercept = kn, linetype = "dotted")
q5

# save plot
ggsave("index/figures/ch3_toy5.png", plot = q5, device = "png", 
       width = 5, height = 3)

# fit natural spline regression
spline_toy_nat <- lm(Y ~ ns(x, knots = kn), data = df)
p_nat <- predict(spline_toy_nat, se = T)
df$Yhats_nat <- p_nat$fit

# plot spline regression estimation
q6 <- ggplot(df) +
  geom_point(aes(x, Y)) +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                linetype = "dashed", color = "darkorange") +
  geom_line(aes(x, Yhats_nat), color = "deepskyblue3") +
  geom_vline(xintercept = c(5, 10, 15, 20), linetype = "dotted")
q6

# save plot
ggsave("index/figures/ch3_toy6.png", plot = q6, device = "png", 
       width = 5, height = 3)

# see what happens outside of the bounds
x_longer <- seq(-5, 30, length.out = 81)
y_longer_cub <- predict(spline_toy_cub, newdata = data.frame(x = x_longer))
y_longer_nat <- predict(spline_toy_nat, newdata = data.frame(x = x_longer))

df_longer <- data.frame(
  x = c(x_longer, x_longer), 
  spline = c(rep("Cubic", 81), rep("Natural", 81)), 
  Yhat = c(y_longer_cub, y_longer_nat)
)

# plot outside of bounds
qbounds <- ggplot(df_longer) + 
  geom_line(aes(x, Yhat), color = "deepskyblue3") +
  geom_function(fun = function(x) exp(x/10) + 2*sin(x/2), 
                linetype = "dashed", color = "darkorange") +
  geom_vline(xintercept = c(0, 25), linetype = "dotted") +
  facet_wrap(~spline) 
qbounds

# save plot
ggsave("index/figures/ch3_toybounds.png", plot = qbounds, device = "png", 
       width = 7, height = 3)
