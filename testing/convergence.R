library(bkmr)
library(NLinteraction)

# bkmr fits
bkmr_sm_am2_203 <- readRDS("testing/bkmr_sm_am2_203.RDS")
bkmr_lg_am2_284 <- readRDS("testing/bkmr_lg_am2_284.RDS")

# bsr fits
bsr_sm_am2_203df2 <- readRDS("testing/bsr_sm_am2_203df2.RDS")
bsr_lgf_am2_284df2 <- readRDS("testing/bsr_lgf_am2_284df2.RDS")

png("index/figures/traceplots/bksm_traceplot.png", width = 8, height = 6, units = "in", res = 360)
par(mfrow = c(2, 2))
TracePlot(fit = bkmr_sm_am2_203, par = "beta") # prior probability
TracePlot(fit = bkmr_sm_am2_203, par = "sigsq.eps") # variance of residuals
TracePlot(fit = bkmr_sm_am2_203, par = "r", comp = 1) # prob of each
TracePlot(fit = bkmr_sm_am2_203, par = "r", comp = 5)
dev.off()

# TracePlot(fit = bkmr_sm_am2_203, par = "lambda") # vertical scale 

png("index/figures/traceplots/bklg_traceplot.png", width = 8, height = 6, units = "in", res = 360)
par(mfrow = c(2, 2))
TracePlot(fit = bkmr_lg_am2_284, par = "beta")
TracePlot(fit = bkmr_lg_am2_284, par = "sigsq.eps")
TracePlot(fit = bkmr_lg_am2_284, par = "r", comp = 1)
TracePlot(fit = bkmr_lg_am2_284, par = "r", comp = 5)
dev.off()

# tau = prior probability of inclusion 
# sigma = variance of residuals
# beta equiv to rho

# bsr sm
h <- bsr_sm_am2_203df2$posterior

png("index/figures/traceplots/bssm_traceplot.png", width = 8, height = 6, units = "in", res = 360)
par(mfrow = c(2, 2))
htau <- t(h[["tau"]][, , 1])
plot(htau[, 1], type = "l", 
     main = paste0("(tau = ", round(mean(htau[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(htau[, 1]), col = "blue", lwd = 2)

hsigma <- t(h[["sigma"]])
plot(hsigma[, 1], type = "l", 
     main = paste0("(sigma = ", round(mean(hsigma[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(hsigma[, 1]), col = "blue", lwd = 2)

hzeta1 <- t(h[["zeta"]][, , 1, 2])
plot(hzeta1[, 1], type = "l", 
     main = paste0("(zeta1 = ", round(mean(hzeta1[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(hzeta1[, 1]), col = "blue", lwd = 2)

hzeta5 <- t(h[["zeta"]][, , 5, 2])
plot(hzeta5[, 1], type = "l", 
     main = paste0("(zeta5 = ", round(mean(hzeta5[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(hzeta5[, 1]), col = "blue", lwd = 2)
dev.off()

# bsr lg
h <- bsr_lgf_am2_284df2$posterior

png("index/figures/traceplots/bslg_traceplot.png", width = 8, height = 6, units = "in", res = 360)
par(mfrow = c(2, 2))
htau <- t(h[["tau"]][, , 1])
plot(htau[, 1], type = "l", 
     main = paste0("(tau = ", round(mean(htau[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(htau[, 1]), col = "blue", lwd = 2)

hsigma <- t(h[["sigma"]])
plot(hsigma[, 1], type = "l", 
     main = paste0("(sigma = ", round(mean(hsigma[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(hsigma[, 1]), col = "blue", lwd = 2)

hzeta1 <- t(h[["zeta"]][, , 1, 2])
plot(hzeta1[, 1], type = "l", 
     main = paste0("(zeta1 = ", round(mean(hzeta1[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(hzeta1[, 1]), col = "blue", lwd = 2)

hzeta5 <- t(h[["zeta"]][, , 5, 2])
plot(hzeta5[, 1], type = "l", 
     main = paste0("(zeta5 = ", round(mean(hzeta5[,1]), 2), ")"), 
     xlab = "scan", ylab = "parameter value")
abline(h = mean(hzeta5[, 1]), col = "blue", lwd = 2)
dev.off()

# put all together

library(magick)

# list of input files
file_list <- c("index/figures/traceplots/bksm_traceplot.png", 
               "index/figures/traceplots/bklg_traceplot.png", 
               "index/figures/traceplots/bssm_traceplot.png", 
               "index/figures/traceplots/bslg_traceplot.png")

# output file name
output_file <- "index/figures/traceplots/bksm_traceplotmerged.png"

# merge images
merge_and_save <- function(file_list, output_file) {
  images <- image_read(file_list)
  
  # two columns
  top_row <- image_append(images[1:2], stack = TRUE)
  bottom_row <- image_append(images[3:4], stack = TRUE)
  
  # combine columns
  combined_image <- image_append(c(top_row, bottom_row), stack = FALSE)
  
  # add labels
  combined_image <- image_annotate(combined_image, "a", location = "+30+20", 
                                   size = 120, color = "black", font = "Roboto", weight = 700)
  combined_image <- image_annotate(combined_image, "b", location = "+2900+20", 
                                   size = 120, color = "black", font = "Roboto", weight = 700)
  combined_image <- image_annotate(combined_image, "c", location = "+30+2110", 
                                   size = 120, color = "black", font = "Roboto", weight = 700)
  combined_image <- image_annotate(combined_image, "d", location = "+2900+2110", 
                                   size = 120, color = "black", font = "Roboto", weight = 700)
  
  image_write(combined_image, path = output_file, format = "png")
}

merge_and_save(file_list, output_file)
