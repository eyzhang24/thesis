library(bkmr)

#generate data
set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X

#look at data
z1 <- seq(min(dat$Z[, 1]), max(dat$Z[, 1]), length = 20)
z2 <- seq(min(dat$Z[, 2]), max(dat$Z[, 2]), length = 20)
hgrid.true <- outer(z1, z2, function(x,y) apply(cbind(x,y), 1, dat$HFun))

res <- persp(z1, z2, hgrid.true, theta = 30, phi = 20, expand = 0.5, 
             col = "lightblue", xlab = "", ylab = "", zlab = "")

#fit bkmr
set.seed(111)
fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 10000, verbose = FALSE, varsel = TRUE)

#investigate model convergence
TracePlot(fit = fitkm, par = "beta")
TracePlot(fit = fitkm, par = "sigsq.eps")
TracePlot(fit = fitkm, par = "r", comp = 1)

#PIPs
ExtractPIPs(fitkm)

#estimating h
med_vals <- apply(Z, 2, median)
Znew <- matrix(med_vals, nrow = 1)
h_true <- dat$HFun(Znew)

h_est1 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "approx")
h_est2 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "exact")
set.seed(111)
samps3 <- SamplePred(fitkm, Znew = Znew, Xnew = cbind(0))

h_est_compare <- data.frame(
  method = c("truth", 1:3),
  post_mean = c(h_true, h_est1$postmean, h_est2$postmean, mean(samps3)),
  post_sd = c(NA, sqrt(h_est1$postvar), sqrt(h_est2$postvar), sd(samps3))
)
h_est_compare

#summarize model output
pred.resp.univar <- PredictorResponseUnivar(fit = fitkm)
library(ggplot2)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")

#tuning parameters
data.frame(fitkm$control.params)

#look at prior fits of r
priorfits <- InvestigatePrior(y = y, Z = Z, X = X, 
                              q.seq = c(2, 1/2, 1/4, 1/16))

PlotPriorFits(y = y, Z = Z, X = X, 
              fits = priorfits)
