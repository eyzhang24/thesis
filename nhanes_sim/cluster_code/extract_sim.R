library(tidyverse)
library(bkmr)
library(NLinteraction)

# set theme for plots
theme_set(theme_light())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_update(
  strip.background = element_rect(color="gray", fill="white"), 
  strip.text = element_text(color = "gray30")
)

########
# test for bkmr
#########

mod1 <- read_rds("_rslurm_bkmr_sm03/mods/bkmr_sm_am2_201.RDS")

#look at trace plots
TracePlot(mod1, par = "r")

#extract PIPs
#returns a dataframe with exposures + PIP's
pip1 <- ExtractPIPs(mod1)
pip1

#extract univariate relationships
#returns data frame with grid of univariate values, estimates, se's
pred.resp.univar <- PredictorResponseUnivar(mod1)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")

#extract bivariate relationships
#returns data frame with grid of bivariate values, estimates, se's
#full bivariate takes a long time to run
# pred.resp.bivar <- PredictorResponseBivar(mod1)
colnames(mod1$Z)
pred.resp.bivar <- PredictorResponseBivar(mod1, 
                                          z.pairs = data.frame(
                                            variable1 = 4, #Hg index
                                            variable2 = 5  #Ni index
                                          ))

#raster
pred.resp.bivar |> 
  ggplot(aes(z1, z2, fill = est)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = c("deepskyblue3", "white", "darkorange"),
    # low = "deepskyblue3", mid = "white", high = "darkorange", 
    na.value = NA) 

#compare relationship at various quantiles
pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar, 
  Z = mod1$Z, qs = c(0.1, 0.5, 0.9))

pred.resp.bivar.levels |> 
  ggplot(aes(z1, est)) + 
  facet_wrap(~variable1) +
  geom_smooth(aes(col = quantile), stat = "identity")

# try interaction estimate
risks.int <- SingVarIntSummaries(mod1, 
                                 which.z = c(4, 5), 
                                 qs.diff = c(0.25, 0.75), 
                                 qs.fixed = c(0.25, 0.75),
                                 method = "approx")

ggplot(risks.int, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip() +
  geom_hline(aes(yintercept = 0), linetype = "dashed")

#######
# try bivariate interaction function
#######
colnames(mod1$Z)
bivarinter <- function(fit, z1, z2, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), 
                       q.rest = 0.5) {
  #' fit = bkmr model
  #' z1 = index of chemical to estimate diff in response for
  #' z2 = index of chemical to fix at diff quantiles
  #' qs.diff = quantiles to calculate response diff for
  #' qs.fixed = quantiles to fix other chemical at
  #' q.rest = quantile to fix rest of chemicals at
  
  # extract fit
  y <- fit$y
  Z <- fit$Z
  X <- fit$X
  
  # fix z2, test z1
  
  # fix all chems at 0.5
  point2 <- point1 <- apply(Z, 2, quantile, q.rest)
  # fix z2 at lower
  point2[z2] <- point1[z2] <- quantile(Z[, z2], qs.fixed[1])
  # fix z1 at lower and upper
  point2[z1] <- quantile(Z[, z1], qs.diff[2])
  point1[z1] <- quantile(Z[, z1], qs.diff[1])
  newz.q1 <- rbind(point1, point2) # has all lower quantiles of z2
  # fix all chems at 0.5
  point2 <- point1 <- apply(Z, 2, quantile, q.rest)
  # fix z2 at higher
  point2[z2] <- point1[z2] <- quantile(Z[, z2], qs.fixed[2])
  # fix z1 at lower and upper
  point2[z1] <- quantile(Z[, z1], qs.diff[2])
  point1[z1] <- quantile(Z[, z1], qs.diff[1])
  newz.q2 <- rbind(point1, point2) # has all upper quantiles of z2
  
  # prepare
  cc <- c(-1 * c(-1, 1), c(-1, 1))
  newz <- rbind(newz.q1, newz.q2)
  
  # default to using approximate calculation
  preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz)
  
  # extract intervals
  int <- drop(cc %*% preds$postmean)
  int.se <- drop(sqrt(cc %*% preds$postvar %*% cc))
  ints <- c(est = int, sd = int.se)
  
  df <- dplyr::tibble(variable = colnames(Z)[z1], fixedat = colnames(Z)[z2], 
                      est = ints["est"], sd = ints["sd"])
  return(df)
}

# test
bivar.test <- bivarinter(mod1, 5, 4, 
                         qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), 
                         q.rest = 0.5)
bivar.test
ggplot(bivar.test, aes(variable, est, ymin = est - 1.96*sd, 
                      ymax = est + 1.96*sd)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip() +
  geom_hline(aes(yintercept = 0), linetype = "dashed")

#SingVarIntSummaries
function (fit, y = NULL, Z = NULL, X = NULL, which.z = 1:ncol(Z), 
          qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), method = "approx", 
          sel = NULL, z.names = colnames(Z), ...) 
{
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) 
      y <- fit$y
    if (is.null(Z)) 
      Z <- fit$Z
    if (is.null(X)) 
      X <- fit$X
  }
  if (is.null(z.names)) 
    z.names <- paste0("z", 1:ncol(Z))
  ints <- sapply(which.z, function(whichz) SingVarIntSummary(whichz = whichz, 
                                                             fit = fit, Z = Z, X = X, y = y, 
                                                             qs.diff = qs.diff, qs.fixed = qs.fixed, 
                                                             method, sel = sel, ...))
  df <- dplyr::tibble(variable = factor(z.names[which.z], levels = z.names), 
                      est = ints["est", ], sd = ints["sd", ])
}

#SingVarIntSummary 
function (whichz = 1, fit, y = NULL, Z = NULL, X = NULL, 
          qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), 
          method = "approx", sel = NULL, 
          ...) 
{
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) 
      y <- fit$y
    if (is.null(Z)) 
      Z <- fit$Z
    if (is.null(X)) 
      X <- fit$X
  }
  q.fixed <- qs.fixed[1]
  point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
  point2[whichz] <- quantile(Z[, whichz], qs.diff[2])
  point1[whichz] <- quantile(Z[, whichz], qs.diff[1])
  newz.q1 <- rbind(point1, point2) #has all lower quantiles
  q.fixed <- qs.fixed[2]
  point2 <- point1 <- apply(Z, 2, quantile, q.fixed)
  point2[whichz] <- quantile(Z[, whichz], qs.diff[2])
  point1[whichz] <- quantile(Z[, whichz], qs.diff[1])
  newz.q2 <- rbind(point1, point2) #has all upper quantiles
  if (method %in% c("approx", "exact")) {
    preds.fun <- function(znew) {
      ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = znew, 
                          sel = sel, method = method) 
    }
    interactionSummary <- interactionSummary.approx
  }
  else {
    stop("method must be one of c('approx', 'exact')")
  }
  interactionSummary(newz.q1, newz.q2, preds.fun, ...)
}

#interactionSummary.approx
function (newz.q1, newz.q2, preds.fun, ...) 
{
  cc <- c(-1 * c(-1, 1), c(-1, 1))
  newz <- rbind(newz.q1, newz.q2)
  preds <- preds.fun(newz, ...)
  int <- drop(cc %*% preds$postmean)
  int.se <- drop(sqrt(cc %*% preds$postvar %*% cc))
  c(est = int, sd = int.se)
}
#interactionSummary.samp
function (newz.q1, newz.q2, preds.fun, ...) 
{
  cc <- c(-1 * c(-1, 1), c(-1, 1))
  newz <- rbind(newz.q1, newz.q2)
  preds <- preds.fun(newz, ...)
  int.preds <- drop(preds %*% cc)
  c(est = mean(int.preds), sd = sd(int.preds))
}



########
# test for bsr
########

cnames <- c("As", "Cd", "Co", "Hg", "Ni", "Tl", "Pb", "Mo", "Sb", "Sn")

mod2 <- read_rds("_rslurm_bsr_sm15/mods/bsr_sm_dm2_1401df4.RDS")
mod2$MainPIP
mod2$InteractionPIP
# plotInt(mod2)
InteractionProb(mod2, Xsub = c(4, 5, 6))

data.frame(variable = cnames, PIP = mod2$MainPIP)

out1_resp1 <- read_rds("sim/sim_resp_sm_a.RDS")

# prepare data
df <- out1_resp1[[1401]]
X <- df |> 
  select(As:Sn) |> 
  as.matrix.data.frame()
C <- df |>
  bind_cols(
    data.frame(model.matrix(~ race-1, data = 
                              mutate(df, race = as.factor(race))))
  ) |> 
  select(race2:race5, smoke:bmi) |> 
  as.matrix.data.frame()
Y <- df$y

source("extract_fxns.R")

test <- trivarsurf_bsr(mod2, X, C, j1 = 4, j2 = 5, j3 = 6, gridLength = 50, 
                           quantile_j23 = c(0.1, 0.5, 0.9), quantile_rest = 0.5)

plotSurface1d(mod2, X = X, C = C, j1 = 4, j2 = 5, 
                gridLength = 50, quantile_j2 = 0.5, quantile_rest = 0.5)


NLmod <- mod2
j1 <- 4
j2 <- 5
gridLength <- 50
quantile_j2 <- 0.5
quantile_rest <- 0.5

getbivarsurf <- function(NLmod, X, C, j1, j2, gridLength = 50, quantile_j2, 
                         quantile_rest = 0.5) {
  
  # define parameters
  n = dim(X)[1]
  ns = NLmod$ns
  k = NLmod$k
  p = dim(X)[2]
  Xstar = array(NA, dim = c(n, p, ns + 1))
  Xstar[, , 1] = 1
  for (j in 1:p) {
    Xstar[, j, 2:(ns + 1)] = scale(splines::ns(X[, j], df = ns))
  }
  
  # define posteriors
  zetaPost = NLmod$posterior$zeta
  betaList = NLmod$posterior$beta
  betaCPost = NLmod$posterior$betaC
  totalScans = dim(NLmod$posterior$betaC)[2]
  nChains = dim(NLmod$posterior$betaC)[1]
  
  # create design of covariates
  pc = dim(C)[2]
  NewDesignC = matrix(NA, gridLength, pc + 1)
  NewDesignC[, 1] = 1
  for (jc in 1:pc) {
    NewDesignC[, jc + 1] = mean(C[, jc])
  }

  # create design of chemicals
  df <- purrr::map_df(quantile_j2, \(quantile_j2) {
    n = dim(X)[1]
    NewDesignMat = matrix(NA, gridLength, p)
    for (j in 1:p) {
      NewDesignMat[, j] = quantile(X[, j], quantile_rest)
    }
    NewDesignMat[, j1] = seq(quantile(X[, j1], 0.025), 
                             quantile(X[, j1], 0.975), length = gridLength)
    NewDesignMat[, j2] = quantile(X[, j2], quantile_j2)
    NewDesign = array(NA, dim = c(gridLength, p, ns + 1))
    NewDesign[, , 1] = 1
    for (j in 1:p) {
      temp_ns_object = splines::ns(X[, j], df = ns)
      temp_sds = apply(temp_ns_object, 2, sd)
      temp_means = apply(temp_ns_object, 2, mean)
      NewDesign[, j, 2:(ns + 1)] = t((t(predict(temp_ns_object, 
                                                NewDesignMat[, j])) - temp_means)/temp_sds)
    }
    
    # generate predictions
    predictions = NLinteraction:::PredictionsMixture(
      XstarOld = Xstar, XstarNew = NewDesign, 
      designC = NewDesignC, totalScans = totalScans, nChains = nChains, 
      zetaPost = zetaPost, betaList = betaList, betaCPost = betaCPost, 
      k = k, ns = ns)
    
    # get surface
    return(data.frame(
      j1val = NewDesignMat[, j1], 
      j2quant = rep(quantile_j2, gridLength), 
      est = apply(predictions$PredictedPost, 3, mean), 
      lower = apply(predictions$PredictedPost, 3, quantile, 0.025), 
      upper = apply(predictions$PredictedPost, 3, quantile, 0.975)
    ))
  })
  return(df)
}

test <- getbivarsurf(mod2, X = X, C = C, j1 = 4, j2 = 5, 
                      gridLength = 50, quantile_j2 = c(0.1, 0.5, 0.9), 
                      quantile_rest = 0.5)

test |> 
  filter(j23quant == 0.5) |> 
  ggplot(aes(x = j1val)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80") +
  geom_line(aes(y = est))

test |> 
  ggplot(aes(x = j1val, y = est, color = as.factor(j23quant))) +
  geom_line()

