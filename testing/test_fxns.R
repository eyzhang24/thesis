library(tidyverse)
library(bkmr)
library(NLinteraction)

df <- read_csv("testing/am2_201df.csv")

bk <- read_rds("testing/bkmr_sm_am2_201.RDS")

fit <- bk
z1 <- 5
z2 <- 4
qs.diff = c(0.25, 0.75)
qs.fixed = c(0.25, 0.75)
q.rest = 0.5


bivarinter <- function(fit, z1, z2, qs.diff = c(0.25, 0.75), qs.fixed = c(0.25, 0.75), 
                       q.rest = 0.5) {
  
  # extract fit
  y <- fit$y
  Z <- fit$Z
  X <- fit$X
  
  # fix all chems at q.rest
  point2 <- point1 <- apply(Z, 2, quantile, q.rest)
  # fix z2 at lower
  point2[z2] <- point1[z2] <- quantile(Z[, z2], qs.fixed[1])
  # fix z1 at lower and upper
  point2[z1] <- quantile(Z[, z1], qs.diff[2])
  point1[z1] <- quantile(Z[, z1], qs.diff[1])
  newz.q1 <- rbind(point1, point2) # has all lower quantiles of z2
  # fix all chems at q.rest
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
  
  return(data.frame(z1 = colnames(Z)[z1], z2 = colnames(Z)[z2], 
                    est = ints["est"], sd = ints["sd"], row.names = NULL))
}

bs <- read_rds("testing/bsr_sm_am2_201df3.RDS")

NLmod <- bs
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
j1 <- 5
j2 <- 4
gridLength <- 50
quantile_j2 <- 0.5
quantile_rest <- 0.5

getbivarsurf <- function(NLmod, X, C, j1, j2, gridLength = 50, quantile_j2, 
                         quantile_rest = 0.5) {
  
  # define parameters
  n  <-  dim(X)[1]
  ns <- NLmod$ns
  k <- NLmod$k
  p <- dim(X)[2]
  Xstar <- array(NA, dim = c(n, p, ns + 1))
  Xstar[, , 1] <- 1
  for (j in 1:p) {
    Xstar[, j, 2:(ns + 1)] <- scale(splines::ns(X[, j], df = ns))
  }
  
  # define posteriors
  zetaPost <- NLmod$posterior$zeta
  betaList <- NLmod$posterior$beta
  betaCPost <- NLmod$posterior$betaC
  totalScans <- dim(NLmod$posterior$betaC)[2]
  nChains <- dim(NLmod$posterior$betaC)[1]
  
  # create design of covariates
  pc <- dim(C)[2]
  NewDesignC <- matrix(NA, gridLength, pc + 1)
  NewDesignC[, 1] <- 1
  for (jc in 1:pc) {
    NewDesignC[, jc + 1] <- mean(C[, jc])
  }
  
  # for each quantile of j2
  df <- purrr::map_df(quantile_j2, \(quantile_j2) {
    # create design of chemicals
    n <- dim(X)[1]
    NewDesignMat <- matrix(NA, gridLength, p)
    for (j in 1:p) {
      NewDesignMat[, j] <- quantile(X[, j], quantile_rest)
    }
    NewDesignMat[, j1] <- seq(quantile(X[, j1], 0.025), 
                              quantile(X[, j1], 0.975), length = gridLength)
    NewDesignMat[, j2] <- quantile(X[, j2], quantile_j2)
    NewDesign <- array(NA, dim = c(gridLength, p, ns + 1))
    NewDesign[, , 1] <- 1
    for (j in 1:p) {
      temp_ns_object <- splines::ns(X[, j], df = ns)
      temp_sds <- apply(temp_ns_object, 2, sd)
      temp_means <- apply(temp_ns_object, 2, mean)
      NewDesign[, j, 2:(ns + 1)] <- t((t(predict(temp_ns_object, 
                                                 NewDesignMat[, j])) - temp_means)/temp_sds)
    }
    
    # generate predictions
    predictions <- NLinteraction:::PredictionsMixture(
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