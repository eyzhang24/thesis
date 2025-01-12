# example with output from a BSR model

# packages
# library(tidyverse)
library(NLinteraction)

# read in a bsr model
# fit on n=252 data w/ 10 chems, 8 covariates
# higher effect, mult. interxn between Hg (index 4) and Ni (index 5)
if(basename(getwd()) != "thesis") setwd(path.expand("~/repo/thesis")) # added by AK

bsrmod <- readRDS("testing/bsr_posteriors/bsr_sm_am2_201df3.RDS")

# in order to get predictions, we need the original data
# read data
df <- readRDS("testing/bsr_posteriors/data_sm_am2_201.RDS")
# re-format
X <- df |> 
  dplyr::select(As:Sn) |> 
  as.matrix.data.frame()
C <- df |>
  cbind(
    data.frame(model.matrix(~ race-1, data = 
                              dplyr::mutate(df, race = as.factor(race))))
  ) |> 
  dplyr::select(race2:race5, smoke:bmi) |> 
  as.matrix.data.frame()
Y <- df$y

# here's code for a larger size (n=1000) dataset
# higher effect, mult. interxn between Hg and Ni

# bsrmod <- readRDS("testing/bsr_posteriors/bsr_lgf_am2_201df3.RDS")
# df <- readRDS("testing/bsr_posteriors/data_lg_am2_201.RDS")
# X <- df |>
#   dplyr::select(As:Sn) |>
#   as.matrix.data.frame()
# C <- df |>
#   cbind(
#     data.frame(model.matrix(~ race-1, data =
#                               dplyr::mutate(df, race = as.factor(race))))
#   ) |>
#   dplyr::select(race2:race5, smoke:bmi) |>
#   as.matrix.data.frame()
# Y <- df$y

# and also code for a larger size (n=1000) stratified model
# effect of Hg is doubled in race cat 5
# listmod <- readRDS("testing/bsr_posteriors/bsr_lg_ep2_301df2.RDS")
# bsrmod <- listmod[[5]] # choose race category 5
# dffull <- readRDS("testing/bsr_posteriors/data_lg_ep2_301.RDS")
# df <- dffull[dffull$race == 5, ] # choose race category 5
# X <- df |>
#   dplyr::select(As:Sn) |>
#   as.matrix.data.frame()
# C <- df |>
#   dplyr::select(smoke:bmi) |>
#   as.matrix.data.frame()
# Y <- df$y

#####################
# built in functions ------------------------------------------------------
#####################

# pip's
univar <- bsrmod$MainPIP
bivar <- bsrmod$InteractionPIP

# access the posteriors
posteriors <- bsrmod$posterior
  # sigma = variance of residuals
  # tau = prior probability of inclusion (i think?)
  # zeta = indicator for inclusion, PIP (i think?)
  # beta = coefficients on chems
  # betac = coefficients on cov's

# built in functions to visualize pred-exp line
par(mfrow=c(1,3), pty='s')
plotSurface1d(NLmod = bsrmod, X=X, C=C, j1=4, j2=5,
              gridLength=30, quantile_j2=0.2, quantile_rest = 0.5,
              xlab="Hg (4)", ylab="Posterior predictive distribution",
              main="20th quantile of Ni (5)")
plotSurface1d(NLmod = bsrmod, X=X, C=C, j1=4, j2=5,
              gridLength=30, quantile_j2=0.5, quantile_rest = 0.5,
              xlab="Hg (4)", ylab="Posterior predictive distribution",
              main="50th quantile of Ni (5)")
plotSurface1d(NLmod = bsrmod, X=X, C=C, j1=4, j2=5,
              gridLength=30, quantile_j2=0.8, quantile_rest = 0.5,
              xlab="Hg (4)", ylab="Posterior predictive distribution",
              main="80th quantile of Ni (5)")

# built in functions to visualize bivariate pred-exp surface
# takes a second to run
par(mfrow=c(1,2), pty='s')
plotSurface2dMean(NLmod = bsrmod, X=X, C=C, j1=4, j2=5,
                  gridLength_j1=20, gridLength_j2 = 20,
                  quantile_rest = 0.5, xlab='Hg (4)', ylab='Ni (5)',
                  main="Posterior mean")
points(X[,4], X[,5], pch=16, cex=0.5)
plotSurface2dSD(NLmod = bsrmod, X=X, C=C, j1=4, j2=5,
                gridLength_j1=20, gridLength_j2 = 20,
                quantile_rest = 0.5, xlab='Hg (4)', ylab='Ni (5)',
                main="Posterior SD")
points(X[,4], X[,5], pch=16, cex=0.5)

#####################
# getting predictions -----------------------------------------------------
#####################

# built in fxns don't give access to the estimates themselves
# i wrote these two functions so i could make some of the plots in my thesis
# based on NLinteraction:::plotSurface1d

# estimate univar exp-resp relationships
univarsurf_bsr <- function(NLmod, X, C, j1, gridLength = 50, 
                           quantile_rest = 0.5) {
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' gridLength = number of points to estimate response at
  #' j1 = index of first chemical
  #' j2 = index of second chemical
  #' quantile_j2 = vector of quantiles to fix second chemical at
  #' quantile_rest = quantile to fix other chemicals at
  #' 
  #' NOTE order of j1 and j2 doesn't matter
  
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
  
  # create design of chemicals
  n <- dim(X)[1]
  NewDesignMat <- matrix(NA, gridLength, p)
  for (j in 1:p) {
    NewDesignMat[, j] <- quantile(X[, j], quantile_rest)
  }
  NewDesignMat[, j1] <- seq(quantile(X[, j1], 0.025),
                            quantile(X[, j1], 0.975), length = gridLength)
  NewDesign <- array(NA, dim = c(gridLength, p, ns + 1))
  NewDesign[, , 1] <- 1
  for (j in 1:p) {
    temp_ns_object <- splines::ns(X[, j], df = ns)
    temp_sds <- apply(temp_ns_object, 2, sd)
    temp_means <- apply(temp_ns_object, 2, mean)
    NewDesign[, j, 2:(ns + 1)] <- t((t(
      predict(temp_ns_object,
              NewDesignMat[, j])
    ) - temp_means) / temp_sds)
  }
  
  # generate predictions
  predictions <- NLinteraction:::PredictionsMixture(
    XstarOld = Xstar,
    XstarNew = NewDesign,
    designC = NewDesignC,
    totalScans = totalScans,
    nChains = nChains,
    zetaPost = zetaPost,
    betaList = betaList,
    betaCPost = betaCPost,
    k = k,
    ns = ns
  )
  
  return(data.frame(
    j1val = NewDesignMat[, j1], 
    est = apply(predictions$PredictedPost, 3, mean), 
    lower = apply(predictions$PredictedPost, 3, quantile, 0.025), 
    upper = apply(predictions$PredictedPost, 3, quantile, 0.975)))
}

# quick example
univarex <- univarsurf_bsr(NLmod = bsrmod, X = X, C = C, 
                           j1 = 4, gridLength = 50, quantile_rest = 0.5)

# estimate bivar exp-resp relationships, fixing at diff quantiles of 2nd chem
bivarsurf_bsr <- function(NLmod, X, C, j1, j2, gridLength = 50, 
                          quantile_j2 = c(0.1, 0.5, 0.9), quantile_rest = 0.5) {
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' gridLength = number of points to estimate response at
  #' j1 = index of first chemical
  #' j2 = index of second chemical
  #' quantile_j2 = vector of quantiles to fix second chemical at
  #' quantile_rest = quantile to fix other chemicals at
  #' 
  #' NOTE order of j1 and j2 does matter
  
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

# quick example
bivarex <- bivarsurf_bsr(NLmod = bsrmod, X = X, C = C, 
                         j1 = 4, j2 = 5, gridLength = 50, 
                         quantile_j2 = c(0.1, 0.5, 0.9), quantile_rest = 0.5)



#### quick examples by AK
# RE: statement in the thesis
# "We were not able to find a method in the literature for combining variances from
#   estimated responses at two different sets of predictor values in the spline regression
#   framework."

# examples below: 
#  1) estimating joint effect of all exposures 
#  2) estimating independent effect of a single exposure at a fixed quantile of all other exposures
#  3) estimating "stratified" effect of one exposure at different quantiles of 
#     another exposure and a fixed quantile of all other exposures


# fit
bsrmod <- readRDS("testing/bsr_posteriors/bsr_sm_am2_201df3.RDS")
# original data
df <- readRDS("testing/bsr_posteriors/data_sm_am2_201.RDS")

NLmod <- bsrmod


X <- df |> 
  dplyr::select(As:Sn) |> 
  as.matrix.data.frame()
C <- df |>
  cbind(
    data.frame(model.matrix(~ race-1, data = 
                              dplyr::mutate(df, race = as.factor(race))))
  ) |> 
  dplyr::select(race2:race5, smoke:bmi) |> 
  as.matrix.data.frame()
Y <- df$y


get_joint_effect_estimate <- function(
    NLmod, X, C, ref_quantile=0.5, index_quantiles=setdiff(seq(0.1,.9,.1), .5)){
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' ref_quantile = quantile to fix chemical of interest as reference
  #' index_quantiles = quantile to fix chemical of interest for comparisons (in calculations of mean difference)
  
  
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
  
  # create design of covariates (arithmetic mean, for simplicity)
  pc <- dim(C)[2]
  NewDesignC <- matrix(NA, length(index_quantiles)+1, pc + 1)
  NewDesignC[, 1] <- 1
  for (jc in 1:pc) {
    NewDesignC[, jc + 1] <- mean(C[, jc])
  }
  
  # create design of chemicals
  allquantiles  = c(ref_quantile, index_quantiles)
  n <- dim(X)[1]

  NewDesignMat <- apply(X,2, quantile, allquantiles)
  
  NewDesign <- array(NA, dim = c(length(index_quantiles)+1, p, ns + 1))
  NewDesign[, , 1] <- 1
  for (j in 1:p) {
    temp_ns_object <- splines::ns(X[, j], df = ns)
    temp_sds <- apply(temp_ns_object, 2, sd)
    temp_means <- apply(temp_ns_object, 2, mean)
    NewDesign[, j, 2:(ns + 1)] <- t((t(                            
      predict(temp_ns_object,
              NewDesignMat[, j])
    ) - temp_means) / temp_sds)
  }
  
  
  

  
  # generate predictions
  predictions <- NLinteraction:::PredictionsMixture(
    XstarOld = Xstar,
    XstarNew = NewDesign,
    designC = NewDesignC,
    totalScans = totalScans,
    nChains = nChains,
    zetaPost = zetaPost,
    betaList = betaList,
    betaCPost = betaCPost,
    k = k,
    ns = ns
  )
  
  # nchains X total scans/iterations X # test data points
  refpreds = as.numeric(predictions$PredictedPost[,,1])
  idxpreds = lapply(1:length(index_quantiles), function(x) as.numeric(predictions$PredictedPost[,,x+1]))
  pred_difference = lapply(1:length(index_quantiles), function(x) idxpreds[[x]]-refpreds)
  
  return(data.frame(
    #j1val = NewDesignMat[, j1], 
    quantile = c(ref_quantile, index_quantiles), 
    est = c(0, as.numeric(lapply(pred_difference, mean))), 
    lower =  c(0, as.numeric(lapply(pred_difference, quantile, 0.025))), 
    upper =  c(0, as.numeric(lapply(pred_difference, quantile, 0.975))))
    )
}

predest = get_joint_effect_estimate(bsrmod, X, C, ref_quantile=0.5, index_quantiles=setdiff(seq(0.05,.9,.05), .5))

# this will look similar to some BKMR plots, and effect estimates are held in the data frame
library(ggplot2)
ggplot(data=predest, aes(x=quantile, y=est)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) + 
  geom_line() + 
  labs(y="Mean difference in Y", x="Joint quantile of all exposures") + theme_bw()



get_independent_effect_estimate <- function(
    NLmod, X, C, j1, quantile_rest=0.5, ref_quantile=0.5, 
    index_quantiles=setdiff(seq(0.1,.9,.1), .5)){
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' j1 = index of  chemical of interest
  #' quantile_rest = quantile to fix other chemicals at
  #' ref_quantile = quantile to fix chemical of interest as reference
  #' index_quantiles = quantile to fix chemical of interest for comparisons (in calculations of mean difference)
  
  
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
  
  # create design of covariates (arithmetic mean, for simplicity)
  pc <- dim(C)[2]
  NewDesignC <- matrix(NA, length(index_quantiles)+1, pc + 1)
  NewDesignC[, 1] <- 1
  for (jc in 1:pc) {
    NewDesignC[, jc + 1] <- mean(C[, jc])
  }
  
  # create design of chemicals
  allquantiles  = c(ref_quantile, index_quantiles)
  n <- dim(X)[1]
  
  NewDesignMat <- apply(X,2, function(x) rep(quantile(x, quantile_rest), length(allquantiles)))
  indexq = quantile(X[,j1], allquantiles)
  NewDesignMat[,j1] <- indexq
  rownames(NewDesignMat) <- names(indexq)
  
  NewDesign <- array(NA, dim = c(length(index_quantiles)+1, p, ns + 1))
  NewDesign[, , 1] <- 1
  for (j in 1:p) {
    temp_ns_object <- splines::ns(X[, j], df = ns)
    temp_sds <- apply(temp_ns_object, 2, sd)
    temp_means <- apply(temp_ns_object, 2, mean)
    NewDesign[, j, 2:(ns + 1)] <- t((t(                            
      predict(temp_ns_object,
              NewDesignMat[, j])
    ) - temp_means) / temp_sds)
  }
  
  
  # generate predictions
  predictions <- NLinteraction:::PredictionsMixture(
    XstarOld = Xstar,
    XstarNew = NewDesign,
    designC = NewDesignC,
    totalScans = totalScans,
    nChains = nChains,
    zetaPost = zetaPost,
    betaList = betaList,
    betaCPost = betaCPost,
    k = k,
    ns = ns
  )
  
  # nchains X total scans/iterations X # test data points
  refpreds = as.numeric(predictions$PredictedPost[,,1])
  idxpreds = lapply(1:length(index_quantiles), function(x) as.numeric(predictions$PredictedPost[,,x+1]))
  pred_difference = lapply(1:length(index_quantiles), function(x) idxpreds[[x]]-refpreds)
  
  return(data.frame(
    #j1val = NewDesignMat[, j1], 
    quantile = c(ref_quantile, index_quantiles), 
    est = c(0, as.numeric(lapply(pred_difference, mean))), 
    lower =  c(0, as.numeric(lapply(pred_difference, quantile, 0.025))), 
    upper =  c(0, as.numeric(lapply(pred_difference, quantile, 0.975))))
  )
}


predest_independent = get_independent_effect_estimate(
  bsrmod, X, C, j1=4, ref_quantile=0.5, index_quantiles=setdiff(seq(0.05,.9,.05), .5))

# this will look similar to some BKMR plots, and effect estimates are held in the data frame
library(ggplot2)
ggplot(data=predest_independent, aes(x=quantile, y=est)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) + 
  geom_line() + 
  labs(y="Mean difference", x=paste("Quantile of", colnames(X)[4])) + theme_bw()



get_commonreferent_effect_estimate <- function(
    NLmod, X, C, j1, j2, quantile_j2 = c(0.5, 0.1, 0.9), quantile_rest=0.5, 
    ref_quantile=0.5, index_quantiles=setdiff(seq(0.1,.9,.1), .5)){
  # common referent estimates for an interaction
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' j1 = index of  chemical of interest
  #' j2 = index of  chemical of at levels of which interactions are assessed
  #' quantile_j2  = levels of j2 over which interaction is considered, 
  #'      the first value will form a common referent, so order is important
  #' quantile_rest = quantile to fix other chemicals at
  #' ref_quantile = quantile to fix chemical of interest as reference
  #' index_quantiles = quantile to fix chemical of interest for comparisons (in calculations of mean difference)
  
  
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
  
  # create design of covariates (arithmetic mean, for simplicity)
  pc <- dim(C)[2]
  NewDesignC <- matrix(NA, (length(index_quantiles)+1)*length(quantile_j2), pc + 1)
  NewDesignC[, 1] <- 1
  for (jc in 1:pc) {
    NewDesignC[, jc + 1] <- mean(C[, jc])
  }
  
  # create design of chemicals
  allquantiles  = c(ref_quantile, index_quantiles)
  n <- dim(X)[1]
  
  NewDesignMat <- apply(X,2, function(x) rep(quantile(x, quantile_rest), length(allquantiles)*length(quantile_j2)))
  indexq = quantile(X[,j1], allquantiles)
  interq = quantile(X[,j2], quantile_j2)
  NewDesignMat[,j1] <- rep(indexq, length(quantile_j2))
  NewDesignMat[,j2] <- rep(interq, each=length(allquantiles))
  rownames(NewDesignMat) <- rep(names(indexq), length(quantile_j2))
  
  NewDesign <- array(NA, dim = c((length(index_quantiles)+1)*length(quantile_j2), p, ns + 1))
  NewDesign[, , 1] <- 1
  for (j in 1:p) {
    temp_ns_object <- splines::ns(X[, j], df = ns)
    temp_sds <- apply(temp_ns_object, 2, sd)
    temp_means <- apply(temp_ns_object, 2, mean)
    NewDesign[, j, 2:(ns + 1)] <- t((t(                            
      predict(temp_ns_object,
              NewDesignMat[, j])
    ) - temp_means) / temp_sds)
  }
  
  
  # generate predictions
  predictions <- NLinteraction:::PredictionsMixture(
    XstarOld = Xstar,
    XstarNew = NewDesign,
    designC = NewDesignC,
    totalScans = totalScans,
    nChains = nChains,
    zetaPost = zetaPost,
    betaList = betaList,
    betaCPost = betaCPost,
    k = k,
    ns = ns
  )
  
  # nchains X total scans/iterations X # test data points
  refpreds = as.numeric(predictions$PredictedPost[,,1])
  idxpreds = lapply(1:(length(allquantiles)*length(quantile_j2)-1), function(x) as.numeric(predictions$PredictedPost[,,x+1]))
  pred_difference = lapply(1:(length(allquantiles)*length(quantile_j2)-1), function(x) idxpreds[[x]]-refpreds)
  
  
  
  
  return(
    data.frame(
    #j1val = NewDesignMat[, j1], 
    quantile = rep(allquantiles, length(quantile_j2)), 
    quantile_j2 = rep(quantile_j2, each=length(allquantiles)), 
    est = c(0, as.numeric(lapply(pred_difference, mean))), 
    lower =  c(0, as.numeric(lapply(pred_difference, quantile, 0.025))), 
    upper =  c(0, as.numeric(lapply(pred_difference, quantile, 0.975))))
  )
}


predest_independent_commonint = get_commonreferent_effect_estimate(
  bsrmod, X, C, j1=4, j2=5, quantile_j2=c(0.1, 0.5, 0.9), ref_quantile=0.5, 
  index_quantiles=setdiff(seq(0.05,.9,.05), .5))

library(ggplot2)
ggplot(data=predest_independent_commonint, aes(x=quantile, y=est, colour=factor(quantile_j2), fill=factor(quantile_j2))) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.5) + # makes it super messy
  geom_line() + 
  labs(y="Mean difference", x=paste("Quantile of", colnames(X)[4])) +
  scale_color_discrete(name=paste("Quantile of", colnames(X)[5]))+
  scale_fill_discrete(name=paste("Quantile of", colnames(X)[5])) + theme_bw()



get_stratified_effect_estimate <- function(
    NLmod, X, C, j1, j2, quantile_j2 = c(0.5, 0.1, 0.9), 
    quantile_rest=0.5, ref_quantile=0.5, 
    index_quantiles=setdiff(seq(0.1,.9,.1), .5)){
  # stratified effect estimates for an interaction (different referent for each level of j2)
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' j1 = index of  chemical of interest
  #' j2 = index of  chemical of at levels of which interactions are assessed
  #' quantile_j2  = levels of j2 over which interaction is considered
  #' quantile_rest = quantile to fix other chemicals at
  #' ref_quantile = quantile to fix chemical of interest as reference
  #' index_quantiles = quantile to fix chemical of interest for comparisons (in calculations of mean difference)
  
  
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
  
  # create design of covariates (arithmetic mean, for simplicity)
  pc <- dim(C)[2]
  NewDesignC <- matrix(NA, (length(index_quantiles)+1)*length(quantile_j2), pc + 1)
  NewDesignC[, 1] <- 1
  for (jc in 1:pc) {
    NewDesignC[, jc + 1] <- mean(C[, jc])
  }
  
  # create design of chemicals
  allquantiles  = c(ref_quantile, index_quantiles)
  n <- dim(X)[1]
  
  NewDesignMat <- apply(X,2, function(x) rep(quantile(x, quantile_rest), length(allquantiles)*length(quantile_j2)))
  indexq = quantile(X[,j1], allquantiles)
  interq = quantile(X[,j2], quantile_j2)
  NewDesignMat[,j1] <- rep(indexq, length(quantile_j2))
  NewDesignMat[,j2] <- rep(interq, each=length(allquantiles))
  rownames(NewDesignMat) <- rep(names(indexq), length(quantile_j2))
  
  NewDesign <- array(NA, dim = c((length(index_quantiles)+1)*length(quantile_j2), p, ns + 1))
  NewDesign[, , 1] <- 1
  for (j in 1:p) {
    temp_ns_object <- splines::ns(X[, j], df = ns)
    temp_sds <- apply(temp_ns_object, 2, sd)
    temp_means <- apply(temp_ns_object, 2, mean)
    NewDesign[, j, 2:(ns + 1)] <- t((t(                            
      predict(temp_ns_object,
              NewDesignMat[, j])
    ) - temp_means) / temp_sds)
  }
  
  
  # generate predictions
  predictions <- NLinteraction:::PredictionsMixture(
    XstarOld = Xstar,
    XstarNew = NewDesign,
    designC = NewDesignC,
    totalScans = totalScans,
    nChains = nChains,
    zetaPost = zetaPost,
    betaList = betaList,
    betaCPost = betaCPost,
    k = k,
    ns = ns
  )
  
  # nchains X total scans/iterations X # test data points
  refidx = seq(1,length(allquantiles)*length(quantile_j2), length(allquantiles))
  allidx = 1:(length(allquantiles)*length(quantile_j2))
  indexidx = setdiff(allidx, refidx)
  pred_difference_list = predictions$PredictedPost*NA # lazy way to create list of NAs of right size
  for(kidx in 1:dim(pred_difference_list)[3]){
    whichref = refidx[max(which(refidx<=kidx))]
    cat(paste(kidx, whichref, "\n"))
    pred_difference_list[,,kidx] = predictions$PredictedPost[,,kidx] - predictions$PredictedPost[,,whichref]
  }
  pred_difference = lapply(1:dim(pred_difference_list)[3], function(x) as.numeric(pred_difference_list[,,x]))

  return(
    data.frame(
      #j1val = NewDesignMat[, j1], 
      quantile = rep(allquantiles, length(quantile_j2)), 
      quantile_j2 = rep(quantile_j2, each=length(allquantiles)), 
      est = as.numeric(lapply(pred_difference, mean)), 
      lower =  as.numeric(lapply(pred_difference, quantile, 0.025)), 
      upper =  as.numeric(lapply(pred_difference, quantile, 0.975)), 
      var = as.numeric(lapply(pred_difference, var)), 
      n = as.numeric(lapply(pred_difference, length)))
  )
}


predest_independent_stratified = get_stratified_effect_estimate(
  bsrmod, X, C, j1=5, j2=4, quantile_j2=c(.25, .75), 
  ref_quantile=0.5, index_quantiles=setdiff(seq(0.05,.9,.05), .5))

bsrmod$InteractionPIP[4,5] # little evidence of interaction

library(ggplot2)
predest_independent_stratified |> 
  # mutate(se = sqrt(var/n), 
  #        upper = est + 1.96*se, 
  #        lower = est  - 1.96*se) |> 
  ggplot(aes(x=quantile, y=est, colour=factor(quantile_j2), fill=factor(quantile_j2))) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.15,linetype=2) + # makes it super messy
  geom_line() + 
  labs(y="Mean difference", x=paste("Quantile of", colnames(X)[4])) +
  scale_color_discrete(name=paste("Quantile of", colnames(X)[5]))+
  scale_fill_discrete(name=paste("Quantile of", colnames(X)[5])) + theme_bw()

# est interaction
get_inter_estimate <- function(
    NLmod, X, C, j1, j2, 
    quantile_j2 = c(0.25, 0.75), 
    quantile_rest=0.5, ref_quantile=0.25, 
    index_quantile=0.75){
  # stratified effect estimates for an interaction (different referent for each level of j2)
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' j1 = index of  chemical of interest
  #' j2 = index of  chemical of at levels of which interactions are assessed
  #' quantile_j2  = levels of j2 over which interaction is considered
  #' quantile_rest = quantile to fix other chemicals at
  #' ref_quantile = quantile to fix chemical of interest as reference
  #' index_quantiles = quantile to fix chemical of interest for comparisons (in calculations of mean difference)
  
  
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
  
  # create design of covariates (arithmetic mean, for simplicity)
  pc <- dim(C)[2]
  NewDesignC <- matrix(NA, (length(index_quantile)+1)*length(quantile_j2), pc + 1)
  NewDesignC[, 1] <- 1
  for (jc in 1:pc) {
    NewDesignC[, jc + 1] <- mean(C[, jc])
  }
  
  # create design of chemicals
  allquantiles  = c(ref_quantile, index_quantile)
  n <- dim(X)[1]
  
  NewDesignMat <- apply(X,2, function(x) rep(quantile(x, quantile_rest), length(allquantiles)*length(quantile_j2)))
  indexq = quantile(X[,j1], allquantiles)
  interq = quantile(X[,j2], quantile_j2)
  NewDesignMat[,j1] <- rep(indexq, length(quantile_j2))
  NewDesignMat[,j2] <- rep(interq, each=length(allquantiles))
  rownames(NewDesignMat) <- rep(names(indexq), length(quantile_j2))
  
  NewDesign <- array(NA, dim = c((length(index_quantile)+1)*length(quantile_j2), p, ns + 1))
  NewDesign[, , 1] <- 1
  for (j in 1:p) {
    temp_ns_object <- splines::ns(X[, j], df = ns)
    temp_sds <- apply(temp_ns_object, 2, sd)
    temp_means <- apply(temp_ns_object, 2, mean)
    NewDesign[, j, 2:(ns + 1)] <- t((t(                            
      predict(temp_ns_object,
              NewDesignMat[, j])
    ) - temp_means) / temp_sds)
  }
  
  
  # generate predictions
  predictions <- NLinteraction:::PredictionsMixture(
    XstarOld = Xstar,
    XstarNew = NewDesign,
    designC = NewDesignC,
    totalScans = totalScans,
    nChains = nChains,
    zetaPost = zetaPost,
    betaList = betaList,
    betaCPost = betaCPost,
    k = k,
    ns = ns
  )
  
  # nchains X total scans/iterations X # test data points
  
  refidx = seq(1,length(allquantiles)*length(quantile_j2), length(allquantiles))
  allidx = 1:(length(allquantiles)*length(quantile_j2))
  indexidx = setdiff(allidx, refidx)
  pred_difference_list = predictions$PredictedPost*NA # lazy way to create list of NAs of right size
  for(kidx in 1:dim(pred_difference_list)[3]){
    whichref = refidx[max(which(refidx<=kidx))]
    # cat(paste(kidx, whichref, "\n"))
    pred_difference_list[,,kidx] = predictions$PredictedPost[,,kidx] - predictions$PredictedPost[,,whichref]
  }
  pred_difference = lapply(1:dim(pred_difference_list)[3], function(x) as.numeric(pred_difference_list[,,x]))
  
  ests <- pred_difference[[indexidx[2]]]-pred_difference[[indexidx[1]]]
  
  return(
    data.frame(
      est = mean(ests), 
      lower =  as.numeric(quantile(ests, 0.025)), 
      upper =  as.numeric(quantile(ests, 0.975)), 
      se = sqrt(var(ests)/length(ests))
      )
  )
}

get_inter_estimate(bsrmod, X, C, 4, 5)
get_inter_estimate(bsrmod, X, C, 5, 4)
