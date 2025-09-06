
bivarinter_bkmr <- function(fit, z1, z2, qs.diff = c(0.25, 0.75), 
                            qs.fixed = c(0.25, 0.75), q.rest = 0.5) {
  #' fit = bkmr model
  #' z1 = index of chemical 1
  #' z2 = index of chemical 2
  #' qs.diff = quantiles to calculate response diff for
  #' qs.fixed = quantiles to fix other chemical at
  #' q.rest = quantile to fix rest of chemicals at
  #' 
  #' note that order of z1 and z2 don't matter
  
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

trivarinter_bkmr <- function(fit, z1, z2, z3, qs.diff = c(0.25, 0.75), 
                             qs.fixed = c(0.25, 0.75), q.rest = 0.5) {
  #' fit = bkmr model
  #' z1 = index of chemical 1
  #' z2 = index of chemical 2
  #' z3 = index of chemical 3
  #' qs.diff = quantiles to calculate response diff for
  #' qs.fixed = quantiles to fix other chemical at
  #' q.rest = quantile to fix rest of chemicals at
  #' 
  #' note that order of z1, z2, z3 don't matter
  
  # extract fit
  y <- fit$y
  Z <- fit$Z
  X <- fit$X
  
  df <- purrr::map_df(1:3, \(x) {
    if(x == 1) {
      c1 <- z1; c2 <- z2; c3 <- z3
    } else if (x == 2) {
      c1 <- z2; c2 <- z3; c3 <- z1
    } else {
      c1 <- z3; c2 <- z1; c3 <- z2
    }
    
    # fix all chems at q.rest
    point2 <- point1 <- apply(Z, 2, quantile, q.rest)
    # fix c2 at lower
    point2[c2] <- point1[c2] <- quantile(Z[, c2], qs.fixed[1])
    # fix c3 at lower
    point2[c3] <- point1[c3] <- quantile(Z[, c3], qs.fixed[1])
    # fix c1 at lower and upper
    point2[c1] <- quantile(Z[, c1], qs.diff[2])
    point1[c1] <- quantile(Z[, c1], qs.diff[1])
    newz.q1 <- rbind(point1, point2) # has all lower quantiles of c2, c3
    # fix all chems at q.rest
    point2 <- point1 <- apply(Z, 2, quantile, q.rest)
    # fix c2 at higher
    point2[c2] <- point1[c2] <- quantile(Z[, c2], qs.fixed[2])
    # fix c3 at higher
    point2[c3] <- point1[c3] <- quantile(Z[, c3], qs.fixed[2])
    # fix c1 at lower and upper
    point2[c1] <- quantile(Z[, c1], qs.diff[2])
    point1[c1] <- quantile(Z[, c1], qs.diff[1])
    newz.q2 <- rbind(point1, point2) # has all upper quantiles of c2, c3
    
    # prepare
    cc <- c(-1 * c(-1, 1), c(-1, 1))
    newz <- rbind(newz.q1, newz.q2)
    
    # default to using approximate calculation
    preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz)
    
    # extract intervals
    int <- drop(cc %*% preds$postmean)
    int.se <- drop(sqrt(cc %*% preds$postvar %*% cc))
    ints <- c(est = int, sd = int.se)
    
    return(data.frame(variable = colnames(Z)[c1], fixedat1 = colnames(Z)[c2], 
                      fixedat2 = colnames(Z)[c3],
                      est = ints["est"], sd = ints["sd"], row.names = NULL))
  })
  
 return(df)
}

trivarsurf_bkmr <- function(fit, z1, z2, z3, qs.diff = c(0.1, 0.5, 0.9), 
                            q.fixed = 0.5, ngrid = 50) {
  #' fit = bkmr model
  #' z1 = index of chemical 1
  #' z2 = index of chemical 2
  #' z3 = index of chemical 3
  #' qs.diff = quantiles to calculate response diff for
  #' q.fixed = quantiles to fix rest of chemicals at
  #' 
  #' note that order of z1, z2, z3 don't matter

  # call from fit
  y <- fit$y
  Z <- fit$Z
  X <- fit$X
  z.names <- colnames(Z)
  
  df <- purrr::map_df(1:3, \(x) {
    if(x == 1) {
      c1 <- z1; c2 <- z2; c3 <- z3
    } else if (x == 2) {
      c1 <- z2; c2 <- z3; c3 <- z1
    } else {
      c1 <- z3; c2 <- z1; c3 <- z2
    }
    
    # create new ordering
    ord <- c(c1, c2, c3, setdiff(1:ncol(Z), c(c1, c2, c3)))
    
    # create grid of z-values to evaluate at
    z1.grid <- seq(min(Z[, ord[1]]), max(Z[, ord[1]]), length = ngrid)
    z2.grid <- quantile(Z[, ord[2]], probs = qs.diff)
    z3.grid <- quantile(Z[, ord[3]], probs = qs.diff)
    z.all <- c(list(z1.grid), list(z2.grid), list(c(-99)))
    if (ncol(Z) > 3) {
      z.others <- lapply(4:ncol(Z), function(x) quantile(Z[, ord[x]], q.fixed))
      z.all <- c(z.all, z.others)
    }
    newz.grid <- expand.grid(z.all)
    newz.grid[, 3] <- rep(z3.grid, each = ngrid)
    z1save <- newz.grid[, 1]
    colnames(newz.grid) <- colnames(Z)[ord]
    newz.grid <- newz.grid[, colnames(Z)]
    
    # evaluate prediction, assume approx fit
    preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid)
    preds.mean <- preds$postmean
    preds.se <- sqrt(diag(preds$postvar))
    
    # return
    return(data.frame(z1_val = z1save, 
                      z23_q = rep(qs.diff, each = ngrid), 
                      est = preds.mean, 
                      se = preds.se, 
                      z1_name = rep(colnames(Z)[c1], length(z1save)), 
                      z2_name = rep(colnames(Z)[c2], length(z1save)), 
                      z3_name = rep(colnames(Z)[c3], length(z1save))))
  })
  
  return(df)
}

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

trivarsurf_bsr <- function(NLmod, X, C, j1, j2, j3, gridLength = 50, 
                           quantile_j23 = c(0.1, 0.5, 0.9), quantile_rest = 0.5) {
  #' NLmod = bsr model
  #' X = matrix or dataframe of chemical values used to fit model
  #' C = matrix or dataframe of covariate values used to fit model
  #' gridLength = number of points to estimate response at
  #' j1 = index of first chemical, chemical used as primary predictor
  #' j2 = index of second chemical
  #' j3 = index of third chemical
  #' quantile_j23 = vector of quantiles to fix second and third chemicals at
  #' quantile_rest = quantile to fix other chemicals at
  #' 
  #' NOTE order of j1, j2, j3 matters

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
  
  # for each quantile of j2, j3
  df <- purrr::map_df(quantile_j23, \(quantile_j23) {
    # create design of chemicals
    n <- dim(X)[1]
    NewDesignMat <- matrix(NA, gridLength, p)
    for (j in 1:p) {
      NewDesignMat[, j] <- quantile(X[, j], quantile_rest)
    }
    NewDesignMat[, j1] <- seq(quantile(X[, j1], 0.025), 
                              quantile(X[, j1], 0.975), length = gridLength)
    NewDesignMat[, j2] <- quantile(X[, j2], quantile_j23)
    NewDesignMat[, j3] <- quantile(X[, j3], quantile_j23)
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
      j23quant = rep(quantile_j23, gridLength), 
      est = apply(predictions$PredictedPost, 3, mean), 
      lower = apply(predictions$PredictedPost, 3, quantile, 0.025), 
      upper = apply(predictions$PredictedPost, 3, quantile, 0.975)
    ))
  })
  return(df)
}

