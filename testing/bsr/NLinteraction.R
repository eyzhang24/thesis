# installling NLinteraction
# library(devtools)
# install_github(repo = "jantonelli111/NLinteraction")
library(NLinteraction)

# simulate data
set.seed(0)
n = 100
p = 10
pc = 1

X = matrix(rnorm(n*p), n, p)

C = matrix(rnorm(n*pc), nrow=n)

TrueH = function(X) {
  return(1.5*(X[,2]*X[,3]) - 1.6*(X[,4]^2 * X[,5]))
}

Y = 5 + C + TrueH(X) + rnorm(n)

# build model
NLmod2 = NLint(Y=Y, X=X, C=C, nIter=10000, nBurn=5000, thin=5, nChains=2, ns=2)
NLmod3 = NLint(Y=Y, X=X, C=C, nIter=10000, nBurn=5000, thin=5, nChains=2, ns=3)

# evaluate WAIC
print(c(NLmod2$waic,NLmod3$waic))

# choose model 2
NLmod = NLmod2

# univariate pip's
NLmod$MainPIP
barplot(NLmod$MainPIP)

# bivariate pip's
NLmod$InteractionPIP
plotInt(NLmod = NLmod)

# visualize exposure response relationship
par(mfrow=c(1,3), pty='s')
plotSurface1d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
              gridLength=30, quantile_j2=0.2, quantile_rest = 0.5,
              xlab="X4", ylab="Posterior predictive distribution",
              main="20th quantile of X5")
plotSurface1d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
              gridLength=30, quantile_j2=0.5, quantile_rest = 0.5,
              xlab="X4", ylab="Posterior predictive distribution",
              main="50th quantile of X5")
plotSurface1d(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
              gridLength=30, quantile_j2=0.8, quantile_rest = 0.5,
              xlab="X4", ylab="Posterior predictive distribution",
              main="80th quantile of X5")

# full bivariate plot for smaller grid lengths
par(mfrow=c(1,2), pty='s')
plotSurface2dMean(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
                  gridLength_j1=20, gridLength_j2 = 20,
                  quantile_rest = 0.5, xlab='X4', ylab='X5',
                  main="Posterior mean")
points(X[,4], X[,5], pch=16, cex=0.5)
plotSurface2dSD(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
                gridLength_j1=20, gridLength_j2 = 20,
                quantile_rest = 0.5, xlab='X4', ylab='X5',
                main="Posterior SD")
points(X[,4], X[,5], pch=16, cex=0.5)

# full bivariate plot for longer grid lengths
plotSurface2dMean(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
                  gridLength_j1=40, gridLength_j2 = 40,
                  quantile_rest = 0.5, xlab='X4', ylab='X5',
                  main="Posterior mean")
points(X[,4], X[,5], pch=16, cex=0.5)
plotSurface2dSD(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
                gridLength_j1=40, gridLength_j2 = 40,
                quantile_rest = 0.5, xlab='X4', ylab='X5',
                main="Posterior SD")
points(X[,4], X[,5], pch=16, cex=0.5)

# prevent extrapolation
par(mfrow=c(1,2), pty='s')
plotSurface2dMean(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
                  gridLength_j1=40, gridLength_j2 = 40,
                  quantile_rest = 0.5, xlab='X4', ylab='X5',
                  main="Posterior mean", minDist=0.5)
points(X[,4], X[,5], pch=16, cex=0.5)
plotSurface2dSD(NLmod = NLmod, X=X, C=C, j1=4, j2=5,
                gridLength_j1=40, gridLength_j2 = 40,
                quantile_rest = 0.5, xlab='X4', ylab='X5',
                main="Posterior SD", minDist=0.5)
points(X[,4], X[,5], pch=16, cex=0.5)