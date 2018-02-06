# CBJ November 2017 (Coryn Bailer-Jones, calj@mpia.de)
# The Metropolis algorithm

library(mvtnorm) # for rmvnorm

# Metropolis (MCMC) algorithm to sample from function func.
# The first argument of func must be a real vector of parameters, 
# the initial values of which are provided by the real vector thetaInit.
# func() returns a two-element vector, the logPrior and logLike 
# (log base 10), the sum of which is taken to be the log of the density 
# function (i.e. unnormalized posterior). If you don't have this separation,
# just set func to return one of them as zero. The MCMC sampling PDF is the 
# multivariate Gaussian with fixed covariance, sampleCov. A total of 
# Nburnin+Nsamp samples are drawn, of which the last Nsamp are kept. As the 
# sampling PDF is symmetric, the Hasting factor cancels, leaving the basic 
# Metropolis algorithm. Diagnostics are printed very verbose^th sample: 
# sample number, acceptance rate so far.
# ... is used to pass data, prior parameters etc. to func().
# If demo=FALSE (default), then
# return a Nsamp * (2+Ntheta) matrix (no names), where the columns are
# 1:  log10 prior PDF
# 2:  log10 likelihood
# 3+: Ntheta parameters
# (The order of the parameters in thetaInit and sampleCov must match.)
# If demo=TRUE, return the above (funcSamp) as well as thetaPropAll, a 
# Nsamp * Ntheta matrix of proposed steps, as a two element named list.
metrop <- function(func, thetaInit, Nburnin, Nsamp, sampleCov, verbose, 
                   demo=FALSE, ...) {

  Ntheta   <- length(thetaInit)
  thetaCur <- thetaInit
  funcCur  <- func(thetaInit, ...) # log10
  funcSamp <- matrix(data=NA, nrow=Nsamp, ncol=2+Ntheta) 
  # funcSamp will be filled and returned
  nAccept  <- 0
  acceptRate <- 0
  if(demo) {
    thetaPropAll <- matrix(data=NA, nrow=Nsamp, ncol=Ntheta)
  }
  
  for(n in 1:(Nburnin+Nsamp)) {

    # Metropolis algorithm. No Hastings factor for symmetric proposal
    if(is.null(dim(sampleCov))) { # theta and sampleCov are scalars
      thetaProp <- rnorm(n=1, mean=thetaCur, sd=sqrt(sampleCov))
    } else {
      thetaProp <- rmvnorm(n=1, mean=thetaCur, sigma=sampleCov, 
                           method="eigen")
    }
    funcProp  <- func(thetaProp, ...) 
    logMR <- sum(funcProp) - sum(funcCur) # log10 of the Metropolis ratio
    #cat(n, thetaCur, funcCur, ":", thetaProp, funcProp, "\n")
    if(logMR>=0 || logMR>log10(runif(1, min=0, max=1))) {
      thetaCur   <- thetaProp
      funcCur    <- funcProp
      nAccept    <- nAccept + 1
      acceptRate <- nAccept/n
    }
    if(n>Nburnin) {
      funcSamp[n-Nburnin,1:2] <- funcCur
      funcSamp[n-Nburnin,3:(2+Ntheta)] <- thetaCur
      if(demo) {
        thetaPropAll[n-Nburnin,1:Ntheta] <- thetaProp
      }
    }

    # Diagnostics
    if( is.finite(verbose) && (n%%verbose==0 || n==Nburnin+Nsamp) ) {
      s1 <- noquote(formatC(n,          format="d", digits=5, flag=""))
      s2 <- noquote(formatC(Nburnin,    format="g", digits=5, flag=""))
      s3 <- noquote(formatC(Nsamp,      format="g", digits=5, flag=""))
      s4 <- noquote(formatC(acceptRate, format="f", digits=4, width=7, 
                            flag=""))
      cat(s1, "of", s2, "+", s3, s4, "\n")
    }

  }

  if(demo) {
    return(list(funcSamp=funcSamp, thetaPropAll=thetaPropAll))
  } else {
    return(funcSamp)
  }
 
}
