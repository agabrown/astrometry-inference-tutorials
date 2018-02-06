# CBJ November 2017 (Coryn Bailer-Jones, calj@mpia.de)
# Emcee algorithm

# affine-invariant ensemble sampler (emcee) to sample from func().
# A set of Nwalker walkers are evolved over Nburnin+Nsamp iterations. (Nwalker>1)
# The initial walkers are passed in as the 2D matrix thetaInit[Nwalker, Ntheta].
# These must all give finite function values.
# boundaries is either NULL, or a list of length Ntheta, each elements of which is either NULL or a
# two element vectors of form c(lo,hi) with hi>lo which indicate periodic boundaries for that parameter.
# If the latter case, the thetaInit must lie in the range [lo,hi] (is not checked).
# This means that the maximum distance between any two values of theta is (hi-lo)/2.
# If the parameters have been transformed, then we need to supply the function jacobian to
# return the Jacobian (a scalar) of the parameters. (Note that the parameters within
# metrop() are the transformed parameters.)
# Return: funcSamp, which is a 3D array funcSamp[iteration, walker, theta], of size
# (Nburnin+Nsamp) * Nwalker * (2+Ntheta) matrix (no names)
# which is the set of all iterations (including the burn-in) for all walkers.
# funcSamp[,walker,] is a 2D array with and one row per sample and columns
# 1:  log10 prior PDF
# 2:  log10 likelihood
# 3+: Ntheta parameters
# Note that the walkers are updated asynchronously, i.e. a walker is updated in an iteration
# immediately. Thus the walker w_j used to update walker w_k may already have been updated in this
# iteration. This is required to maintain detailed balance. Note that Algorithm 2 in arXiv:1202.3665v3
# is actually wrong in this respect, as it claims that w_j is only drawn from those at iteration t.
emcee <- function(func, jacobian=NULL, thetaInit, boundaries=NULL, Nburnin, Nsamp, verbose, ...) {
  
  Nwalker <- nrow(thetaInit)
  Ntheta  <- ncol(thetaInit)
  if(is.finite(verbose)) {
    cat("Sampling the posterior with emcee algorithm with", Nwalker, "walkers\n")
    cat("iteration # of Nburnin + Nsamp acceptRate\n")
  }
  
  # setup
  acon <- 2 # if acceptance rate is too low, increase it by lowering acon (and vice versa)
  walker <- thetaInit
  walkerFunc <- matrix(data=NA, nrow=Nwalker, ncol=2)  # function values of the walkers
  walkerLJac <- vector(mode="numeric", length=Nwalker) # log10(Jacobian) of the walkers
  funcSamp <- array(data=NA, dim=c(Nsamp+Nburnin, Nwalker, 2+Ntheta)) # this will be filled and returned
  if(is.null(jacobian)) {
    jacobian <- function(...) {
      return(1)
    }
  }
  if(!is.null(boundaries)) {
    if(!is.list(boundaries) || length(boundaries)!=Ntheta) {
      stop("emcee(): boundaries must be NULL or a list of length equal to number of parameters")
    } # if boundaries is NULL, then so are all its elements
  }
  
  # initialize
  for(w in 1:Nwalker) {
    walkerFunc[w,] <- func(walker[w,], ...) # log10
    walkerLJac[w]  <- log10(jacobian(walker[w,], ...))
  }
  if(!(all(is.finite(walkerFunc)))) {
    stop("emcee(): Initialization gives non-finite function values")
  }
  if(!(all(is.finite(walkerLJac)))) {
    stop("emcee(): Initialization gives non-finite Jacobian")
  }
  thetaProp <- rep(NA, Ntheta)
  nAccept  <- 0
  acceptRate <- 0
  startTime=Sys.time()

    # walk 
  for(n in 1:(Nburnin+Nsamp)) {
    
    for(w in 1:Nwalker) { 
      #   for(w in sample.int(Nwalker)) { # do this instead to randomly permute order of updates
      z <- (1/acon)*(1 + (acon-1)*runif(min=0, max=1, n=1))^2
      sel <- sample(x=setdiff(1:Nwalker, w), size=1) # select another walker
      for(j in 1:Ntheta) { # Compute proposed step, taking care of periodic boundary conditions per parameter
                           # Don't need to separately fold distance and final parameter.
                           # Just compute proposal, then fold.
        thetaProp[j] <- walker[sel,j] + z*(walker[w,j] - walker[sel,j])
        if(!is.null(boundaries[[j]])) {
          thetaProp[j] <- (thetaProp[j] - boundaries[[j]][1]) %% diff(boundaries[[j]]) + boundaries[[j]][1]
        } 
      }
      if(all(is.finite(thetaProp))) {
        funcProp <- func(thetaProp, ...)
        ljacProp <- log10(jacobian(thetaProp, ...))
        logSR <- (Ntheta-1)*log10(z) + sum(funcProp) - sum(walkerFunc[w,]) + ljacProp - walkerLJac[w] # log10 of selection ratio
      } else {
        warning("emcee(): thetaProp = ", thetaProp, "\n  Not all thetaProp are finite so the proposal is being rejected.\n  Probably the theta values are very large (positive or negative), which is not good behaviour")
        logSR <- -Inf
      }
      if(is.nan(logSR)) { # should never occur given other traps
        stop("emcee(): logSR=NaN at n=", n, "w=", w, "funcProp=", funcProp, "walkerFunc=", walkerFunc)
      }
      if(logSR>=0 || logSR>log10(runif(1, min=0, max=1))) { # runif will not generate 0 or 1 exactly
        walker[w,] <- thetaProp
        walkerFunc[w,] <- funcProp
        walkerLJac[w]  <- ljacProp
        nAccept <- nAccept + 1
      } # otherwise walker is unchanged
      funcSamp[n,w,1:2] <- walkerFunc[w,]
      funcSamp[n,w,3:(2+Ntheta)] <- walker[w,]
    } # End of walker loop
    acceptRate <- nAccept/(n*Nwalker)
    
    # diagnostics
    if( is.finite(verbose) && (n%%verbose==0 || n==Nburnin+Nsamp) ) {
      sdump1 <- noquote( formatC(n,          format="d", digits=5, flag="") )
      sdump2 <- noquote( formatC(Nburnin,    format="g", digits=5, flag="") )
      sdump3 <- noquote( formatC(Nsamp,      format="g", digits=5, flag="") )
      sdump4 <- noquote( formatC(acceptRate, format="f", width=7, digits=4, flag="") )
      cat(sdump1, "of", sdump2, "+", sdump3, sdump4, "\n")
    }
    
  }
  
  return(funcSamp)
  
}
