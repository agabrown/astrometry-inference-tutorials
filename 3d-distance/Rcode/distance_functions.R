# CBJ November 2017 (Coryn Bailer-Jones, calj@mpia.de)
# Functions for distance and inference

source("metropolis.R") # for quantiles.distpost3
library(PolynomF)    # for mode.distpost3
library(mvtnorm)     # for d.likemulti

##### Naming conventions

# w      - measured parallax, as
# wsd    - uncertainty in w
# r      - inferred distance, pc
# rTrue  - true distance, pc (if doing simulations)
# fTrue  - wsd/w (fractional parallax uncertainty) 
# rMean/rMedian/rMode - mean/median/mode of posterior
#  d.xxx - density of PDF xxx (normalized)
# ud.xxx - density of PDF xxx (unnormalized)
#  r.xxx - random sample from PDF xxx
# like   - likelihood PDF
# prior  - prior PDF
# post   - posterior PDF
# {}     - in comments indicates a vector

##### General functions

# Return 1 if rlo <= r <= rhi, otherwise zero. Vectorized in any one parameter.
lim <- function(r, rlo, rhi) ifelse(r<rlo | r>rhi, 0, 1)

# Normalized Gaussian likelihood in w. Vectorized in any one parameter.
d.like <- function(w, r, wsd) dnorm(x=w, mean=1/r, sd=wsd)

##### (1) uniform distance prior (between rlo and rhi)

# Vectorized in any one parameter
d.distprior1 <- function(r, rmax) lim(r, 0, rmax)*1/rmax 
r.distprior1 <- function(Nsamp, rmax) runif(Nsamp, min=0, max=rmax)

# Vectorized in any one parameter
ud.distpost1 <- function(r, w, wsd, rlo, rhi) {
  d.like(w, r, wsd)*lim(r, rlo, rhi)
} 


##### (2) uniform space density prior (between rlo and rhi)

d.distprior2 <- function(r, rmax) lim(r, 0, rmax)*(3/rmax^3)*r^2 # Vectorized in r or rmax
r.distprior2 <- function(Nsamp, rmax) {  # Vectorized in Nsamp
  # Use rejection sampling. Expect 1/3 * Nsamp samples to be returned
  r1 <- runif(Nsamp, min=0, max=1)
  r2 <- runif(Nsamp, min=0, max=1)
  r  <- rmax*r1[which(r2<=r1^2)]
  return(r)
}

# Vectorized in any one parameter
ud.distpost2 <- function(r, w, wsd, rlo, rhi) {
  d.like(w, r, wsd)*lim(r, rlo, rhi)*r^2
}

# Define posterior function func() as required by metrop()
func.distpost2 <- function(r, w, wsd, rmax) {
  return( c(log10(d.distprior2(r, rmax)), log10(d.like(w, r, wsd))) )
}

# The mode of the posterior PDF is its smaller root; undefined when f>1/sqrt(8) (gives NaN)
mode.distpost2 <- function(w, wsd) {(w/(4*wsd^2))*(1 - sqrt(1 - 8*(wsd/w)^2))}


##### (3) exponentially decreasing space density prior (with length scale rlen)

# Vectorized in r or rlen
d.distprior3 <- function(r, rlen) lim(r,0,Inf)*(1/(2*rlen^3))*r^2*exp(-r/rlen) 

# Cannot use rejection sampling as support of r is semi-infinite. So use MCMC.
# Acceptance rate is around 0.7. I don't deal with correlations, so small Nsamp can and does
# return some identical values.
r.distprior3 <- function(Nsamp, rlen) {
  return(metrop(func=func.distprior3, thetaInit=rlen, Nburnin=50, 
                Nsamp=Nsamp, sampleCov=2*rlen^2, verbose=Inf, rlen=rlen)[,3])
}
# Test prior sampling by comparing the following
# rlen <- 1e3
# x=seq(from=1,to=(10*rlen),length.out=1e3); plot(x, x^2*exp(-x/rlen), type="l")
# truehist(r.distprior3(1e5, rlen))

# Vectorized in any one parameter
ud.distpost3 <- function(r, w, wsd, rlen) {
  d.like(w, r, wsd)*exp(-r/rlen)*r^2
}

# Define M=ud.distpost3(r=rMode) and F(r)=log(ud.distpost3(r) - log(M/2)
# Function returns F(r) which we can use in a root finding algorithm to find r,
# which will gives bound of FWHM. (If wanted x*Max instead of half*Max, 
# use x*M instead of M/2 in expression.)
rootfunc.ud.distpost3 <- function(r, w, wsd, rlen, M) {
  #return(2*log(r) - r/rlen - (1/(2*wsd^2))*(w-1/r)^2 - log(wsd*M/2) - 0.5*log(2*pi))
  return(log(ud.distpost3(r, w, wsd, rlen))-log(M/2))
  }

# Define func() as required by metrop() for prior and posterior sampling
func.distprior3 <- function(r, rlen) {
  return( c(log10(d.distprior3(r, rlen)), 0) )
}
func.distpost3 <- function(r, w, wsd, rlen) {
  return( c(log10(d.distprior3(r, rlen)), log10(d.like(w, r, wsd))) )
}

# Return posterior mode - a single real number (or NA if something is wrong).
# If retall=TRUE, return all three roots of derivative of PDF (default is FALSE).
# Inputs:
# w    - parallax,             unrestricted
# wsd  - parallax uncertainty, must be >0
# rlen - prior length scale,   must be >0 and <Inf
# There are then only two types of solutions of the cubic equation: 
# 1. only one root is real => one maxima. Take this as the mode.
# 2. all three roots are real => two maxima. 
#    if w>=0, take the smallest.
#    if w<0   take the positive one (there should be only one).
mode.distpost3 <- function(w, wsd, rlen, retall=FALSE) {
  # special cases:
  # w<=0 works normally, except for w=-Inf
  if(w==-Inf)  return(Inf)
  if(w==Inf)   return(0)
  if(wsd==0)   return(1/w)
  if(wsd==Inf) return(2*rlen)
  if(wsd<0)    return(NA) # makes no sense
  if(rlen<=0)  return(NA) # makes no sense
  r <- polynom()
  p <- r^3/rlen - 2*r^2 + (w/wsd^2)*r - 1/wsd^2
  roots <- solve(p) # complex vector of length 3, sorted by increasing size of real part
  rMode <- switch(EXPR = toString(length(which(Im(roots)==0))),
                  "0" = NA,
                  "1" = Re(roots[which(Im(roots)==0)]),
                  "2" = NA,
                  "3" = ifelse(w>0, min(roots), roots[roots>0]) # should be real and unique
  )
  if(retall) {
    return(roots)
  } else {
    return(rMode)
  }
}

# Compute bounds of FWHM by root finding. Providing rMode avoids recalculating this.
# rmax is maximum value of upper bound (needed for root finding)
fwhm.distpost3 <- function(w, wsd, rMode=NA, rlen, rmax=1e6) {
  if(is.na(rMode)) {
    rMode <- mode.distpost3(w, wsd, rlen)
  }
  M <- ud.distpost3(r=rMode, w, wsd, rlen)
  lo <- uniroot(rootfunc.ud.distpost3, interval=c(0, rMode), 
                w=w, wsd=wsd, rlen=rlen, M=M)$root
  hi <- uniroot(rootfunc.ud.distpost3, interval=c(rMode, rmax), 
                w=w, wsd=wsd, rlen=rlen, M=M)$root
  return(c(lo=lo,hi=hi))
}

# Return value(s) at quantile(s) specified by probs (default is just median)
# Computed via MCMC, as integral appears to be non-analytic.
# Both Metropolis and HybridMC are slow and sometimes inaccurate.
# If diagnose=TRUE and MCMC is performed, plot histogram of samples (as well as 
# posterior PDF Z!=NA), report mid and final acceptance rates and step size.
quantiles.distpost3 <- function(w, wsd, rlen, probs=0.5, Nsamp=1e5, diagnose=FALSE) {
  if(diagnose) {cat("Inputs:", w, wsd, rlen, probs, Nsamp, "\n")}
  if(wsd<=0)  return(NA) # makes no sense
  if(rlen<=0) return(NA) # makes no sense
  rMode <- mode.distpost3(w, wsd, rlen)
  stepSize <- rMode # originally used rSD if it existed, otherwise rMode
  verbose <- ifelse(diagnose, Nsamp/2, Inf)
  samp  <- metrop(func=func.distpost3, thetaInit=rMode, Nburnin=1e-3*Nsamp, Nsamp=Nsamp,
                  sampleCov=stepSize^2, verbose=verbose, w=w, wsd=wsd, rlen=rlen)[,3]
  #samp <- as.vector(hybridMC(y.start=rMode, n.samp=Nsamp, logDens=logDens.post3, 
  #                           dLogDens=dLogDens.post3, epsilon=0.1, LFsteps=1, 
  #                           w=w, wsd=wsd, rlen=rlen))
  if(diagnose) {
    cat("stepSize =", stepSize, "\n")
    truehist(samp)
    # need to compute normalization constant Z for the following
    # if may not exist (see original code)
    #Z <- Z.post3(w, wsd, rlen)
    #r <- seq(from=min(samp), to=max(samp), length.out=1e4)
    #lines(fac*r, d.post3(r, w, wsd, rlen, Z))
  }
  return(quantile(samp, probs=probs))
}

##### (3b) - as (3), but with multiple stars and covariance

# Return correlation coefficient between parallaxes for two
# sources with specified sky positions
# and amplitude amp (0-1) and length scale len (degrees)
# using an exponential correlation model.
# Note that this uses a simplified angular separation computation
# which does not work for large separations or near the poles.
parcor <- function(ra1, dec1, ra2, dec2, amp=0.5, len=0.25) {
  sep <- sqrt(((ra1-ra2)*cos(conv*0.5*(dec1+dec2)))^2 + (dec1-dec2)^2)
  return( amp*exp(-sep/len))
}

# Return covariance matrix between parallaxes for all sources in 
# data frame dat, using correlation function parcor.
parcovmat <- function(dat, amp, len) {
  V <- matrix(data=0, nrow=nrow(dat), ncol=nrow(dat))
  for(i in 2:nrow(dat)) {
    for(j in 1:(i-1)) {
      V[i,j] <- dat$parallax_error[i]*dat$parallax_error[j] * 
        parcor(ra1=dat$ra[i], dec1=dat$dec[i], 
               ra2=dat$ra[j], dec2=dat$dec[j], amp=amp, len=len)
      V[j,i] <- V[i,j]
    }
  }
  return(V + diag(dat$parallax_error^2))
}

# Return density of a normalized dim(w) dimensional Gaussian likelihood 
# in {w} with mean (1/r)*{1} and covariance matrix parcovmat.
# w and wsd must be vectors of same size, and r a scalar.
# If parcovmat=NULL, then assume product of dim(w) 1D Gaussians.
d.likemulti <- function(w, r, wsd, parcovmat=NULL) {
  if(length(w)!=length(wsd)) stop("w and wsd must be same length")
  if(length(r)==1) {
    if(is.null(parcovmat)) {
      return(prod(dnorm(x=w, mean=1/r, sd=wsd)))
    } else {
      return(dmvnorm(x=w, mean=rep.int(1/r, length(w)), sigma=parcovmat))
    }
  } else {
    like <- double(length(r))
    for(i in 1:length(r)) {
      if(is.null(parcovmat)) {
        like[i] <- prod(dnorm(x=w, mean=1/r[i], sd=wsd))
      } else {
        like[i] <- dmvnorm(x=w, mean=rep.int(1/r[i], length(w)), sigma=parcovmat)
      }
    }
    return(like)
  }
}

# Requires vector w and wsd, scalar rlen, and scalar or vector r.
# Returns nonormalized posterior values for each r.
ud.distpost3multi <- function(r, w, wsd, parcovmat=NULL, rlen) {
  d.likemulti(w, r, wsd, parcovmat)*exp(-r/rlen)*r^2
}

##### Likelihood for multiple stars and finite cluster size (no correlations)

# The functions below compute the product of N 1-dimensional integrals,
# each integral being over the true, unknown distance for each star (likelihood for that star).
# The functions differ in cluster geometry assumptions and how the integral is done, but they
# all compute P({w} | rc, sc, ...)

# 1. Assume true distances from cluster center are described by a 1D Gaussian
# with mean 0 and stdev sc, along l.o.s: cluster has negligible extent 
# transverse to l.o.s. Uses binomial approximation in quadratic expansion 
# to make integral a standard one (with error function).
# (It turns out this is quite a bad approximation.)
# Given vectors or scalars {w}, {wsd}, and scalar rc, sc, 
# compute (normalized) likelihood P({w} | {wsd}, rc, sc), 
# where rc and sc are cluster distance.
d.likecluster1 <- function(w, wsd, rc, sc) {
  if(any(wsd<=0) || rc<=0 || sc<=0) stop("wsd, rc, sc must all be positive")
  like <- double(length(w))
  for(i in 1:length(w)) {
    gg <- (w[i]*rc-1)/(wsd[i]^2*rc^3)
    bb <- wsd[i]^2*rc^4/(2*(1+wsd[i]^2*rc^4/sc^2))
    integ <- sqrt(pi*bb)*exp(bb*gg^2)*(1+sqrt(pi)/2-sqrt(pi)*pnorm(gg*sqrt(2*bb)))
    like[i] <- integ*exp(-(w[i]*rc-1)^2/(2*wsd[i]^2*rc^2))/(2*pi*wsd[i]*sc)
  }
  return(prod(like))
}

# 2. Assume cluster is spherical, i.e. true distances from cluster centre
# are described by an isotropic 3D Gaussian with mean 0 and stdev sc.
# No other geometric approximation, so each integral is done numerically.
# costheta = cos(theta) where theta is angular separation between assumed 
# cluster centre and star (as seen by observer).
# However, if theta<10 deg or so, cos(theta)~1 and can reduce to 1D case
# again, but this time with numerical integral, not binomial approximation.
# In that case no need to pass costheta (as NULL will trigger the 1D case).
# Given vectors or scalars {w}, {wsd}, {costheta}, and scalar rc, sc, 
# compute (normalized) likelihood P({w} | {wsd}, {costheta}, rc, sc), 
# where rc and sc are cluster distance.
# If retlog=TRUE then return the (natural) log likelihood
d.likecluster2 <- function(w, wsd, costheta=NULL, rc, sc, retlog=TRUE) {
  if(any(wsd<=0) || rc<=0 || sc<=0) stop("wsd, rc, sc must all be positive")
  # integrand takes a vector r (others all scalar) and returns vector of same length
  intFail <- 0
  integrand <- function(r, w, wsd, costheta, rc, sc) { 
    if(is.null(costheta)) { 
      return( d.like(w=w, r=r, wsd=wsd)*dnorm(x=r, mean=rc, sd=sc) )
    } else {
    #stop("unsolved problem for 3D case (costheta=/=NULL): numerical integration often doesn't converge") 
    x   <- sqrt(r^2 + rc^2 - 2*r*rc*costheta) # term in sqrt always non-negative
    tl  <- double(length(x))
    sel <- which(x/rc<1e-6) # special case to avoid divide by small x (i.e. neglect costheta).
                            # keep limit very small, or else it introduces artefacts.
    tl[sel] <- d.like(w=w, r=r[sel], wsd=wsd)*dnorm(x=r[sel], mean=rc, sd=sc)
    sel <- which(x/rc>=1e-6)
    tl[sel] <- d.like(w=w, r=r[sel], wsd=wsd)*dnorm(x=x[sel], mean=0, sd=sc)*(r[sel]-rc*costheta)/x[sel]
    return(tl)
    # The following is the principle, but has no trap for small x.
    #return( d.like(w=w, r=r, wsd=wsd)*dnorm(x=x, mean=0, sd=sc)*(r-rc*costheta)/x )
    }
  }
  like <- double(length(w))
  for(i in 1:length(w)) {
    # In principle limits are (0,Inf), but finite ones are faster.
    # Also works okay if costheta is NULL.
    like[i] <- integrate.func(integrand, lower=max(0,rc-5*sc), upper=rc+5*sc,
                              w=w[i], wsd=wsd[i], costheta=costheta[i], rc=rc, sc=sc)
    if(is.na(like[i])) { # if numerical integral fails, revert to 1D case.
                         # this may not be valid, but what else can we do?
      intFail <- intFail+1
      like[i] <- integrate.func(integrand, lower=max(0,rc-5*sc), upper=rc+5*sc,
                                w=w[i], wsd=wsd[i], costheta=NULL, rc=rc, sc=sc)
    }
    #cat(like[i], "\n")
  }
  if(intFail>0) cat("Numerical integration failures:", intFail, "\n")
  #cat(intFail, prod(like), "\n")
  if(retlog) {
    return(sum(log(like))) 
  } else {
    return(prod(like)) 
  }
}

