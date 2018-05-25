# CBJ November 2017 (Coryn Bailer-Jones, calj@mpia.de)
# Functions for converting astrometry to phase space

# source("./metropolis.R")

library(mvtnorm)
source("constants.R") # only need kf

# Data model:
# 3D astrometry  (astro3d) = c(parallax, pmra, pmdec) in (mas, mas/yr, mas/yr)
# 3D phase space (phase3d) = c(r, v, phi) in (pc, km/s, radians)
# r is distance, v is tangential speed, phi is position angle from North to East
# astro3dCov is 3x3 covariance matrix of astro3d in same units (e.g. made by cov.astro3d)
# rlen is length scale of distance prior; vtmax is maximum tangential speed (parameter of speed prior)
# astromdat: matrix with named columns, from which we require:
# - astro3d
# - standard deviations in astro3d as well as their corresponding correlations. 
#   (name of which are given in cov.astro3d() below).

# Transform given astro3d to phase3d
astro3d.to.phase3d <- function(astro3d) {
  dist  <- 1e3/astro3d[1]
  speed <- kf*sqrt(astro3d[2]^2 + astro3d[3]^2)/astro3d[1]
  phi   <- atan(astro3d[2]/astro3d[3]) # to range -pi/2 to +pi/2
  if(astro3d[3]<0) { # quadrant correction, and put in range -pi to +pi
    phi <- ifelse(phi<0, phi+pi, phi-pi)
  }
  return(as.vector(c(dist, speed, phi))) # as.vector to get rid of names
}

# Evaluate Jacobian matrix |d_phase3d / d_astro3d| at a given astro3d
jac.phase3d.astro3d <- function (astro3d) {
  musq <- astro3d[2]^2 + astro3d[3]^2
  mu <- sqrt(musq)
  J <- matrix(data=c(-1e3/astro3d[1], 0, 0,
                     -kf*mu/astro3d[1]^2, kf*astro3d[2]/(mu*astro3d[1]), kf*astro3d[3]/(mu*astro3d[1]),
                     0, astro3d[3]/musq, -astro3d[2]/musq),
              nrow=3, ncol=3, byrow=TRUE)
  return(J)
}

# Make covariance matrix for astro3d from astromdat
cov.astro3d <- function(astromdat) {
  # V is upper diagonal of covariance matrix  
  V <- matrix(data=c(0, 
                     astromdat["parallax_error"]*astromdat["pmra_error"]*astromdat["parallax_pmra_corr"],
                     astromdat["parallax_error"]*astromdat["pmdec_error"]*astromdat["parallax_pmdec_corr"],
                     0, 0,
                     astromdat["pmra_error"]*astromdat["pmdec_error"]*astromdat["pmra_pmdec_corr"],
                     0, 0, 0),
              nrow=3, ncol=3, byrow=TRUE)
  S <- diag(astromdat[c("parallax_error", "pmra_error", "pmdec_error")])^2 # variances!
  return(t(V) + S + V)
}

# Return log10 (normalized) likelihood: P(3D astrometry | 3D phase space, Covariance)
loglike.astro3d <- function(phase3d, astro3d, astro3dCov) {
  a <- phase3d[2]/(kf*phase3d[1])
  pred <- 1e3*c(1/phase3d[1], a*sin(phase3d[3]), a*cos(phase3d[3]))
  return( (1/log(10))*dmvnorm(x=astro3d, mean=pred, sigma=astro3dCov, log=TRUE) )
}

# Return log10 (unnormalized) prior: P(3D phase space | rlen, vtmax)
logprior.phase3d <- function(phase3d, rlen, vtmax) {
    distPrior  <- ifelse(phase3d[1]>0, (1/(2*rlen^3))*phase3d[1]^2*exp(-phase3d[1]/rlen), 0)
    speedPrior <- dbeta(phase3d[2]/vtmax, shape1=2, shape2=3)
    anglePrior <- 1
    logPrior <- sum( log10(distPrior), log10(speedPrior), log10(anglePrior) )
    return(logPrior)
}

# Return log10 (unnormalized) posterior: P(3D phase space | 3D astrometry, Covariance, rlen, vtmax)
logpost.astro3d <- function(phase3d, astro3d, astro3dCov, rlen, vtmax) {
  logprior <- logprior.phase3d(phase3d, rlen, vtmax)
  if(is.finite(logprior)) { # only evaluate model if parameters are sensible
    return( c(logprior, loglike.astro3d(phase3d, astro3d, astro3dCov)) ) 
  } else {
    return( c(-Inf, -Inf) )
  }
}

# Initialize walkers for emcee by drawing Nwalker numbers independently for the
# three parameters from narrow distributions which respect boundary conditions of priors.
# Return as Nwalker x 3 matrix
init.phase3d <- function(Nwalker, phase3d, vtmax) {
  distwidth <- 10
  if(phase3d[1]<distwidth/2) {
    dist <- runif(Nwalker, min=0, max=distwidth)
  } else {
    dist <- runif(Nwalker, min=phase3d[1]-distwidth/2, max=phase3d[1]+distwidth/2)
  }
  speedwidth <- 10
  if(phase3d[2]<speedwidth/2) {
    speed <- runif(Nwalker, min=0, max=distwidth)
  }
  if(phase3d[2]>vtmax-speedwidth/2) {
    speed <- runif(Nwalker, min=vtmax-speedwidth, max=vtmax)
  }
  if(phase3d[2]>=speedwidth/2 && phase3d[2]<=vtmax-speedwidth/2) {
    speed <- runif(Nwalker, min=phase3d[2]-speedwidth/2, max=phase3d[2]+speedwidth/2)
  }
  anglewidth <- 10
  angle <- runif(Nwalker, min=phase3d[3]-anglewidth/2, max=phase3d[3]+anglewidth/2)
  return( cbind(dist, speed, angle) )
}

# Given a vector of angles (radian), return the same folded to either:
# range 0 to 2pi, or
# range -pi to +pi if this gives a range less than pi (i.e. is a smaller range)
fold.angular.vector <- function(inphi) {
  alphi <- inphi %% (2*pi)      # alphi has range 0 to 2pi
  bephi <- alphi
  sel   <- alphi>pi
  bephi[sel] <- alphi[sel]-2*pi # bephi has range -pi to +pi
  if(diff(range(bephi))<pi) {return(bephi)} else {return(alphi)}
}

# Plot chains, 1D marginal and 2D marginal posteriors, by writing
# to currently open device. 1D marginals are KDE for distance and speed,
# but histogram for angle.
plot.mcmc <- function(postSamp, phase3dNom, phase3dSig, phase3dnames, sourcename) {
  par(mfrow=c(3,3), mar=c(3.0,3.5,0.5,0.5), oma=0.5*c(1,1,4,2), mgp=c(1.8,0.6,0), cex=0.9)
  phi <- fold.angular.vector(postSamp[,5])
    for(j in 3:4) { # first two columns of postSamp (distance, speed)
    plot(1:nrow(postSamp), postSamp[,j], type="l",
        xlab="iteration", ylab=phase3dnames[j-2])
    postDen <- density(postSamp[,j], n=2^10)
    plot(postDen$x, postDen$y, type="l", lwd=1.5, yaxs="i", 
        ylim=1.05*c(0,max(postDen$y)), xlab=phase3dnames[j-2], ylab="density")
    abline(v=phase3dNom[j-2], col="blue")
    abline(v=quantile(postSamp[,j], probs=c(0.05,0.95)), col="red", lty=2)
    abline(v=median(postSamp[,j]), col="red")
    if(j==3) plot(postSamp[,3], postSamp[,4], xlab=phase3dnames[1], ylab=phase3dnames[2], pch=".")
    if(j==4) plot(postSamp[,3], phi,          xlab=phase3dnames[1], ylab=phase3dnames[3], pch=".")
  }
  # Plot for angle
  plot(1:nrow(postSamp), postSamp[,5], type="l", xlab="iteration", ylab=phase3dnames[3]) # intentionally not folded
  f  <- hist(phi, breaks=50, plot=FALSE)
  bs <- f$mids[2]-f$mids[1]
  plot(c(f$breaks[1]-bs,f$breaks), c(0,f$density,0), type="s", ylim=c(0,1.02*max(f$density)),
       xlab=phase3dnames[3], ylab="density")
  # Plot nominal value for both ranges (-pi/+pi, 0/2pi) - only one will show
  abline(v=c(0,2*pi)+phase3dNom[3], col="blue")
  abline(v=quantile(phi, probs=c(0.05,0.95)), col="red", lty=2)
  abline(v=median(phi), col="red")
  plot(postSamp[,4], phi, xlab=phase3dnames[2], ylab=phase3dnames[3], pch=".")
  #
  mtext(paste(sourcename, formatC(n, width=5)), adj=1, outer=TRUE)
}

# Given a covariance matrix, return a matrix with standard deviations
# on the leading diagonal and correlation coefficients on the off-diagonals
cov.to.sdcor <- function(covmat) {
  N <- nrow(covmat)
  if(N!=ncol(covmat)) {stop("cov.to.sdcor: matrix not square")}
  sdcor <- matrix(data=NA, nrow=N, ncol=N)
  for(i in 1:N) {
    sdcor[i,i] <- sqrt(covmat[i,i])
    for(j in setdiff(1:N, i)) {
      sdcor[i,j] <- covmat[i,j]/sqrt(covmat[i,i]*covmat[j,j])
    }
  }
  return(sdcor)
}

