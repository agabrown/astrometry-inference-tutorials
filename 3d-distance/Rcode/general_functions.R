# CBJ November 2017 (Coryn Bailer-Jones, calj@mpia.de). Corrected 21 June 2018.
# General functions

# Return density of a beta function shifted and scaled to lie in range
# specified by xrange (two element vector).
d.betagen <- function(x=NA, xrange=c(0,1), shape1=NA, shape2=NA) {
  if(xrange[2]<=xrange[1]) {stop("Need xrange[2]>xrange[1]")}
    return( dbeta((x-xrange[1])/(xrange[2]-xrange[1]), shape1=shape1, shape2=shape2)/(xrange[2]-xrange[1]) )
}

# Return integral of f over its first argument from lower to upper.
# Returns NA if it cannot be calculated.
integrate.func <- function(f, lower, upper, ..., subdivisions=1e3, rel.tol=1e-6, abs.tol=1e-6) {
  integ <- integrate(f=f, lower=lower, upper=upper, ..., subdivisions=subdivisions, 
                     rel.tol=rel.tol, abs.tol=abs.tol, stop.on.error=FALSE)
  return( ifelse(integ$message=="OK" && integ$value>0 && integ$abs.error>0, integ$value, NA) )
}

# Given a 1D unnormalized posteior PDF "dense" computed on a grid "x", 
# return 3-element named list of normalization constant, mean, SD (essentially moments 0,1,2).
# The grid must have sufficient density and span to cover essentially all of the density,
# but it need not be uniform (CHECK)
# See PBI section 5.1.1 for details.
pdfmom <- function(dense, x) {
  if(length(x)!=length(x)) stop("dense and x must have same size")
  if(length(x)<=1) stop("dense and x must have at least two elements")
  deltax <- mean(diff(x))
  xNorm <- deltax*sum(dense)
  dense <- dense/xNorm # dense is now normalized
  xMean <- deltax*sum(x*dense)
  xSD   <- sqrt( deltax * sum((x-xMean)^2 * dense) )
  return(list(Z=xNorm, mean=xMean, sd=xSD))
}
