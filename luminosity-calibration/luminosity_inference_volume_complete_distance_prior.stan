/*
 * Stan model for the problem of inferring the mean absolute magnitude (and its variance) from a
 * parallax survey in which apparent magnitudes and parallaxes of stars are observed and the errors
 * on the measurements are reported. The model assumes a uniform distribution of stars in space
 * (constant density), with their distances limited by an upper and lower bound.
 *
 * The survey is assumed to be volume complete (i.e. there is no magnitude limit).
 *
 * See chapter 16 in 'Astrometry for Astrophysics', 2013, edited by W.F. van Altena (Cambridge
 * University Press).
 *
 * Anthony Brown Nov 2017 - Jan 2018
 */

functions {
    /**
     * Provide the PDF of the distance distribution of stars around the observer which varies as
     * r^alpha, where alpha>=0. For stars distributed with uniform space density out to a certain
     * distance from the Sun alpha=2.
     * 
     * The distribution is given by:
     *
     * p(r) = (alpha+1) * r^alpha / (H^(alpha+1) - L^(alpha+1)),
     *
     * where the distribution is defined for L <= r <= H.
     *
     * @param r Vector of values where to evaluate the PDF.
     * @param rmin Lower bound L of the distribution (L >= 0).
     * @param rmax Upper bound H of the distribution (H > L).
     * @param alpha Shape parameter of the distribution (alpha >= 0).
     *
     * @return Natural logarithm of the probability density at r.
     * @throws If the parameters of the PDF do not satisfy the conditions above.
     */
    real distance_prior_lpdf(vector r, real rmin, real rmax, real alpha) {
        if (rmin < 0 || rmax <= rmin || alpha < 0) {
            reject("The following conditions must be satisfied: rmin>=0, rmax>rmin, alpha>=0; found rmin = ", rmin, ", rmax = ", rmax, ", alpha = ", alpha);
        }
        if (min(r) < rmin || max(r) > rmax) {
            return negative_infinity();
        }
        return num_elements(r) * ( log(alpha+1) - log(pow(rmax, alpha+1) - pow(rmin, alpha+1)) ) + alpha*sum(log(r));
    }
}

data {
    real minDist;                    // Assumed minimum possible distance
    real maxDist;                    // Assumed maximum possible distance
    int<lower=0> N;                  // Number of stars in survey
    vector[N] obsPlx;                // List of observed parallaxes
    vector<lower=0>[N] errObsPlx;    // List of parallax errors
    vector[N] obsMag;                // List of observed apparent magnitudes
    vector<lower=0>[N] errObsMag;    // List of apparent magnitude errors
}

transformed data {
}

parameters {
    real meanAbsMag;                            // Mean of absolute magnitude distribution
    real<lower=0.01> sigmaAbsMag;               // Standard deviation of absolute magnitude distribution (lower bound of 0.01 to avoid sampling sigma=0)
    vector<lower=minDist, upper=maxDist>[N] r;  // True distance, bounded by assumed minimum and maximum values
    vector[N] absMag;                           // True absolute magnitudes
}

transformed parameters {
    vector[N] predictedMag = absMag + 5.0*log10(r) - 5.0;  // Predicted true apparent magnitude for each star (distance in pc)
    vector<lower=0>[N] plx;                                // True parallax (in units of mas)
    for (k in 1:N) {
        plx[k] = 1000.0/r[k];
    }
}

model {
    meanAbsMag ~ normal(5.5, 10.5);            // Hyperprior on mean absolute magnitude (broadly between -5 and 16)
    sigmaAbsMag ~ normal(0, 2) T[0,];          // Hyperprior on standard deviation of absolute magnitude (half-normal positive)
    absMag ~ normal(meanAbsMag, sigmaAbsMag);  // Model absolute magnitude distribution for single class of stars
    r ~ distance_prior(minDist, maxDist, 2.0); // Prior on distance distribution
    obsPlx ~ normal(plx, errObsPlx);           // Likelihood observed parallax
    obsMag ~ normal(predictedMag, errObsMag);  // Likelihood observed apparent magnitude
}

generated quantities {
}
