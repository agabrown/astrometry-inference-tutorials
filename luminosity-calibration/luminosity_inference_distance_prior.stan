/*
 * Stan model for the problem of inferring the mean absolute magnitude (and its variance) from a
 * parallax survey in which apparent magnitudes and parallaxes of stars are observed and the errors
 * on the measurements are reported. The model assumes a uniform distribution of stars in space
 * (constant density), with their distances limited by an upper and lower bound.
 *
 * The survey is modelled as magnitude limited with a the limit applied to the *true* apparent
 * magnitudes. This is only an approximation of what happens in the corresponding simulated parallax
 * survey (where the limit is applied to the *observed* apparent magnitudes).
 *
 * See chapter 16 in 'Astrometry for Astrophysics', 2013, edited by W.F. van Altena (Cambridge
 * University Press) for a discussion of the volume limited case.
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
     * @param r Value or r where to evaluate the PDF.
     * @param rmin Lower bound L of the distribution (L >= 0).
     * @param rmax Upper bound H of the distribution (H > L).
     * @param alpha Shape parameter of the distribution (alpha >= 0).
     *
     * @return Natural logarithm of the probability density at r.
     * @throws If the parameters of the PDF do not satisfy the conditions above.
     */
    real distance_prior_lpdf(real r, real rmin, real rmax, real alpha) {
        if (rmin < 0 || rmax <= rmin || alpha < 0) {
            reject("The following conditions must be satisfied: rmin>=0, rmax>rmin, alpha>=0; found rmin = ", rmin, ", rmax = ", rmax, ", alpha = ", alpha);
        }
        if (r < rmin || r > rmax) {
            return negative_infinity();
        }
        return log(alpha+1) - log(pow(rmax, alpha+1) - pow(rmin, alpha+1)) + alpha*log(r);
    }

    /**
     * Provide the PDF of the joint prior on distance and true absolute magnitude. This distribution
     * is the product of the distance prior above and the normal distribution for the true absolute
     * magnitude
     *
     * @param rM Vector representing distance and true absolute magnitude (in that order!).
     * @param rmin Lower bound of the distance distribution
     * @param rmax Upper bound of the distance distribution
     * @param alpha Shape parameter of the distance distribution (alpha >=0)
     * @param mu Mean of the normal distribution for the true absolute magnitude
     * @param sigma Standard deviation of the normal distribution for the true absolute magnitude
     * @param mlim The suvey limit
     */
    real distance_absmag_joint_prior_lpdf(vector rM, real rmin, real rmax, real alpha, real mu, real sigma, real mlim) {
        if (rM[2]+5*log10(rM[1])-5 > mlim)
            return negative_infinity();
        return distance_prior_lpdf(rM[1] | rmin, rmax, alpha) + normal_lpdf(rM[2] | mu, sigma);
    }

    /**
     * Provides the integrand for the approximate normalization of the posterior in the presence of
     * a magnitude limit. This is the marginal distribution for the distance:
     *
     * p(r) = 3*r^2/A * Phi((mlim-(mu+5*log10(r)-5))/sigma)
     *
     * where Phi is the cumulative standard normal distribution and the error in the observed
     * apparent magnitudes has been ignored.
     *
     * @param r Distance for which to evaluate the function above.
     * @param sp Array of survey parameters in following order:
     *      sp[1] -> rmin : Lower bound of the distance distribution
     *      sp[2] -> rmax : Upper bound of the distance distribution
     *      sp[3] -> alpha : Shape parameter of the distance distribution (alpha >=0)
     *      sp[4] -> mu : Mean of the normal distribution for the true absolute magnitude
     *      sp[5] -> sigma : Standard deviation of the normal distribution for the true absolute magnitude
     *      sp[6] -> mlim : The suvey limit
     */
    real marginal_distance_pdf(real r, real[] sp) {
        real A = pow(sp[2], sp[3]+1) - pow(sp[1], sp[3]+1);
        return exp( log(sp[3]+1) - log(A) + sp[3]*log(r) + normal_lcdf(sp[6] - sp[4] - 5*log10(r) + 5| 0, sp[5]) );
    }

    /**
     * Provides the ODE to solve for the integral of the marginal distribution in distance which
     * results in the normalization of the posterior in the presence of a magnitude limit.
     *
     * NOTE: this way of calculating the numerical integral needed for the normalization is a bit of
     * a hack. Stan has no numerical quadrature tools, but does provide a means of integrating ODEs.
     * Comparison of the results obtained with Wolfram Alpha and Python/scipy/numpy shows that the
     * result is accurate to 5 decimal places for the default tolerance and maximum-steps parameters
     * for the Stan integrate_ode_rk45() function.
     *
     * @param r The distance for which to evaluate the ODE equations ('time' variable in Stan docs).
     * @param y The state variable for the ODE (will contain the integral in the end).
     * @param theta Parameter array which has to be provided in this order:
     *      theta[1] -> rmin : Lower bound of the distance distribution
     *      theta[2] -> rmax : Upper bound of the distance distribution
     *      theta[3] -> alpha : Shape parameter of the distance distribution (alpha >=0)
     *      theta[4] -> mu : Mean of the normal distribution for the true absolute magnitude
     *      theta[5] -> sigma : Standard deviation of the normal distribution for the true absolute magnitude
     *      theta[6] -> mlim : The suvey limit
     * @param x_r Not used here, only present to satisfy interface requirements (see Stan docs).
     * @param x_i Not used here, only present to satisfy interface requirements (see Stan docs).
     */
    real[] normalization_ode(real r, real[] y, real[] theta, real[] x_r, int[] x_i) {
        real dydr[1];
        dydr[1] = marginal_distance_pdf(r, theta);
        return dydr;
    }

}

data {
    real minDist;                         // Assumed minimum possible distance
    real maxDist;                         // Assumed maximum possible distance
    real surveyLimit;                     // Apparent magnitude limit of survey
    int<lower=0> N;                       // Number of stars in survey
    vector[N] obsPlx;                     // List of observed parallaxes
    vector<lower=0>[N] errObsPlx;         // List of parallax errors
    vector<upper=surveyLimit>[N] obsMag;  // List of observed apparent magnitudes
    vector<lower=0>[N] errObsMag;         // List of apparent magnitude errors
}

transformed data {
    real x_r[0];                          // For ODE integration below, not used so leave empty
    int x_i[0];
    real rs[1];
    rs[1] = maxDist;
}

parameters {
    real meanAbsMag;                            // Mean of absolute magnitude distribution
    real<lower=0.01> sigmaAbsMag;               // Standard deviation of absolute magnitude distribution (lower bound of 0.01 to avoid sampling sigma=0)
    vector<lower=minDist, upper=maxDist>[N] r;  // True distance, bounded by assumed minimum and maximum values
    vector[N] absMag;                           // True absolute magnitudes. NOTE: this needs to be initialized by the user explicitly to avoid
                                                // that Stan chooses initial values of the absolute magnitudes that lead to apparent magnitudes beyond
                                                // the survey limit (without more information Stan will intialize with values between -2 and 2 and in
                                                // distant samples this can lead to choices that cannot occur in the survey). As initial values choose
                                                // something below surveyLimit-5*log10(maxDist)+5.
}

transformed parameters {
    vector[N] predictedMag = absMag + 5.0*log10(r) - 5.0;  // Predicted true apparent magnitude for each star (distance in pc)
    vector<lower=0>[N] plx;                                // True parallax (in units of mas)
    matrix[2,N] r_absMag;                                  // Distance and true absolute magnitude latent variables.
    for (k in 1:N) {
        r_absMag[1,k] = r[k];
        r_absMag[2,k] = absMag[k];
        plx[k] = 1000.0/r[k];
    }
}

model {
    real theta[6];
    real y0[1];
    real yhat[1,1];

    theta[1] = minDist;
    theta[2] = maxDist;
    theta[3] = 2.0;
    theta[6] = surveyLimit;

    meanAbsMag ~ normal(5.5, 10.5);            // Hyperprior on mean absolute magnitude (broadly between -5 and 16)
    sigmaAbsMag ~ normal(0, 2);          // Hyperprior on standard deviation of absolute magnitude (truncated below zero)

    // Joint distribution of r and M_true
    for (k in 1:N) {
        r_absMag[:,k] ~ distance_absmag_joint_prior(minDist, maxDist, 2.0, meanAbsMag, sigmaAbsMag, surveyLimit);
    }
    
    obsPlx ~ normal(plx, errObsPlx);           // Likelihood observed parallax
    obsMag ~ normal(predictedMag, errObsMag);  // Likelihood for observed apparent magnitude
    
    // Normalization of posterior to account for magnitude limit
    theta[4] = meanAbsMag;
    theta[5] = sigmaAbsMag;
    y0[1] = marginal_distance_pdf(minDist, theta);
    yhat = integrate_ode_rk45(normalization_ode, y0, minDist, rs, theta, x_r, x_i);//, 1.0e-10, 1.0e-10, 1e6);
    target += -N*log(yhat[1,1]);
}

generated quantities {
}
