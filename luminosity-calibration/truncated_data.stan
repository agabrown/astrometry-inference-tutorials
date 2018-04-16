/*
 * Stan model for inferring the mean and variance of a normal distribution for which only a truncated set
 * of samples is available. That is samples x_i were generated from the normal distribution and then all
 * samples with x_i > xlim were discarded (without reporting the number of discarded samples).
 *
 * Anthony G.A. Brown Nov 2017 - Jan 2018
 * <brown@strw.leidenuniv.nl>
 */

functions {
}

data {
    real xlimit;    // Sample values are truncated above this (known) limit.
    int<lower=0> K; // Number of samples available after truncation.
    vector[K] x;    // List of samples.
}

transformed data {
}

parameters {
    real mu;             // Mean of normal distribution.
    real<lower=0> sigma; // Standard deviation of normal distribution.
}

transformed parameters {
}

model {
    mu ~ normal(0, 5);                     // Prior on mu
    sigma ~ normal(0, 3);                  // Prior on sigma
    for (k in 1:K)
        x[k] ~ normal(mu, sigma) T[,xlimit]; // Likelihood for x[k], truncated at xlimit.
}

generated quantities {
}
