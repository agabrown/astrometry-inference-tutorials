"""
Provides plotting utilities for MCMC related data.

Anthony Brown Nov 2017 - Nov 2017
"""
import numpy as np

def convert_to_stdev_nan(logL):
    """
    Given a grid of log-likelihood values, convert them to cumulative
    standard deviation. This is useful for drawing contours from a
    grid of likelihoods.

    Code from astroML: https://github.com/astroML/astroML

    THIS version is robust to a logL array that contains NaNs.

    Parameters
    ----------

    logL : float array
        Input array of log-likelihood (or log-probability) values. Can contain NaNs.

    Returns
    -------

    Array of values that represent the cumulative probability for logL.
    """
    sigma = np.exp(logL)

    shape = sigma.shape
    sigma = sigma.ravel()

    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(sigma)[::-1]
    i_unsort = np.argsort(i_sort)

    sigma_sorted = sigma[i_sort]
    notnan = np.logical_not(np.isnan(sigma_sorted))
    sigma_cumsum = np.empty(sigma.size)
    sigma_cumsum[:] = np.nan
    sigma_cumsum[notnan] = sigma_sorted[notnan].cumsum()
    sigma_cumsum /= sigma_cumsum[-1]

    return sigma_cumsum[i_unsort].reshape(shape)
