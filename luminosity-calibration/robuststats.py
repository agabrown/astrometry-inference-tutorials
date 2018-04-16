"""
Provides various methods for robust estimates of simple statistics such as the mean and variance, which
in this case are estimated through the median and RSE.

Anthony G.A. Brown May 2015
<brown@strw.leidenuniv.nl>
"""

from numpy import median
from scipy.stats import scoreatpercentile

def rse(x):
  """
  Calculate the Robust Scatter Estimate for an array of values (see GAIA-C3-TN-ARI-HL-007).

  Parameters
  ----------

  x - Array of input values (can be of any dimension)

  Returns
  -------

  The Robust Scatter Estimate (RSE), defined as 0.390152 * (P90-P10), where P10 and P90 are the 10th and
  90th percentile of the distribution of x.
  """
  return 0.390152*(scoreatpercentile(x,90)-scoreatpercentile(x,10))

def robust_stats(x):
  """
  Provide robust statistics of the values in array x (which can be of any dimension).

  Parameters
  ----------

  x - input array (numpy array is assumed)

  Returns
  -------

  Dictionary {'median':median, 'rse':RSE, 'lowerq':lower quartile, 'upperq':upper quartile, 'min':minimum
  value, 'max':maximum value}
  """

  med=median(x)
  therse=rse(x)
  lowerq=scoreatpercentile(x,25)
  upperq=scoreatpercentile(x,75)

  return {'median':med, 'rse':therse, 'lowerq':lowerq, 'upperq':upperq, 'min':x.min(), 'max':x.max(),
      'ndata':x.size}
