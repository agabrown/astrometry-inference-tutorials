# luminosity-calibration-tutorial
Tutorial on use of parallaxes in luminosity calibration; accompanies the Gaia DR2 paper on the use
of parallaxes (Luri et al. 2018, A&A Special Issue for Gaia DR2).

_Code and notebooks in this tutorial by Anthony G.A. Brown <brown@strw.leidenuniv.nl>_

This tutorial expands on the luminosity inference example discussed in chapter 16 in 'Astrometry for
Astrophysics', 2013, edited by W.F. van Altena (Cambridge University Press). It can also be seen as
a simplified version of [Hakwins et al.
2017](https://ui.adsabs.harvard.edu/#abs/2017MNRAS.471..722H/abstract).

__Any errors in this tutorial are of course mine. Comments are welcome__

## Notebooks

The following python notebooks are included with this tutorial:
* [Luminosity calibration with a distance prior.](./Luminosity_Inference_DistPrior.ipynb) This is
  the main tutorial explaining how to infer from a set of measured parallaxes and apparent
  magnitudes of stars, their mean absolute magnitude and the standard deviation around the mean. In
  particular this tutorial is intended to demonstrate that negative and low-quality parallaxes can
  be used without any problem in the analysis. The four notebooks below provide supporting material.
* [Simulated parallax surveys.](./Parallax_survey_simulation.ipynb) This notebook explains in more
  detail how the simulated parallax survey was created which is used in the luminosity calibration
  tutorial. It also demonstrates aspects of the effect of selecting stars with high-quality
  (introducing a bias on the inferred absolute magnitude) and of the effect of the survey magnitude
  limit (including the mathematics of prediction the apparent magnitude distribution).
* [Handling data truncation.](./Handling_Data_Truncation.ipynb) A simplified demonstration of how to
  handle data truncation (i.e. selection functions) in a Bayesian analysis.  It is meant to provide
  more insight into the treatment of the survey magnitude limit in the luminosity calibration
  problem.
* [What's with the negative parallaxes?](./DemoNegativeParallax.ipynb) A simplified explanation on
  how negative parallaxes arise in astrometric surveys. The aim is to show why negative parallaxes
  are perfectly legitimate outcomes of the astrometric measurement process.
* [Distribution of quantities calculated from parallax data.](./Parallax_related_quantities.ipynb)
  Provides complementary demonstrations to section 3 of Luri et al. (2018) as to why calculating
  distances or tangential velocities of stars from parallaxes through the naive inversion of the
  latter leads to problems. This is done by contrasting the naive approach that taken in
  [Bailer-Jones (2015)](https://ui.adsabs.harvard.edu/#abs/2015PASP..127..994B/abstract).

## Prerequisites

This tutorial depends on the following Python packages (besides [scipy](https://www.scipy.org/),
[numpy](http://numpy.org/), [matplotlib](https://matplotlib.org/)):
* [PyStan](http://mc-stan.org/users/interfaces/pystan.html)
* [scikit-learn](http://scikit-learn.org)
* [astropy](http://www.astropy.org/index.html)
* [corner.py](https://github.com/dfm/corner.py)
* [PyGaia](https://github.com/agabrown/PyGaia)
* [daft](https://github.com/dfm/daft) (Optional)

## TGAS parallax error model

The [parallax survey simulation module](./parallaxsurveys.py) contains a Hipparcos and a TGAS error
model. For TGAS the distribution of parallax errors is assumed to be independent of magnitude. The
probability density of the parallax error distribution is obtained as a kernel density estimate from all
parallax errors in the TGAS catalogue. One can use the
[fit-tgas-parallax-errors.py](./fit-tgas-parallax-errors.py) Python code to try to fit this distribution
with an inverse Gamma distribution. The necessary data can be extracted from the [Gaia
archive](https://archives.esac.esa.int/gaia) with the following query (store as
"TGAS-allPlxErrorsVsGmag.fits", _note the FITS format_):

`select parallax_error, phot_g_mean_mag from gaiadr1.tgas_source`

The survey simulation code makes use directly of the KDE mentioned above.

## Subfolder contents

pgm : Probabilistic graphical models for the tutorial
