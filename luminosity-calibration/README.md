# luminosity-calibration-tutorial
Tutorial on use of parallaxes in luminosity calibration; accompanies the Gaia DR2 paper on the use
of parallaxes (Luri et al. 2018 in prep).

_Code and notebooks by Anthony Brown <brown@strw.leidenuniv.nl>_

This tutorial expands on the luminosity inference examples discussed in chapter 16 in 'Astrometry for
Astrophysics', 2013, edited by W.F. van Altena (Cambridge University Press). It can also be seen as a
simplified version of [Hakwins et al.
2017](https://ui.adsabs.harvard.edu/#abs/2017MNRAS.471..722H/abstract).

__Any errors in this tutorial are of course mine. Comments are welcome__

## Prerequisites

This tutorial depends on the following Python packages (besides numpy, scipy, matplotlib):
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
